package starbeast3.operators;



import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Executors;


import beast.app.BeastMCMC;
import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.ParallelMCMC;
import beast.core.State;
import beast.core.StateNode;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;
import beast.core.util.Log;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.operators.AdaptableVarianceMultivariateNormalOperator;
import beast.evolution.sitemodel.SiteModelInterface;
import beast.evolution.substitutionmodel.Frequencies;
import beast.evolution.substitutionmodel.SubstitutionModel;
import beast.util.Transform;
import starbeast3.evolution.branchratemodel.BranchRateModelSB3;

@Description("Run MCMC on different treelikelihood parts of the model in parallel before combining them in a single Gibbs move")
public class ParallelMCMCRealParameterOperator extends MultiStepOperator {
	

	final public Input<CompoundDistribution> distributionInput = new Input<>("distribution", 
			"compound distribution of all likelihoods",
			Validate.REQUIRED);
    
    
	@Override
	public void initAndValidate() {
	    otherState = otherStateInput.get();
		List<Distribution> distributions = distributionInput.get().pDistributions.get();
		 
		nrOfThreads = maxNrOfThreadsInput.get() > 0 ? 
				Math.min(BeastMCMC.m_nThreads, maxNrOfThreadsInput.get()) : 
				BeastMCMC.m_nThreads;
		if (nrOfThreads > distributions.size()) {
			nrOfThreads = distributions.size();
		}
	    exec = Executors.newFixedThreadPool(nrOfThreads);
	    mcmcs = new ArrayList<>();
	    


	    // determine statenode tabu list: all stateNodes that are shared among threads
	    // 1. determine potential state nodes per thread first
	    Set<StateNode>[] stateNodes = new Set[nrOfThreads];
	    int start = 0;
	    for (int i = 0; i < nrOfThreads; i++) {
	    	int end = (i + 1) * distributions.size() / nrOfThreads;
	    	stateNodes[i] = new HashSet<>();
			for (Distribution d : distributions.subList(start, end)) {
				getRealParameterStateNodes(d, otherState.stateNodeInput.get(), stateNodes[i]);
			}
	    	start = end;
	    }
	    
	    // 2. determine overlaps in potential state nodes per thread, and count dimensions
	    Set<StateNode> tabu = new HashSet<>();
	    int totalDim = 0;
	    if (nrOfThreads == 1) {
	    	for (StateNode s : stateNodes[0]) {
	    		totalDim += s.getDimension();
	    	}
	    }else {
		    for (int i = 0; i < nrOfThreads; i++) {
		    	for (int j = i + 1; j < nrOfThreads; j++) {
		    		for (StateNode s : stateNodes[i]) {
		    			if (stateNodes[j].contains(s)) {
		    				tabu.add(s);
		    			}else {
		    				totalDim += s.getDimension();
		    			}
		    		}
		    	}
		    }
	    }

	    
	    
	    // 3. Determine chain length
	    if (chainLengthInput.get() != null) {
	    	chainLength = chainLengthInput.get();
	    }else if(coverageInput.get() != null) {
	    	chainLength = (long) (coverageInput.get() * totalDim);
	    	Log.warning(this.getID() + ": dimensional chain length: " + chainLength);
	    }else {
	    	throw new IllegalArgumentException("Please provide either 'chainLength' or 'coverageInput' but not both");
	    }
	    
	    
	    
	    // Create mcmc objects
	    start = 0;
	    for (int i = 0; i < nrOfThreads; i++) {
	    	int end = (i + 1) * distributions.size() / nrOfThreads;
	    	mcmcs.add(createParallelMCMC(distributions.subList(start, end), (int)(1.0*chainLength / this.nrOfThreads), tabu));
	    	start = end;
	    }
	    
	    for (ParallelMCMC pMCMC : mcmcs) {
	    	pMCMC.setOtherState(otherState);
	    }
	    
	    
	  
	    
	    // 1 thread and 1 chain mode
	    if (learningInput.get() || (this.mcmcs.size() == 1 && chainLength == 1)) {
	    	ParallelMCMC mcmc;
	    	if (this.mcmcs.size() == 1 && chainLength == 1) mcmc = this.mcmcs.get(0);
	    	else mcmc = createParallelMCMC(distributions, chainLength, new HashSet<>());
	    	this.singleStepOperators = mcmc.operatorsInput.get();
	    }
	    
	    
	    super.initAndValidate();
	}

	@Override
	public int stepCount() {
		return mcmcs.size();
	}

	
	private ParallelMCMC createParallelMCMC(List<Distribution> distributions, long chainLength, Set<StateNode> tabu) {
		List<Distribution> distrs = new ArrayList<>();
		List<StateNode> stateNodes = new ArrayList<>();
		List<Operator> operators = new ArrayList<>();
		
		
		for (Distribution d : distributions) {
			distrs.add(d);
			
			Set<StateNode> stateNodeSet = new HashSet<>();
			getRealParameterStateNodes(d, otherState.stateNodeInput.get(), stateNodeSet);
			stateNodeSet.removeAll(tabu);
			if (stateNodeSet.size() > 0) {
				stateNodes.addAll(stateNodeSet);
				
				
				
				List<Distribution> priorsList = new ArrayList<>();
				getRealParameterPriors(stateNodeSet, priorsList);
				distrs.addAll(priorsList);			
				
				int dim = 0;
				for (StateNode s : stateNodeSet) {
					dim += s.getDimension();
				}
				
				List<Transform> transformations = new ArrayList<>();
				for (StateNode s : stateNodeSet) {
					Log.warning("Adding " + s.getID());
					Transform f;
					// TODO: check priors instead of ID to determine whether it is a
					// scale parameter
					// location parameter
					// simplex parameter
					if (s.getID().startsWith("freq")) {
						f = new Transform.LogConstrainedSumTransform(s, 1.0);
					} else {
						f = new Transform.LogTransform(s);
					}
					transformations.add(f);
				}
				
				AdaptableVarianceMultivariateNormalOperator AVMNOperator = new AdaptableVarianceMultivariateNormalOperator();
				AVMNOperator.initByName("weight", 1.0, 
						"coefficient", 1.0, 
						"scaleFactor", 1.0, 
						"beta", 0.05, 
						"every", 1,
						"initial", 200 * dim, 
						"burnin", 100 * dim, 
						"transformations", transformations);
				operators.add(AVMNOperator);
			}
		}
		
		
		CompoundDistribution sampleDistr = new CompoundDistribution();
		sampleDistr.initByName("distribution", distrs);
		
		State state = new State();
		state.initByName("stateNode", stateNodes);
		
		
		// Learn the chain length
		int nregression = this.nrOfThreads > 1 ? nregressionInput.get() : 0;

		
		ParallelMCMC mcmc = new ParallelMCMC();
		mcmc.initByName("state", state, "operator", operators, "distribution", sampleDistr, "chainLength", chainLength, "robust", false, "nregression", nregression);
		return mcmc;
		
	}

	private void getRealParameterPriors(Set<StateNode> stateNodeList, List<Distribution> priorsList) {
		for (StateNode sn : stateNodeList) {
			for (BEASTInterface o : sn.getOutputs()) {
				if (o instanceof Distribution) {
					for (BEASTInterface o2 : o.getOutputs()) {
						if (o2.getID() != null && o2.getID().equals("prior")) {
							priorsList.add((Distribution) o);
						}
					}
				}				
			}
		}		
	}

	public static void getRealParameterStateNodes(BEASTInterface d, List<StateNode> otherStateNodes, Set<StateNode> stateNodes) {
		for (Object o : d.listActiveBEASTObjects()) {
			if (o instanceof StateNode && otherStateNodes.contains(o) && o instanceof RealParameter) {
				// must have subst, or site model in its outputs, or be mean clock rate
				for (BEASTInterface o2 : ((BEASTInterface)o).getOutputs()) {
					if (o2 instanceof SubstitutionModel ||
						o2 instanceof SiteModelInterface ||
						o2 instanceof Frequencies ||
						(o2 instanceof BranchRateModel && !(o2 instanceof BranchRateModelSB3) && o2.getInput("clock.rate").get() == o)) {
						stateNodes.add((StateNode) o);
						break;
					}
				}
			} else if (o instanceof BEASTInterface) {
				getRealParameterStateNodes((BEASTInterface) o, otherStateNodes, stateNodes);				
			}			
 		}
	}

	@Override
	public double proposal() {
		return super.proposal();
	}







}
