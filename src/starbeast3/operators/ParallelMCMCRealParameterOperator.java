package starbeast3.operators;



import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Executors;


import beastfx.app.beast.BeastMCMC;
import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.CompoundRealParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.CompoundDistribution;
import beast.base.core.Log;
import beast.base.core.ProgramStatus;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.operator.kernel.AdaptableVarianceMultivariateNormalOperator;
import beast.base.inference.operator.DeltaExchangeOperator;
import beast.base.evolution.operator.ScaleOperator;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.substitutionmodel.Frequencies;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.inference.operator.kernel.Transform;
import starbeast3.core.ParallelMCMC;
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
				Math.min(ProgramStatus.m_nThreads, maxNrOfThreadsInput.get()) : 
				ProgramStatus.m_nThreads;
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
	    
	    
		 // Learn the chain length
		int nregression = this.nrOfThreads > 1 ? nregressionInput.get() : 0;
	    
	    
	    // Create mcmc objects
	    start = 0;
	    for (int i = 0; i < nrOfThreads; i++) {
	    	int end = (i + 1) * distributions.size() / nrOfThreads;
	    	mcmcs.add(createParallelMCMC(distributions.subList(start, end), (int)(1.0*chainLength / this.nrOfThreads), tabu, otherState, nregression));
	    	start = end;
	    }
	    
	    for (ParallelMCMC pMCMC : mcmcs) {
	    	pMCMC.setOtherState(otherState);
	    }
	    
	    
	  
	    
	    // 1 thread and 1 chain mode
	    if (learningInput.get() || (this.mcmcs.size() == 1 && chainLength == 1)) {
	    	ParallelMCMC mcmc;
	    	if (this.mcmcs.size() == 1 && chainLength == 1) mcmc = this.mcmcs.get(0);
	    	else mcmc = createParallelMCMC(distributions, chainLength, new HashSet<>(), otherState, nregression);
	    	this.singleStepOperators = mcmc.operatorsInput.get();
	    }
	    
	    
	    super.initAndValidate();
	}

	@Override
	public int stepCount() {
		return mcmcs.size();
	}

	
	public static ParallelMCMC createParallelMCMC(List<Distribution> distributions, long chainLength, Set<StateNode> tabu, State otherState, int nregression) {
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
				
				
				Set<Distribution> priorsSet = new HashSet<>();
				getParameterPriors(stateNodeSet, priorsSet);
				distrs.addAll(priorsSet);			
				
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
					
					//if (s.getID().endsWith("nuclear_genome")) {
						//continue;
					//}
					
					if (s.getID().startsWith("freq")) {
						
						DeltaExchangeOperator op = new DeltaExchangeOperator();
						op.initByName("parameter", s, "delta", 0.2, "weight", 2.0);
						operators.add(op);
						f = new Transform.LogConstrainedSumTransform(s, 1.0);
						
					} else {
						
	
						ScaleOperator op = new ScaleOperator();
						op.initByName("parameter", s, "scaleFactor", 0.5, "weight", 1.0);
						operators.add(op);
						f = new Transform.LogTransform(s);
					}
					
					transformations.add(f);
					
				}
				
				
				
				
				if (!transformations.isEmpty()) {
					AdaptableVarianceMultivariateNormalOperator AVMNOperator = new AdaptableVarianceMultivariateNormalOperator();
					AVMNOperator.initByName(
							"weight", 2.0, 
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
		}
		

		
		CompoundDistribution sampleDistr = new CompoundDistribution();
		sampleDistr.initByName("distribution", distrs);
		
		State state = new State();
		state.initByName("stateNode", stateNodes);
		
		
		

		
		ParallelMCMC mcmc = new ParallelMCMC();
		mcmc.initByName("state", state, "operator", operators, "distribution", sampleDistr, "chainLength", chainLength, "robust", false, "nregression", nregression);
		return mcmc;
		
	}




	/**
	 * Get real parameter state nodes, if they are part of the likelihood
	 * @param d
	 * @param otherStateNodes
	 * @param stateNodes
	 */
	public static void getRealParameterStateNodes(BEASTInterface d, List<StateNode> otherStateNodes, Set<StateNode> stateNodes) {
		for (Object o : d.listActiveBEASTObjects()) {
			if (o instanceof StateNode && otherStateNodes.contains(o) && o instanceof RealParameter) {
				
				StateNode state = (StateNode)o;
				if (!stateNodes.contains(state) && parameterIsPartOfLikelihood((RealParameter)o)) {
					stateNodes.add(state);
				}
				
			} else if (o instanceof BEASTInterface) {
				getRealParameterStateNodes((BEASTInterface) o, otherStateNodes, stateNodes);				
			}			
 		}
	}
	
	
	
	public static boolean parameterIsPartOfLikelihood(RealParameter o) {
		
		// Part of likelihood?
		for (BEASTInterface o2 : ((BEASTInterface)o).getOutputs()) {
			if (o2 instanceof SubstitutionModel ||
				o2 instanceof SiteModelInterface ||
				o2 instanceof Frequencies ||
				(o2 instanceof BranchRateModel && !(o2 instanceof BranchRateModelSB3) && o2.getInput("clock.rate").get() == o)) {
				
				return true;
			}
		}
		
		
		// Is it part of a compound parameter that is part of the likelihood?
		for (BEASTInterface o2 : ((BEASTInterface)o).getOutputs()) {
			if (o2 instanceof CompoundRealParameter) {
				if (parameterIsPartOfLikelihood((RealParameter)o2)) {
					return true;
				}
			}
		}
		
		return false;
	}
	

	@Override
	public double proposal() {
		return super.proposal();
	}



}




