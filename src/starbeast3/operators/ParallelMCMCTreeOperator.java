package starbeast3.operators;




import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.concurrent.Executors;

import beast.base.core.BEASTInterface;
import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.inference.Distribution;
import beast.base.inference.MCMC;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.CompoundDistribution;
import beast.base.core.Log;
import beast.base.core.ProgramStatus;
import beast.base.evolution.operator.ScaleOperator;
import beast.base.evolution.operator.kernel.AdaptableVarianceMultivariateNormalOperator;
import beast.base.evolution.operator.kernel.BactrianScaleOperator;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.tree.Tree;
import beast.base.inference.operator.kernel.BactrianDeltaExchangeOperator;
import beast.base.inference.operator.kernel.BactrianIntervalOperator;
import beast.base.inference.operator.kernel.BactrianRandomWalkOperator;
import beast.base.inference.operator.kernel.Transform;
import starbeast3.core.ParallelMCMC;
import starbeast3.evolution.speciation.GeneTreeForSpeciesTreeDistribution;

@Description("Run MCMC on different gene tree parts of the model in parallel before combining them in a single Gibbs move")
public class ParallelMCMCTreeOperator extends MultiStepOperator {
	
    final public Input<Boolean> useBactrianOperatorsInput = new Input<>("bactrian", "flag to indicate that bactrian operators should be used where possible", true);
    final public Input<Boolean> includeRealParametersInput = new Input<>("includeRealParameters", "flag to include Real Parameters for each of the partitions in the analysis", true);
    final public Input<List<StateNode>> excludeInput = new Input<>("exclude", "parameters to ensure are not operated on", new ArrayList<>());
    final public Input<List<StateNode>> includeInput = new Input<>("include", "parameters to add to the state, in case they are not automatically detected", new ArrayList<>());
    
    final public Input<List<StateNode>> logsumInput = new Input<>("logsum", "these parameters, if operated on, will have their total sum respected by the avmn operator", new ArrayList<>());
    
    
    
	final public Input<List<ParallelDistSet>> distributionInput = new Input<>("distribution", 
			"Distribution on a tree conditinionally independent from all other distributions given the state of the rest"
			+ "of parameter space. ",
			new ArrayList<>());
    
	
	final public Input<Boolean> unthreadInput = new Input<>("unthread", "flag to convert ThreadedTreeLikelihood back into TreeLikelihood when this is called", false);
	  

	boolean sampleFromPrior = false;
    
    List<ParallelDistSet> distributions;
    

    
	@Override
	public void initAndValidate() {
		this.distributions = distributionInput.get();
		mcmcs = new ArrayList<>();
		
		

		
		// Tidy the distributions
		tidyDistributions(this.distributions);
		
		
		
		if (distributions.isEmpty()) {
			Log.warning("ParallelMCMCTreeOperator: Please provide at least one 'distribution'");
			return;
			
		}
		
	    otherState = otherStateInput.get();
	    otherState.initialise();
	    
	    
		
	  
		
		 
		nrOfThreads = maxNrOfThreadsInput.get() > 0 ?
				Math.min(ProgramStatus.m_nThreads, maxNrOfThreadsInput.get()) : 
				ProgramStatus.m_nThreads;
		nrOfThreads = Math.min(nrOfThreads, distributions.size());
		Log.warning("Running " + this.getID() + " with " + this.nrOfThreads + " threads");
		//System.exit(1);
	    exec = Executors.newFixedThreadPool(nrOfThreads);
	    
	    
	    
	    
	    // Sample from prior?
	    this.sampleFromPrior = this.sampleFromPriorEnabled();
	    
	    
	    
	    // Load balancing. Ensure a roughly equal distribution of site patterns across all threads
	    int totalDim = 0;
	    Collections.sort(distributions);
	    List<List<ParallelMCMCTreeOperatorTreeDistribution>> balancedDistributions = new ArrayList<>();
	    for (int i = 0; i < nrOfThreads; i++) balancedDistributions.add(new ArrayList<>());
	    int threadNum = 0;
	    for (int i = 0; i < distributions.size(); i ++) {
	    	ParallelDistSet d = distributions.get(i);
	    	totalDim += d.getTaxonCount();
	    	//System.out.println("patterns " + d.getNumberPatterns());
	    	balancedDistributions.get(threadNum).addAll(d.getDists());
	    	threadNum++;
	    	if (threadNum >= nrOfThreads) threadNum = 0;
	    }
	    
	    
	    
	    // Determine chain length
	    if (chainLengthInput.get() != null) {
	    	chainLength = chainLengthInput.get();
	    }else if(coverageInput.get() != null) {
	    	chainLength = (long) (coverageInput.get() * totalDim);
	    	Log.warning(this.getID() + ": dimensional chain length: " + chainLength);
	    }else {
	    	throw new IllegalArgumentException("Please provide either 'chainLength' or 'coverageInput' but not both");
	    }
	    
	    if (this.runtimeInput.get() <= 0 && this.nrOfThreads == 1) chainLength = 1;
	    
	    
	    // Ensure that there are no fixed or duplicated state nodes across the threads
	    List<StateNode> doNotInclude = getTabooStateNodes(balancedDistributions, otherState);
	    if (!doNotInclude.isEmpty()) {
	    	String tabooStr = "";
	    	for (StateNode state : doNotInclude) {
	    		tabooStr += "'" + state.getID() + "' ";
	    	}
	    	Log.warning("The following stateNodes will NOT be operated on by " + this.getClass().getSimpleName() + " because they are not part of the state, or they appear in more than one thread: " + tabooStr);
	    }
	    for (StateNode state : excludeInput.get()) {
	    	
	    	if (includeInput.get().contains(state)) {
	    		throw new IllegalArgumentException(state.getID() + " must not be on the 'exclude' and 'include' lists!");
	    	}
	    	
	    	if (!doNotInclude.contains(state)) {
	    		Log.warning("Excluding " + state.getID() + " from operator by request");
	    		doNotInclude.add(state);
	    	}
	    }
	    
	    
	    
	    // Create parallel MCMCs
	    for (int i = 0; i < nrOfThreads; i++) {
	    	mcmcs.add(createParallelMCMC(balancedDistributions.get(i), (int)(1.0*chainLength / this.nrOfThreads), doNotInclude));
	    }
	    Log.warning(this.getClass().getCanonicalName() + ": total chain length is " + chainLength);

	    
	    for (ParallelMCMC pMCMC : mcmcs) {
	    	pMCMC.setOtherState(otherState);
	    }
	    
	    
	    // 1 thread and 1 chain mode
	    if (learningInput.get() || (this.mcmcs.size() == 1 && chainLength == 1)) {
	    	ParallelMCMC mcmc;
	    	if (this.mcmcs.size() == 1 && chainLength == 1) mcmc = this.mcmcs.get(0);
	    	else mcmc = createParallelMCMC(balancedDistributions.get(0), chainLength, doNotInclude);
	    	this.singleStepOperators = mcmc.operatorsInput.get();
	    }
	    
	  
	    super.initAndValidate();
	  
	    
	    
	}
	




	/*
	 * If there is one tree which appears in more than 1 likelihood or prior, then merge the threads together
	 */
	private void tidyDistributions(List<ParallelDistSet> distributions) {
		
		
		// Dist 1
		for (int i = 0; i < distributions.size(); i++) {
			ParallelDistSet dist1 = distributions.get(i);
			List<Tree> treeSet1 = dist1.getTrees();
			
			// Dist 2
			for (int j = i+1; j < distributions.size(); j++) {
				ParallelDistSet dist2 = distributions.get(j);
				List<Tree> treeSet2 = dist2.getTrees();
				
				
				// If any of the trees are the same, then join dist1 with dist2
				for (Tree tree1 : treeSet1) {
					if (treeSet2.contains(tree1)) {
						
						
						// Update list and repeat
						ParallelDistSet dist3 = new ParallelDistSet(dist1, dist2);
						distributions.remove(j);
						distributions.remove(i);
						distributions.add(dist3);
						dist3.setID(dist1.getID() + "_" + dist2.getID());
						Log.warning("Merging " + dist1.getID() + " with " + dist2.getID() + " because they have the same tree (" + tree1.getID() + ")" );
						tidyDistributions(distributions);
						return;
						
					}
				}
				
				
				
			}
			
			
		}
		
	}



	@Override
	public int stepCount() {
		return mcmcs.size();
	}

	private ParallelMCMC createParallelMCMC(List<ParallelMCMCTreeOperatorTreeDistribution> distributions, long chainLength,  List<StateNode> doNotInclude) {
		List<Distribution> distrs = new ArrayList<>();
		List<StateNode> stateNodes = new ArrayList<>();
		List<Operator> operators = new ArrayList<>();
		for (ParallelMCMCTreeOperatorTreeDistribution d : distributions) {
			
			
			
			
			if (!distrs.contains(d.geneprior)) {
				Log.warning("Adding dist " + d.getGeneprior().getID());
				distrs.add(d.getGeneprior());
			}
			for (Distribution d2 : d.getOtherDists()) {
				if (!distrs.contains(d2)) {
					Log.warning("Adding dist " + d2.getID());
					distrs.add(d2);
				}
			}
			if (!sampleFromPrior && !distrs.contains(d.treelikelihood)) {
				Log.warning("Adding dist " + d.getTreelikelihood().getID());
				distrs.add(d.getTreelikelihood()); 
			}
			
			
			// Any other tree prior distributions?
			Set<Distribution> priorsSet = new HashSet<>();
			getParameterPriors(d.tree, priorsSet);
			for (Distribution dist : priorsSet) {
				if (!distrs.contains(dist)) {
					Log.warning("Adding dist " + dist.getID());
					distrs.add(dist);	
				}
			}
			
			
			//System.out.println("Adding: " + d.geneprior.getID() + " / " + d.treelikelihood.getID());
			
			
			// Include the tree?
			if (!doNotInclude.contains(d.tree)) {
				
				if (!stateNodes.contains(d.tree)) stateNodes.add(d.tree);
				
				
				// Uniform operator
				beast.base.evolution.operator.Uniform UniformOperator = new beast.base.evolution.operator.Uniform();
				UniformOperator.initByName("tree", d.tree, "weight", 30.0);
				operators.add(UniformOperator);
				
				if (useBactrianOperatorsInput.get()) {
					
					
					// Root scaler
					beast.base.evolution.operator.kernel.BactrianScaleOperator treeRootScaler = new beast.base.evolution.operator.kernel.BactrianScaleOperator();
					treeRootScaler.initByName("scaleFactor", 0.5, "tree", d.tree, "weight", 10.0, "rootOnly", true);
					operators.add(treeRootScaler);
					
					// Bactrian interval
					beast.base.evolution.operator.kernel.BactrianNodeOperator intervalOperator = new beast.base.evolution.operator.kernel.BactrianNodeOperator();
					intervalOperator.initByName("tree", d.tree, "weight", 10.0);
					operators.add(intervalOperator);
					
					// Subtree slide
					beast.base.evolution.operator.kernel.BactrianSubtreeSlide SubtreeSlide = new beast.base.evolution.operator.kernel.BactrianSubtreeSlide();
					SubtreeSlide.initByName("tree", d.tree, "weight", 10.0);
					operators.add(SubtreeSlide);
					
					
					// Adaptable operators for tree scaling
					beast.base.evolution.operator.AdaptableOperatorSampler adaptable = new beast.base.evolution.operator.AdaptableOperatorSampler();
					List<Operator> adaptOperators = new ArrayList<>();
					
					// Tree scaler
					beast.base.inference.operator.kernel.BactrianUpDownOperator treeScaler = new beast.base.inference.operator.kernel.BactrianUpDownOperator();
					treeScaler.initByName("scaleFactor", 0.5, "up", d.tree, "weight", 1.0);
					adaptOperators.add(treeScaler);

					
					// Subtree slide
					adaptOperators.add(SubtreeSlide);
					
					// Uniform
					adaptOperators.add(UniformOperator);
					
					adaptable.initByName("tree", d.tree, "weight", 100.0, "operator", adaptOperators);
					operators.add(adaptable);
					
					
					
				} else {
	//		        <operator id="treeScaler.t:$(n)" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:$(n)" weight="3.0"/>
					BactrianScaleOperator treeScaler = new BactrianScaleOperator();
					treeScaler.initByName("scaleFactor", 0.5, "tree", d.tree, "weight", 3.0);
					operators.add(treeScaler);
	//		    	<operator id="treeRootScaler.t:$(n)" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:$(n)" weight="3.0"/>
					ScaleOperator treeRootScaler = new ScaleOperator();
					treeRootScaler.initByName("scaleFactor", 0.5, "tree", d.tree, "weight", 3.0, "rootOnly", true);
					operators.add(treeRootScaler);
	
	//		    	<operator id="SubtreeSlide.t:$(n)" spec="SubtreeSlide" tree="@Tree.t:$(n)" weight="15.0"/>
					BactrianSubtreeSlide SubtreeSlide = new BactrianSubtreeSlide();
					SubtreeSlide.initByName("tree", d.tree, "weight", 15.0);
					operators.add(SubtreeSlide);
				}
				
				
	
	//		    <operator id="narrow.t:$(n)" spec="Exchange" tree="@Tree.t:$(n)" weight="15.0"/>
				beast.base.evolution.operator.Exchange narrow = new beast.base.evolution.operator.Exchange();
				narrow.initByName("tree", d.tree, "weight", 15.0);
				operators.add(narrow);
	//	    	<operator id="wide.t:$(n)" spec="Exchange" isNarrow="false" tree="@Tree.t:$(n)" weight="3.0"/>
				beast.base.evolution.operator.Exchange wide = new beast.base.evolution.operator.Exchange();
				wide.initByName("tree", d.tree, "isNarrow", false, "weight", 15.0);
				operators.add(wide);
	//		    <operator id="WilsonBalding.t:$(n)" spec="WilsonBalding" tree="@Tree.t:$(n)" weight="3.0"/>
				beast.base.evolution.operator.WilsonBalding WilsonBalding = new beast.base.evolution.operator.WilsonBalding();
				WilsonBalding.initByName("tree", d.tree, "weight", 15.0);
				operators.add(WilsonBalding);
			}
				
				
			if (includeRealParametersInput.get()) {
				
				// Real parameters
				Set<StateNode> stateNodeList = new HashSet<>();
				stateNodeList.addAll(includeInput.get());
				stateNodeList.addAll(d.getInclude());
				ParallelMCMCRealParameterOperator.getRealParameterStateNodes(d.treelikelihood, otherState.stateNodeInput.get(), stateNodeList);
				
				
				
				// Remove forbidden statenodes
				Set<StateNode> stateNodeList2 = new HashSet<>();
				for (StateNode s : stateNodeList) {
					if (!doNotInclude.contains(s)) {
						stateNodeList2.add(s);
					}
				}
				stateNodeList = stateNodeList2;
				
				// Add non-duplicates
				for (StateNode s : stateNodeList) {
					if (!stateNodes.contains(s)) stateNodes.add(s);
				}
				
				
				/**
				 * Add priors to dist
				 */
				priorsSet = new HashSet<>();
				getParameterPriors(stateNodeList, priorsSet);
				for (Distribution dist : priorsSet) {
					if (!distrs.contains(dist)) {
						Log.warning("Adding dist " + dist.getID());
						distrs.add(dist);	
					}
				}
				
	
				
				int dim = 0;
				for (StateNode s : stateNodeList) {
					dim += s.getDimension();
				}
				
				List<Transform> transformations = new ArrayList<>();
				Transform f = null;
				
				// Add the tree?
				if (!doNotInclude.contains(d.tree)) {
					f = new Transform.LogTransform(d.tree);
					transformations.add(f);
					Log.warning("Adding " + d.tree.getID());
				}
				
				for (StateNode s : stateNodeList) {
					
					
					if (!(s instanceof RealParameter)) continue;
					RealParameter rp = (RealParameter)s;
					
					
					// TODO: check priors instead of ID to determine whether it is a
					// scale parameter
					// location parameter
					// simplex parameter
					if (s.getID().startsWith("freq") || logsumInput.get().contains(s)) {
						
						BactrianDeltaExchangeOperator op = new BactrianDeltaExchangeOperator();
						op.initByName("parameter", s, "delta", 0.2, "weight", 0.5);
						operators.add(op);
						double sum = 0;
						for (int i = 0; i < s.getDimension(); i++) {
							sum += s.getArrayValue(i);
						}
						f = new Transform.LogConstrainedSumTransform(s, sum);
						
						
						Log.warning("Adding logsum(" + sum + ") for " + s.getID());
						
					} 
					
					// Logit
					else if (rp.getLower() != Double.NEGATIVE_INFINITY && rp.getUpper() != Double.POSITIVE_INFINITY) {
						
						BactrianIntervalOperator op = new BactrianIntervalOperator();
						op.initByName("parameter", s, "scaleFactor", 0.5, "weight", 0.5);
						operators.add(op);
						
						List<Function> l = new ArrayList<>();
						l.add(s);
						f = new Transform.LogitTransform(l);
						
						Log.warning("Adding logit(" + rp.getLower() + "," + rp.getUpper() + ") " + s.getID());
						
						
						
						
					}
					
					
					
					
					// Scale
					else if (rp.getLower() >= 0)  {
						
	
						BactrianScaleOperator op = new BactrianScaleOperator();
						op.initByName("parameter", s, "scaleFactor", 0.5, "weight", 0.5);
						operators.add(op);
						f = new Transform.LogTransform(s);
						
						Log.warning("Adding scale " + s.getID());
					}
					
					
					// No transform
					else {
						
						
						BactrianRandomWalkOperator op = new BactrianRandomWalkOperator();
						op.initByName("parameter", s, "scaleFactor", 0.5, "weight", 0.5);
						operators.add(op);
						
						
						List<Function> l = new ArrayList<>();
						l.add(s);
						f = new Transform.NoTransform(l);
						Log.warning("Adding notransform " + s.getID());
						
						
					}
					
					if (f != null) transformations.add(f);
					
				}
				
				if (!transformations.isEmpty()) {
					AdaptableVarianceMultivariateNormalOperator AVMNOperator = new AdaptableVarianceMultivariateNormalOperator();
					AVMNOperator.initByName(
							"weight", 5.0, 
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
		
		
		List<String> stateNodeIDs = new ArrayList<>();
		for (StateNode stateNode : stateNodes) {
			if (stateNodeIDs.contains(stateNode.getID())) {
				Log.warning("Duplicate statenode : " + stateNode.getID());
			}
			
			stateNodeIDs.add(stateNode.getID());
			
		}
		Collections.sort(stateNodeIDs);
		Log.info("ParallelMCMC State: " + stateNodeIDs);
		
		State state = new State();
		state.initByName("stateNode", stateNodes);
		
		
		// Learn the chain length
		int nregression = this.nrOfThreads > 1 ? nregressionInput.get() : 0;

		ParallelMCMC mcmc = new ParallelMCMC();
		mcmc.initByName("state", state, "operator", operators, "distribution", sampleDistr, "chainLength", chainLength, "robust", false, "nregression", nregression);
		return mcmc;
	}
	
	
	
	
	
	/**
	 * Get list of state nodes which appear in more than 1 family, or are not part of the state
	 * These should not be operated on
	 * @param dists
	 * @return
	 */
	public List<StateNode> getTabooStateNodes(List<List<ParallelMCMCTreeOperatorTreeDistribution>> dists, State mainState){
		
		List<StateNode> taboo = new ArrayList<>();
		List<String> stateNodeIds = new ArrayList<>();
		//Log.warning(" States : " + mainState.toString());
		for (int i = 0; i < mainState.getNrOfStateNodes(); i ++) {
			
			
			
			StateNode state = mainState.getStateNode(i);
			stateNodeIds.add(state.getID());
			
			
			// bModelTest: do not add the rates
			if (state instanceof RealParameter) {
				for (BEASTInterface o : state.getOutputs()) {
					if (o.getClass().getCanonicalName().equals("bmodeltest.evolution.substitutionmodel.NucleotideRevJumpSubstModel")) {
						GeneralSubstitutionModel sm = (GeneralSubstitutionModel)o;
						if (sm.ratesInput.get() == state) {
							if (!taboo.contains(state)) taboo.add(state);
							Log.warning("NucleotideRevJumpSubstModel" + sm.getID() + " contains bModeltest rates " + state.getID() + " so it will not be operated on by " + this.getID());
							continue;
						}
					}
				}
			}
			
		}
		
		
		
		
		
		for (int i = 0; i < dists.size(); i ++) {
			List<ParallelMCMCTreeOperatorTreeDistribution> dist1 = dists.get(i);
		
			
			// Get list of state nodes in dist1
			Set<StateNode> dist1StateNodes = new HashSet<>();
			for (ParallelMCMCTreeOperatorTreeDistribution d : dist1) {
				dist1StateNodes.add(d.tree);
				ParallelMCMCRealParameterOperator.getRealParameterStateNodes(d.treelikelihood, otherState.stateNodeInput.get(), dist1StateNodes);
			}
			
			
			// Check if they are part of the state
			for (StateNode state : dist1StateNodes) {
				if (!stateNodeIds.contains(state.getID())){
					if (!taboo.contains(state)) {
						taboo.add(state);
					}
				}
				
			}
			
			
			
			
			for (int j = i+1; j < dists.size(); j ++) {
				List<ParallelMCMCTreeOperatorTreeDistribution> dist2 = dists.get(j);
				
				
				// Get list of state nodes in dist2
				for (ParallelMCMCTreeOperatorTreeDistribution d : dist2) {
					
					// Same tree appears in multiple threads? Do not operate on it
					if (dist1StateNodes.contains(d.tree)) {
						if (!taboo.contains(d.tree)) taboo.add(d.tree);
					}
					
					Set<StateNode> dist2StateNodes = new HashSet<>();
					ParallelMCMCRealParameterOperator.getRealParameterStateNodes(d.treelikelihood, otherState.stateNodeInput.get(), dist2StateNodes);
					for (StateNode s : dist2StateNodes) {
					
						// Same parameter appears in multiple threads? Do not operate on it
						if (dist1StateNodes.contains(s)) {
							if (!taboo.contains(s)) taboo.add(s);
						}
						
					}
				
					
				}
				
				
			}
			
			
		}
		
		
		
		
		return taboo;
		
	}
	
	

	@Override
	public double proposal() {
		
		
		/*
		// Stop threading the tree likelihood
		//boolean wasThreading = this.likelihoodInput.get().useThreads();
		if ((this.appliedRegression || !this.doRegression) && unthreadInput.get()) {
			//this.likelihoodInput.get().useThreads(false);
			for (ParallelMCMCTreeOperatorTreeDistribution distr : this.distributions) {
				//Log.warning("Threading off");
				distr.stopThreading();
				
			}
		}
		*/
		double logHR = super.proposal();
		
		
		/*
		// Start threading the tree likelihood again
		if ((this.appliedRegression || !this.doRegression) && unthreadInput.get()) {
			//this.likelihoodInput.get().useThreads(wasThreading);
			for (ParallelMCMCTreeOperatorTreeDistribution distr : this.distributions) {
				//Log.warning("Threading on");
				distr.startThreading();
			}
		}
		*/
		return logHR;
		
	}


	/**
	 * Check if we are sampling from the prior
	 * @param otherState
	 * @return
	 */
	private boolean sampleFromPriorEnabled() {
		

		
		for (BEASTInterface obj : this.getOutputs()) {
			Log.warning(obj.getID() + "," + obj.getClass().getCanonicalName());
            if (obj instanceof MCMC) {
            	MCMC mcmc = (MCMC)obj;
            	if (mcmc.sampleFromPriorInput.get()) {
            		Log.warning(this.getID() + ": sampling from prior is enabled");
            		return true;
            	}
            }
            
		}
		
		Log.warning(this.getID() + ": sampling from prior is disabled");
		return false;
	}

    
    

}
