package starbeast3.operators;



import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Executors;

import beast.app.BeastMCMC;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Function;
import beast.core.Input;
import beast.core.Operator;
import beast.core.ParallelMCMC;
import beast.core.State;
import beast.core.StateNode;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;
import beast.core.util.Log;
import beast.evolution.operators.*;
import beast.util.Transform;
import starbeast3.GeneTreeForSpeciesTreeDistribution;

@Description("Run MCMC on different gene tree parts of the model in parallel before combining them in a single Gibbs move")
public class ParallelMCMCTreeOperator extends MultiStepOperator {
	
    final public Input<Boolean> useBactrianOperatorsInput = new Input<>("bactrian", "flag to indicate that bactrian operators should be used where possible", true);
    final public Input<Boolean> includeRealParametersInput = new Input<>("includeRealParameters", "flag to include Real Parameters for each of the partitions in the analysis", true);

	final public Input<List<ParallelMCMCTreeOperatorTreeDistribution>> distributionInput = new Input<>("distribution", 
			"Distribution on a tree conditinionally independent from all other distributions given the state of the rest"
			+ "of parameter space. ",
			new ArrayList<>());
    
	
	final public Input<Boolean> unthreadInput = new Input<>("unthread", "flag to convert ThreadedTreeLikelihood back into TreeLikelihood when this is called", false);
	  

    
    List<ParallelMCMCTreeOperatorTreeDistribution> distributions;
    

    

    
    
	@Override
	public void initAndValidate() {
		this.distributions = distributionInput.get();
		mcmcs = new ArrayList<>();
		
		if (distributions.isEmpty()) {
			Log.warning("ParallelMCMCTreeOperator: Please provide at least one 'distribution'");
			return;
			
		}
		
	    otherState = otherStateInput.get();
		 
		nrOfThreads = maxNrOfThreadsInput.get() > 0 ?
				Math.min(BeastMCMC.m_nThreads, maxNrOfThreadsInput.get()) : 
				BeastMCMC.m_nThreads;
		nrOfThreads = Math.min(nrOfThreads, distributions.size());
		Log.warning("Running " + this.getID() + " with " + this.nrOfThreads + " threads");
		//System.exit(1);
	    exec = Executors.newFixedThreadPool(nrOfThreads);
	    
	    
	    
	    
	    // Load balancing. Ensure a roughly equal distribution of site patterns across all threads
	    int totalDim = 0;
	    Collections.sort(distributions);
	    List<List<ParallelMCMCTreeOperatorTreeDistribution>> balancedDistributions = new ArrayList<>();
	    for (int i = 0; i < nrOfThreads; i++) balancedDistributions.add(new ArrayList<>());
	    int threadNum = 0;
	    for (int i = 0; i < distributions.size(); i ++) {
	    	ParallelMCMCTreeOperatorTreeDistribution d = distributions.get(i);
	    	totalDim += d.getTree().getTaxaNames().length;
	    	//System.out.println("patterns " + d.getNumberPatterns());
	    	balancedDistributions.get(threadNum).add(d);
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
	    
	    
	    // Create parallel MCMCs
	    for (int i = 0; i < nrOfThreads; i++) {
	    	mcmcs.add(createParallelMCMC(balancedDistributions.get(i), (int)(1.0*chainLength / this.nrOfThreads)));
	    }
	    Log.warning(this.getClass().getCanonicalName() + ": total chain length is " + chainLength);

	    
	    for (ParallelMCMC pMCMC : mcmcs) {
	    	pMCMC.setOtherState(otherState);
	    }
	    
	    
	    // 1 thread and 1 chain mode
	    if (learningInput.get() || (this.mcmcs.size() == 1 && chainLength == 1)) {
	    	ParallelMCMC mcmc;
	    	if (this.mcmcs.size() == 1 && chainLength == 1) mcmc = this.mcmcs.get(0);
	    	else mcmc = createParallelMCMC(distributions, chainLength);
	    	this.singleStepOperators = mcmc.operatorsInput.get();
	    }
	    
	  
	    super.initAndValidate();
	  
	    
	    
	}

	@Override
	public int stepCount() {
		return mcmcs.size();
	}

	private ParallelMCMC createParallelMCMC(List<ParallelMCMCTreeOperatorTreeDistribution> distributions, long chainLength) {
		List<Distribution> distrs = new ArrayList<>();
		List<StateNode> stateNodes = new ArrayList<>();
		List<Operator> operators = new ArrayList<>();
		for (ParallelMCMCTreeOperatorTreeDistribution d : distributions) {
			distrs.add(d.geneprior);
			distrs.add(d.treelikelihood); 
			
			
			System.out.println("Adding: " + d.geneprior.getID() + " / " + d.treelikelihood.getID());
			
			stateNodes.add(d.tree);
			
			// Uniform operator
			beast.evolution.operators.Uniform UniformOperator = new beast.evolution.operators.Uniform();
			UniformOperator.initByName("tree", d.tree, "weight", 30.0);
			operators.add(UniformOperator);
			
			if (useBactrianOperatorsInput.get()) {
				
				
				// Root scaler
				beast.evolution.operators.BactrianScaleOperator treeRootScaler = new beast.evolution.operators.BactrianScaleOperator();
				treeRootScaler.initByName("scaleFactor", 0.5, "tree", d.tree, "weight", 10.0, "rootOnly", true);
				operators.add(treeRootScaler);
				
				// Bactrian interval
				beast.evolution.operators.BactrianNodeOperator intervalOperator = new beast.evolution.operators.BactrianNodeOperator();
				intervalOperator.initByName("tree", d.tree, "weight", 10.0);
				operators.add(intervalOperator);
				
				// Subtree slide
				beast.evolution.operators.BactrianSubtreeSlide SubtreeSlide = new beast.evolution.operators.BactrianSubtreeSlide();
				SubtreeSlide.initByName("tree", d.tree, "weight", 10.0);
				operators.add(SubtreeSlide);
				
				
				// Adaptable operators for tree scaling
				orc.operators.AdaptableOperatorSampler adaptable = new orc.operators.AdaptableOperatorSampler();
				List<Operator> adaptOperators = new ArrayList<>();
				
				// Tree scaler
				beast.evolution.operators.BactrianUpDownOperator treeScaler = new beast.evolution.operators.BactrianUpDownOperator();
				treeScaler.initByName("scaleFactor", 0.5, "up", d.tree, "weight", 1.0);
				adaptOperators.add(treeScaler);
				
				// Epoch scaler
				EpochOperator epochOperator = new EpochOperator();
				List<GeneTreeForSpeciesTreeDistribution> trees = new ArrayList<>();
				trees.add(d.geneprior);
				epochOperator.initByName("scaleFactor", 0.5, "gene", trees, "weight", 1.0);
				adaptOperators.add(epochOperator);
				
				// Subtree slide
				adaptOperators.add(SubtreeSlide);
				
				// Uniform
				adaptOperators.add(UniformOperator);
				
				adaptable.initByName("tree", d.tree, "weight", 100.0, "operator", adaptOperators);
				operators.add(adaptable);
				
				
				
			} else {
//		        <operator id="treeScaler.t:$(n)" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:$(n)" weight="3.0"/>
				ScaleOperator treeScaler = new ScaleOperator();
				treeScaler.initByName("scaleFactor", 0.5, "tree", d.tree, "weight", 3.0);
				operators.add(treeScaler);
//		    	<operator id="treeRootScaler.t:$(n)" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:$(n)" weight="3.0"/>
				ScaleOperator treeRootScaler = new ScaleOperator();
				treeRootScaler.initByName("scaleFactor", 0.5, "tree", d.tree, "weight", 3.0, "rootOnly", true);
				operators.add(treeRootScaler);

//		    	<operator id="SubtreeSlide.t:$(n)" spec="SubtreeSlide" tree="@Tree.t:$(n)" weight="15.0"/>
				beast.evolution.operators.SubtreeSlide SubtreeSlide = new beast.evolution.operators.SubtreeSlide();
				SubtreeSlide.initByName("tree", d.tree, "weight", 15.0);
				operators.add(SubtreeSlide);
			}
			
			

//		    <operator id="narrow.t:$(n)" spec="Exchange" tree="@Tree.t:$(n)" weight="15.0"/>
			beast.evolution.operators.Exchange narrow = new beast.evolution.operators.Exchange();
			narrow.initByName("tree", d.tree, "weight", 15.0);
			operators.add(narrow);
//	    	<operator id="wide.t:$(n)" spec="Exchange" isNarrow="false" tree="@Tree.t:$(n)" weight="3.0"/>
			beast.evolution.operators.Exchange wide = new beast.evolution.operators.Exchange();
			wide.initByName("tree", d.tree, "isNarrow", false, "weight", 15.0);
			operators.add(wide);
//		    <operator id="WilsonBalding.t:$(n)" spec="WilsonBalding" tree="@Tree.t:$(n)" weight="3.0"/>
			beast.evolution.operators.WilsonBalding WilsonBalding = new beast.evolution.operators.WilsonBalding();
			WilsonBalding.initByName("tree", d.tree, "weight", 15.0);
			operators.add(WilsonBalding);
			if (includeRealParametersInput.get()) {
				Set<StateNode> stateNodeList = new HashSet<>();
				ParallelMCMCRealParameterOperator.getRealParameterStateNodes(d.treelikelihood, otherState.stateNodeInput.get(), stateNodeList);
				stateNodes.addAll(stateNodeList);
				
				
				/**
				 * Add priors to dist
				 */
				Set<Distribution> priorsSet = new HashSet<>();
				getRealParameterPriors(stateNodeList, priorsSet);
				distrs.addAll(priorsSet);	
				for (Distribution d2 : priorsSet) {
					Log.warning("Adding dist " + d2.getID());
				}
				
				
				int dim = 0;
				for (StateNode s : stateNodeList) {
					dim += s.getDimension();
				}
				
				List<Transform> transformations = new ArrayList<>();
				
				// Add the tree
				Transform f = new Transform.LogTransform(d.tree);
				transformations.add(f);
				Log.warning("Adding " + d.tree.getID());
				
				for (StateNode s : stateNodeList) {
					
					
					if (!(s instanceof RealParameter)) continue;
					RealParameter rp = (RealParameter)s;
					
					
					// TODO: check priors instead of ID to determine whether it is a
					// scale parameter
					// location parameter
					// simplex parameter
					if (s.getID().startsWith("freq")) {
						
						DeltaExchangeOperator op = new DeltaExchangeOperator();
						op.initByName("parameter", s, "delta", 0.2, "weight", 0.5);
						operators.add(op);
						f = new Transform.LogConstrainedSumTransform(s, 1.0);
						
						
						Log.warning("Adding logsum " + s.getID());
						
					} 
					
					// Logit
					else if (rp.getLower() != Double.NEGATIVE_INFINITY && rp.getUpper() != Double.POSITIVE_INFINITY) {
						
						BactrianIntervalOperator op = new BactrianIntervalOperator();
						op.initByName("parameter", s, "scaleFactor", 0.5, "weight", 0.5);
						operators.add(op);
						
						List<Function> l = new ArrayList<>();
						l.add(s);
						f = new Transform.LogitTransform(l);
						
						Log.warning("Adding logit " + s.getID());
						
						
						
						
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
					
					transformations.add(f);
					
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



    
    

}
