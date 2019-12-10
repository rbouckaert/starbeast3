package starbeast3.operators;



import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;

import beast.app.BeastMCMC;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.ParallelMCMC;
import beast.core.Param;
import beast.core.State;
import beast.core.StateNode;
import beast.core.util.CompoundDistribution;
import beast.core.util.Log;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.operators.*;
import beast.evolution.tree.Tree;
import starbeast3.GeneTreeForSpeciesTreeDistribution;
import beast.util.Transform;

@Description("Run MCMC on different gene tree parts of the model in parallel before combining them in a single Gibbs move")
public class ParallelMCMCTreeOperator extends Operator implements MultiStepOperator {
	
	@Description("Distribution n a tree conditinionally independent from all other distributions given the state of the rest of parameter space")
	public class TreeDistribution {
		Tree tree;
		GenericTreeLikelihood treelikelihood;
		GeneTreeForSpeciesTreeDistribution geneprior;
		
		public Tree getTree() {return tree;}
		public void setTree(Tree tree) {this.tree = tree;}
		public GenericTreeLikelihood getTreelikelihood() {return treelikelihood;}
		public void setTreelikelihood(GenericTreeLikelihood treelikelihood) {this.treelikelihood = treelikelihood;}
		public GeneTreeForSpeciesTreeDistribution getGeneprior() {return geneprior;}
		public void setGeneprior(GeneTreeForSpeciesTreeDistribution geneprior) {this.geneprior = geneprior;}

		public TreeDistribution(@Param(name="tree", description="tree for which") Tree tree,
				@Param(name="treelikelihood", description="treelikelihood part of the distribution") GenericTreeLikelihood treelikelihood,
				@Param(name="geneprior", description="prior on the gene tree") GeneTreeForSpeciesTreeDistribution geneprior) {
			this.tree = tree;
			this.treelikelihood = treelikelihood;
			this.geneprior = geneprior;
		}
	}
	
    final public Input<Long> chainLengthInput =
            new Input<>("chainLength", "Length of the MCMC chain: each individual ParallelMCMC performs chainLength/nrOfThreads samples",
                    Input.Validate.REQUIRED);

	final public Input<List<TreeDistribution>> distributionInput = new Input<>("distribution", 
			"Distribution on a tree conditinionally independent from all other distributions given the state of the rest"
			+ "of parameter space. ",
			new ArrayList<>());
	
    final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads","maximum number of threads to use, if "
    		+ "less than 1 the number of threads in BeastMCMC is used (default -1)", -1);

    final public Input<State> otherStateInput = new Input<>("otherState", "main state containing all statenodes for this analysis");
    final public Input<Boolean> includeRealParametersInput = new Input<>("includeRealParameters", "flag to include Real Parameters for each of the partitions in the analysis", true);

    private ExecutorService exec;
    private CountDownLatch countDown;
    private List<ParallelMCMC> mcmcs;
    private State otherState;
    
	@Override
	public void initAndValidate() {
		List<TreeDistribution> distributions = distributionInput.get();
	    otherState = otherStateInput.get();
		 
		int nrOfThreads = maxNrOfThreadsInput.get() > 0 ? 
				Math.min(BeastMCMC.m_nThreads, maxNrOfThreadsInput.get()) : 
				BeastMCMC.m_nThreads;
	    exec = Executors.newFixedThreadPool(nrOfThreads);
	    mcmcs = new ArrayList<>();
	    
	    int start = 0;
	    for (int i = 0; i < nrOfThreads; i++) {
	    	int end = (i + 1) * distributions.size() / nrOfThreads;
	    	mcmcs.add(createParallelMCMC(distributions.subList(start, end), chainLengthInput.get()/nrOfThreads));
	    	start = end;
	    }
	    
	    for (ParallelMCMC pMCMC : mcmcs) {
	    	pMCMC.setOtherState(otherState);
	    }
	}

	@Override
	public int stepCount() {
		return mcmcs.size();
	}

	private ParallelMCMC createParallelMCMC(List<TreeDistribution> distributions, long chainLength) {
		List<Distribution> distrs = new ArrayList<>();
		List<StateNode> stateNodes = new ArrayList<>();
		List<Operator> operators = new ArrayList<>();
		for (TreeDistribution d : distributions) {
			distrs.add(d.geneprior);
			distrs.add(d.treelikelihood);
			
			stateNodes.add(d.tree);
			
//		    <operator id="treeScaler.t:$(n)" spec="ScaleOperator" scaleFactor="0.5" tree="@Tree.t:$(n)" weight="3.0"/>
			ScaleOperator treeScaler = new ScaleOperator();
			treeScaler.initByName("scaleFactor", 0.5, "tree", d.tree, "weight", 3.0);
			operators.add(treeScaler);
//	    	<operator id="treeRootScaler.t:$(n)" spec="ScaleOperator" rootOnly="true" scaleFactor="0.5" tree="@Tree.t:$(n)" weight="3.0"/>
			ScaleOperator treeRootScaler = new ScaleOperator();
			treeRootScaler.initByName("scaleFactor", 0.5, "tree", d.tree, "weight", 3.0, "rootOnly", true);
			operators.add(treeRootScaler);
//		    <operator id="UniformOperator.t:$(n)" spec="Uniform" tree="@Tree.t:$(n)" weight="30.0"/>
			Uniform UniformOperator = new Uniform();
			UniformOperator.initByName("tree", d.tree, "weight", 30.0);
			operators.add(UniformOperator);
//	    	<operator id="SubtreeSlide.t:$(n)" spec="SubtreeSlide" tree="@Tree.t:$(n)" weight="15.0"/>
			SubtreeSlide SubtreeSlide = new SubtreeSlide();
			SubtreeSlide.initByName("tree", d.tree, "weight", 15.0);
			operators.add(SubtreeSlide);
//		    <operator id="narrow.t:$(n)" spec="Exchange" tree="@Tree.t:$(n)" weight="15.0"/>
			Exchange narrow = new Exchange();
			narrow.initByName("tree", d.tree, "weight", 15.0);
			operators.add(narrow);
//	    	<operator id="wide.t:$(n)" spec="Exchange" isNarrow="false" tree="@Tree.t:$(n)" weight="3.0"/>
			Exchange wide = new Exchange();
			wide.initByName("tree", d.tree, "isNarrow", false, "weight", 3.0);
			operators.add(wide);
//		    <operator id="WilsonBalding.t:$(n)" spec="WilsonBalding" tree="@Tree.t:$(n)" weight="3.0"/>
			WilsonBalding WilsonBalding = new WilsonBalding();
			WilsonBalding.initByName("tree", d.tree, "weight", 3.0);
			operators.add(WilsonBalding);
			
			if (includeRealParametersInput.get()) {
				List<StateNode> stateNodeList = new ArrayList<>();
				ParallelMCMCRealParameterOperator.getRealParameterStateNodes(d.treelikelihood, otherState.stateNodeInput.get(), stateNodeList);
				stateNodes.addAll(stateNodeList);
				
				int dim = 0;
				for (StateNode s : stateNodeList) {
					dim += s.getDimension();
				}
				
				List<Transform> transformations = new ArrayList<>();
				for (StateNode s : stateNodeList) {
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
				AVMNOperator.initByName("weight", 1.0, "coefficient", 1.0, "scaleFactor", 1.0, "beta", 0.05, "every", 1,
						"initial", 200 * dim, "burnin", 100 * dim, "transformations", transformations);
				operators.add(AVMNOperator);

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

		ParallelMCMC mcmc = new ParallelMCMC();
		mcmc.initByName("state", state, "operator", operators, "distribution", sampleDistr, "chainLength", chainLength);
		return mcmc;
	}

	@Override
	public double proposal() {
		proposeUsingThreads();
		otherState.setEverythingDirty(true);
		return Double.POSITIVE_INFINITY;
	}

    class CoreRunnable implements Runnable {
        MCMC mcmc;

        CoreRunnable(MCMC core) {
        	mcmc = core;
        }

        @Override
		public void run() {
            try {
            	mcmc.run();
            } catch (Exception e) {
                Log.err.println("Something went wrong in a calculation of " + mcmc.getID());
                e.printStackTrace();
                System.exit(1);
            }
            countDown.countDown();
        }

    } // CoreRunnable

    private void proposeUsingThreads() {
        try {

            countDown = new CountDownLatch(mcmcs.size());
            // kick off the threads
            for (MCMC mcmc : mcmcs) {
                CoreRunnable coreRunnable = new CoreRunnable(mcmc);
                exec.execute(coreRunnable);
            }
            countDown.await();
        } catch (RejectedExecutionException | InterruptedException e) {
            Log.err.println("Stop using threads: " + e.getMessage());
        }
    }
    
    @Override
    public List<StateNode> listStateNodes() {
    	List<StateNode> stateNodes = new ArrayList<>();
    	for (ParallelMCMC mcmc : mcmcs) {
    		stateNodes.addAll(mcmc.startStateInput.get().stateNodeInput.get());
    	}
    	return stateNodes;
    }
}
