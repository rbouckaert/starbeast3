package starbeast3.operators;



import java.util.ArrayList;
import java.util.HashSet;
import java.util.List;
import java.util.Set;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;

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

@Description("Run MCMC on different treelikelihood parts of the model in parallel before combining them in a single Gibbs move")
public class ParallelMCMCRealParameterOperator extends Operator implements MultiStepOperator {
	
    final public Input<Long> chainLengthInput =
            new Input<>("chainLength", "Length of the MCMC chain: each individual ParallelMCMC performs chainLength/nrOfThreads samples",
                    Input.Validate.REQUIRED);

	final public Input<CompoundDistribution> distributionInput = new Input<>("distribution", 
			"compound distribution of all likelihoods",
			Validate.REQUIRED);
	
    final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads","maximum number of threads to use, if "
    		+ "less than 1 the number of threads in BeastMCMC is used (default -1)", -1);

    final public Input<State> otherStateInput = new Input<>("otherState", "");

    private ExecutorService exec;
    private CountDownLatch countDown;
    private List<ParallelMCMC> mcmcs;
    private State otherState;
    
	@Override
	public void initAndValidate() {
	    otherState = otherStateInput.get();
		List<Distribution> distributions = distributionInput.get().pDistributions.get();
		 
		int nrOfThreads = maxNrOfThreadsInput.get() > 0 ? 
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
	    
	    // 2. determine overlaps in potential state nodes per thread	    
	    Set<StateNode> tabu = new HashSet<>();
	    for (int i = 0; i < nrOfThreads; i++) {
	    	for (int j = i + 1; j < nrOfThreads; j++) {
	    		for (StateNode s : stateNodes[i]) {
	    			if (stateNodes[j].contains(s)) {
	    				tabu.add(s);
	    			}
	    		}
	    	}
	    }

	    
	    
	    start = 0;
	    for (int i = 0; i < nrOfThreads; i++) {
	    	int end = (i + 1) * distributions.size() / nrOfThreads;
	    	mcmcs.add(createParallelMCMC(distributions.subList(start, end), chainLengthInput.get()/nrOfThreads, tabu));
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
		
		ParallelMCMC mcmc = new ParallelMCMC();
		mcmc.initByName("state", state, "operator", operators, "distribution", sampleDistr, "chainLength", chainLength);
		return mcmc;
	}

	private void getRealParameterPriors(Set<StateNode> stateNodeList, List<Distribution> priorsList) {
		for (StateNode sn : stateNodeList) {
			for (BEASTInterface o : sn.getOutputs()) {
				if (o instanceof Distribution) {
					for (BEASTInterface o2 : o.getOutputs()) {
						if (o2.getID().equals("prior")) {
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
						(o2 instanceof BranchRateModel && o2.getInput("clock.rate").get() == o)) {
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
		proposeUsingThreads();
        for (int i = 0; i < otherState.stateNodeInput.get().size(); i++) {
        	otherState.getStateNode(i).index = i;
        }
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
