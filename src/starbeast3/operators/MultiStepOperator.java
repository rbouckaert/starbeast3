package starbeast3.operators;

import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.RejectedExecutionException;

import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.OperatorScheduleRecalculator;
import beast.core.ParallelMCMC;
import beast.core.State;
import beast.core.StateNode;
import beast.core.util.Log;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

@Description("Operator that does proposals that count for one step or more steps in the MCMC")
public abstract class MultiStepOperator extends Operator {
	

	
	
	
    final public Input<Long> chainLengthInput =
            new Input<>("chainLength", "Length of the MCMC chain: each individual ParallelMCMC performs chainLength/nrOfThreads samples");
    
    final public Input<Double> coverageInput =
            new Input<>("chainCoverage", "The MCMC chain length is the coverage times the number of parameters",
                    Input.Validate.XOR, chainLengthInput);


	
    final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads", "maximum number of threads to use, if "
    		+ "less than 1 the number of threads in BeastMCMC is used (default -1)", -1);

    final public Input<State> otherStateInput = new Input<>("otherState", "main state containing all statenodes for this analysis");

    
    final public Input<Boolean> learningInput =  new Input<>("learning", "Learn whether to parallelise (n threads) or not (1 thread 1 operator)", true);
    final public Input<Integer> burninInput =  new Input<>("burnin", "How many operator calls before thread learning kicks in. Learning will begin after chainLength regression.", 10000);
    final public Input<OperatorScheduleRecalculator> scheduleInput =  new Input<>("schedule", "Operator schedule (if learning is applied)");
    
    final public Input<Integer> nregressionInput =  new Input<>("nregression", "Number of MCMC chainLengths vs runtimes to sample in order to learn chainLengths, for load"
    		+ " balancing. Set to <5 to skip the training.", 200);
    final public Input<Double> targetCPUInput =  new Input<>("targetCPU", "Proportion of threads allocated (if > 1) that should be spent on MCMC, as opposed to overhead."
    		+ "Larger targetCPU will mean longer chains and lower operator weight. If targetCPU=0, then the load balancing will match the slowest thread. Set nregression=0 to omit this step.", 0.8);
    final public Input<Double> targetWeightInput =  new Input<>("targetWeight", "Target effective weight of this operator, to be learned if regression is applied."
    		+ " The effective weight is the operator weight * chainLength sum. Set this to 0 (or nregression=0) to omit this step.", 0.0);
    final public Input<Double> runtimeInput =  new Input<>("runtime", "Max runtime of MCMC chains during training (only applicable if load balancing is being trained). If this is"
    		+ "set to -1, then chain lengths are sampled instead of runtimes.", -1.0);
    
    
    final public Input<Tree> dummyInput = new Input<>("speciesTree", "an optional dummy input so that beauti can load the template (hack)", Input.Validate.OPTIONAL);
    
      //final public Input<CompoundDistribution> likelihoodInput = new Input<>("likelihood", "the likelihood", Input.Validate.REQUIRED);

	
    protected long chainLength;
    protected int nrOfThreads;
    
    protected State otherState;
    protected boolean useMCMC;
	protected ExecutorService exec;
    protected CountDownLatch countDown;
    protected List<ParallelMCMC> mcmcs;
    protected List<CoreRunnable> runnables;
    
    
    
    // For when the chain length is 1 and the number of threads is 1
    protected ParallelMCMCThreadLearner learner;
    protected List<Operator> singleStepOperators;
    protected double[] operatorProbs;
    
    
    // Thread learning
    boolean learnThreads;
    
    // Regression
    boolean doRegression;
    boolean appliedRegression;
    
    long nproposals;
	
	/** number of steps to be performed by operator **/
	public int stepCount() {
		return 1;
	}
	
	
	public void useMCMC(boolean useMCMC) {
		this.useMCMC = useMCMC;
	}
	
	@Override
	public void initAndValidate() {
		
		
		this.nproposals = 0;
		
		// Doing regression on chainlengths?
		if (nregressionInput.get() >= 5 && this.mcmcs.size() > 1) {
			this.doRegression = true;
			if (targetCPUInput.get() >= 1) {
				throw new IllegalArgumentException("targetCPU must be less than 1");
			}
		}else {
			this.doRegression = false;
		}
		this.appliedRegression = false;
		
		
		this.useMCMC = true;
		
		// 1 thread and 1 chain mode
	    if (learningInput.get()) {
	    	
	    	// Learn after regression
	    	this.learnThreads = !this.doRegression;
	    	
	    	
	    	if (scheduleInput.get() == null) {
	    		throw new IllegalArgumentException("Please provide an operator schedule (or set learning=false)");
	    	}
	    	
	    	this.learner = new ParallelMCMCThreadLearner(this, chainLength, burninInput.get(), this.scheduleInput.get(), this.targetWeightInput.get());
	    	
			
	    	this.learnThreads = false;
	    	this.learner = null;
	    }
	    
	    
	    // Single operator sampling
	    if (learningInput.get() || !this.useMCMC || (this.mcmcs.size() == 1 && chainLength == 1)) {
	    	this.operatorProbs = new double[this.singleStepOperators.size()];
	    	double weightSum = 0;
	    	double cumSum = 0;
	    	for (Operator op : this.singleStepOperators) weightSum += op.getWeight();
	    	for (int i = 0; i < this.singleStepOperators.size(); i ++) {
	    		Operator op = this.singleStepOperators.get(i);
	    		cumSum += op.getWeight() / weightSum;
	    		this.operatorProbs[i] = cumSum;
	    	}
	    }else {
	    }
	    
	    
	    this.runnables = new ArrayList<>();
	    for (MCMC mcmc : this.mcmcs) {
	    	CoreRunnable runnable = new CoreRunnable(mcmc);
	    	this.runnables.add(runnable);
	    }
		  
	}
	
	
	@Override 
	public void accept() {
		train();
	}
	
	@Override 
	public void reject() {
		train();
	}
	
	@Override 
	public void reject(int reasonNr) {
		train();
	}
	
	
	/**
	 * Update thread learning and chainLength learning
	 */
	private void train() {
		
		// The thread learner will cache what it has learned
		if (this.learnThreads) this.learner.stop();
		
		// Regression model
		if (!this.appliedRegression && this.doRegression && this.mcmcs.get(0).finishedRegression()) {
			
			// Find the slowest thread, by runtime vs chainlength slope
			int slowestThread = 0;
			double slowestSlope = this.mcmcs.get(0).getRuntimeSlope();
			for (int i = 1; i < this.mcmcs.size(); i ++) {
				double slope = this.mcmcs.get(i).getRuntimeSlope();
				if (slope > slowestSlope) {
					slowestThread = i;
					slowestSlope = slope;
				}
			}
			
			
			double targetRuntime;
			
			double meanIntercept = 0;
			double meanSlope = 0;
			for (ParallelMCMC mcmc : this.mcmcs) {
				meanIntercept += mcmc.getRuntimeIntercept() / this.mcmcs.size();
				meanSlope += mcmc.getRuntimeSlope() / this.mcmcs.size();
			}
			
			
			// Set the runtime of the slowest thread to match the user-specified CPU use
			if (targetCPUInput.get() > 0) {
				
				// The overhead proportion must be less than this value to justify using this many threads
				// eg. if there are 2 threads, then we must make sure that less than 50% of the thread time goes to overhead
				double maximumOverhead = (this.nrOfThreads-1.0) / this.nrOfThreads;
				
				// If we are to use targetCPU of the additional threads (beyond the first one), then this should be the overhead proportion
				// eg. if there are 3 threads, then targetCPU=0.5 would mean that overhead should be 1/3
				double targetOverhead = maximumOverhead * (1 - targetCPUInput.get());
				
				
				
				// Set the target runtime such that the target overhead is attained for this thread (the slowest thread)
				double slope = meanSlope;
				double intercept = meanIntercept;
				int targetChainLength = (int)((intercept/targetOverhead - intercept) / slope);
				targetChainLength = Math.max(targetChainLength, 1);
				targetRuntime = this.mcmcs.get(slowestThread).predict(targetChainLength);
				Log.warning(this.getID() + ": want all mcmcs to have a runtime of " + targetRuntime + "ms. This will give the average thread "
						+ "an overhead of " + (targetOverhead*100) + "%. On average, runtime(ms) = " + meanIntercept + " + " + meanSlope + "*chainLength.");
				

				
			}
			
			// Set the runtime of the slowest thread to match the user-specified chain length
			else {
				targetRuntime = this.mcmcs.get(slowestThread).predict(this.chainLength / this.nrOfThreads);
				
				
				
				Log.warning(this.getID() + ": want all mcmcs to have a runtime of " + targetRuntime + "ms. This is the average runtime of thread " + (1+slowestThread) + ". "
						+ "On average, runtime(ms) = " + meanIntercept + " + " + meanSlope + "*chainLength.");
			
			
			}
			
			
			// Set the runtime of all chains to that of the slowest
			for (ParallelMCMC mcmc : this.mcmcs) {
				if (targetRuntime > 0) mcmc.setChainlengthToTargetRuntime(targetRuntime);
				else mcmc.setChainlengthToTargetRuntime(chainLength*1.0 / this.nrOfThreads);
			}
			
			
			
			
			
			
			// Adjust the weight of this operator, such that the new weight is equal to the effective weight divided by the chainLength sum
			if (targetWeightInput.get() > 0) {
				
				long chainLengthSum = 0;
				for (ParallelMCMC mcmc : this.mcmcs) {
					chainLengthSum += mcmc.getChainLength();
				}
				if (this.learnThreads) this.learner.setChainLength(chainLengthSum);
				double newWeight = targetWeightInput.get() / chainLengthSum;
				Log.warning(this.getID() + ": setting operator weight to " + newWeight + " to attain an effective weight of " + targetWeightInput.get());
				this.m_pWeight.set(newWeight);
				this.scheduleInput.get().reweight();
				
				
				
			}

		
			this.appliedRegression = true;
			
			
			// Can start learning threads now that the chain lengths have optimised
			if (this.learner != null) this.learnThreads = true;
			
		}
		
	}
	
	
	
	@Override
	public double proposal() {
		
		double logHR;
	
		// The learner will sample whether MCMC is run or not
		if (this.learnThreads) this.learner.start();
		
		// Do not waste time creating a whole mcmc chain if it has only 1 step
		if (!this.useMCMC || (this.mcmcs.size() == 1 && chainLength == 1)) {
			
			// Sample operator
			int opIndex = Randomizer.randomChoice(this.operatorProbs);
			//Log.warning("sampling operator at index " + opIndex + " out of " + this.singleStepOperators.size());
			logHR = this.singleStepOperators.get(opIndex).proposal();
			
			
			nproposals++;
			
			
		// Run 1 or more MCMCs
		} else {
		
			
			// Do not waste time creating threads if there is only 1 thread
			if (this.mcmcs.size() == 1) {
				try {
					ParallelMCMC mcmc = this.mcmcs.get(0);
					mcmc.run();
					nproposals += mcmc.getChainLength();
				} catch (Exception e) {
					e.printStackTrace();
				}
			}else {
				proposeUsingThreads();
			}
			otherState.setEverythingDirty(true);
			
			logHR = Double.POSITIVE_INFINITY;
		
		
		}
		
		
		

		return logHR;
		
	}
	
    private void proposeUsingThreads() {
        try {
        	
            countDown = new CountDownLatch(mcmcs.size());
            
            // Sample runtime
            if (!this.appliedRegression && this.doRegression && runtimeInput.get() > 0) {
            	double runtime = Randomizer.nextDouble() * runtimeInput.get();
            	//Log.warning("Runtime " + runtime);
	            for (ParallelMCMC mcmc : mcmcs) {
	            	if (runtimeInput.get() > 0) {
	            		mcmc.setRuntime((long)runtime);
	            	}
	            }
            }
            
            // Kick off the threads
            for (CoreRunnable runnable : this.runnables) {
                exec.execute(runnable);
            }
            
            
            countDown.await();
            
            
            for (ParallelMCMC mcmc : mcmcs) {
            	//Log.warning("there were " + mcmc.getChainLength() + " proposals");
            	nproposals += mcmc.getChainLength();
            }
            
        } catch (RejectedExecutionException | InterruptedException e) {
            Log.err.println("Stop using threads: " + e.getMessage());
        }
    }
    
    
    public class CoreRunnable implements Runnable {
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
    
    
    
    
    public long getNrProposals() {
    	return nproposals;
    }


    
    @Override
    public List<StateNode> listStateNodes() {
    	List<StateNode> stateNodes = new ArrayList<>();
    	for (ParallelMCMC mcmc : mcmcs) {
    		stateNodes.addAll(mcmc.startStateInput.get().stateNodeInput.get());
    	}
    	if (stateNodes.isEmpty()) stateNodes.add(dummyInput.get());
    	return stateNodes;
    }
    
    
	protected static void getRealParameterPriors(Set<StateNode> stateNodeList, Set<Distribution> priorsList) {
		for (StateNode sn : stateNodeList) {
			for (BEASTInterface o : sn.getOutputs()) {
				if (o instanceof Distribution) {
					getRealParameterPriors(o, (Distribution) o, priorsList);
				}
			}
		}		
	}
	
	
	

	protected static void getRealParameterPriors(BEASTInterface o, Distribution distr, Set<Distribution> priorsList) {
		for (BEASTInterface o2 : o.getOutputs()) {
			if (o2.getID() != null && o2.getID().equals("prior")) {
				priorsList.add(distr);
				break;
			}
			getRealParameterPriors(o2, distr, priorsList);
		}
	}
    
	
	
}
