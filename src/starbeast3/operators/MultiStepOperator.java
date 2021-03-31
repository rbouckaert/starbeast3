package starbeast3.operators;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.RejectedExecutionException;

import beast.core.Description;
import beast.core.Input;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.OperatorScheduleRecalculator;
import beast.core.ParallelMCMC;
import beast.core.State;
import beast.core.StateNode;
import beast.core.util.Log;
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
    final public Input<Boolean> includeRealParametersInput = new Input<>("includeRealParameters", "flag to include Real Parameters for each of the partitions in the analysis", true);

    
    final public Input<Boolean> learningInput =  new Input<>("learning", "Learn whether to parallelise (n threads) or not (1 thread 1 operator)", true);
    final public Input<Integer> burninInput =  new Input<>("burnin", "How many operator calls before thread learning kicks in", 10000);
    final public Input<OperatorScheduleRecalculator> scheduleInput =  new Input<>("schedule", "Operator schedule (if learning is applied)");
    
    final public Input<Integer> nregressionInput =  new Input<>("nregression", "Number of MCMC chainLengths vs runtimes to sample in order to learn chainLengths, for load"
    		+ " balancing. Set to <5 to skip the training.", 2000);
    
      //final public Input<CompoundDistribution> likelihoodInput = new Input<>("likelihood", "the likelihood", Input.Validate.REQUIRED);

	
    protected long chainLength;
    protected int nrOfThreads;
    
    protected State otherState;
    protected boolean useMCMC;
	protected ExecutorService exec;
    protected CountDownLatch countDown;
    protected List<ParallelMCMC> mcmcs;
    
    
    
    // For when the chain length is 1 and the number of threads is 1
    protected ParallelMCMCThreadLearner learner;
    protected List<Operator> singleStepOperators;
    protected double[] operatorProbs;
    
    
    // Regression
    boolean doRegression;
    boolean appliedRegression;
	
	/** number of steps to be performed by operator **/
	public int stepCount() {
		return 1;
	}
	
	
	public void useMCMC(boolean useMCMC) {
		this.useMCMC = useMCMC;
	}
	
	@Override
	public void initAndValidate() {
		
		// Doing regression on chainlengths?
		if (nregressionInput.get() >= 5 && this.mcmcs.size() > 1) {
			this.doRegression = true;
		}else {
			this.doRegression = false;
		}
		this.appliedRegression = false;
		
		
		this.useMCMC = true;
		
		// 1 thread and 1 chain mode
	    if (learningInput.get()) {
	    	
	    	
	    	if (scheduleInput.get() == null) {
	    		throw new IllegalArgumentException("Please provide an operator schedule (or set learning=false)");
	    	}
	    	
	    	this.learner = new ParallelMCMCThreadLearner(this, chainLength, burninInput.get(), this.scheduleInput.get());
	    	
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
	    	this.learner = null;
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
		if (this.learner != null) this.learner.stop();
		
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
			double targetRuntime = this.mcmcs.get(slowestThread).predict(this.chainLength / this.nrOfThreads);
			Log.warning("Want all mcmcs to have a runtime of " + targetRuntime + ". This is the average runtime of thread " + slowestThread);
			
			// Set the chainLength of all chains to match it
			for (ParallelMCMC mcmc : this.mcmcs) {
				mcmc.applyRegression(targetRuntime);
			}
			
		
			this.appliedRegression = true;
			
		}
		
	}
	
	
	
	@Override
	public double proposal() {
		
		double logHR;
	
		// The learner will sample whether MCMC is run or not
		if (this.learner != null) this.learner.start();
		
		// Do not waste time creating a whole mcmc chain if it has only 1 step
		if (!this.useMCMC) {
			
			// Sample operator
			int opIndex = Randomizer.randomChoice(this.operatorProbs);
			//Log.warning("sampling operator at index " + opIndex + " out of " + this.singleStepOperators.size());
			logHR = this.singleStepOperators.get(opIndex).proposal();
			
			
		// Run 1 or more MCMCs
		} else {
		
			
			// Do not waste time creating threads if there is only 1 thread
			if (this.mcmcs.size() == 1) {
				try {
					ParallelMCMC mcmc = this.mcmcs.get(0);
					mcmc.run();
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
    
    
    


    
    @Override
    public List<StateNode> listStateNodes() {
    	List<StateNode> stateNodes = new ArrayList<>();
    	for (ParallelMCMC mcmc : mcmcs) {
    		stateNodes.addAll(mcmc.startStateInput.get().stateNodeInput.get());
    	}
    	return stateNodes;
    }
    
	
	
}
