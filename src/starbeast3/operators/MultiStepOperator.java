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
	
	/** number of steps to be performed by operator **/
	public int stepCount() {
		return 1;
	}
	
	
	public void useMCMC(boolean useMCMC) {
		this.useMCMC = useMCMC;
	}
	
	@Override
	public void initAndValidate() {
		
		this.useMCMC = true;
		
		// 1 thread and 1 chain mode
	    if (learningInput.get()) {
	    	
	    	this.learner = new ParallelMCMCThreadLearner(this, chainLength, burninInput.get());
	    	
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
		
		
		// The learner will cache what it has learned
		if (this.learner != null) this.learner.stop();

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
