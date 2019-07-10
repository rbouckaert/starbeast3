package starbeast3.operators;

import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;

import beast.app.BeastMCMC;
import beast.core.Description;
import beast.core.Input;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.ParallelMCMC;
import beast.core.State;
import beast.core.StateNode;
import beast.core.util.Log;

@Description("Run MCMC on different parts of the model in parallel before combining them in a single Gibbs move")
public class ParallelMCMCOperator extends Operator {
	final public Input<List<ParallelMCMC>> mcmcInput = new Input<>("mcmc", "mcmc specification on a small subset of parameters. "
			+ "Each mcmc has to be on an independent part of parameter space, for instance, mutation rate moves "
			+ "for individual partitions.", new ArrayList<>());
    final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads","maximum number of threads to use, if "
    		+ "less than 1 the number of threads in BeastMCMC is used (default -1)", -1);

    final public Input<State> otherStateInput = new Input<>("otherState", "");

    private ExecutorService exec;
    private CountDownLatch countDown;
    private List<ParallelMCMC> mcmcs;
    private State otherState;
    
	@Override
	public void initAndValidate() {
		mcmcs = mcmcInput.get();
		int nrOfThreads = maxNrOfThreadsInput.get() > 0 ? 
				Math.min(BeastMCMC.m_nThreads, maxNrOfThreadsInput.get()) : 
				BeastMCMC.m_nThreads;
	    exec = Executors.newFixedThreadPool(nrOfThreads);
	    otherState = otherStateInput.get();
	    for (ParallelMCMC pMCMC : mcmcs) {
	    	pMCMC.setOtherState(otherState);
	    }
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
