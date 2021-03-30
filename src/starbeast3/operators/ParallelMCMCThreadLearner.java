package starbeast3.operators;

import beast.core.Description;
import beast.core.util.Log;
import beast.util.Randomizer;


@Description("Learns the number of threads to alot to a parallel operator by maximising samples/hr.")
public class ParallelMCMCThreadLearner {
	
	long numCalls = 0;
	long burnin;
	boolean burnedin;
	
	MultiStepOperator operator;
	
	long numStatesWhenThreading;
	double meanRuntimeWhenThreading;
	long nSamplesThreading;
	
	
	long numStatesWhenNotThreading;
	double meanRuntimeWhenNotThreading;
	long nSamplesNotThreading;
	
	
	long startTime;
	boolean sampledThreads;
	
	
	
	/**
	 * Initialises the learning for this operator
	 * @param operator
	 * @param maxNrThreads
	 * @param chainLength
	 */
	public ParallelMCMCThreadLearner(MultiStepOperator operator, long chainLength, int burnin) {
		this.operator = operator;
		this.numCalls = 0;
		this.burnin = burnin;
		this.burnedin = false;
		this.numStatesWhenThreading = chainLength;
		this.numStatesWhenNotThreading = 1;
		this.meanRuntimeWhenThreading = 0;
		this.meanRuntimeWhenNotThreading = 0;
		this.nSamplesThreading = 0;
		this.nSamplesNotThreading = 0;
		this.sampledThreads = false;
	}
	
	
	
	/**
	 * Begin the timer and sample a thread setting
	 */
	public void start() {
		
		if (this.burnedin) {
			return;
		}
		
		this.startTime = System.currentTimeMillis();
		this.sampledThreads = Randomizer.nextBoolean();
		this.operator.useMCMC(this.sampledThreads);
	}
	
	
	/**
	 * Stop the timer and record results
	 */
	public void stop() {
		
		
		if (this.burnedin) {
			return;
		}
		
		
		long endTime = System.currentTimeMillis();
		long runtime = endTime - this.startTime;
		this.numCalls ++;
		
		if (this.sampledThreads) {
			this.meanRuntimeWhenThreading = (this.meanRuntimeWhenThreading*this.nSamplesThreading + runtime) / (this.nSamplesThreading + 1);
			this.nSamplesThreading ++;
		}else {
			this.meanRuntimeWhenNotThreading = (this.meanRuntimeWhenNotThreading*this.nSamplesNotThreading + runtime) / (this.nSamplesNotThreading + 1);
			this.nSamplesNotThreading ++;
		}
		
		
		// Burnin reached?
		if (this.nSamplesThreading >= this.burnin/2.0 && this.nSamplesNotThreading >= this.burnin/2.0) {
			
			// Make a decision on what to use now
			double threading = this.numStatesWhenThreading / this.meanRuntimeWhenThreading * 1000 * 3600;
			double notThreading = this.numStatesWhenNotThreading / this.meanRuntimeWhenNotThreading * 1000 * 3600;
			if (threading > notThreading) {
				Log.warning(operator.getID() + " will use threading (" + threading + "states/hr  >  " + notThreading + "states/hr.");
				this.sampledThreads = true;
			}else {
				Log.warning(operator.getID() + " will NOT use threading (" + notThreading + " states/hr  >  " + threading + " states/hr.");
				this.sampledThreads = false;
			}
			
			this.operator.useMCMC(this.sampledThreads);
			this.burnedin = true;
		}
	}
	
	

}
