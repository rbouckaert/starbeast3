package beast.core;

import java.io.IOException;
import java.math.BigDecimal;
import java.math.MathContext;
import java.util.HashSet;
import java.util.List;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.math3.stat.regression.SimpleRegression;
import org.xml.sax.SAXException;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input.Validate;
import beast.core.MCMC;
import beast.core.Operator;
import beast.core.State;
import beast.core.StateNode;
import beast.core.util.CompoundDistribution;
import beast.core.util.Evaluator;
import beast.core.util.Log;
import beast.util.Randomizer;

@Description("Implements MCMC without logging, or resume and suppresses all screen output. Used for ParallelMCMCOperator.")
public class ParallelMCMC extends MCMC {
	
	
	final public Input<Boolean> robustInput = new Input<>("robust", "whether to periodically robustly recalculate posterior", true);
	final public Input<Integer> nregressionInput = new Input<>("nregression", "number of chainLength/runtime samples to make, or 0 for no learning", 0);
	
	
	private final int MIN_CHAIN_LENGTH = 50;
    private long sampleCount = 0;
    private boolean robust;
    
    // Regression
    private SimpleRegression regression = null;
    private long initChainLength;
    private int regressionNr;
    private boolean doRegression;
    private boolean finishedRegression;
    private boolean appliedRegression;
    private long[] chainLengths;
    private double[] runtimes;
    
    private long targetRuntime = -1;
    private long numStepsInChain;
    private long startTime;

	public void setOtherState(State otherState) {
		this.otherState = otherState;
	}
	
	public ParallelMCMC() {
		loggersInput.setRule(Validate.OPTIONAL);
		startStateInput.setRule(Validate.REQUIRED);
	}
	
	public void setLogPosterior(double logP) {
		this.oldLogLikelihood = logP;
	}
	

    private static boolean printDebugInfo = false;

    private State otherState;
    int [] otherStateNr;
    
    @Override
    public void initAndValidate() {
    	
    	
    	// Regression for learning chain length
    	this.doRegression = nregressionInput.get() >= 5;
    	this.finishedRegression = false;
    	this.appliedRegression = false;
    	if (this.doRegression) {
    		this.initChainLength = Math.max(1, chainLengthInput.get());
    		this.regressionNr = 0;
    		this.chainLengths = new long[nregressionInput.get()];
    		this.runtimes = new double[nregressionInput.get()];
    	}
    	
    	this.robust = robustInput.get();
    	state = startStateInput.get();
    	otherStateNr = new int[state.stateNodeInput.get().size()];
    	
        operatorSchedule = operatorScheduleInput.get();
        for (final Operator op : operatorsInput.get()) {
            operatorSchedule.addOperator(op);
        }

        if (sampleFromPriorInput.get()) {
            // remove beastObject with id likelihood from posterior, if it is a CompoundDistribution
            if (posteriorInput.get() instanceof CompoundDistribution) {
                final CompoundDistribution posterior = (CompoundDistribution) posteriorInput.get();
                final List<Distribution> distrs = posterior.pDistributions.get();
                final int distrCount = distrs.size();
                for (int i = 0; i < distrCount; i++) {
                    final Distribution distr = distrs.get(i);
                    final String id = distr.getID();
                    if (id != null && id.equals("likelihood")) {
                        distrs.remove(distr);
                        break;
                    }
                }
                if (distrs.size() == distrCount) {
                    throw new RuntimeException("Sample from prior flag is set, but distribution with id 'likelihood' is " +
                            "not an input to posterior.");
                }
            } else {
                throw new RuntimeException("Don't know how to sample from prior since posterior is not a compound distribution. " +
                        "Suggestion: set sampleFromPrior flag to false.");
            }
        }

        // State initialisation
        final HashSet<StateNode> operatorStateNodes = new HashSet<>();
        for (final Operator op : operatorsInput.get()) {
            for (final StateNode stateNode : op.listStateNodes()) {
                operatorStateNodes.add(stateNode);
            }
        }
        if (startStateInput.get() != null) {
            this.state = startStateInput.get();
            if (storeEveryInput.get() > 0) {
                this.state.m_storeEvery.setValue(storeEveryInput.get(), this.state);
            }
        } else {
            // create state from scratch by collecting StateNode inputs from Operators
            this.state = new State();
            for (final StateNode stateNode : operatorStateNodes) {
                this.state.stateNodeInput.setValue(stateNode, this.state);
            }
            this.state.m_storeEvery.setValue(storeEveryInput.get(), this.state);
        }

        // grab the interval for storing the state to file
        if (storeEveryInput.get() > 0) {
            storeEvery = storeEveryInput.get();
        } else {
            storeEvery = state.m_storeEvery.get();
        }

        this.state.initialise();
        this.state.setPosterior(posteriorInput.get());

        // sanity check: all operator state nodes should be in the state
        final List<StateNode> stateNodes = this.state.stateNodeInput.get();
        for (final Operator op : operatorsInput.get()) {
            List<StateNode> nodes = op.listStateNodes();
            if (nodes.size() == 0) {
                    throw new RuntimeException("Operator " + op.getID() + " has no state nodes in the state. "
                                    + "Each operator should operate on at least one estimated state node in the state. "
                                    + "Remove the operator or add its statenode(s) to the state and/or set estimate='true'.");
                    // otherwise the chain may hang without obvious reason
            }
	        for (final StateNode stateNode : op.listStateNodes()) {
	            if (!stateNodes.contains(stateNode)) {
	                throw new RuntimeException("Operator " + op.getID() + " has a statenode " + stateNode.getID() + " in its inputs that is missing from the state.");
	            }
	        }
	    }
    
        // sanity check: at least one operator required to run MCMC
        if (operatorsInput.get().size() == 0) {
        	Log.warning.println("Warning: at least one operator required to run the MCMC properly, but none found.");
        }
        
        // sanity check: all state nodes should be operated on
        for (final StateNode stateNode : stateNodes) {
            if (!operatorStateNodes.contains(stateNode)) {
                Log.warning.println("Warning: state contains a node " + stateNode.getID() + " for which there is no operator.");
            }
        }
        state.m_storeEvery.setValue(Integer.MAX_VALUE, state);

        burnIn = 0;//burnInInput.get();
        chainLength = chainLengthInput.get();
        state.setEverythingDirty(true);
        posterior = posteriorInput.get();
    } // init
    
    
    
    /**
     * Run the next chain for this long (or -1 for auto)
     * @param runtime
     */
    public void setRuntime(long runtime) {
    	this.targetRuntime = runtime;
    }
    
    
    @Override
    public void run() throws IOException, SAXException, ParserConfigurationException {
    	

    	
    	// Regression
    	startTime = 0;
    	if (this.doRegression && !this.finishedRegression) {
    		
    		startTime = System.currentTimeMillis();
    		
    		// Sample chain length, measure runtime
    		if (targetRuntime <= 0) {
	    		chainLength = Randomizer.nextInt((int)(5*this.initChainLength));
	    		this.chainLengths[this.regressionNr] = chainLength;
    		}
    		

    		else {
    			chainLength = (long) 1e10;
    		}
    		
    	}
    	
    	
        // set up state (again). Other beastObjects may have manipulated the
        // StateNodes, e.g. set up bounds or dimensions
        state.initAndValidate();
        // also, initialise state with the file name to store and set-up whether to resume from file
//        state.setStateFileName(stateFileName);
//        operatorSchedule.setStateFileName(stateFileName);

  
        state.storeCalculationNodes();
        
        // do the sampling
        logAlpha = 0;
        debugFlag = Boolean.valueOf(System.getProperty("beast.debug"));

        doLoop();
        
        // Log the runtime
        if (this.doRegression && !this.finishedRegression) {
        	

        	long runTime = System.currentTimeMillis() - startTime;
        	this.runtimes[this.regressionNr] = runTime;
        	
        	
    		// Sample runtime, measure chain length
    		if (targetRuntime > 0) {
    			this.chainLengths[this.regressionNr] = numStepsInChain;
    			chainLength = numStepsInChain;
    			//Log.warning("Observed: " + runTime + "," + numStepsInChain + "(" + this.targetRuntime + ")");
    		}
        	
    		
        	this.regressionNr++;
        	if (this.regressionNr >= this.runtimes.length) {
        		this.finishedRegression = true;
        		this.trainRegression();
        	}
        }
        
    } // run;

    
    
    

    /**
     * main MCMC loop 
     * @throws IOException *
     */
    protected void doLoop() throws IOException {
        int corrections = 0;
        final boolean isStochastic = posterior.isStochastic();
        
        // make local state the current state
        int index = 0;
        for (StateNode stateNode : state.stateNodeInput.get()) {
        	stateNode.state = state;
    		otherStateNr[index] = stateNode.index;
        	stateNode.index = index++;
    	}
        oldLogLikelihood =  posterior.calculateLogP();
        //oldLogLikelihood = robustlyCalcPosterior(posterior);
        
        if (burnIn > 0) {
        	throw new IllegalArgumentException("Burnin should be 0");
        	///Log.warning.println("Please wait while BEAST takes " + burnIn + " pre-burnin samples");
        }
        
        numStepsInChain = 0;
        for (long sampleNr = sampleCount; sampleNr <= chainLength + sampleCount; sampleNr++) {
        	
        	
        	// Time over?
        	if (this.doRegression && !this.finishedRegression && this.targetRuntime > 0) {
        		long runTime = System.currentTimeMillis() - startTime;
        		if (runTime > this.targetRuntime) break;
        	}
        	
        	
            final Operator operator = propagateState(sampleNr);

            if (this.robust && (debugFlag && sampleNr % 3 == 0 || sampleNr % 10000 == 0)) {
                // check that the posterior is correctly calculated at every third
                // sample, as long as we are in debug mode
            	final double originalLogP = isStochastic ? posterior.getNonStochasticLogP() : oldLogLikelihood;
                final double logLikelihood = isStochastic ? state.robustlyCalcNonStochasticPosterior(posterior) : state.robustlyCalcPosterior(posterior);
                if (isTooDifferent(logLikelihood, originalLogP)) {
                    reportLogLikelihoods(posterior, "");
                    Log.err.println("At sample " + sampleNr + "\nLikelihood incorrectly calculated: " + originalLogP + " != " + logLikelihood
                    		+ "(" + (originalLogP - logLikelihood) + ")"
                            + " Operator: " + operator.getName());
                }
                if (sampleNr > NR_OF_DEBUG_SAMPLES * 3) {
                    // switch off debug mode once a sufficient large sample is checked
                    debugFlag = false;
                    if (isTooDifferent(logLikelihood, originalLogP)) {
                        // incorrect calculation outside debug period.
                        // This happens infrequently enough that it should repair itself after a robust posterior calculation
                        corrections++;
                        if (corrections > 100) {
                            // after 100 repairs, there must be something seriously wrong with the implementation
                        	Log.err.println("Too many corrections. There is something seriously wrong that cannot be corrected");
                            state.storeToFile(sampleNr);
                            operatorSchedule.storeToFile();
                            System.exit(1);
                        }
                        oldLogLikelihood = state.robustlyCalcPosterior(posterior);;
                    }
                } else {
                    if (isTooDifferent(logLikelihood, originalLogP)) {
                        // halt due to incorrect posterior during intial debug period
                        System.exit(1);
                    }
                }
            } else {
                if (sampleNr >= 0) {
                	operator.optimize(logAlpha);
                }
            }
            callUserFunction(sampleNr);

            if (posterior.getCurrentLogP() == Double.POSITIVE_INFINITY) {
            	throw new RuntimeException("Encountered a positive infinite posterior. This is a sign there may be numeric instability in the model.");
            }
            
            
            numStepsInChain++;
            
        }
        if (corrections > 0) {
        	Log.err.println("\n\nNB: " + corrections + " posterior calculation corrections were required. This analysis may not be valid!\n\n");
        }

        // restore to original State
        index = 0;
        for (StateNode stateNode : state.stateNodeInput.get()) {
        	stateNode.state = otherState;
    		stateNode.index = otherStateNr[index++];
        }
        sampleCount += chainLength;
    }

    /**
     * Perform a single MCMC propose+accept/reject step.
     *
     * @param sampleNr the index of the current MCMC step
     * @return the selected {@link beast.core.Operator}
     */
    protected Operator propagateState(final long sampleNr) {
        state.store(sampleNr);
//            if (m_nStoreEvery > 0 && sample % m_nStoreEvery == 0 && sample > 0) {
//                state.storeToFile(sample);
//            	operatorSchedule.storeToFile();
//            }

        final Operator operator = operatorSchedule.selectOperator();

        if (printDebugInfo) System.err.print("\n" + sampleNr + " " + operator.getName()+ ":");

        final Distribution evaluatorDistribution = operator.getEvaluatorDistribution();
        Evaluator evaluator = null;

        if (evaluatorDistribution != null) {
            evaluator = new Evaluator() {
                @Override
                public double evaluate() {
                    double logP = 0.0;

                    state.storeCalculationNodes();
                    state.checkCalculationNodesDirtiness();

                    try {
                        logP = evaluatorDistribution.calculateLogP();
                    } catch (Exception e) {
                        e.printStackTrace();
                        System.exit(1);
                    }

                    state.restore();
                    state.store(sampleNr);

                    return logP;
                }
            };
        }
        final double logHastingsRatio = operator.proposal(evaluator);

        if (logHastingsRatio != Double.NEGATIVE_INFINITY) {

            if (operator.requiresStateInitialisation()) {
                state.storeCalculationNodes();
                state.checkCalculationNodesDirtiness();
            }

            newLogLikelihood = posterior.calculateLogP();

            logAlpha = newLogLikelihood - oldLogLikelihood + logHastingsRatio; //CHECK HASTINGS
            if (printDebugInfo) System.err.print(logAlpha + " " + newLogLikelihood + " " + oldLogLikelihood);

            if (logAlpha >= 0 || Randomizer.nextDouble() < Math.exp(logAlpha)) {
                // accept
                oldLogLikelihood = newLogLikelihood;
                state.acceptCalculationNodes();

                if (sampleNr >= 0) {
                    operator.accept();
                }
                if (printDebugInfo) System.err.print(" accept");
            } else {
                // reject
                if (sampleNr >= 0) {
                    operator.reject(newLogLikelihood == Double.NEGATIVE_INFINITY ? -1 : 0);
                }
                state.restore();
                state.restoreCalculationNodes();
                if (printDebugInfo) System.err.print(" reject");
            }
            state.setEverythingDirty(false);
        } else {
            // operation failed
            if (sampleNr >= 0) {
                operator.reject(-2);
            }
            state.restore();
            if (!operator.requiresStateInitialisation()) {
                state.setEverythingDirty(false);
                state.restoreCalculationNodes();
            }
            if (printDebugInfo) System.err.print(" direct reject");
        }
        return operator;
    }

    private boolean isTooDifferent(double logLikelihood, double originalLogP) {
    	//return Math.abs((logLikelihood - originalLogP)/originalLogP) > 1e-6;
    	return Math.abs(logLikelihood - originalLogP) > 1e-6;
	}

    
    /**
     * Has regression finished?
     * @return
     */
    public boolean finishedRegression() {
    	return this.finishedRegression;
    }
    
    
    /**
     * Is regression being applied?
     * @return
     */
    public boolean isDoingRegression() {
    	return this.doRegression;
    }
    
    
    public double getRuntimeSlope() {
    	if (regression == null) return 0;
    	return regression.getSlope();
    }
    
    
    public double getRuntimeIntercept() {
    	if (regression == null) return 0;
    	return regression.getIntercept();
    }
    
    
    /**
     * Train the runtime vs chainLength model
     */
    private void trainRegression() {
    	
    	// Train the model (log space)
    	regression = new SimpleRegression();
    	for (int i = 0; i < Math.min(this.regressionNr, runtimes.length); i ++) {
    		double len = this.chainLengths[i];
    		double time = this.runtimes[i];
    		//if (len <= 0 || time <= 0) continue;
    		
    		//Log.warning("Training regression " + len + ", " + time);
    		
    		
    		regression.addData(len, time);
    		//regression.addData(Math.log(len), Math.log(time));
    		//Log.warning(this.chainLengths[i] + "," + this.runtimes[i]);
    	}
    	regression.regress();
    }
    
    /**
     * Use the regression model to predict runtime from chainLength
     * @return
     */
    public double predict(long chainLen) {
    	if (regression == null) return 0;
    	//return Math.exp(regression.predict(1.0*chainLen));
    	double result = regression.predict(1.0*chainLen);
    	if (Double.isNaN(result)) {
    		
    		Log.warning("Warning: regression predicted NaN for chain length " + chainLen + " (slope = " + regression.getSlope() + ", intercept = " + regression.getIntercept() + ")");
    		result = 1000;
    	}
    	return result;
    }
    
    
    
    /**
     * Apply the learned linear model between runtime and number of states
     * By setting the chain length such that the mean runtime is 'targetRuntime'
     */
    public void setChainlengthToTargetRuntime(double targetRuntime) {
    	
		// Avoid numerical errors
		final double minSlope = 1e-6;
		
    	double targetRunTime_log = targetRuntime; // Math.log(targetRuntime);
    	
    	if (this.appliedRegression) return;
    	if (!this.finishedRegression) return;
    	if (!this.doRegression) return;
    	
    	if (this.regression == null) this.trainRegression();
    	
    	
    	double slope = regression.getSlope();
		if (slope <= minSlope || slope == Double.NaN) slope = minSlope;
    	double intercept = regression.getIntercept();
    	if (intercept == Double.NaN) intercept = 0;
    	
    	// Apply the model
    	long targetChainlength_log = (long)((targetRunTime_log - intercept) / slope);
    	//this.chainLength = (long) Math.max(1, Math.exp(targetChainlength_log));
    	this.chainLength = (long) Math.max(MIN_CHAIN_LENGTH, targetChainlength_log);
    	Log.warning("Setting chain length to " + this.chainLength);
    	
    	
    	
    	// Print the model
    	double r2 = regression.getR();
    	Log.warning("Trained model: runtime(ms) = " + sf(intercept) + " + " + sf(slope) + "*chainLength\t\t(R2=" + sf(r2) + ")");
    	
    	
    	
    	this.appliedRegression = true;
    	
    }
    
    
    /**
     * Round to to 3sf
     * @param d
     * @return
     */
    private double sf(double d) {
    	BigDecimal bd = new BigDecimal(d);
    	bd = bd.round(new MathContext(3));
    	return bd.doubleValue();
    }

    
    /**
     * Get the chain length
     * @return
     */
	public long getChainLength() {
		return this.chainLength;
	}



    
    
    
    
    
}
