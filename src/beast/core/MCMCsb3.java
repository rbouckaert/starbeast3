package beast.core;

import java.io.IOException;

import beast.core.util.Evaluator;
import beast.core.util.Log;
import beast.util.Randomizer;
import starbeast3.operators.MultiStepOperator;
import starbeast3.operators.MultiStepOperatorSchedule;

@Description("Variant of MCMC that takes MultiStepOperators in account")
public class MCMCsb3 extends MCMC {
	
	@Override
    protected void doLoop() throws IOException {
		if (!(operatorSchedule instanceof MultiStepOperatorSchedule)) {
			throw new IllegalArgumentException(getClass().getName() + " should be used with " + MultiStepOperatorSchedule.class.getName() + " as operator schedule");
		}
		
		
        int corrections = 0;
        final boolean isStochastic = posterior.isStochastic();
        
        if (burnIn > 0) {
        	Log.warning.println("Please wait while BEAST takes " + burnIn + " pre-burnin samples");
        }
        for (long sampleNr = -burnIn; sampleNr <= chainLength; sampleNr++) {
        	
            final Operator operator = propagateState(sampleNr);

            if (debugFlag && sampleNr % 2 == 0 || sampleNr % 10000 == 0) {
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
                        state.storeToFile(sampleNr);
                        operatorSchedule.storeToFile();
                        System.exit(1);
                    }
                }
            } else {
                if (sampleNr >= 0) {
                	operator.optimize(logAlpha);
                }
            }
            callUserFunction(sampleNr);

            // make sure we always save just before exiting
            if (storeEvery > 0 && (sampleNr + 1) % storeEvery == 0 || sampleNr == chainLength) {
                /*final double logLikelihood = */
                state.robustlyCalcNonStochasticPosterior(posterior);
                state.storeToFile(sampleNr);
                operatorSchedule.storeToFile();
            }
            
            if (posterior.getCurrentLogP() == Double.POSITIVE_INFINITY) {
            	throw new RuntimeException("Encountered a positive infinite posterior. This is a sign there may be numeric instability in the model.");
            }
            
            if (operator instanceof MultiStepOperator) {
            	sampleNr += ((MultiStepOperator)operator).stepCount() - 1;
            }            
        }
        if (corrections > 0) {
        	Log.err.println("\n\nNB: " + corrections + " posterior calculation corrections were required. This analysis may not be valid!\n\n");
        }
    }
	
	
    protected Operator propagateState(final long sampleNr) {
        state.store(sampleNr);
                
        final Operator operator = ((MultiStepOperatorSchedule) operatorSchedule).selectOperator(sampleNr);

        // if (printDebugInfo) System.err.print("\n" + sampleNr + " " + operator.getName()+ ":");

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
            // reportLogLikelihoods(posterior, "");

            logAlpha = newLogLikelihood - oldLogLikelihood + logHastingsRatio; //CHECK HASTINGS
            // if (printDebugInfo) System.err.print(logAlpha + " " + newLogLikelihood + " " + oldLogLikelihood);

            if (logAlpha >= 0 || Randomizer.nextDouble() < Math.exp(logAlpha)) {
                // accept
                oldLogLikelihood = newLogLikelihood;
                state.acceptCalculationNodes();

                if (sampleNr >= 0) {
                    operator.accept();
                }
                // if (printDebugInfo) System.err.print(" accept");
            } else {
                // reject
                if (sampleNr >= 0) {
                    operator.reject(newLogLikelihood == Double.NEGATIVE_INFINITY ? -1 : 0);
                }
                state.restore();
                state.restoreCalculationNodes();
                // if (printDebugInfo) System.err.print(" reject");
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
            // if (printDebugInfo) System.err.print(" direct reject");
        }
        log(sampleNr);
        return operator;
    }
    
    private boolean isTooDifferent(double logLikelihood, double originalLogP) {
    	return Math.abs((logLikelihood - originalLogP)/originalLogP) > 1e-6;
    	//return Math.abs(logLikelihood - originalLogP) > 1e-6;
	}

}
