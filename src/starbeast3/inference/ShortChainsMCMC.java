package starbeast3.inference;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;

import org.apache.commons.math3.stat.inference.*;

import beast.app.BeastMCMC;
import beast.app.tools.LogCombiner;
import beast.core.Description;
import beast.core.Input;
import beast.core.Logger;
import beast.core.Logger.LogFileMode;
import beast.core.MCMC;
import beast.core.util.Log;
import beast.util.XMLParser;
import beast.util.XMLParserException;
import beast.util.XMLProducer;


@Description("Runs short MCMC chains -- one per posterior sample -- "
		+ "from fixed start position till convergence is reached")
public class ShortChainsMCMC extends MCMC {
    final public Input<Integer> maxNrOfThreadsInput = new Input<>("shortThreads","maximum number of threads to use, if less than 1 the number of threads in BeastMCMC is used (default -1)", -1);
    final public Input<Integer> sampleCountInput = new Input<>("sampleCount","number of samples to take to represent posterior distribution", 100);
    final public Input<Double> pLevelInput = new Input<>("pLevel","significance level for testing posterior distribution has not changed", 0.05);

	/** plugins representing MCMC with model, loggers, etc **/
	private MCMC [] mcmcs;

	private int nrOfThreads;
	private static ExecutorService exec;
	private double pLevel;
	private int round;

	
	@Override
	public void initAndValidate() {
		round = 0;
		pLevel = pLevelInput.get();
		nrOfThreads = Math.min(maxNrOfThreadsInput.get(), BeastMCMC.m_nThreads);
		exec = Executors.newFixedThreadPool(nrOfThreads);
		chainLength = chainLengthInput.get();

		mcmcs = new MCMC[nrOfThreads];

		// the difference between the various chains is
		// 1. it runs an MCMC, not a ShortChainsMCMC
		// 2. remove shortThreads and sampleCount attributes
		// 3. output logs change for every chain
		// 4. log to stdout is removed to prevent clutter on stdout
		String sXML = new XMLProducer().toXML(this);
		sXML = sXML.replaceAll("shortThreads=[^ >]*", "");		
		sXML = sXML.replaceAll("pLevel=[^ >]*", "");
		sXML = sXML.replaceAll("sampleCount=[^ >]*", "");
		String ShortMCMCLogger = ShortMCMCLogger.class.getName();
		sXML = sXML.replaceAll("spec=\"Logger\"", "spec=\"" + ShortMCMCLogger + "\"");
		String ShortChainsMCMC = this.getClass().getName();
		while (ShortChainsMCMC.length() > 0) {
			sXML = sXML.replaceAll("\\b"+ShortChainsMCMC+"\\b", ShortMCMC.class.getName());
			if (ShortChainsMCMC.indexOf('.') >= 0) {
				ShortChainsMCMC = ShortChainsMCMC.substring(ShortChainsMCMC.indexOf('.')+1);
			} else {
				ShortChainsMCMC = "";
			}
		}

		
		// create new chains
		XMLParser parser = new XMLParser();
		for (int i = 0; i < mcmcs.length; i++) {
			String sXML2 = sXML;
			sXML2 = sXML2.replaceAll("fileName=\"","fileName=\"" + i);
			if (sXML2.equals(sXML)) {
				// Uh oh, no seed in log name => logs will overwrite
				throw new IllegalArgumentException("Use $(seed) in log file name to guarantee log files do not overwrite");
			}
			try {
				mcmcs[i] = (MCMC) parser.parseFragment(sXML2, true);
			} catch (XMLParserException e) {
				throw new IllegalArgumentException(e);
			}
			// remove log to stdout, if any
			for (int iLogger = mcmcs[i].loggersInput.get().size()-1; iLogger >= 0; iLogger--) {
				if (mcmcs[i].loggersInput.get().get(iLogger).fileNameInput.get() == null) {
					mcmcs[i].loggersInput.get().remove(iLogger);
				} else {
					// only log after chainLengthInput.get() samples
					Logger l = mcmcs[i].loggersInput.get().get(iLogger);
					l.everyInput.setValue((int) chainLength, l);
				}
			}
		}
	
	} // initAndValidate
	
	
    class CoreRunnable implements Runnable {
        MCMC mcmc;
        double [] posterior;
        int start, end;
        
        CoreRunnable(MCMC core, int start, int end, double [] posterior) {
            mcmc = core;
            this.start = start;
            this.end = end;
            this.posterior = posterior;
        }

        @Override
		public void run() {
            try {
            	for (int i = start; i < end; i++) {
            		mcmc.setStateFile(stateFileName + i, true);
            		mcmc.run();
            		posterior[i] = mcmc.robustlyCalcPosterior(mcmc.posteriorInput.get());
            		//Log.info("Round " + round + " sample " + i);
            	}
            } catch (Exception e) {
                Log.err.println("Something went wrong in a calculation of " + mcmc.getID());
                e.printStackTrace();
                System.exit(1);
            }
            countDown.countDown();
        }
        
    } // CoreRunnable

    CountDownLatch countDown;

    private double [] calculateLogPUsingThreads(long chainLength) {
        try {

            countDown = new CountDownLatch(nrOfThreads);
            // kick off the threads
            int start = 0;
            CoreRunnable [] coreRunnable = new CoreRunnable[nrOfThreads];
            double [] posterior = new double[sampleCountInput.get()];
            for (int i = 0; i < nrOfThreads; i++) {
    			// need this to keep regression testing time reasonable
                int end = sampleCountInput.get() * (i+1) / nrOfThreads;
                coreRunnable[i] = new CoreRunnable(mcmcs[i], start, end, posterior);
                start = end;
                exec.execute(coreRunnable[i]);
            }
            countDown.await();
            Log.info("End of round " + round);
            
            return posterior;
        } catch (RejectedExecutionException | InterruptedException e) {
        	e.printStackTrace();
        	throw new RuntimeException(e);
        }
    }
	
	
	@Override 
	public void run() {
		long startTime = System.currentTimeMillis();
		// copy starting state
		for (int i = 0; i < sampleCountInput.get(); i++) {
			try {
		    	Files.copy(new File(stateFileName).toPath(), new File(stateFileName+i).toPath(), StandardCopyOption.REPLACE_EXISTING);
			} catch (IOException e) {
				e.printStackTrace();
				throw new IllegalArgumentException(e);
			}
		}

		Logger.FILE_MODE = LogFileMode.resume;
		
		long chainLength = chainLengthInput.get();
		double [] newSample = calculateLogPUsingThreads(chainLength);
		double [] oldSample;
		
		do {
			round++;
			oldSample = newSample;
			newSample = calculateLogPUsingThreads(chainLength);
			//chainLength *= 2;
		} while (!stop(oldSample, newSample));
		
		long endTime = System.currentTimeMillis();
        Log.info.println("Total calculation time: " + (endTime - startTime) / 1000.0 + " seconds");

        // combine logs
        for (Logger log : loggersInput.get()) {
        	if (log.fileNameInput.get() != null) {
        		List<String> args = new ArrayList<>();
        		args.add("-renumber");
        		args.add("-burnin");
        		args.add(100 * round / (round + 1) + "");        		
        		args.add("-log");
        		for (int i = 0; i < nrOfThreads; i++) {
        			args.add(i + log.fileNameInput.get());
        		}
        		args.add("-o");
        		args.add(log.fileNameInput.get());
            	LogCombiner.main(args.toArray(new String[]{}));
        	}
        }
	} // run
	
	
	
	boolean stop(double [] sample1, double [] sample2) {
		KolmogorovSmirnovTest test = new KolmogorovSmirnovTest();
		double pValue = test.kolmogorovSmirnovTest(sample1, sample2);
		Log.info("pValue = " + pValue);
		return (pValue < pLevel);
		
//		MannWhitneyUTest test = new MannWhitneyUTest();
//		double pValue = test.mannWhitneyUTest(sample1, sample2);
//		return (pValue < 0.05);
	}
	
}
