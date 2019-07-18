package starbeast3.inference;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.*;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;
import java.util.stream.Stream;

import org.apache.commons.math3.stat.inference.*;

import beast.app.BeastMCMC;
import beast.core.Description;
import beast.core.Input;
import beast.core.Logger;
import beast.core.MCMC;
import beast.core.util.Log;
import beast.util.XMLParser;
import beast.util.XMLParserException;
import beast.util.XMLProducer;


@Description("Runs short MCMC chains -- one per posterior sample -- "
		+ "from fixed start position till convergence is reached")
public class ShortChainsMCMC extends MCMC {
    final public Input<Integer> maxNrOfThreadsInput = new Input<>("threads","maximum number of threads to use, if less than 1 the number of threads in BeastMCMC is used (default -1)", -1);
    final public Input<Integer> sampleCountInput = new Input<>("sampleCount","number of samples to take to represent posterior distribution", 100);

	/** plugins representing MCMC with model, loggers, etc **/
	private MCMC [] m_chains;

	private int nrOfThreads;
	private static ExecutorService exec;

	
	@Override
	public void initAndValidate() {
		
		nrOfThreads = Math.min(maxNrOfThreadsInput.get(), BeastMCMC.m_nThreads);
		exec = Executors.newFixedThreadPool(nrOfThreads);

		m_chains = new MCMC[nrOfThreads];

		// the difference between the various chains is
		// 1. it runs an MCMC, not a  MultiplMCMC
		// 2. remove chains attribute
		// 3. output logs change for every chain
		// 4. log to stdout is removed to prevent clutter on stdout
		String sXML = new XMLProducer().toXML(this);
		sXML = sXML.replaceAll("chains=[^ /]*", "");
		String sMultiMCMC = this.getClass().getName();
		while (sMultiMCMC.length() > 0) {
			sXML = sXML.replaceAll("\\b"+sMultiMCMC+"\\b", MCMC.class.getName());
			if (sMultiMCMC.indexOf('.') >= 0) {
				sMultiMCMC = sMultiMCMC.substring(sMultiMCMC.indexOf('.')+1);
			} else {
				sMultiMCMC = "";
			}
		}

		// create new chains
		XMLParser parser = new XMLParser();
		for (int i = 0; i < m_chains.length; i++) {
			String sXML2 = sXML;
			sXML2 = sXML2.replaceAll("\\$\\(seed\\)", i+"");
			if (sXML2.equals(sXML)) {
				// Uh oh, no seed in log name => logs will overwrite
				throw new IllegalArgumentException("Use $(seed) in log file name to guarantee log files do not overwrite");
			}
			try {
				m_chains[i] = (MCMC) parser.parseFragment(sXML2, true);
            	m_chains[i].setStateFile(stateFileName + "." + i, true);
			} catch (XMLParserException e) {
				throw new IllegalArgumentException(e);
			}
			// remove log to stdout, if any
			for (int iLogger = m_chains[i].loggersInput.get().size()-1; iLogger >= 0; iLogger--) {
				if (m_chains[i].loggersInput.get().get(iLogger).fileNameInput.get() == null) {
					m_chains[i].loggersInput.get().remove(iLogger);
				} else {
					Logger l = m_chains[i].loggersInput.get().get(iLogger);
					l.everyInput.setValue(chainLengthInput.get(), l);
				}
			}
		}
	
	} // initAndValidate
	
	
    class CoreRunnable implements Runnable {
        MCMC mcmc;
        double [] posterior;
        
        CoreRunnable(MCMC core, int sampleCount) {
            mcmc = core;
            posterior = new double[sampleCount];
        }

        @Override
		public void run() {
            try {
            	for (int i = 0; i < posterior.length; i++) {
            		mcmc.run();
            		posterior[i] = mcmc.robustlyCalcPosterior(mcmc.posteriorInput.get());
            	}
            } catch (Exception e) {
                Log.err.println("Something went wrong in a calculation of " + mcmc.getID());
                e.printStackTrace();
                System.exit(1);
            }
            countDown.countDown();
        }

        double [] getPosterior() {
        	return posterior;
        }
        
    } // CoreRunnable

    CountDownLatch countDown;

    private double [] calculateLogPUsingThreads(long chainLength) {
        try {

            countDown = new CountDownLatch(nrOfThreads);
            // kick off the threads
            int start = 0;
            CoreRunnable [] coreRunnable = new CoreRunnable[nrOfThreads];
            for (int i = 0; i < nrOfThreads; i++) {
            	Files.copy(new File(stateFileName).toPath(), new File(stateFileName+i).toPath());
    			// need this to keep regression testing time reasonable
                int end = sampleCountInput.get() * (i+1) / nrOfThreads;
                coreRunnable[i] = new CoreRunnable(m_chains[i], end - start);
                start = end;
                exec.execute(coreRunnable[i]);
            }
            countDown.await();
            List<Double> posterior = new ArrayList<>(sampleCountInput.get());
            for (CoreRunnable core : coreRunnable) {
                double [] sample = core.getPosterior();
                for (double d : sample) {
                	posterior.add(d);
                }
            }
            
            return Stream.of(posterior.toArray(new Double[]{})).mapToDouble(Double::doubleValue).toArray();
        } catch (RejectedExecutionException | InterruptedException | IOException e) {
        	e.printStackTrace();
        	throw new RuntimeException(e);
        }
    }
	
	
	@Override 
	public void run() {
		long chainLength = chainLengthInput.get();
		double [] newSample = calculateLogPUsingThreads(chainLength);
		double [] oldSample;
		
		do {
			oldSample = newSample;
			newSample = calculateLogPUsingThreads(chainLength);
			chainLength *= 2;
		} while (!stop(oldSample, newSample));
	} // run
	
	
	
	boolean stop(double [] sample1, double [] sample2) {
		KolmogorovSmirnovTest test = new KolmogorovSmirnovTest();
		double pValue = test.kolmogorovSmirnovTest(sample1, sample2);
		return (pValue < 0.05);
		
//		MannWhitneyUTest test = new MannWhitneyUTest();
//		double pValue = test.mannWhitneyUTest(sample1, sample2);
//		return (pValue < 0.05);
	}
	
}
