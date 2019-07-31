package starbeast3.inference;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintStream;
import java.nio.file.Files;
import java.nio.file.StandardCopyOption;
import java.util.concurrent.CountDownLatch;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.RejectedExecutionException;

import javax.xml.parsers.ParserConfigurationException;

import org.apache.commons.math3.stat.inference.*;
import org.xml.sax.SAXException;

import beast.app.BeastMCMC;
import beast.core.Description;
import beast.core.Input;
import beast.core.Logger;
import beast.core.Logger.LogFileMode;
import beast.core.MCMC;
import beast.core.StateNodeInitialiser;
import beast.core.util.Log;
import beast.util.XMLParser;
import beast.util.XMLParserException;
import beast.util.XMLProducer;


@Description("Runs short MCMC chains -- one per posterior sample -- "
		+ "from fixed start position till convergence is reached")
public class ShortChainsMCMC extends MCMC {
    final public Input<Integer> maxNrOfThreadsInput = new Input<>("shortThreads","maximum number of threads to use, if less than 1 the number of threads in BeastMCMC is used (default -1)", -1);
    final public Input<Integer> sampleCountInput = new Input<>("sampleCount","number of samples to take to represent posterior distribution", 100);
    final public Input<Double> pLevelInput = new Input<>("pLevel","significance level for testing posterior distribution has not changed", 0.5);
    final public Input<Double> multiplierInput = new Input<>("chainLengthMultiplier", "chain length multiplier, used to increase chain length after each iteration", 1.0);

	/** plugins representing MCMC with model, loggers, etc **/
	private ShortMCMC [] mcmcs;

	private int nrOfThreads;
	private static ExecutorService exec;
	private double pLevel;
	private int round;
	private double multiplier;

	
	@Override
	public void initAndValidate() {
		round = 0;
		pLevel = pLevelInput.get();
		nrOfThreads = Math.min(maxNrOfThreadsInput.get(), BeastMCMC.m_nThreads);
		// prevent any of the short chain MCMSs threading:
		BeastMCMC.m_nThreads = 1;
		exec = Executors.newFixedThreadPool(nrOfThreads);
		chainLength = chainLengthInput.get();

		mcmcs = new ShortMCMC[nrOfThreads];

		// the difference between the various chains is
		// 1. it runs an MCMC, not a ShortChainsMCMC
		// 2. remove shortThreads and sampleCount attributes
		// 3. output logs change for every chain
		// 4. log to stdout is removed to prevent clutter on stdout
		String sXML = new XMLProducer().toXML(this);
		sXML = sXML.replaceAll("shortThreads=[^ >]*", "");		
		sXML = sXML.replaceAll("pLevel=[^ >]*", "");
		sXML = sXML.replaceAll("sampleCount=[^ >]*", "");
		sXML = sXML.replaceAll("chainLengthMultiplier=[^ >]*", "");		
		sXML = sXML.replaceAll("logEvery=\"[0-9]*\"", "logEvery=\"" + Integer.MAX_VALUE + "\"");
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
		
		multiplier = multiplierInput.get();
		
		// create new chains
        countDown = new CountDownLatch(nrOfThreads);
		for (int i = 0; i < mcmcs.length; i++) {
            exec.execute(new CoreInit(sXML,i));
		}
		try {
			countDown.await();
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	
	} // initAndValidate
	
	
    class CoreInit implements Runnable {
    	String sXML;
    	int i;
    	
        CoreInit(String xml, int index) {
        	this.sXML = xml;
        	this.i = index;
        }

        @Override
		public void run() {
            try {
    			String sXML2 = sXML;
    			sXML2 = sXML2.replaceAll("fileName=\"","fileName=\"" + i);
    			if (sXML2.equals(sXML)) {
    				// Uh oh, no seed in log name => logs will overwrite
    				throw new IllegalArgumentException("Use $(seed) in log file name to guarantee log files do not overwrite");
    			}
    			try {
    				XMLParser parser = new XMLParser();
    				mcmcs[i] = (ShortMCMC) parser.parseFragment(sXML2, true);
    			} catch (XMLParserException e) {
    				throw new IllegalArgumentException(e);
    			}
    			// remove log to stdout, if any
    			for (int iLogger = mcmcs[i].loggersInput.get().size()-1; iLogger >= 0; iLogger--) {
    				if (mcmcs[i].loggersInput.get().get(iLogger).fileNameInput.get() == null) {
    					mcmcs[i].loggersInput.get().remove(iLogger);
    				}
    			}
            } catch (Exception e) {
                Log.err.println("Something went wrong initialising MCMC[" + i + "]");
                e.printStackTrace();
                System.exit(1);
            }
            countDown.countDown();
        }
        
    } // CoreInit
	
    class CoreRunnable implements Runnable {
    	ShortMCMC mcmc;
        double [] posterior;
        int start, end;
        
        CoreRunnable(ShortMCMC core, int start, int end, double [] posterior) {
            mcmc = core;
            this.start = start;
            this.end = end;
            this.posterior = posterior;
        }

        @Override
		public void run() {
            try {
            	for (int i = start; i < end; i++) {
            		mcmc.setChainLength(chainLength);
            		mcmc.setStateFile(stateFileName + i, true);
            		mcmc.run();
            		posterior[i] = mcmc.robustlyCalcPosterior(mcmc.posteriorInput.get());
            		Log.warning.print(i + " ");
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

    private double [] calculateUsingThreads(long chainLength) {
        try {

            int start = 0;
            double [] posterior = new double[sampleCountInput.get()];

            if (nrOfThreads > 1) {
                // kick off nrOfThreads-1 threads
	            countDown = new CountDownLatch(nrOfThreads - 1);
	            CoreRunnable [] coreRunnable = new CoreRunnable[nrOfThreads];
	            for (int i = 0; i < nrOfThreads-1; i++) {
	    			// need this to keep regression testing time reasonable
	                int end = sampleCountInput.get() * (i+1) / nrOfThreads;
	                coreRunnable[i] = new CoreRunnable(mcmcs[i], start, end, posterior);
	                start = end;
	                exec.execute(coreRunnable[i]);
	            }
            }
            
            
            // do the work of thread [nrOfThreads] in main thread
            ShortMCMC mcmc = mcmcs[nrOfThreads - 1];
        	for (int i = start; i < posterior.length; i++) {
        		mcmc.setStateFile(stateFileName + i, true);
        		mcmc.setChainLength(chainLength);
        		try {
					mcmc.run();
				} catch (IOException | SAXException | ParserConfigurationException e) {
					e.printStackTrace();
				}
        		posterior[i] = mcmc.robustlyCalcPosterior(mcmc.posteriorInput.get());
        		Log.warning.print(i + " ");
        		//Log.info("Round " + round + " sample " + i);
        	}

            
            if (nrOfThreads > 1) {
            	countDown.await();
            }
            Log.info("\nEnd of round " + round);
            
            return posterior;
        } catch (RejectedExecutionException | InterruptedException e) {
        	e.printStackTrace();
        	throw new RuntimeException(e);
        }
    }
	
	
	@Override 
	public void run() {
		long startTime = System.currentTimeMillis();
		
		
		if (Logger.FILE_MODE != Logger.LogFileMode.resume) {
			
			startStateInput.get().initialise();
			startStateInput.get().setStateFileName(stateFileName);	
            for (final StateNodeInitialiser initialiser : initialisersInput.get()) {
                initialiser.initStateNodes();
            }
            startStateInput.get().storeToFile(0);
		}
				
		
		
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
		
		chainLength = chainLengthInput.get();
		double [] newSample = calculateUsingThreads(chainLength);
		double [] oldSample;
		
		do {
			round++;
			oldSample = newSample;
			newSample = calculateUsingThreads(chainLength);
			chainLength *= multiplier;
		} while (!stop(oldSample, newSample));
		
		long endTime = System.currentTimeMillis();
        Log.info.println("Total calculation time: " + (endTime - startTime) / 1000.0 + " seconds");

        
        try {
			combineLogs();
		} catch (IOException e) {
			e.printStackTrace();
		}

	} // run
	
	private void combineLogs() throws IOException {
        // combine logs
        for (Logger log : loggersInput.get()) {
        	if (log.fileNameInput.get() != null) {
        		PrintStream out = new PrintStream(log.fileNameInput.get());
        		boolean isTreeLog = processHeader(log.fileNameInput.get(), out);
        		        		
        		int start = 0;
        		int sampleCount = 0;
        		for (int i = 0; i < nrOfThreads; i++) {
        			int lineCount = 0;
        	        BufferedReader fin = new BufferedReader(new FileReader(i + log.fileNameInput.get()));
        	        while (fin.ready()) {
        	        	lineCount++;
        	        	fin.readLine();
        	        }
        	        fin.close();        			

        	        fin = new BufferedReader(new FileReader(i + log.fileNameInput.get()));
	                int end = sampleCountInput.get() * (i+1) / nrOfThreads;
	                int samples = (end - start) * 2;
	                int skip = lineCount - samples - (isTreeLog ? 1 : 0);
	                for (int k = 0; k < skip; k++) {
	                	fin.readLine();
	                }
	                for (int k = 0; k < samples; k++) {
	                	String str = fin.readLine();
	                	// renumber
	                	if (isTreeLog) {
	                		str = "tree STATE_" + sampleCount + str.substring(str.indexOf('='));
	                	} else {
	                		str = sampleCount + str.substring(str.indexOf('\t'));
	                	}
	                	sampleCount++;
	                	out.println(str);
	                }
        			start = end;
        	        fin.close();
        		}        		
        		
        		out.close();
        		
//        		List<String> args = new ArrayList<>();
//        		args.add("-renumber");
//        		args.add("-burnin");
//        		args.add(100 * round / (round + 1) + "");        		
//        		args.add("-log");
//        		for (int i = 0; i < nrOfThreads; i++) {
//        			args.add(i + log.fileNameInput.get());
//        		}
//        		args.add("-o");
//        		args.add(log.fileNameInput.get());
//            	LogCombiner.main(args.toArray(new String[]{}));
        	}
        }		
	}


	private boolean processHeader(String fileName, PrintStream out) throws IOException {
        BufferedReader fin = new BufferedReader(new FileReader("0" + fileName));
        String str = fin.readLine();
        if (str.toUpperCase().startsWith("#NEXUS")) {
    		out.println(str);
        	while (fin.ready()) {
        		str = fin.readLine();
        		if (str.startsWith("tree STATE")) {
        			return true;
        		}
        		out.println(str);
        	}
        	fin.close();
        	return true;
        }
        if (!str.startsWith("#")) {
    		out.println(str);
        	fin.close();
    		return false;
        }
    	while (fin.ready()) {
    		str = fin.readLine();
    		if (!str.startsWith("#")) {
        		out.println(str);
    	    	fin.close();
    			return false;
    		}
    	}		
    	fin.close();
		return false;
	}


	final static boolean debug = true;
	
	boolean stop(double [] sample1, double [] sample2) {
		if (debug) {
			saveFile(sample1, "/tmp/sample1.log");
			saveFile(sample2, "/tmp/sample2.log");
		}
		
		
		KolmogorovSmirnovTest test = new KolmogorovSmirnovTest();
		double pValue = test.kolmogorovSmirnovTest(sample1, sample2);
		Log.info("pValue = " + pValue);
		return (pValue > pLevel);
		
//		MannWhitneyUTest test = new MannWhitneyUTest();
//		double pValue = test.mannWhitneyUTest(sample1, sample2);
//		return (pValue < 0.05);
	}


	private void saveFile(double[] sample2, String file) {
		try {
			StringBuilder log = new StringBuilder();
			log.append("sampleNr\tvalue\n");
			for (int i = 0; i < sample2.length; i++) {
				log.append(i + "\t" + sample2[i] + "\n");
			}
			FileWriter outfile = new FileWriter(new File(file));
			outfile.write(log.toString());
			outfile.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
		
	}
	
}
