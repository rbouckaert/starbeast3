package starbeast3.inference;

import java.io.IOException;
import java.io.OutputStream;
import java.io.PrintStream;
import java.util.Collections;

import javax.xml.parsers.ParserConfigurationException;

import org.xml.sax.SAXException;

import beast.base.core.Description;
import beast.base.inference.Logger;
import beast.base.inference.MCMC;
import beast.base.core.Log;

@Description("Runs MCMC silently inside ShortChainMCMC -- do not use directly")
public class ShortMCMC extends MCMC {
	private PrintStream nullstream = null;
	
    @Override
    public void run() throws IOException, SAXException, ParserConfigurationException {
        // set up state (again). Other beastObjects may have manipulated the
        // StateNodes, e.g. set up bounds or dimensions
        state.initAndValidate();
        // also, initialise state with the file name to store and set-up whether to resume from file
        state.setStateFileName(stateFileName);
        operatorSchedule.setStateFileName(stateFileName);

        burnIn = burnInInput.get();
        // chainLength = chainLengthInput.get();
        state.setEverythingDirty(true);
        posterior = posteriorInput.get();
        
        Log.setLevel(Log.Level.error);
        state.restoreFromFile();
        Log.setLevel(Log.Level.info);


        if (nullstream == null) {
        	nullstream = new PrintStream(new OutputStream() {
        		
        		@Override
        		public void write(int b) throws IOException {
        		}
        	});

        }
        
        PrintStream orgErr = System.err;
//        System.setErr(nullstream);
//        operatorSchedule.restoreFromFile();
//        System.setErr(orgErr);
        
        burnIn = 0;
        oldLogLikelihood = state.robustlyCalcPosterior(posterior);

        state.storeCalculationNodes();

        
        // do the sampling
        logAlpha = 0;
        debugFlag = Boolean.valueOf(System.getProperty("beast.debug"));

        if (Double.isInfinite(oldLogLikelihood) || Double.isNaN(oldLogLikelihood)) {
            reportLogLikelihoods(posterior, "");
            throw new RuntimeException("Could not find a proper state to initialise. Perhaps try another seed.\nSee http://www.beast2.org/2018/07/04/fatal-errors.html for other possible solutions.");
        }

        loggers = loggersInput.get();

        // put the loggers logging to stdout at the bottom of the logger list so that screen output is tidier.
        Collections.sort(loggers, (o1, o2) -> {
            if (o1.isLoggingToStdout()) {
                return o2.isLoggingToStdout() ? 0 : 1;
            } else {
                return o2.isLoggingToStdout() ? -1 : 0;
            }
        });

        // initialises log so that log file headers are written, etc.
        for (final Logger log : loggers) {
            log.init();
            ((ShortMCMCLogger) log).setEvery((int) chainLength);
        }

        doLoop();

        close();

        state.storeToFile(chainLength);
        operatorSchedule.storeToFile();
    } // run;

	public void setChainLength(long chainLength) {
		this.chainLength = chainLength;
	}

}
