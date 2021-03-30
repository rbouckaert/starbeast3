package starbeast3.evolution.likelihood;

import java.io.PrintStream;

import beast.app.BeastMCMC;
import beast.core.Description;
import beast.core.util.Log;
import beast.evolution.likelihood.ThreadedTreeLikelihood;
import beast.evolution.likelihood.TreeLikelihood;


@Description("A tree likelihood which can act as either threaded or unthreaded")
public class MetaTreeLikelihood extends TreeLikelihood {
	
	protected boolean usingThreads;
	protected ThreadedTreeLikelihood threadedLikelihood;
	
	 
	 
	@Override
    public void initAndValidate() {
		 
		 super.initAndValidate();
		 
		 // Initialise the threaded tree likelihood
		 /*
		 this.threadedLikelihood = new ThreadedTreeLikelihood();
		 this.threadedLikelihood.initByName( "data", this.dataInput.get(),
				 						"tree", this.treeInput.get(),
				 						"siteModel", this.siteModelInput.get(),
				 						"branchRateModel", this.branchRateModelInput.get());
		 */
		 this.usingThreads = false;
		 if (BeastMCMC.m_nThreads == 1) this.usingThreads = false;
		 
	}
	
	
	
	/**
	 * Start using threads
	 */
	public void startThreading() {
		if (BeastMCMC.m_nThreads == 1) {
			this.usingThreads = false;
			return;
		}
		if (this.usingThreads) return;
		this.usingThreads = true;
		
		//this.threadedLikelihood.has
		super.store();
		
	}
	
	
	/**
	 * Stop using threads
	 */
	public void stopThreading() {
		if (!this.usingThreads) return;
		this.usingThreads = false;
		
		//this.hasDirt = 1;
		//this.calculateLogP();
		
		this.threadedLikelihood.store();
		
	}
	
	
	@Override
    public double calculateLogP() {
		
		// Log.warning("Initial L1 = " + super.calculateLogP());
		// Log.warning("Initial L2 = " + this.threadedLikelihood.calculateLogP());
		// System.exit(1);
		 
		
		if (this.usingThreads) return this.threadedLikelihood.calculateLogP();
		else return super.calculateLogP();
	}
	
	@Override
    public double getCurrentLogP() {
		if (this.usingThreads) return this.threadedLikelihood.getCurrentLogP();
		else return super.getCurrentLogP();
    }

	@Override
    public double getStoredLogP() {
		if (this.usingThreads) return this.threadedLikelihood.getStoredLogP();
		else return super.getStoredLogP();
    }
	

	@Override
    public void store() {
		if (this.usingThreads) this.threadedLikelihood.store();
		else super.store();
	}
	
	
	@Override
    public void restore() {
		if (this.usingThreads) this.threadedLikelihood.restore();
		else super.restore();
	}
	
	
	
	
	
	
	
	@Override
    public void init(final PrintStream out) {
		if (this.usingThreads) this.threadedLikelihood.init(out);
		else super.init(out);
    }

    @Override
    public void log(final long sample, final PrintStream out) {
    	if (this.usingThreads) this.threadedLikelihood.log(sample, out);
		else super.log(sample, out);
    }




    @Override
    public int getDimension() {
        return 1;
    }

    @Override
    public double getArrayValue() {
    	if (this.usingThreads) return this.threadedLikelihood.getArrayValue();
		else return super.getArrayValue();
    }

    @Override
    public double getArrayValue(final int dim) {
    	if (this.usingThreads) return this.threadedLikelihood.getArrayValue(dim);
		else return super.getArrayValue(dim);
    }

    @Override
    public boolean isStochastic() {
    	if (this.usingThreads) return this.threadedLikelihood.isStochastic();
		else return super.isStochastic();
    }


    @Override
    public double getNonStochasticLogP() {
    	if (this.usingThreads) return this.threadedLikelihood.getNonStochasticLogP();
		else return super.getNonStochasticLogP();
    }
	
	
}

