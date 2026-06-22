/**
 * @author Huw A. Ogilvie
 * @author Jordan Douglas
 */



package starbeast3.evolution.branchratemodel;


import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.util.InputUtil;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.evolution.branchratemodel.Base;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.inference.StateNode;
import beast.base.util.Randomizer;

public class UCRelaxedClockModelSB3 extends Base implements BranchRateModelSB3 {
    
	
	final public Input<TreeInterface> treeInput = new Input<>("tree", "(Species) tree to apply per-branch rates to.", Input.Validate.REQUIRED);
    final public Input<Boolean> estimateRootInput = new Input<>("estimateRoot", "Estimate rate of the root branch.", false);
    final public Input<Boolean> noCacheInput = new Input<>("noCache", "Always recalculate branch rates.", false);
    final public Input<RealScalarParam<PositiveReal>> stdevInput = new Input<>("stdev", "Standard deviation of the log-normal distribution for branch rates. If not supplied uses exponential.");
    final public Input<RealVectorParam<NonNegativeReal>> realRatesInput = new Input<>("realRates", "The real rates associated with nodes in the species tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    
    
    final private double MEAN_CLOCK_RATE = 1.0; // Mean clock rate. Equal to 1/lambda for exponential, or exp(M + S^2/2) for lognormal
    
    public enum Mode {
        categories,
        rates,
        quantiles
    }
    
    public enum RateDistribution {
        exponential,
        lognormal,
        parametric
    }
    
    Mode mode = Mode.categories;
    RateDistribution rateDistribution = RateDistribution.exponential;
    
    private int nBins;
    private double currentLogNormalStdev;
    private double storedLogNormalStdev;
    private double[] binRates;
    private double[] storedBinRates;
    private double[] ratesArray;
    private double[] storedRatesArray;

    private int nEstimatedRates;
    private int rootNodeNumber;
    private boolean estimateRoot;
    private boolean noCache;
    private boolean needsUpdate;
    private boolean binRatesNeedsUpdate;
    
    private RealVectorParam<NonNegativeReal> realRates;
    
    
    

    @Override
    public void initAndValidate() {
        realRates = realRatesInput.get();
        final TreeInterface speciesTree = treeInput.get();
        final Node[] speciesNodes = speciesTree.getNodesAsArray();
        estimateRoot = estimateRootInput.get().booleanValue();
        noCache = noCacheInput.get().booleanValue();
        rootNodeNumber = speciesTree.getRoot().getNr();

        
        // Get the rates and their parameterisation model
        mode = Mode.rates;
        
        
        // Exponential, lognormal, or parametric rate distribution
    	rateDistribution = stdevInput.get() == null ? RateDistribution.exponential : RateDistribution.lognormal;
       

        if (estimateRoot) {
            nEstimatedRates = speciesNodes.length;
        } else {
            nEstimatedRates = speciesNodes.length - 1;
        }


        currentLogNormalStdev = -1.0; 
        storedLogNormalStdev = -1.0;


        ratesArray = new double[speciesNodes.length];
        storedRatesArray = new double[speciesNodes.length];
        
        switch(mode) {
        
        
	        
	        // Initialise real branch rates
	        case rates: {
	        	
	        	if (realRates.size() != nEstimatedRates) {
	        		realRates.setDimension(nEstimatedRates);
	        		
	        		// Randomly draw rates from the appropriate distribution
	        		switch(rateDistribution) {
	        	        
	        			// Initialise continuous exponential rates
	             		case exponential: {
	             			for (int i = 0; i < nEstimatedRates; i++) {
	             				double val = Randomizer.nextExponential(1 / MEAN_CLOCK_RATE);
	             				realRates.set(i, val);
	             			}
	             			break;
	             		}
	             		
	             		// Initialise continuous lognormal rates
	             		case lognormal: {
	             			
	             			currentLogNormalStdev = stdevInput.get().get();
	    	            	storedLogNormalStdev = currentLogNormalStdev;
	             			
	             			// Mean in log space
	                        final double M = Math.log(MEAN_CLOCK_RATE) - (0.5 * currentLogNormalStdev * currentLogNormalStdev);
	    	        		for (int i = 0; i < nEstimatedRates; i++) {
	    	        			double val = Randomizer.nextLogNormal(M, currentLogNormalStdev, false);
	    	        			realRates.set(i, val);
	    	        		}
	             			
	             			break;
	             		}
	        		
	        		 }
				    
				}
	        	
	        	binRatesNeedsUpdate = false;
	        	break;
	        }
	        
	        
        }
        
        
      
        needsUpdate = true;
        
    }
    
    
    
    

    @Override
    public boolean requiresRecalculation() {
    	
    	// If lognormal S changes then branch rates require recalculation from bins
        if (rateDistribution == RateDistribution.lognormal && mode == Mode.categories) {
            final double proposedLogNormalStdev = stdevInput.get().get();
            if (proposedLogNormalStdev != currentLogNormalStdev) {
                binRatesNeedsUpdate = true;
            } 
        }
       
        needsUpdate = binRatesNeedsUpdate || InputUtil.isDirty(meanRateInput) || InputUtil.isDirty(realRatesInput);
        
        return needsUpdate;
    }

    @Override
    public void store() {
        storedLogNormalStdev = currentLogNormalStdev;
        System.arraycopy(ratesArray, 0, storedRatesArray, 0, ratesArray.length);
        
        // Store the bin rates only if using categories/quantiles
        if (mode == Mode.categories || mode == Mode.quantiles) System.arraycopy(binRates, 0, storedBinRates, 0, binRates.length);
        
        super.store();
    }

    @Override
    public void restore() {
        double tmpLogNormalStdev = currentLogNormalStdev;
        double[] tmpRatesArray = ratesArray;

        currentLogNormalStdev = storedLogNormalStdev;
        ratesArray = storedRatesArray;
        
        storedLogNormalStdev = tmpLogNormalStdev;
        storedRatesArray = tmpRatesArray;
        
        

        super.restore();
    }


    private void update() {
    	


        Double estimatedMean;
        if (meanRateInput.get() == null) {
            estimatedMean = 1.0;
        } else {
            estimatedMean = meanRateInput.get().get();
        }

        
        
        // Multiply the raw rate by the clock rate
        switch(mode){
        
	        // Real numbers
	        case rates: {
	        	
	            for (int i = 0; i < nEstimatedRates; i++) {
	                ratesArray[i] = estimatedMean * realRates.get(i);
	            }
	        	
	        	break;
	        }
        
        
        }
        
       
        binRatesNeedsUpdate = false;
        if (!estimateRoot) ratesArray[rootNodeNumber] = estimatedMean;

    }

    @Override
    public double[] getRatesArray() {
        if (needsUpdate || noCache) {
            synchronized (this) {
                update();
                needsUpdate = false;
            }
        }
        


        return ratesArray;
    }
    
    
    @Override
    public double getRateForBranch(Node node) {
        if (needsUpdate || noCache) {
            synchronized (this) {
                update();
                needsUpdate = false;
            }
        }

        assert ratesArray[node.getNr()] > 0.0;
        return ratesArray[node.getNr()];
        
    }


	public Mode getRateMode() {
		return this.mode;
	}
	

    
    
}




