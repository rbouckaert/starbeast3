/**
 * @author Huw A. Ogilvie
 * @author Jordan Douglas
 */



package starbeast3.evolution.branchratemodel;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ExponentialDistribution;
import org.apache.commons.math.distribution.ExponentialDistributionImpl;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.util.Randomizer;
import beast.evolution.branchratemodel.BranchRateModel;

public class UCRelaxedClockModelSB3 extends BranchRateModel.Base implements BranchRateModelSB3 {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "(Species) tree to apply per-branch rates to.", Input.Validate.REQUIRED);
    final public Input<Boolean> estimateRootInput = new Input<>("estimateRoot", "Estimate rate of the root branch.", false);
    final public Input<Boolean> noCacheInput = new Input<>("noCache", "Always recalculate branch rates.", false);
    final public Input<RealParameter> stdevInput = new Input<>("stdev", "Standard deviation of the log-normal distribution for branch rates. If not supplied uses exponential.");
    final public Input<IntegerParameter> discreteRatesInput = new Input<>("discreteRates", "The rate categories associated with nodes in the species tree for sampling of individual rates among branches.");
    final public Input<RealParameter> realRatesInput = new Input<>("realRates", "The real rates associated with nodes in the species tree for sampling of individual rates among branches.", Input.Validate.XOR, discreteRatesInput);
    final public Input<Integer> nBinsInput = new Input<>("nBins", "Number of discrete branch rate bins (default is equal to the number of estimated branch rates). Only used when branch rates are catrgories.", -1);
    
    final private double MEAN_CLOCK_RATE = 1.0; // Mean clock rate. Equal to 1/lambda for exponential, or exp(M + S^2/2) for lognormal
    
    enum Mode {
        categories,
        rates
    }
    
    enum RateDistribution {
        exponential,
        lognormal
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
    
    
    private IntegerParameter categories;
    private RealParameter realRates;

    @Override
    public boolean requiresRecalculation() {
    	
    	// If lognormal S changes then branch rates require recalculation from bins
        if (rateDistribution == RateDistribution.lognormal && mode == Mode.categories) {
            final double proposedLogNormalStdev = stdevInput.get().getValue();
            if (proposedLogNormalStdev != currentLogNormalStdev) {
                binRatesNeedsUpdate = true;
            } else {
                //binRatesNeedsUpdate = false;
            }
        }

        needsUpdate = binRatesNeedsUpdate || meanRateInput.isDirty() || 
        				(mode == Mode.categories && discreteRatesInput.isDirty()) ||
        				(mode == Mode.rates && realRatesInput.isDirty());
        
        return needsUpdate;
    }

    @Override
    public void store() {
        storedLogNormalStdev = currentLogNormalStdev;
        System.arraycopy(ratesArray, 0, storedRatesArray, 0, ratesArray.length);
        
        // Store the bin rates only if using categories
        if (mode == Mode.categories) System.arraycopy(binRates, 0, storedBinRates, 0, binRates.length);
        
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
        
        
        // Restore the bin rates only if using categories
        if (mode == Mode.categories) {
        	double[] tmpBinRates = binRates;
        	binRates = storedBinRates;
        	storedBinRates = tmpBinRates;
        }
        

        super.restore();
    }

    @Override
    public void initAndValidate() {
        categories = discreteRatesInput.get();
        realRates = realRatesInput.get();
        final TreeInterface speciesTree = treeInput.get();
        final Node[] speciesNodes = speciesTree.getNodesAsArray();
        estimateRoot = estimateRootInput.get().booleanValue();
        noCache = noCacheInput.get().booleanValue();
        rootNodeNumber = speciesTree.getRoot().getNr();

        
        mode = categories == null ?  Mode.rates : Mode.categories;
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
        
        
        	// Initialise discrete branch rates
	        case categories: {
	        	
	            final int nBinsSupplied = nBinsInput.get().intValue();
	            nBins = (nBinsSupplied <= 0) ? nEstimatedRates : nBinsSupplied;
	            
	            binRates = new double[nBins];
	            storedBinRates = new double[nBins];
	        	
	            categories.setDimension(nEstimatedRates);
	            Integer[] initialCategories = new Integer[nEstimatedRates];
	            for (int i = 0; i < nEstimatedRates; i++) {
	                initialCategories[i] = Randomizer.nextInt(nBins);
	            }
	            
	            // Set initial values of rate categories
	            IntegerParameter other = new IntegerParameter(initialCategories);
	            categories.assignFromWithoutID(other);
	            categories.setLower(0);
	            categories.setUpper(nBins - 1);
	            
	            switch(rateDistribution) {
	            
		            // Initialise categorical exponential rates
		            case exponential: {
		            	
		            	// Special case: if using an exponential prior then set the bin rates now because they will never change
		            	final ExponentialDistribution exponentialDistr = new ExponentialDistributionImpl(MEAN_CLOCK_RATE);
		            	try {
	                        for (int i = 0; i < nBins; i++) {
	                        	binRates[i] = exponentialDistr.inverseCumulativeProbability((i + 0.5) / nBins);
	                        }
	                    } catch (MathException e) {
	                        throw new RuntimeException("Failed to compute exponential distribution inverse cumulative probability!");
	                    }
		            	
		            	binRatesNeedsUpdate = false;
		            	break;
		            }
		            
		            // Initialise categorical lognormal rates
		            case lognormal: {
		            	
		            	currentLogNormalStdev = stdevInput.get().getValue();
		            	storedLogNormalStdev = currentLogNormalStdev;
		            	binRatesNeedsUpdate = true;
		            	break;
		            }
	            
	            }
	            

	            
	            break;
	        }
	        
	        // Initialise real branch rates
	        case rates: {
	        	
	        	if (realRates.getDimension() != nEstimatedRates) {
	        		realRates.setDimension(nEstimatedRates);
	        		
	        		// Randomly draw rates from the appropriate distribution
	        		Double [] initialRates = new Double[nEstimatedRates];
	        		switch(rateDistribution) {
	        	        
	        			// Initialise continuous exponential rates
	             		case exponential: {
	             			for (int i = 0; i < nEstimatedRates; i++) {
	             				initialRates[i] = Randomizer.nextExponential(1 / MEAN_CLOCK_RATE);
	             			}
	             			break;
	             		}
	             		
	             		// Initialise continuous lognormal rates
	             		case lognormal: {
	             			
	             			currentLogNormalStdev = stdevInput.get().getValue();
	    	            	storedLogNormalStdev = currentLogNormalStdev;
	             			
	             			// Mean in log space
	                        final double M = Math.log(MEAN_CLOCK_RATE) - (0.5 * currentLogNormalStdev * currentLogNormalStdev);
	    	        		for (int i = 0; i < nEstimatedRates; i++) {
	    	        			initialRates[i] = Randomizer.nextLogNormal(M, currentLogNormalStdev, false);
	    	        		}
	             			
	             			break;
	             		}
	        		
	        		 }

				    RealParameter other = new RealParameter(initialRates);
				    realRates.assignFromWithoutID(other);
				    
				}
	        	
	        	realRates.setLower(0.0);
	        	binRatesNeedsUpdate = false;
	        	break;
	        }
	        
        }
        
        
      
        needsUpdate = true;
        
    }
    
    

    private void update() {
    	

        switch(mode) {
        
        	// Update discrete branch rates
	        case categories: {
	        
        		switch(rateDistribution) {
        	        
        			// Update continuous exponential rates
             		case exponential: {
             			
             			// Special case: no updating required
             			
             			break;
             		}
             		
             		// Update continuous lognormal rates
             		case lognormal: {
             			
             			if (binRatesNeedsUpdate || noCache) {
             	            // set the mean in real space to equal 1
             	            currentLogNormalStdev = stdevInput.get().getValue();
             	            
             	            // Mean in log space
             	            final double M = Math.log(MEAN_CLOCK_RATE) - (0.5 * currentLogNormalStdev * currentLogNormalStdev);
             	            final NormalDistribution normalDistr = new NormalDistributionImpl(M, currentLogNormalStdev);

             	            try {
             	                for (int i = 0; i < nBins; i++) {
             	                	
             	                	// Discrete LogNormal distributed rates
             	                    binRates[i] = Math.exp(normalDistr.inverseCumulativeProbability((i + 0.5) / nBins));
             	                }
             	            } catch (MathException e) {
             	                throw new RuntimeException("Failed to compute inverse cumulative probability!");
             	            }
             	        }
             			
             			break;
             		}
        		
        		 }
	        	break;
	        }
	        
	        // Update real branch rates
	        case rates: {
	        
        		// No updates required
	        	break;
	        }
	        
        }
    	
    	
    	

        Double estimatedMean;
        final RealParameter estimatedMeanParameter = meanRateInput.get();
        if (estimatedMeanParameter == null) {
            estimatedMean = 1.0;
        } else {
            estimatedMean = estimatedMeanParameter.getValue();
        }

        
        
        // Multiply the raw rate by the clock rate
        switch(mode){
        
	        case categories: {
	        	
	        	final Integer[] branchRatePointers = discreteRatesInput.get().getValues();
	            for (int i = 0; i < nEstimatedRates; i++) {
	                int b = branchRatePointers[i];
	                ratesArray[i] = estimatedMean * binRates[b];
	            }
	        	break;
	        }
	        
	        case rates: {
	        	
	            for (int i = 0; i < nEstimatedRates; i++) {
	                ratesArray[i] = estimatedMean * realRates.getValue(i);
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
    


    
    
}




