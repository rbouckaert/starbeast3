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
import beast.math.distributions.ParametricDistribution;
import beast.math.distributions.PiecewiseLinearDistribution;
import beast.util.Randomizer;
import beast.evolution.branchratemodel.BranchRateModel;

public class UCRelaxedClockModelSB3 extends BranchRateModel.Base implements BranchRateModelSB3 {
    
    final public Input<Boolean> estimateRootInput = new Input<>("estimateRoot", "Estimate rate of the root branch.", false);
    final public Input<Boolean> noCacheInput = new Input<>("noCache", "Always recalculate branch rates.", false);
    final public Input<RealParameter> stdevInput = new Input<>("stdev", "Standard deviation of the log-normal distribution for branch rates. If not supplied uses exponential.");
    final public Input<IntegerParameter> discreteRatesInput = new Input<>("discreteRates", "The rate categories associated with nodes in the species tree for sampling of individual rates among branches.");
    final public Input<RealParameter> realRatesInput = new Input<>("realRates", "The real rates associated with nodes in the species tree for sampling of individual rates among branches.");
    final public Input<RealParameter> quantilesInput = new Input<>("rateQuantiles", "The rate quantiles for sampling of individual rates among branches.");
    final public Input<Integer> nBinsInput = new Input<>("nBins", "Number of discrete branch rate bins (default is equal to the number of estimated branch rates). Only used when branch rates are catrgories.", -1);
    final public Input<ParametricDistribution> rateDistInput = new Input<>("distr", "the distribution governing the rates among branches, for quantiles only. Must have mean of 1. The clock.rate parameter can be used to change the mean rate.");
    
    
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
    
    
    private IntegerParameter categories;
    private RealParameter realRates;
    private RealParameter quantiles;
    
    final int LATTICE_SIZE_FOR_DISCRETIZED_RATES = 100;
    
    // Piecewise linear approximation for quantiles
    ParametricDistribution quantileDistribution = null; // Distribution of rates for quantile parameterisation
    //PiecewiseLinearDistribution piecewiseLinearQuantiles = null;
    
    
    
    

    @Override
    public void initAndValidate() {
        categories = discreteRatesInput.get();
        realRates = realRatesInput.get();
        quantiles = quantilesInput.get();
        final TreeInterface speciesTree = treeInput.get();
        final Node[] speciesNodes = speciesTree.getNodesAsArray();
        estimateRoot = estimateRootInput.get().booleanValue();
        noCache = noCacheInput.get().booleanValue();
        rootNodeNumber = speciesTree.getRoot().getNr();

        
        // Get the rates and their parameterisation model
        if ((realRates != null && categories != null) ||
            	(realRates != null && quantiles != null) ||
            	(quantiles != null && categories != null)) {
            	throw new IllegalArgumentException("Only one of rateCategories, rateQuantiles or rates should be specified");
            }
        if (categories != null) mode = Mode.categories;
        else if (realRates != null) mode = Mode.rates;
        else if (quantiles != null) mode = Mode.quantiles;
        
        
        // Exponential, lognormal, or parametric rate distribution
        if (mode == Mode.quantiles) {
        	if (rateDistInput.get() == null) throw new IllegalArgumentException("distr must be provided when using quantiles parameterisation");
        	quantileDistribution = rateDistInput.get();
        	rateDistribution = RateDistribution.parametric;
        }
        else {
        	rateDistribution = stdevInput.get() == null ? RateDistribution.exponential : RateDistribution.lognormal;
        }
       

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
	            
	        	if (categories.getDimension() != nEstimatedRates) {
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
	        
	        
	        // Initialise quantiles
	        case quantiles: {
	        	
	        	
	            binRates = new double[LATTICE_SIZE_FOR_DISCRETIZED_RATES];
	            storedBinRates = new double[LATTICE_SIZE_FOR_DISCRETIZED_RATES];
	            
	            //currentLogNormalStdev = stdevInput.get().getValue();
            	//storedLogNormalStdev = currentLogNormalStdev;
            	binRatesNeedsUpdate = true;
	            
	            if (quantiles.getDimension() != nEstimatedRates) {
	        	
	            	quantiles.setDimension(nEstimatedRates);
	                Double[] initialQuantiles = new Double[nEstimatedRates];
	                for (int i = 0; i < nEstimatedRates; i++) {
	                    initialQuantiles[i] = Randomizer.nextDouble();
	                }
	                RealParameter other = new RealParameter(initialQuantiles);
	                quantiles.assignFromWithoutID(other);
	                quantiles.setLower(0.0);
	                quantiles.setUpper(1.0);
	                
	            }
	        	break;
	        }
	        
        }
        
        
      
        needsUpdate = true;
        
    }
    
    
    
    

    @Override
    public boolean requiresRecalculation() {
    	
    	// If lognormal S changes then branch rates require recalculation from bins
        if (rateDistribution == RateDistribution.lognormal && mode == Mode.categories) {
            final double proposedLogNormalStdev = stdevInput.get().getValue();
            if (proposedLogNormalStdev != currentLogNormalStdev) {
                binRatesNeedsUpdate = true;
            } 
        }
        
        if (mode == Mode.quantiles && quantileDistribution.isDirtyCalculation()) {
        	binRatesNeedsUpdate = true;
        }

        needsUpdate = binRatesNeedsUpdate || meanRateInput.isDirty() || 
        				(mode == Mode.categories && discreteRatesInput.isDirty()) ||
        				(mode == Mode.rates && realRatesInput.isDirty()) ||
        				(mode == Mode.quantiles && quantilesInput.isDirty());
        
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
        
        
        // Restore the bin rates only if using categories/quantiles
        if (mode == Mode.categories || mode == Mode.quantiles) {
        	double[] tmpBinRates = binRates;
        	binRates = storedBinRates;
        	storedBinRates = tmpBinRates;
        }
        

        super.restore();
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
	        
	        
	        // Update quantiles using a linear approximation (but only the quantiles which are being used)
	        case quantiles: {
	        	
	        	for (int i = 0; i < binRates.length; i++) {
                    binRates[i] = 0;
                }
	        	
	        	for (int nodeNumber = 0; nodeNumber < nEstimatedRates; nodeNumber++) {
					
				    // Use cached rates
				    double q = quantiles.getValue(nodeNumber);
				    double v = q * (binRates.length - 1);
				    int i = (int) v;
				    
				    // Make sure cached rates are calculated
				    if (binRates[i] == 0.0) {
				        try {
				        	if (i > 0) {
				        		binRates[i] = quantileDistribution.inverseCumulativeProbability(((double)i) / (binRates.length-1));
				        	} else {
				        		binRates[i] = quantileDistribution.inverseCumulativeProbability(0.1 / (binRates.length-1));
				        	}
				        } catch (MathException e) {
				            throw new RuntimeException("Failed to compute inverse cumulative probability!");
				        }
				    }
				    if (i < binRates.length - 1 && binRates[i + 1] == 0.0) {
				        try {
				        	if (i < binRates.length - 2) {
				        		binRates[i + 1] = quantileDistribution.inverseCumulativeProbability(((double)(i + 1)) / (binRates.length-1));
				        	} else {
				        		binRates[i + 1] = quantileDistribution.inverseCumulativeProbability((binRates.length - 1 - 0.1) / (binRates.length-1));
				        	}
				        } catch (MathException e) {
				            throw new RuntimeException("Failed to compute inverse cumulative probability!");
				        }
				    }
				    
				    // Calculate piecewise linear approximation
				    double r = binRates[i];
				    if (i < binRates.length - 1) {
				    	r += (binRates[i+1] - binRates[i]) * (v - i);
				    }
				    //System.out.println(nodeNumber + ": q = " + q + " rate = " + r);
				    ratesArray[nodeNumber] = r;
				    
				}
		        	
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
        
        
        	// Categories
	        case categories: {
	        	
	        	final Integer[] branchRatePointers = categories.getValues();
	            for (int i = 0; i < nEstimatedRates; i++) {
	                int b = branchRatePointers[i];
	                ratesArray[i] = estimatedMean * binRates[b];
	            }
	        	break;
	        }
	        
	        
	        // Real numbers
	        case rates: {
	        	
	            for (int i = 0; i < nEstimatedRates; i++) {
	                ratesArray[i] = estimatedMean * realRates.getValue(i);
	            }
	        	
	        	break;
	        }
	        
	        
	        // Quantiles
	        case quantiles: {
	        	
	        	for (int i = 0; i < nEstimatedRates; i++) {
	                ratesArray[i] = estimatedMean * ratesArray[i];
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
        
        //final double M = Math.log(MEAN_CLOCK_RATE) - (0.5 * currentLogNormalStdev * currentLogNormalStdev);
        //System.out.println(quantileDistribution + " q1 = " + quantiles.getValue(0) +  " r1 = " + ratesArray[0]);
        //for (int i = 0; i < ratesArray.length; i ++) {
        	//System.out.print(ratesArray[i] + "\t");
        //}
        //System.out.println("");

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
	
	
	public PiecewiseLinearDistribution getPiecewiseQuantileApproximation() {
		if (this.quantileDistribution instanceof PiecewiseLinearDistribution) {
			return (PiecewiseLinearDistribution)this.quantileDistribution;
		}
		else return null;
	}
    
    
    
}




