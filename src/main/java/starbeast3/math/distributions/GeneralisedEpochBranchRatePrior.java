package starbeast3.math.distributions;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.domain.Real;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.inference.distribution.LogNormal;
import beast.base.spec.inference.distribution.ScalarDistribution;
import beast.base.spec.inference.distribution.TensorDistribution;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.RealScalar;
import beast.base.spec.inference.parameter.RealScalarParam;
import beast.base.spec.inference.parameter.BoolScalarParam;
import beast.base.spec.type.RealVector;
import beast.base.util.Randomizer;

import org.apache.commons.statistics.distribution.LogNormalDistribution;
import org.apache.commons.statistics.distribution.NormalDistribution;


@Description("Log normal prior distribution on branch rates, with a mean branch rate determiend by epochs")
public class GeneralisedEpochBranchRatePrior extends TensorDistribution<RealVector<NonNegativeReal>, Double>  {
	
	final public Input<Tree> treeInput = new Input<>("tree", "the tree that contains the branch rates", Input.Validate.REQUIRED); 
	final public Input<RealScalar<PositiveReal>> sigmaInput = new Input<>("sigma", "relaxed clock standard deviation (log normal).", Input.Validate.REQUIRED); 
	final public Input<RealScalar<PositiveReal>> betaInput = new Input<>("beta", "dependency on trait (log space).", Input.Validate.REQUIRED); 
	final public Input<BoolScalarParam> indicatorInput = new Input<>("indicator", "use relaxed clock if this is false", Input.Validate.REQUIRED);
	final public Input<BoolScalarParam> logTemperaturesInput = new Input<>("logTemperatures", "take the logarithm of temperaturess?", Input.Validate.REQUIRED);
	
	
	final public Input<RealVector<PositiveReal>> epochTimesInput = new Input<>("epoch", "ages of temp values", Input.Validate.REQUIRED); 
	final public Input<RealVector<Real>> temperatureInput = new Input<>("temperature", "temp values, with one element less than the number of epochs", Input.Validate.REQUIRED); 
	
	
	int nepochs;
	RealVector<PositiveReal> epochs;
	RealVector<Real> temperatures;
	
	
	@Override
    public void initAndValidate() {
        super.initAndValidate();
        
        
        epochs = epochTimesInput.get();
        temperatures = temperatureInput.get();
        nepochs = epochs.size()+1;
        
        if (epochs.size() != temperatures.size()-1) {
        	throw new IllegalArgumentException("epochs must have one more dimension than temperature");
        }
        
        if (temperatures.size() < 2) {
        	throw new IllegalArgumentException("please give at least 2 temperatures");
        }
        
        if (epochs.get(0) <= 0) {
        	throw new IllegalArgumentException("the first epoch element must be older than 0");
        }
        
        double e = epochs.get(0);
        for (int i = 1; i < nepochs-1;  i++) {
        	double e2 = epochs.get(i);
        	if (e2 <= e) {
        		throw new IllegalArgumentException("epoch ages must proceed backwards in time (e.g. 1, 3, 5, 10)");
        	}
        	
        	e=e2;
        }
        
    }
	
	@Override
	public double calculateLogP() {
	
		RealVector<NonNegativeReal> branchRates = paramInput.get();
		logP = calcLogP(branchRates.getElements());
		return logP;
	
	}



	@Override
	public void refresh() {
		// TODO Auto-generated method stub
		
	}



	@Override
	public Double getLowerBoundOfParameter() {
		return 0.0;
	}



	@Override
	public Double getUpperBoundOfParameter() {
		// TODO Auto-generated method stub
		return Double.POSITIVE_INFINITY;
	}
	
	public double getTemperature(int epochNr) {
		
		
		double epochTemp = temperatures.get(epochNr);
		if (logTemperaturesInput.get().get()) {
			epochTemp = Math.log(epochTemp);
		}

		return epochTemp;
	}



	@Override
	protected double calcLogP(Double... value) {
		return this.calcLogP(Arrays.asList(value));
	}
	
	
	private double getMeanOfBranch(int nodeNr) {
		
		if (!indicatorInput.get().get()) {
			return 0;
		}
		
		
		Node node = treeInput.get().getNode(nodeNr);
		double beta = betaInput.get().get();
		if (node.isRoot()) {
			return beta * getTemperature(this.nepochs-1);
		}
		
		
		
		double t1 = node.getHeight();
		double t2 = node.getParent().getHeight();
		
		
		// Use the final temperature only
		if (t2 > epochs.get(epochs.size()-1)) {
			return beta * getTemperature(this.nepochs-1);
		}
		
		double tempMean = 0;
		double totalOverlap = 0;
		for (int i = 0; i < this.nepochs; i ++) {
			double epochTemp = getTemperature(i);
			double epochYoung = i == 0 ? 0 : epochs.get(i-1);
			double epochOld = i == this.nepochs-1 ? Double.POSITIVE_INFINITY : epochs.get(i);
			
			
			// Branch is not in this epoch
			if (t1 > epochOld || t2 < epochYoung) {
				continue;
			}
			
			
			double overlap = 0;
			
			// Case 1: epoch is inside the branch
			if (t1 < epochYoung & t2 > epochOld) {
				overlap = epochOld-epochYoung;
			}
			
			// Case 2: epoch crosses the bottom of the branch
			else if (t1 < epochYoung & t2 < epochOld) {
				overlap = epochOld - t1;
			}
			
			// Case 3: epoch crosses the top of the branch
			else if (t1 > epochYoung & t2 > epochOld) {
				overlap = t2 - epochYoung;
			}
			
			// Case 4: branch is inside the epoch
			else {
				overlap = t2-t1;
			}
			
			tempMean += epochTemp*overlap;
			totalOverlap += overlap;
			
			
			//System.out.println("\tbranch " + nodeNr + " has time " + t1 + " - " + t2 + " and temp " + epochTemp + " overlap " + overlap);
			
			
		}
		
		tempMean = tempMean / totalOverlap;
		
		
		//System.out.println("branch " + nodeNr + " has time " + t1 + " - " + t2 + " and mean temp " + tempMean);
		
		
		return beta * tempMean;
		
	}
	
    private double calcLogP(List<Double> branchRates) {
    	
    	// Check sigma is positive
        double s = sigmaInput.get().get();
        if (s <= 0) {
    		return Double.NEGATIVE_INFINITY;
    	}
        
        
        Tree tree = (Tree) treeInput.get();
		int dimension = tree.getNodeCount()-1;
        
		
		double logp = 0;
        for (int nodeNr = 0; nodeNr < dimension; nodeNr ++) {
        	
        	double branchMean = getMeanOfBranch(nodeNr);
        	double rate = branchRates.get(nodeNr);
        	LogNormalDistribution dist = LogNormalDistribution.of(branchMean, s);
        	if (rate < 0) {
        		logp = Double.NEGATIVE_INFINITY;
        		return logp;
        	}

        	
        	logp += dist.logDensity(rate);
        	
        	//System.out.println(rate + " " + branchMean + " " + s);
        	
        }
    	
        return logp;
        
    }
    
    
    @Override
	public List<Double> sample() {
    	
    	return null;
    }

//    
//    public ScalarDistribution<?,?> getDist(){
//    	
//    	final double mean = 1;
//    	double s = sigmaInput.get().get();;
//    	
//    	
//    	RealScalar<Real> mParam = new RealScalarParam<>(mean, Real.INSTANCE);
//    	RealScalar<PositiveReal> sParam = new RealScalarParam<>(s, PositiveReal.INSTANCE);
//    	
//    	
//    	LogNormal dist = new LogNormal();
//    	dist.initByName("M", mParam, "S", sParam, "meanInRealSpace", true);
//    	
//    	
//    	return dist;
//    	
//    }


	


}
