package starbeast3.operators;

import java.text.DecimalFormat;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.parameter.BooleanParameter;
import beast.base.inference.parameter.RealParameter;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.evolution.operator.ScaleOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import starbeast3.genekernel.GTKOperator;



/**
 * Adapted from BEASTLabs's BactrianScaleOperator and beast2's ScaleOperator
 * Operator is compatible with gene tree kernels
 * 
 * Jordan Douglas, March 2020
 */
@Description("Scale operator that finds scale factor according to a Bactrian distribution (Yang & Rodriguez, 2013), "
		+ "which is a mixture of two Gaussians: p(x) = 1/2*N(x;-m,1-m^2) + 1/2*N(x;+m,1-m^2) and more efficient than RealRandomWalkOperator")
public class BactrianTreeScaleOperator extends GTKOperator {
	
	

    final public Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor: larger means more bold proposals", 1.0);
    final public Input<Boolean> rootOnlyInput = new Input<>("rootOnly", "scale root of a tree only, ignored if tree is not specified (default false)", false);
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);

    
    
    private double m_fScaleFactor;


	@Override
	public void initAndValidate() {
        m_fScaleFactor = scaleFactorInput.get();
        super.initAndValidate();
	}
    
	protected double getScaler(int i, double value) {
    	return kernelDistribution.getScaler(i, value, getCoercableParameterValue());
	}
    

    @Override
    public double proposal() {
        try {

            final Tree tree = this.sampleTree(this);
            
            if (rootOnlyInput.get()) {
                final Node root = tree.getRoot();                    
                final double scale = getScaler(root.getNr(), root.getHeight());
                final double newHeight = root.getHeight() * scale;

                if (newHeight < Math.max(root.getLeft().getHeight(), root.getRight().getHeight())) {
                    return Double.NEGATIVE_INFINITY;
                }
                root.setHeight(newHeight);
                return Math.log(scale);
            } else {
                // scale the beast.tree
                final double scale = getScaler();
                final int internalNodes = tree.scale(scale);
                return Math.log(scale) * internalNodes;
            }

        } catch (Exception e) {
            // whatever went wrong, we want to abort this operation...
            return Double.NEGATIVE_INFINITY;
        }
    }

    @Override
    public void optimize(double logAlpha) {
        // must be overridden by operator implementation to have an effect
    	if (optimiseInput.get()) {
	        double delta = calcDelta(logAlpha);
	        double scaleFactor = getCoercableParameterValue();
	        delta += Math.log(scaleFactor);
	        scaleFactor = Math.exp(delta);
	        setCoercableParameterValue(scaleFactor);
    	}
    }
    
    @Override
    public double getTargetAcceptanceProbability() {
    	return 0.3;
    }
    
    
    protected double getScaler() {
        return (m_fScaleFactor + (Randomizer.nextDouble() * ((1.0 / m_fScaleFactor) - m_fScaleFactor)));
    }
    
    protected boolean outsideBounds(final double value, final RealParameter param) {
        final Double l = param.getLower();
        final Double h = param.getUpper();

        return (value < l || value > h);
        //return (l != null && value < l || h != null && value > h);
    }


    @Override
    public double getCoercableParameterValue() {
        return m_fScaleFactor;
    }
    
    
    @Override
    public String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = getCoercableParameterValue() * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10 || prob > 0.40) {
            return "Try setting scale factor to about " + formatter.format(newWindowSize);
        } else return "";
    }

}
