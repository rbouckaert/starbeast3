package starbeast3.operators;

import java.text.DecimalFormat;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.operators.KernelDistribution;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import genekernel.GTKTreeOperator;


@Description("Scales a parameter or a complete beast.tree (depending on which of the two is specified.")
public class TreeScaleOperator extends GTKTreeOperator {


    public final Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor: range from 0 to 1. Close to zero is very large jumps, close to 1.0 is very small jumps.", 0.75);

    public final Input<Boolean> scaleAllIndependentlyInput =
            new Input<>("scaleAllIndependently", "if true, all elements of a parameter (not beast.tree) are scaled with " +
                    "a different factor, otherwise a single factor is used", false);

    final public Input<Integer> degreesOfFreedomInput = new Input<>("degreesOfFreedom", "Degrees of freedom used when " +
            "scaleAllIndependently=false and scaleAll=true to override default in calculation of Hasting ratio. " +
            "Ignored when less than 1, default 0.", 0);
   
    final public Input<Boolean> rootOnlyInput = new Input<>("rootOnly", "scale root of a tree only, ignored if tree is not specified (default false)", false);

    final public Input<Double> scaleUpperLimit = new Input<>("upper", "Upper Limit of scale factor", 1.0 - 1e-8);
    final public Input<Double> scaleLowerLimit = new Input<>("lower", "Lower limit of scale factor", 1e-8);


    /**
     * shadows input *
     */
    private double scaleFactor;
    private double upper, lower;

    @Override
    public void initAndValidate() {
        scaleFactor = scaleFactorInput.get();
        upper = scaleUpperLimit.get();
        lower = scaleLowerLimit.get();
        super.initAndValidate();
    }



    protected boolean outsideBounds(final double value, final RealParameter param) {
        final Double l = param.getLower();
        final Double h = param.getUpper();

        return (value < l || value > h);
        //return (l != null && value < l || h != null && value > h);
    }

    
    /**
     * Following the kernel distribution
     * @param i
     * @param value
     * @return
     */
	protected double getScaler(int i, double value) {
    	return kernelDistribution.getScaler(i, value, getCoercableParameterValue());
	}
	
	
	/**
	 * Standard scaling
	 * @return
	 */
    protected double getScaler() {
        return (scaleFactor + (Randomizer.nextDouble() * ((1.0 / scaleFactor) - scaleFactor)));
    }

    /**
     * override this for proposals,
     *
     * @return log of Hastings Ratio, or Double.NEGATIVE_INFINITY if proposal should not be accepted *
     */
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


    /**
     * automatic parameter tuning *
     */
    @Override
    public void optimize(final double logAlpha) {
        double delta = calcDelta(logAlpha);
        delta += Math.log(1.0 / scaleFactor - 1.0);
        setCoercableParameterValue(1.0 / (Math.exp(delta) + 1.0));
    }

    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        scaleFactor = Math.max(Math.min(value, upper), lower);
    }

    @Override
    public String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double sf = Math.pow(scaleFactor, ratio);

        final DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else if (prob > 0.40) {
            return "Try setting scaleFactor to about " + formatter.format(sf);
        } else return "";
    }

} // class ScaleOperator
