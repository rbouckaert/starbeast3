package starbeast3.operators;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.operators.KernelDistribution;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;

@Description("Scale operator that scales random epoch in a tree")
public class EpochOperator extends Operator {
    final public Input<List<Tree>> treeInput = new Input<>("tree", "beast.tree on which this operation is performed", new ArrayList<>());
    final public Input<List<RealParameter>> upInput = new Input<>("up", "list of parameters to increase when tree increases (optional)", new ArrayList<>());
    final public Input<List<RealParameter>> downInput = new Input<>("down", "list of parameters to decrease when tree increases (optional)", new ArrayList<>());
    final public Input<KernelDistribution> kernelDistributionInput = new Input<>("kernelDistribution", "provides sample distribution for proposals", 
    		KernelDistribution.newDefaultKernelDistribution());
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor -- positive number that determines size of the jump: higher means bigger jumps.", 0.1);

    KernelDistribution kernelDistribution;
    double scaleFactor;
    final double updownFactor = 1.0/3.0;
    
    @Override
	public void initAndValidate() {
    	kernelDistribution = kernelDistributionInput.get();
    	scaleFactor = scaleFactorInput.get();
    	if (treeInput.get().isEmpty()) throw new IllegalArgumentException("Please provide at least 1 'tree'");
    	
	}	
	
	public EpochOperator(){}
	
	public EpochOperator(Tree tree, double weight) {
		initByName("tree", tree, "weight", weight);
	}
	
	
	

    @Override
    public double proposal() {
    	
    	List<Tree> trees = treeInput.get();
    	List<RealParameter> ups = upInput.get();
    	List<RealParameter> downs = downInput.get();

    	// Get upper and lower
    	double upperFwd = 0;
    	double lowerFwd = Double.NEGATIVE_INFINITY;
    	for (Tree tree : trees) {
	    	double upper_t = tree.getRoot().getHeight();
			double lower_t = 0;
			Node [] nodes = tree.getNodesAsArray();
			for (int i = 0; i < tree.getLeafNodeCount(); i++) {
				lower_t = Math.max(nodes[i].getHeight(), lower_t);
			}
			
			upperFwd = Math.max(upperFwd, upper_t);
			lowerFwd = Math.max(lowerFwd, lower_t);
    	}

		
	
    	// Lower bound of scaling
		double l = lowerFwd + Randomizer.nextDouble() * (upperFwd - lowerFwd);
		double pfromFwd = -2*Math.log(upperFwd - lowerFwd);
		
		
		// Upper bound of scaling
		//double u = lowerFwd + scale * (l - lowerFwd);
		double u = lowerFwd + Randomizer.nextDouble() * (upperFwd - lowerFwd);
		
		
		// Sample a scale factor
		double scale = kernelDistribution.getScaler(1, scaleFactor);
		
		
		// Ensure 'l' is lower than 'u'
		if (u < l) {
			double tmp = l;
			l = u;
			u = tmp;
		}
		//double delta = u-l;
		double oldRange = upperFwd - lowerFwd;
		double newRange = (l-lowerFwd) + (u-l)*scale + (upperFwd-u);
		double delta = newRange - oldRange;
		
		
		
		
		int scaled=0, totalNodes=0, goingUp=0, goingDown=0;
		double oldLength = 0, newLength = 0;
		
		
		// Scale the trees
		for (Tree tree : trees) {
			
			Node [] nodes = tree.getNodesAsArray();

			for (int i = tree.getLeafNodeCount(); i < nodes.length; i++) {
				
				oldLength += 2*nodes[i].getHeight() - nodes[i].getChild(0).getHeight() - nodes[i].getChild(1).getHeight();
				
				Node node = nodes[i];
				double h = node.getHeight();
				totalNodes ++;
				
				// If above u, then sum by constant amount
				if (h > u) {
					
					h = h + delta;
					node.setHeight(h);
					
				}
				
				
				// If between l and u, scale it
				else if (h > l && h < u) {
					h = l + scale*(h - l);
					node.setHeight(h);
					scaled++;
				}
				
				
				newLength += 2*nodes[i].getHeight() - nodes[i].getChild(0).getHeight() - nodes[i].getChild(1).getHeight();
				
				// If under l then leave it
				
				/*
				if (h > lowerFwd && h < l) {
					h = lowerFwd + scale * (h-lowerFwd);
					node.setHeight(h);
					scaled++;
				} else if (h > l) {				
					h += delta;
					node.setHeight(h);
				}
				*/
			}
	
			for (Node node0 : nodes) {
				if (node0.getLength() < 0) {
					return Double.NEGATIVE_INFINITY;
				}
			}
		}
		
		
		// Hastings ratio for reverse
    	double upperBck = 0;
    	double lowerBck = Double.NEGATIVE_INFINITY;
    	for (Tree tree : trees) {
	    	double upper_t = tree.getRoot().getHeight();
			double lower_t = 0;
			Node [] nodes = tree.getNodesAsArray();
			for (int i = 0; i < tree.getLeafNodeCount(); i++) {
				lower_t = Math.max(nodes[i].getHeight(), lower_t);
			}
			
			upperBck = Math.max(upperBck, upper_t);
			lowerBck = Math.max(lowerBck, lower_t);
    	}
    	double pfromBck = -2*Math.log(upperBck - lowerBck);
    	

		
		double scale2 = upperBck / upperFwd;
		//double scale2 = newLength / oldLength;
		double proportionScaled = 1.0 * scaled / totalNodes;
		if (scale2 > 1) scale2 = 1 + proportionScaled*(scale2-1);
		else scale2 = 1 - proportionScaled*(1-scale2);
		
		//Log.warning("scale " + scale + " scale2 " + scale2);
    	
    	
		// Scale parameters up
		for (RealParameter up : ups) {
			for (int i = 0; i < up.getDimension(); i++) {
				 up.setValue(i, up.getValue(i) * scale2);
				 goingUp++;
			}
		}
		
		
		// Scale parameters down
		for (RealParameter down : downs) {
			for (int i = 0; i < down.getDimension(); i++) {
				 down.setValue(i, down.getValue(i) / scale2);
				 goingDown++;
			}
		}
		
		
		return pfromBck-pfromFwd + (scaled)*Math.log(scale) + (goingUp-goingDown)*Math.log(scale2);
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
    	return 0.4;
    }
    
    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
    	scaleFactor = value; // Math.max(Math.min(value, upper), lower);
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




