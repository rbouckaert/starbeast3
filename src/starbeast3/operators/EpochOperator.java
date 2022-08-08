package starbeast3.operators;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import starbeast3.evolution.speciation.GeneTreeForSpeciesTreeDistribution;

@Description("Scale operator that scales random epoch in a tree")
public class EpochOperator extends Operator {
    final public Input<List<GeneTreeForSpeciesTreeDistribution>> genesInput = new Input<>("gene", "list of gene trees", new ArrayList<>());
    final public Input<Boolean> moveSpeciesTreeInput = new Input<>("moveSpeciesTree", "whether to move the species tree or not", false);
    
    final public Input<List<RealParameter>> upInput = new Input<>("up", "list of parameters to increase when tree increases (optional)", new ArrayList<>());
    final public Input<List<RealParameter>> downInput = new Input<>("down", "list of parameters to decrease when tree increases (optional)", new ArrayList<>());
    final public Input<KernelDistribution> kernelDistributionInput = new Input<>("kernelDistribution", "provides sample distribution for proposals", 
    		KernelDistribution.newDefaultKernelDistribution());
    
    
    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);
    final public Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor -- positive number that determines size of the jump: higher means bigger jumps.", 0.1);

    
    Tree speciesTree;
    KernelDistribution kernelDistribution;
    double scaleFactor;
    final double updownFactor = 1.0/3.0;
    
    @Override
	public void initAndValidate() {
    	kernelDistribution = kernelDistributionInput.get();
    	//kernelDistribution = new KernelDistribution.Bactrian(KernelDistribution.Bactrian.mode.uniform); // Uniform only
    	scaleFactor = scaleFactorInput.get();
    	
    	//if (genesInput.get().isEmpty()) throw new IllegalArgumentException("Please provide at least 1 gene tree distribution");
    	
    	// Get species tree
    	speciesTree = null;
    	for (GeneTreeForSpeciesTreeDistribution gene : genesInput.get()) {
    		
    		Tree sp = gene.speciesTreeInput.get();
    		if (speciesTree == null) {
    			speciesTree = sp;
    		}else {
    			if (speciesTree != sp) throw new IllegalArgumentException("Please ensure that all genes share the same species tree");
    		}
    				
    		
    	}
    	
    	
	}	
	
	public EpochOperator(){}
	
	public EpochOperator(Tree tree, double weight) {
		initByName("tree", tree, "weight", weight);
	}
	
	
	

    @Override
    public double proposal() {
    	
    	
    	// All trees
    	List<Tree> trees = new ArrayList<>();
    	if (speciesTree != null && moveSpeciesTreeInput.get()) trees.add(speciesTree);
    	for (GeneTreeForSpeciesTreeDistribution gene : genesInput.get()) {
    		trees.add((Tree) gene.getGeneTree());
    	}
    	if (trees.isEmpty()) return Double.NEGATIVE_INFINITY;
    	
    	List<RealParameter> ups = upInput.get();
    	List<RealParameter> downs = downInput.get();

    	
    	
    	// Get upper and lower
    	double upperFwd = Double.POSITIVE_INFINITY;
    	double lowerFwd = Double.NEGATIVE_INFINITY;
    	for (Tree tree : trees) {
	    	double upper_t = tree.getRoot().getHeight();
			double lower_t = 0;
			Node [] nodes = tree.getNodesAsArray();
			for (int i = 0; i < tree.getLeafNodeCount(); i++) {
				lower_t = Math.max(nodes[i].getHeight(), lower_t);
			}
			
			upperFwd = Math.min(upperFwd, upper_t);
			lowerFwd = Math.max(lowerFwd, lower_t);
    	}
    	
    	
    	
    	// Get lower and upper ranks
    	int lowerRank = Randomizer.nextInt(speciesTree.getNodeCount()); // From leaves to root
    	int upperRank = speciesTree.getLeafNodeCount() + Randomizer.nextInt(speciesTree.getInternalNodeCount()); // From internal nodes to above root
    	if (upperRank < lowerRank) {
    		int tmp = lowerRank;
    		lowerRank = upperRank;
    		upperRank = tmp;
    	}
    	if (upperRank == lowerRank) return Double.NEGATIVE_INFINITY;
    	
    	double l = speciesTree.getNode(lowerRank).getHeight();
    	double u = upperRank == speciesTree.getNodeCount() ? upperFwd : speciesTree.getNode(upperRank).getHeight();
    	
    	
    	
    	
    	// Mid point
    	//double midFwd = (upperFwd - lowerFwd) / 2;
	
    	// Lower bound of scaling
		//double l = lowerFwd + Randomizer.nextDouble() * (upperFwd - lowerFwd);
		//double pfromFwd = -2*Math.log(upperFwd - lowerFwd);
		
		
		
		// Upper bound of scaling
		//double u = lowerFwd + Randomizer.nextDouble() * (upperFwd - lowerFwd);
		
		
		
		// Sample a scale factor. Ensure that the window (log space) does not force genes below their species, and is symmetric
		double lower_s_fwd = this.getLowerScaleLimit(l, u);
		double upper_s_fwd;
		if (lower_s_fwd == 0) upper_s_fwd = Double.POSITIVE_INFINITY;
		else upper_s_fwd = 1.0 / lower_s_fwd;
		
		
		double scale = kernelDistribution.getScaler(1, scaleFactor);
		if (scale <= lower_s_fwd || scale >= upper_s_fwd) {
			return Double.NEGATIVE_INFINITY;
		}
		
		//double pScaleFwd = 
		
		//double pScaleFwd
		
		
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
				//if (node.isRoot()) continue; //tmp
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
    	double upperBck = Double.POSITIVE_INFINITY;
    	double lowerBck = Double.NEGATIVE_INFINITY;
    	for (Tree tree : trees) {
	    	double upper_t = tree.getRoot().getHeight();
			double lower_t = 0;
			Node [] nodes = tree.getNodesAsArray();
			for (int i = 0; i < tree.getLeafNodeCount(); i++) {
				lower_t = Math.max(nodes[i].getHeight(), lower_t);
			}
			
			upperBck = Math.min(upperBck, upper_t);
			lowerBck = Math.max(lowerBck, lower_t);
    	}
    	//double midBck = (upperBck - lowerBck) / 2;
    	//double pfromBck = -2*Math.log(midBck);
    	
    	//double pfromBck = -2*Math.log(upperBck - lowerBck);
    	

		
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
		
		
		// Reverse scale sampling
		double lower_s_bck = this.getLowerScaleLimit(l, u);
		if (lower_s_bck != 0) {
			//double upper_s_bck = 1.0/lower_s_bck;
			
			
			// Log.warning("forward " + lower_s_fwd + " / reverse " + lower_s_bck);
			
		}
		
		//Log.warning("upperFwd " + upperFwd + " lowerFwd " + lowerFwd);
		//Log.warning("upperBck " + upperBck + " lowerBck " + lowerBck);
		
		//return pfromBck-pfromFwd + (scaled)*Math.log(scale) + (goingUp-goingDown)*Math.log(scale2);
		return (scaled)*Math.log(scale) + (goingUp-goingDown)*Math.log(scale2);
		
		
		
    }
    
    
    /**
     * Get lower scale limit 
     * @param l
     * @param u
     * @return
     */
    private double getLowerScaleLimit(double l, double u) {
    	
    	double lower_s = 0;
    	
    	if (false && !moveSpeciesTreeInput.get()) {
			for (Node speciesNode : speciesTree.getNodesAsArray()) {
				
				
				// If the species node is outside of the window then there are no issues
				if (speciesNode.getHeight() <= l || speciesNode.getHeight() >= u) continue;
				
				// Find the maximum amount each node can be scaled down by without breaking the species tree
				for (GeneTreeForSpeciesTreeDistribution gene : genesInput.get()) {
					for (Node geneNode : gene.mapSpeciesNodeToGeneTreeNodes(speciesNode)) {
						if (geneNode.getHeight() > l && geneNode.getHeight() < u) {
							
							// This is the maximum amount the gene can be scaled (down) by without breaking the species tree
							double diff = (speciesNode.getHeight()-l) / (geneNode.getHeight()-l);
							lower_s = Math.max(diff, lower_s);
							
						}
					}
				}
				
			}
		}
    	
    	return lower_s;
    	
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
    
    
    
    @Override
    public List<StateNode> listStateNodes() {
    	List<StateNode> nodes = super.listStateNodes();
    	if (moveSpeciesTreeInput.get()) nodes.add(speciesTree);
    	for (GeneTreeForSpeciesTreeDistribution gene : genesInput.get()) {
    		nodes.add((StateNode)gene.getGeneTree());
    	}
    	return nodes;
    	
    }
    
    
}




