package starbeast3.operators;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.operator.kernel.BactrianNodeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.util.InputUtil;
import beast.base.util.Randomizer;
import starbeast3.evolution.speciation.GeneTreeForSpeciesTreeDistribution;


@Description("Moves a node in a gene tree so that it does not go below the species node")
public class GeneNodeOperator extends BactrianNodeOperator {
	
	
	final public Input<GeneTreeForSpeciesTreeDistribution> geneTreeDistributionsInput = new Input<>("gene", "gene tree for species tree distribution", Input.Validate.REQUIRED);
	   
	
    @Override
    public double proposal() {
    	
    	GeneTreeForSpeciesTreeDistribution gene = geneTreeDistributionsInput.get();
        Tree tree = (Tree) InputUtil.get(treeInput, this);

        // randomly select internal node
        int nodeCount = tree.getNodeCount();
        
        // Abort if no non-root internal nodes
        if (tree.getInternalNodeCount()==1)
            return Double.NEGATIVE_INFINITY;
        
        Node node;
        do {
            int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(nodeNr);
        } while (node.isRoot() || node.isLeaf());
        
        
        // Get matching species tree node
        Node speciesNode = gene.mapGeneNodeToSpeciesNode(node.getNr());
        double speciesLower = speciesNode.getHeight();
        
        
        double upper = node.getParent().getHeight();
        double lower = Math.max(Math.max(node.getLeft().getHeight(), node.getRight().getHeight()), speciesLower);
        
        double scale = kernelDistribution.getScaler(0, Double.NaN, scaleFactor);

        // transform value
        double value = node.getHeight();
        double y = (upper - value) / (value - lower);
        y *= scale;
        double newValue = (upper + lower * y) / (y + 1.0);
        
        if (newValue < lower || newValue > upper) {
        	return Double.NEGATIVE_INFINITY;
        	//throw new RuntimeException("programmer error: new value proposed outside range");
        }
        
        node.setHeight(newValue);

        double logHR = Math.log(scale) + 2.0 * Math.log((newValue - lower)/(value - lower));
        return logHR;
    }


}
