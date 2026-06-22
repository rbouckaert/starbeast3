package starbeast3.operators;

import java.util.*;


import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import starbeast3.evolution.speciation.GeneTreeForSpeciesTreeDistribution;
import starbeast3.tree.SpeciesTree;


public abstract class CoordinatedOperator extends Operator {
	
	public final Input<List<GeneTreeForSpeciesTreeDistribution>> geneTreesInput = new Input<>("gene", "list of gene trees that constrain species tree movement", new ArrayList<>());
    public Input<SpeciesTree> speciesTreeInput = new Input<>("speciesTree", "The species tree state node.", Validate.REQUIRED);

    protected int nGeneTrees;
    protected List<Tree> geneTrees;

    @Override
    public void initAndValidate() {
    	this.geneTrees = new ArrayList<Tree>();
		for (GeneTreeForSpeciesTreeDistribution t : this.geneTreesInput.get() ) {
			this.geneTrees.add((Tree) t.getGeneTree());
		}
        nGeneTrees = geneTrees.size();
    }

    protected Set<String> findDescendants(Node speciesTreeNode, int speciesTreeNodeNumber) {
        final Map<Integer, Set<String>> numberTipMap = speciesTreeInput.get().getNumberTipMap();
        final Set<String> descendantNames = new HashSet<>();

        if (speciesTreeNode.isLeaf()) {
            descendantNames.addAll(numberTipMap.get(speciesTreeNodeNumber));
        } else {
            final Node leftChild = speciesTreeNode.getLeft();
            final Node rightChild = speciesTreeNode.getRight();
            final int leftChildNumber = leftChild.getNr();
            final int rightChildNumber = rightChild.getNr();

            descendantNames.addAll(findDescendants(leftChild, leftChildNumber));
            descendantNames.addAll(findDescendants(rightChild, rightChildNumber));
        }

        return descendantNames;
    }
    
    
 // store and return a single double value
 // value if never set() is positive infinity
 // if set() is called multiple times, the smallest value will be stored
 final class MinimumDouble {
     private double storedDouble;

     public MinimumDouble() {
         storedDouble = Double.POSITIVE_INFINITY;
     }

     public void set(double inputDouble) {
         if (inputDouble < storedDouble) {
             storedDouble = inputDouble;
         }
     }

     public double get() {
         return storedDouble;
     }
 }

}
