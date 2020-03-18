package starbeast3.operators;

import java.util.*;


import beast.core.Input;
import beast.core.Input.Validate;
import beast.evolution.tree.Node;
import genekernel.GTKOperator;
import starbeast3.SpeciesTree;


public abstract class CoordinatedOperator extends GTKOperator {
    public Input<SpeciesTree> speciesTreeInput = new Input<>("speciesTree", "The species tree state node.", Validate.REQUIRED);

    protected int nGeneTrees;

    @Override
    public void initAndValidate() {
    	geneTrees = this.getTrees(this);
        nGeneTrees = this.getGeneTreeCount(this);
        super.initAndValidate();
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
