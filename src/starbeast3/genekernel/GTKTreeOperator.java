package starbeast3.genekernel;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;


/**
 * A class for Tree Operators which accept a kernel / list of gene trees, and also a species tree
 * 'GeneTreeKernelTreeOperator' is to 'GeneTreeKernelOperator' what 'TreeOperator' is to 'Operator'
 * @author Jordan Douglas
 *
 */
public abstract class GTKTreeOperator extends GTKOperator {

	
	final public Input<Tree> treeInput = new Input<>("tree", "beast.tree on which this operation is performed", Validate.OPTIONAL);
    final public Input<Boolean> markCladesInput = new Input<>("markclades", "Mark all ancestors of nodes changed by the operator as changed," +
            " up to the MRCA of all nodes changed by the operator.", false);
	
    
    Tree inputTree;
    
	
	@Override
	public void initAndValidate() {
		
		// This tree is most likely a species tree.
		// If there are any gene tree, then it will most likely be retrieved using sampleTree()
		inputTree = treeInput.get();
		
		super.initAndValidate();
	}
    
    
    
	/**
     * @param parent the parent
     * @param child  the child that you want the sister of
     * @return the other child of the given parent.
     */
    protected Node getOtherChild(final Node parent, final Node child) {
        if (parent.getLeft().getNr() == child.getNr()) {
            return parent.getRight();
        } else {
            return parent.getLeft();
        }
    }

    /**
     * replace child with another node
     *
     * @param node
     * @param child
     * @param replacement
     */
    public void replace(final Node node, final Node child, final Node replacement) {
    	node.removeChild(child);
    	node.addChild(replacement);
        node.makeDirty(Tree.IS_FILTHY);
        replacement.makeDirty(Tree.IS_FILTHY);
    }
	
}
