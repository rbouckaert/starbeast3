package starbeast3.tree;

import beast.evolution.tree.Node;


/** node in a binary tree **/
public class BinaryNode extends Node {
	BinaryTree tree;
	
	public BinaryNode(int labelNr, BinaryTree tree) {
		this.labelNr = labelNr;
		this.tree = tree;
	}
	
	@Override
	public String getID() {
		if (labelNr < tree.taxa.length) {
			return tree.taxa[labelNr];
		}
		return null;
	}
	
	@Override
	public double getHeight() {
		return tree.height[labelNr];
	}
	
    @Override
    public void setHeight(final double height) {
        tree.startEditing();
        tree.height[labelNr] = height;
    }

    @Override
	public Node getLeft() {		
		return tree.nodes[tree.left[labelNr]];
	}

	@Override
	public Node getRight() {		
		return tree.nodes[tree.right[labelNr]];
	}

	@Override
	public Node getParent() {		
		return tree.nodes[tree.parent[labelNr]];
	}

	@Override
	public int getChildCount() {
		if (isLeaf()) {
			return 0;
		}
		return 2;
	}
	
	@Override
	public boolean isRoot() {
		return labelNr == tree.nodes.length - 1;
	}
	
	@Override
	public boolean isLeaf() {
		return labelNr < tree.taxa.length;
	}

    public String toSortedNewick(final int[] maxNodeInClade) {
        return toSortedNewick(maxNodeInClade, false);
    }

    @Override
    public void setMetaData(String pattern, Object value) {
    	// TODO Auto-generated method stub
    	super.setMetaData(pattern, value);
    }
    
    @Override
    public Object getMetaData(String pattern) {
    	// TODO Auto-generated method stub
    	return super.getMetaData(pattern);
    }
}
