package starbeast3.tree;

import java.util.ArrayList;
import java.util.List;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;


/** node in a binary tree **/
public class BinaryNode extends Node {
	private BinaryTree tree;
    private List<Node> children;
	
	public BinaryNode(int labelNr, BinaryTree tree) {
		this.labelNr = labelNr;
		this.tree = tree;
	}
	
	@Override
	public String getID() {
		if (labelNr < tree.getTaxaNames().length) {
			return tree.getTaxaNames()[labelNr];
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
    	int k = tree.left[labelNr];
    	if (k >= 0) {
    		return tree.nodes[k];
    	}
    	return null;
	}

	@Override
	public Node getRight() {		
    	int k = tree.right[labelNr];
    	if (k >= 0) {
    		return tree.nodes[k];
    	}
    	return null;
	}

	@Override
	public Node getParent() {		
    	if (tree.parent[labelNr] >= 0) {
    		return tree.nodes[tree.parent[labelNr]];
    	}
    	return null;
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
		return getParent() == null;
	}
	
	@Override
	public boolean isLeaf() {
		return labelNr < tree.getLeafNodeCount();
	}

    public String toSortedNewick(final int[] maxNodeInClade) {
        return toSortedNewick(maxNodeInClade, false);
    }

    @Override
    public Node getChild(final int childIndex) {
    	if (childIndex == 0) {
    		return getLeft();
    	}
        return getRight();
    }

    @Override
    public List<Node> getChildren() {
    	if (children == null) {
    		children = new ArrayList<Node>(); 
    		children.add(getLeft());
    		children.add(getRight());
    		return children;
    	}
    	children.set(0, getLeft());
    	children.set(1, getRight());
    	return children;
    }

    @Override
    public int getNodeCount() {
		if (isLeaf()) {
			return 1;
		}
        int nodes = 1;
        nodes += getLeft().getNodeCount();
        nodes += getRight().getNodeCount();
        return nodes;
    }
    
    @Override
    public void removeChild(Node child) {
        startEditing();
        if (tree.left[labelNr] == child.getNr()) {
        	tree.left[labelNr] = -1;
        } else {
        	tree.right[labelNr] = -1;
        }
    }
    
    @Override
    public void removeAllChildren(final boolean inOperator) {
    	throw new RuntimeException("Not implemented yet");
    }

    @Override
    public void addChild(Node child) {
        startEditing();
        tree.parent[child.getNr()] = labelNr;
        if (tree.left[labelNr] < 0) {
    		tree.left[labelNr] = child.getNr();
    	} else {
    		tree.right[labelNr] = child.getNr();
    	}
    }
    
    @Override
    public void setParent(Node parent) {
    	if (parent == null) {
    		tree.parent[labelNr] = -1;
    	} else {
    		tree.parent[labelNr] = parent.getNr();
    	}
    }
    
    @Override
    public int scale(final double scale) {
        startEditing();

        int dof = 0;

        makeDirty(Tree.IS_DIRTY);
        if (!isLeaf() && !isFake()) {
            setHeight(getHeight() * scale);

            if (isRoot() || getParent().getHeight() != getHeight())
                dof += 1;
        }
        if (!isLeaf()) {
            dof += getLeft().scale(scale);
            if (getRight() != null) {
                dof += getRight().scale(scale);
            }
            if (height < getLeft().getHeight() || height < getRight().getHeight()) {
                throw new IllegalArgumentException("Scale gives negative branch length");
            }
        }

        return dof;
    }
}
