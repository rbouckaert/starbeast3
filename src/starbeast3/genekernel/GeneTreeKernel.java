package starbeast3.genekernel;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Iterator;
import java.util.List;
import java.util.ListIterator;

import org.w3c.dom.Node;

import beast.base.core.BEASTObject;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.core.Input.Validate;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.tree.Tree;

public class GeneTreeKernel extends StateNode {
	
	
	List<GTKGeneTree> trees;
	List<GTKGeneTree> storedTrees;
	
	
	public GeneTreeKernel() {
		
	}
	
	public GeneTreeKernel(List<GTKGeneTree> trees) {
		this.trees = trees;
		this.storedTrees = new ArrayList<GTKGeneTree>();
	}
	
	
	
	@Override
	public void initAndValidate() {
		this.trees = new ArrayList<GTKGeneTree>();
		this.storedTrees = new ArrayList<GTKGeneTree>();
	}
	
	/**
	 * @return The list of trees contained in the kernel
	 */
	public List<GTKGeneTree> getTrees() {
		return this.trees;
	}
	
	
	/**
	 * @return The list of trees contained in the kernel
	 */
	public List<Tree> getCoercedTrees() {
		List<Tree> coercedTrees = new ArrayList<Tree>();
		for (GTKGeneTree tree : this.trees) coercedTrees.add((Tree) tree);
		return coercedTrees;
	}
	
	
	/**
	 * Deletes the specified tree
	 * Please ensure that all elements pointing to a tree after this index have their indices decremented
	 * @param index
	 */
	public void deleteTree(int index) {
		this.startEditing();
		this.trees.remove(index);
	}

	
	
	/**
	 * @return The the tree in the kernel indexed by 'index'
	 */
	public Tree getTree(int index) {
		return this.trees.get(index);
	}
	
	
	/**
	 * @param newTrees - the list of trees for the kernel
	 */
	public void setTrees(List<GTKGeneTree> newTrees) {
		this.startEditing();
		this.trees = newTrees;
	}
	
	/**
	 * @param tree - a new tree to add to the kernel
	 */
	public void addTree(GTKGeneTree tree) {
		this.startEditing();
		this.trees.add(tree);
	}
	

	@Override
	public void init(PrintStream out) {
		
		
	}

	@Override
	public void close(PrintStream out) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public int getDimension() {
		return trees.size();
	}

	@Override
	public double getArrayValue() {
        return trees.get(0).getArrayValue();
    }

    @Override
	public double getArrayValue(int value) {
        return trees.get(value).getArrayValue();
    }

	@Override
	public StateNode copy() {
		List<GTKGeneTree> copiedList = new ArrayList<GTKGeneTree>();
		for (GTKGeneTree tree : this.trees) {
			GTKGeneTree copy = new GTKGeneTree(tree.getRoot());
			copiedList.add(copy);
		}
		return new GeneTreeKernel(copiedList);
	}

	@Override
	public void assignTo(StateNode other) {
		GeneTreeKernel kernel = (GeneTreeKernel) other;
		kernel.setTrees(this.getTrees());
	}

	@Override
	public void assignFrom(StateNode other) {
		GeneTreeKernel kernel = (GeneTreeKernel) other;
		this.trees = kernel.getTrees();
	}

	@Override
	public void assignFromFragile(StateNode other) {
		// TODO Optimise
		this.assignFrom(other);
	}

	@Override
	public void fromXML(Node node) {
		// TODO Auto-generated method stub
	}

	@Override
	public int scale(double scale) {
		int dof = 0;
		for (Tree tree : this.trees) {
			dof += tree.scale(scale);
		}
		return dof;
	}
	
	
	public void startEditing() {
        this.startEditing(null);
    }
	
	@Override
    public void startEditing(final Operator operator) {
        super.startEditing(operator);
    }


	@Override
	protected void store() {
		this.storedTrees.clear();
		for (GTKGeneTree tree : this.trees) {
			this.storedTrees.add(tree);
			tree.store();
		}
		
		
	}

	@Override
	public void restore() {

		// Point 'temp' to the stored trees
		List<GTKGeneTree> temp = new ArrayList<GTKGeneTree>();
		for (GTKGeneTree tree : this.storedTrees) {
			temp.add(tree);
			tree.restore();
		}
		
		// Point 'storedTrees' to the changed trees
		this.storedTrees = this.trees;
		
		// Point 'trees' to temp
		this.trees = temp;

		hasStartedEditing = false;
		
		
	}

	
	
	@Override
	public void setEverythingDirty(boolean isDirty) {
		setSomethingIsDirty(isDirty);
		for (Tree tree : this.trees) {
			tree.setEverythingDirty(isDirty);
		}
	}

	@Override
	public void log(long sample, PrintStream out) {
		// TODO Auto-generated method stub
		
	}


	
	
	

}
