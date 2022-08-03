package starbeast3.genekernel;

import java.util.List;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;


/**
 * Points to a tree in the gene tree kernel
 * @author Jordan Douglas
 *
 */
public class GTKPointerTree extends Tree {

	final public Input<GeneTreeKernel> kernelInput = 
			new Input<>("kernel", "The gene tree kernel which contains a set of trees", Input.Validate.REQUIRED);
	
	final public Input<IntegerParameter> indicatorInput = 
			new Input<>("indicator", "A parameter which points each observed gene tree to one in the kernel", Input.Validate.REQUIRED);
	
	final public Input<Integer> treeIndexInput = 
			new Input<>("index", "Indexes this tree into the 'indicator'. An integer from 0 to g-1, where g is the number of gene trees.", Input.Validate.REQUIRED);
	
	
	int treeIndex;
	int indicator;
	IntegerParameter indicatorVector;
	GeneTreeKernel kernel;
	
	@Override
    public void initAndValidate() {
		 super.initAndValidate();
		 
		 this.treeIndex = treeIndexInput.get();
		 this.indicatorVector = indicatorInput.get();
		 this.indicator = -1; 
		 this.kernel = kernelInput.get();
		 
	}
	
	
	/**
	 * @return The index of the tree in the gene tree kernel which this tree is pointing to
	 */
	public int getIndicatorValue() {
		return this.indicatorVector.getValue(this.treeIndex);
	}
	
	/**
	 * @return The index of this tree, ie. the index of the element in 'indicatorVector' which this tree corresponds to
	 */
	public int getTreeIndex() {
		return this.treeIndex;
	}
	
	
	/**
	 * Updates the tree this class is pointing to and starts editing
	 */
	public void updateTree() {
		this.startEditing(null);
		this.update();
	}
	
	/**
	 * Updates the tree this class is pointing to without editing
	 */
	public void update() {
		int newIndicator = this.getIndicatorValue();
		// TODO: only update when required
		//if (newIndicator != this.indicator || this.hasStartedEditing) { 
			this.indicator = newIndicator;
			Tree newTree = this.kernel.getTree(this.indicator);
			this.setTree(newTree);
		//}
	}
	
	
	/**
	 * Sets the topology/branch lengths of this tree by pointing to a reference tree
	 * But the data at the tips do not change
	 */
	private void setTree(Tree tree) {
		this.setRoot(tree.getRoot());
		
		this.root = tree.getRoot();
        nodeCount = tree.getNodeCount();
        this.m_nodes = tree.getNodesAsArray();
		//this.initArrays();
		this.setEverythingDirty(true);

        // ensure root is the last node
        if (m_nodes != null && root.getNr() != m_nodes.length - 1) {
            final int rootPos = m_nodes.length - 1;
            Node tmp = m_nodes[rootPos];
            m_nodes[rootPos] = root;
            m_nodes[root.getNr()] = tmp;
            tmp.setNr(root.getNr());
            m_nodes[rootPos].setNr(rootPos);
        }
        
        
        
	}

	

	@Override
	public void store() {
		//super.store();
	}
	
	
	@Override
	public void restore() {
		//super.restore();
		update();
	}

	
	@Override
	public GTKPointerTree copy() {
		GTKPointerTree tree = new GTKPointerTree();
        tree.setID(getID());
        tree.index = index;
        tree.root = root.copy();
        tree.nodeCount = nodeCount;
        tree.internalNodeCount = internalNodeCount;
        tree.leafNodeCount = leafNodeCount;
        tree.treeIndex = treeIndex;
        tree.indicator = indicator;
        tree.indicatorVector =  indicatorVector;
        tree.kernel= kernel;
		return tree;
	}
	
	
	@Override
	public List<BEASTInterface> listActiveBEASTObjects() {
		List<BEASTInterface> beastObjects = super.listActiveBEASTObjects();
		
		// Must remove state nodes, otherwise this class will 
		// be treated as a CalculationNode and stored after the proposal
		for (int i = 0; i < beastObjects.size(); i ++) {
			BEASTInterface beastObject = beastObjects.get(i);
			if (beastObject == this.kernel) {
				beastObjects.remove(i);
				break;
			}
		}
		
		for (int i = 0; i < beastObjects.size(); i ++) {
			BEASTInterface beastObject = beastObjects.get(i);
			if (beastObject == this.indicatorVector) {
				beastObjects.remove(i);
				break;
			}
		}
		
		return beastObjects;
		
	}

	@Override
    public void setEverythingDirty(final boolean isDirty) {
		super.setEverythingDirty(isDirty);
	}
	
}








