package genekernel;

import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;

public class GTKGeneTree extends Tree {
	
	
	public GTKGeneTree() {
		
	}
	
	public GTKGeneTree(final Node rootNode) {
        setRoot(rootNode);
        initArrays();
    }
	
	
	 @Override
	 public void store() {
		 super.store();
	 }

}
