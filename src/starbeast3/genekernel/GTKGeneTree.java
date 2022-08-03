package starbeast3.genekernel;

import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;

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
