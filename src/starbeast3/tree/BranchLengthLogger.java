package starbeast3.tree;



import java.io.PrintStream;

import beast.base.inference.CalculationNode;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.evolution.tree.Tree;


@Description("Logger to report branch lengths of a tree")
public class BranchLengthLogger extends CalculationNode implements Loggable, Function {
    final public Input<Tree> treeInput = new Input<>("tree", "tree to report branch lengths for.", Validate.REQUIRED);
    
    @Override
    public void initAndValidate() {
        // nothing to do
    }

    @Override
    public void init(PrintStream out) {
        final Tree tree = treeInput.get();
        String id = "";
        if (getID() == null || getID().matches("\\s*")) {
           id = tree.getID() + ".length";
        } else {
           id = getID();
        }
        
        for (int nodeNr = 0; nodeNr < getDimension(); nodeNr++) {
        	out.print(id + "." + (nodeNr+1) + "\t");
        }
        
    }

    @Override
    public void log(long sample, PrintStream out) {
        final Tree tree = treeInput.get();
        for (int nodeNr = 0; nodeNr < getDimension(); nodeNr++) {
        	out.print(tree.getNode(nodeNr).getLength() + "\t");
        }
        
    }

	@Override
    public void close(PrintStream out) {
        // nothing to do
    }

    @Override
    public int getDimension() {
    	final Tree tree = treeInput.get();
        return tree.getNodeCount();
    }
    
    @Override
    public double getArrayValue(int dim) {
    	final Tree tree = treeInput.get();
        return tree.getNode(dim).getLength();
    }
}



