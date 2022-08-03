package starbeast3;


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.Parameter;
import starbeast3.genekernel.GeneTreeKernel;
import beast.base.evolution.speciation.TreeTopFinder;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;



@Description("Logs gene tree annotated with metadata in StarBeast format")
public class GeneTreeLogger extends BEASTObject implements Loggable {
    final public Input<Tree> treeInput = new Input<>("tree", "Tree to be logged", Validate.REQUIRED);
    final public Input<StarBeast3Clock> clockModelInput = new Input<>("clockModel", "Clock model used to compute branch rates", Validate.REQUIRED);
    final public Input<List<Function>> metadataInput = new Input<>("metadata", "meta data to be logged with the tree nodes",new ArrayList<>());

    // TreePopSizeFunction popSizeFunction;
    String metaDataLabel;

    static final String rate_label = "rate";

    @Override
    public void initAndValidate() {
        metaDataLabel = "[&" + rate_label + "=";
    }

    @Override
    public void init(final PrintStream out) {
        treeInput.get().init(out);
    }

    @Override
    public void log(final long sample, final PrintStream out) {
    	
        // Make sure we get the current version of the inputs
        final Tree tree = (Tree) treeInput.get().getCurrent();

        List<Function> metadataList = metadataInput.get();
        for (int i = 0; i < metadataList.size(); i++) {
        	if (metadataList.get(i) instanceof StateNode) {
        		metadataList.set(i, ((StateNode) metadataList.get(i)).getCurrent());
        	}
        }
        
        
        // write out the log tree with meta data
        out.print("tree STATE_" + sample + " = ");
        tree.getRoot().sort();
        out.print(toNewick(tree.getRoot(), metadataList));
        //out.print(tree.getRoot().toShortNewick(false));
        out.print(";");
    }


    String toNewick(final Node node,  List<Function> metadataList) {
    	
       
    	
        final StringBuilder buf = new StringBuilder();

        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft(), metadataList));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight(), metadataList));
            }
            buf.append(")");
        } else {
            buf.append(node.getNr()+Tree.taxaTranslationOffset);
        }
        buf.append("[&");
        
        
        // Get the branch rate of this node
        double branchRate = clockModelInput.get().getRateForBranch(node);
        buf.append(rate_label + "=").append(branchRate);
        
        
        
        if (metadataList.size() > 0) {
        	for (Function metadata2 : metadataList) {
	            if (metadataList.indexOf(metadata2) > 0 || buf.length() > 1) {
	            	buf.append(",");
	            }
	            buf.append(((BEASTObject)metadata2).getID());
	            buf.append('=');
	            if (metadata2 instanceof Parameter<?>) {
	            	Parameter<?> p = (Parameter<?>) metadata2;
	            	int dim = p.getMinorDimension1();
	            	if (dim > 1) {
		            	buf.append('{');
		            	for (int i = 0; i < dim; i++) {
			            	buf.append(p.getMatrixValue(node.getNr(), i));
			            	if (i < dim - 1) {
				            	buf.append(',');
			            	}
		            	}
		            	buf.append('}');
	            	} else {
		            	buf.append(metadata2.getArrayValue(node.getNr()));
	            	}
	            } else {
	            	buf.append(metadata2.getArrayValue(node.getNr()));
	            }
        	}
        }
        buf.append(']');
        if (!node.isRoot()) {
            buf.append(":").append(node.getLength());
        }
        return buf.toString();
    }

    double getMetaDataTopValue(final Node node, final Function metadataTop) {
        int nr = node.getNr();
        if (nr >= metadataTop.getDimension()) {
            nr = node.getTree().getRoot().getNr();
        }
        return metadataTop.getArrayValue(nr);
    }

    @Override
    public void close(final PrintStream out) {
        treeInput.get().close(out);
    }

    @Override
    public boolean notCloneable() {
        return true;
    }
}
