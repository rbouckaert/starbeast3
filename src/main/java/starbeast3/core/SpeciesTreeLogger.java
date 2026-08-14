package starbeast3.core;


import java.io.PrintStream;
import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Loggable;
import beast.base.inference.StateNode;
import beast.base.spec.domain.NonNegativeReal;
import beast.base.spec.inference.parameter.RealVectorParam;
import beast.base.spec.type.Tensor;
import beast.base.spec.type.Vector;
import starbeast3.evolution.speciation.SpeciesTreePrior;
import beast.base.evolution.speciation.TreeTopFinder;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;



@Description("Logs tree annotated with metadata in StarBeast format")
public class SpeciesTreeLogger extends BEASTObject implements Loggable {
    final public Input<Tree> treeInput = new Input<>("tree", "tree to be logged", Validate.REQUIRED);
    final public Input<RealVectorParam<NonNegativeReal>> parameterInput = new Input<>("popSize", "population size parameter associated with tree nodes", Validate.REQUIRED);
    final public Input<RealVectorParam<NonNegativeReal>> parameterTopInput = new Input<>("popSizeTop", "population size parameter associated with top of tree branches, only used for non-constant *beast analysis");
    // final public Input<Function> branchRatesInput = new Input<>("realRates", "Real branch rates associated with tree nodes", Validate.REQUIRED);
    final public Input<SpeciesTreePrior> speciesTreePriorInput = new Input<>("speciesTreePrior", "species tree prior, used to find which Population Size Function is used. If not specified, assumes 'constant'");
    final public Input<TreeTopFinder> treeTopFinderInput = new Input<>("treetop", "calculates height of species tree", Validate.REQUIRED);
    final public Input<List<Tensor<?, ?>>> metadataInput = new Input<>("metadata", "meta data to be logged with the tree nodes", new ArrayList<>());
    final public Input<Boolean> sortTreeInput = new Input<>("sort", "whether to sort the tree before logging.", true);

    
    // TreePopSizeFunction popSizeFunction;
    String metaDataLabel;

    static final String pop_label = "pop";
    static final String rate_label = "rate";

    @Override
    public void initAndValidate() {
        metaDataLabel = "[&" + pop_label + "=";
//        if (speciesTreePriorInput.get() != null) {
//            popSizeFunction = speciesTreePriorInput.get().popFunctionInput.get();
//        } else {
//            popSizeFunction = TreePopSizeFunction.constant;
//        }
    }
    


    @Override
    public void init(final PrintStream out) {
        treeInput.get().init(out);
    }

    @Override
    public void log(final long sample, final PrintStream out) {

    	
        final Tree tree = (Tree) treeInput.get().getCurrent();
        RealVectorParam<NonNegativeReal> metadata = parameterInput.get();
        RealVectorParam<NonNegativeReal> metadataTop = parameterTopInput.get();
        List<Tensor<?, ?>> metadataList = metadataInput.get();     

        // write out the log tree with meta data
        out.print("tree STATE_" + sample + " = ");
        if (sortTreeInput.get()) tree.getRoot().sort();
        out.print(toNewick(tree.getRoot(), metadata, metadataTop, /*branchrateMetadata,*/ metadataList));
        //out.print(tree.getRoot().toShortNewick(false));
        out.print(";");
    }


    String toNewick(final Node node, final RealVectorParam<NonNegativeReal> popBtm, final RealVectorParam<NonNegativeReal> popTop, List<Tensor<?, ?>> metadataList) {
        final StringBuilder buf = new StringBuilder();

        if (node.getLeft() != null) {
            buf.append("(");
            buf.append(toNewick(node.getLeft(), popBtm, popTop, metadataList));
            if (node.getRight() != null) {
                buf.append(',');
                buf.append(toNewick(node.getRight(), popBtm, popTop, metadataList));
            }
            buf.append(")");
        } else {
            buf.append(node.getNr()+Tree.taxaTranslationOffset);
        }
        buf.append("[&");
        final double popStart = popBtm.get(node.getNr());
        buf.append(pop_label + "=").append(popStart);

        
        StringBuffer buf2 = new StringBuffer();
		if (metadataList.size() > 0) {
			buf2.append(",");
			boolean needsComma = false;
			for (Tensor<?,?> metadata : metadataList) {
				if (metadata instanceof Vector) {
					
					// Skip node if there is no metadata
					if (node.getNr() >= metadata.size()) {
						continue;
					}
					
					if (needsComma) {
						buf2.append(",");
					}
					
					buf2.append(((BEASTObject) metadata).getID());
					buf2.append('=');
					buf2.append(metadata.get(node.getNr()));
					needsComma = true;

				} else {
					
					// It is a scalar, so just add metadata to node 0
					if (metadata.size() > node.getNr()) {
						if (needsComma) {
							buf2.append(",");
						}
						buf2.append(((BEASTObject) metadata).getID());
						buf2.append('=');
						buf2.append(metadata.get(node.getNr()));
						needsComma = true;
					}
				}
			}
			
		}
		
        if (buf2.length() > 3) {
			buf.append(buf2.toString());
		}
        
        buf.append(']');
        if (!node.isRoot()) {
            buf.append(":").append(node.getLength());
        }
        return buf.toString();
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
