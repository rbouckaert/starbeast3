package starbeast3.tree;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.core.Description;
import beast.core.Input;
import beast.core.StateNode;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;

@Description("Binary tree with efficient store/restore")
public class BinaryTree extends StateNode implements TreeInterface {
    final public Input<TaxonSet> taxonsetInput = new Input<>("taxonset",
            "set of taxa that correspond to the leafs in the tree");

	int [] parent;
	int [] left;
	int [] right;
	double [] height;
	
	int [] storedParent;
	int [] storedLeft;
	int [] storedRight;
	double [] storedHeight;
	
	String [] taxa;

	protected BinaryNode[] nodes = null;

	@Override
	public void initAndValidate() {
		List<String> taxaList = taxonsetInput.get().asStringList();
		taxa = taxaList.toArray(new String[] {});
		int leafCount = taxa.length;
		int nodeCount = leafCount * 2 -1;
				
				
		nodes = new BinaryNode[nodeCount];
		for (int i = 0; i < nodeCount; i++) {
			nodes[i] = new BinaryNode(i, this);
		}
	}

	
    /**
     * print translate block for NEXUS beast.tree file
     */
    public static void printTranslate(final Node node, final PrintStream out, final int nodeCount) {
        final List<String> translateLines = new ArrayList<>();
        printTranslate(node, translateLines, nodeCount);
        Collections.sort(translateLines);
        for (final String line : translateLines) {
            out.println(line);
        }
    }

    static public int taxaTranslationOffset = 1;

    /**
     * need this helper so that we can sort list of entries *
     */
    static void printTranslate(Node node, List<String> translateLines, int nodeCount) {
        if (node.isLeaf()) {
            final String nr = (node.getNr() + taxaTranslationOffset) + "";
            String line = "\t\t" + "    ".substring(nr.length()) + nr + " ";
            if (node.getID().indexOf(' ') > 0) {
            	char c = node.getID().charAt(0);
            	if (c == '\"' || c == '\'') {
                	line += node.getID();
            	} else {
            		line += '\"' + node.getID() + "\"";
            	}
            } else {
            	line += node.getID();
            }
            if (node.getNr() < nodeCount) {
                line += ",";
            }
            translateLines.add(line);
        } else {
            printTranslate(node.getLeft(), translateLines, nodeCount);
            if (node.getRight() != null) {
                printTranslate(node.getRight(), translateLines, nodeCount);
            }
        }
    }

    public static void printTaxa(final Node node, final PrintStream out, final int nodeCount) {
        final List<String> translateLines = new ArrayList<>();
        printTranslate(node, translateLines, nodeCount);
        Collections.sort(translateLines);
        for (String line : translateLines) {
            line = line.substring(line.indexOf(" ", 5)).replace(',', ' ').trim();
            out.println("\t\t\t" + line);
        }
    }
    
    @Override
	public void init(PrintStream out) {
        Node node = getRoot();
        out.println("#NEXUS\n");
        out.println("Begin taxa;");
        out.println("\tDimensions ntax=" + getLeafNodeCount() + ";");
        out.println("\t\tTaxlabels");
        printTaxa(node, out, getNodeCount() / 2);
        out.println("\t\t\t;");
        out.println("End;");

        out.println("Begin trees;");
        out.println("\tTranslate");
        printTranslate(node, out, getNodeCount() / 2);
        out.print(";");
	}

    @Override
	public void log(long sample, PrintStream out) {
        out.print("tree STATE_" + sample + " = ");
        // Don't sort, this can confuse CalculationNodes relying on the tree
        //tree.getRoot().sort();
        final int[] dummy = new int[1];
        final String newick = ((BinaryNode)getRoot()).toSortedNewick(dummy);
        out.print(newick);
        out.print(";");
    }

    @Override
	public void close(PrintStream out) {
        out.print("End;");
	}

	@Override
	public int getDimension() {
        return getNodeCount();
	}

	@Override
	public double getArrayValue() {
        return getRoot().getHeight();
    }

    @Override
	public double getArrayValue(int value) {
        return nodes[value].getHeight();
	}

	@Override
	public int getLeafNodeCount() {
		return taxa.length;
	}

	@Override
	public int getInternalNodeCount() {
		return taxa.length - 1;
	}

	@Override
	public int getNodeCount() {		
		return nodes.length;
	}

	@Override
	public Node getRoot() {
		return nodes[nodes.length - 1];
	}

	@Override
	public Node getNode(int i) {
		return nodes[i];
	}

	@Override
	public Node[] getNodesAsArray() {
		return nodes.clone();
	}

	@Override
	public List<Node> getExternalNodes() {
		List<Node> externalNodes = new ArrayList<>();
		for (int i = 0; i < taxa.length; i++) {
			externalNodes.add(nodes[i]);
		}
		return externalNodes;
	}

	@Override
	public List<Node> getInternalNodes() {
		List<Node> internalNodes = new ArrayList<>();
		for (int i = taxa.length; i < taxa.length*2 - 1; i++) {
			internalNodes.add(nodes[i]);
		}
		return internalNodes;
	}

	@Override
	public TaxonSet getTaxonset() {
        return taxonsetInput.get();
	}

	@Override
	public void getMetaData(Node node, Double[] t, String pattern) {
		// TODO Auto-generated method stub

	}

	@Override
	public void setMetaData(Node node, Double[] t, String pattern) {
		// TODO Auto-generated method stub

	}

	@Override
	public void setEverythingDirty(boolean isDirty) {
        setSomethingIsDirty(isDirty);
	}

	@Override
	public StateNode copy() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void assignTo(StateNode other) {
		// TODO Auto-generated method stub

	}

	@Override
	public void assignFrom(StateNode other) {
		// TODO Auto-generated method stub

	}

	@Override
	public void assignFromFragile(StateNode other) {
		// TODO Auto-generated method stub

	}

	@Override
	public void fromXML(org.w3c.dom.Node node) {
		// TODO Auto-generated method stub

	}

	@Override
	public int scale(double scale) {
        startEditing();
        for (int i = taxa.length; i < taxa.length * 2 -1 ; i++) {
        	height[i] *= scale;
        }
		return taxa.length - 1;
	}

	@Override
	protected void store() {
		System.arraycopy(height, 0, storedHeight, 0, height.length);
		System.arraycopy(left, 0, storedLeft, 0, left.length);
		System.arraycopy(right, 0, storedRight, 0, right.length);
		System.arraycopy(parent, 0, storedParent, 0, parent.length);
	}

	@Override
	public void restore() {
		double [] tmp = height; height = storedHeight; storedHeight = tmp;
		int [] tmp2 = left; left = storedLeft; storedLeft = tmp2;
		int [] tmp3 = right; right = storedRight; storedRight = tmp3; 
		int [] tmp4 = parent; parent = storedParent; storedParent = tmp4; 
	}

	
	void startEditing() {
		startEditing(null);
	}
}
