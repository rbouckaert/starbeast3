package starbeast3.tree;

import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;

import beast.core.BEASTInterface;
import beast.core.Description;
import beast.core.StateNode;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.util.TreeParser;

@Description("Binary tree with efficient store/restore")
public class BinaryTree extends Tree implements TreeInterface {
//    final public Input<TaxonSet> taxonsetInput = new Input<>("taxonset",
//            "set of taxa that correspond to the leafs in the tree");

	int [] parent;
	int [] left;
	int [] right;
	double [] height;
	
	int [] storedParent;
	int [] storedLeft;
	int [] storedRight;
	double [] storedHeight;
	
	protected BinaryNode[] nodes = null;

	@Override
	public void initAndValidate() {
		List<String> taxaList = m_taxonset.get().asStringList();
		m_sTaxaNames = taxaList.toArray(new String[] {});
		leafNodeCount = m_sTaxaNames.length;
		nodeCount = leafNodeCount * 2 - 1;
		internalNodeCount = leafNodeCount - 1;
		
		parent = new int[nodeCount];
		left = new int[nodeCount];
		right = new int[nodeCount];
		height = new double[nodeCount];
				
		storedParent = new int[nodeCount];
		storedLeft = new int[nodeCount];
		storedRight = new int[nodeCount];
		storedHeight = new double[nodeCount];
				
		nodes = new BinaryNode[nodeCount];
		for (int i = 0; i < nodeCount; i++) {
			nodes[i] = new BinaryNode(i, this);
		}
		
		root = nodes[nodeCount - 1];
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
		return m_sTaxaNames.length;
	}

	@Override
	public int getInternalNodeCount() {
		return m_sTaxaNames.length - 1;
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
		for (int i = 0; i < m_sTaxaNames.length; i++) {
			externalNodes.add(nodes[i]);
		}
		return externalNodes;
	}

	@Override
	public List<Node> getInternalNodes() {
		List<Node> internalNodes = new ArrayList<>();
		for (int i = m_sTaxaNames.length; i < m_sTaxaNames.length*2 - 1; i++) {
			internalNodes.add(nodes[i]);
		}
		return internalNodes;
	}

	@Override
	public TaxonSet getTaxonset() {
        return m_taxonset.get();
	}

	   /**
     * copy meta data matching pattern to double array
     *
     * @param node     the node
     * @param t       the integer array to be filled with meta data
     * @param pattern the name of the meta data
     */
    public void getMetaData(final Node node, final Integer[] t, final String pattern) {
        t[Math.abs(node.getNr())] = (Integer) node.getMetaData(pattern);
        if (!node.isLeaf()) {
            getMetaData(node.getLeft(), t, pattern);
            if (node.getRight() != null) {
                getMetaData(node.getRight(), t, pattern);
            }
        }
    }

    /**
     * copy meta data matching pattern to double array
     *
     * @param node     the node
     * @param t       the double array to be filled with meta data
     * @param pattern the name of the meta data
     */
    @Override
	public void getMetaData(final Node node, final Double[] t, final String pattern) {
        t[Math.abs(node.getNr())] = (Double) node.getMetaData(pattern);
        if (!node.isLeaf()) {
            getMetaData(node.getLeft(), t, pattern);
            if (node.getRight() != null) {
                getMetaData(node.getRight(), t, pattern);
            }
        }
    }

    /**
     * traverse tree and assign meta-data values in t to nodes in the
     * tree to the meta-data field represented by the given pattern.
     * This only has an effect when setMetadata() in a subclass
     * of Node know how to process such value.
     */
    @Override
	public void setMetaData(Node node, Double[] t, String pattern) {
        node.setMetaData(pattern, t[Math.abs(node.getNr())]);
        if (!node.isLeaf()) {
            setMetaData(node.getLeft(), t, pattern);
            if (node.getRight() != null) {
                setMetaData(node.getRight(), t, pattern);
            }
        }
	}

	@Override
	public void setEverythingDirty(boolean isDirty) {
        setSomethingIsDirty(isDirty);
        if (!isDirty) {
            for( Node n : nodes ) {
                n.makeDirty(IS_CLEAN);
            }
        } else {
            for( Node n : nodes ) {
                n.makeDirty(IS_FILTHY);
            }
        }
	}

	@Override
	public Tree copy() {
        BinaryTree tree = new BinaryTree();
        tree.initByName("taxonset", m_taxonset.get());
        tree.setID(getID());
        tree.index = index;
        System.arraycopy(height, 0, tree.height, 0, height.length);
		System.arraycopy(left, 0, tree.left, 0, left.length);
		System.arraycopy(right, 0, tree.right, 0, right.length);
		System.arraycopy(parent, 0, tree.parent, 0, parent.length);
        return tree;
	}

	@Override
	public void assignTo(StateNode other) {
		throw new RuntimeException("Not implemented yet");
	}

	@Override
	public void assignFrom(StateNode other) {
		if (other instanceof BinaryTree) {
	        final BinaryTree tree = (BinaryTree) other;        
	        initByName("taxonset", tree.getTaxonset());
	        assignFromFragile(other);
	        setID(tree.getID());
		} else if (other instanceof TreeInterface) {
	        final TreeInterface tree = (TreeInterface) other;
	        if (tree.getTaxonset() != null) {
	        	initByName("taxonset", tree.getTaxonset());
	        } else {
	        	initByName("taxa", ((BEASTInterface)tree).getInput("taxa").get());
	        }
	        setID(tree.getID());
	        assignFromFragile(other);
		} else {
			throw new IllegalArgumentException("Expected state node of type tree");
		}
	}

	@Override
	public void assignFromFragile(StateNode other) {
		if (other instanceof BinaryTree) {
	        final BinaryTree tree = (BinaryTree) other;        
			System.arraycopy(tree.height, 0, height, 0, height.length);
			System.arraycopy(tree.left, 0, left, 0, left.length);
			System.arraycopy(tree.right, 0, right, 0, right.length);
			System.arraycopy(tree.parent, 0, parent, 0, parent.length);
		} else if (other instanceof TreeInterface) {
	        final TreeInterface tree = (TreeInterface) other;	        
	        Node [] otherNodes = tree.getNodesAsArray();
	        for (int i = 0; i < otherNodes.length; i++) {
	        	Node node = otherNodes[i];
	        	height[i] = node.getHeight();
	        	left[i] = node.getLeft() == null ? -1 : node.getLeft().getNr();
	        	right[i] = node.getRight() == null ? -1 : node.getRight().getNr();
	        	parent[i] = node.getParent() == null ? -1 : node.getParent().getNr();
	        }        
		} else {
			throw new IllegalArgumentException("Expected state node of type tree");
		}	
	}

	@Override
	public String toString() {
	    return getRoot().toString();
	}

	@Override
	public void fromXML(org.w3c.dom.Node node) {
        final String newick = node.getTextContent();
        final TreeParser parser = new TreeParser();
        try {
            parser.thresholdInput.setValue(1e-10, parser);
        } catch (Exception e1) {
            e1.printStackTrace();
        }
        try {
            parser.offsetInput.setValue(0, parser);
            parser.parseNewick(newick);
            assignFromFragile(parser);
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

	@Override
	public int scale(double scale) {
        startEditing();
        for (int i = m_sTaxaNames.length; i < m_sTaxaNames.length * 2 -1 ; i++) {
        	height[i] *= scale;
        }
		return m_sTaxaNames.length - 1;
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

    public void setRoot(final Node newroot) {
        final int n = nodes.length - 1;
        final int m = newroot.getNr();
        if (m != n) {
            // ensure root is the last node
        	
        	int leftRoot = left[m];
        	if (leftRoot == n) {leftRoot = m;}
        	int rightRoot = right[m];
        	if (rightRoot == n) {rightRoot = m;}

        	int leftI = left[n];
        	int rightI = right[n];
        	int parentI = parent[n];
        	if (parentI == m) {
        		parentI = n;
        	}
        	
        	left[n] = leftRoot;
        	parent[leftRoot] = n;
        	right[n] = rightRoot;
        	parent[rightRoot] = n;
        	parent[n] = -1;
        	
        	left[m] = leftI;
        	parent[leftI] = m;
        	right[m] = rightI;
        	parent[rightI] = m;
        	parent[m] = parentI;

        	if (left[parentI] == n) {
        		left[parentI] = m;
        	}
            if (right[parentI] == n) {
        		right[parentI] = m;
        	}

        	double tmp2 = height[n];
        	height[n] = height[m];
        	height[m] = tmp2;
        	
        	BinaryNode tmp3 = nodes[n];
        	nodes[n] = nodes[m];
        	nodes[m] = tmp3;
        	nodes[n].setNr(n);
        	nodes[m].setNr(m);
        }
    }
}
