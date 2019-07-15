package starbeast3.operators;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;



@Description("Tree operator which randomly changes the height of a node, " +
        "then reconstructs the tree from node heights.")
public class NodeReheight2 extends Operator {
    public final Input<Tree> treeInput = new Input<>("tree", "the species tree", Validate.REQUIRED);
    public final Input<TaxonSet> taxonSetInput = new Input<>("taxonset", "taxon set describing species tree taxa and their gene trees", Validate.REQUIRED);
    public final Input<List<Tree>> geneTreesInput = new Input<>("geneTree", "list of gene trees that constrain species tree movement", new ArrayList<>());
    Node[] m_nodes;


    /**
     * map node number of leafs in gene trees to leaf nr in species tree *
     */
    //List<Map<Integer, Integer>> m_taxonMap;
    int[][] m_taxonMap;
    int nrOfGeneTrees;
    int nrOfSpecies;

    @Override
    public void initAndValidate() {
        /** build the taxon map for each gene tree **/
        m_taxonMap = new int[geneTreesInput.get().size()][];
        
        final Map<String, Integer> tipNumberMap = new LinkedHashMap<>();

        // generate map of gene tree tip node names to species tree tip node numbers
        Map<String, TaxonSet> speciesMap = new LinkedHashMap<>();
        for (Taxon ts : taxonSetInput.get().taxonsetInput.get()) {
        	speciesMap.put(ts.getID(), (TaxonSet)ts);
        }
        
        for (Node speciesNode : treeInput.get().getExternalNodes()) {
            final String speciesName = speciesNode.getID();
            int speciesNumber = speciesNode.getNr();

            final Set<Taxon> tipSet = new LinkedHashSet<>(speciesMap.get(speciesName).taxonsetInput.get());
            for (Taxon tip: tipSet) {
                final String tipName = tip.getID();
                tipNumberMap.put(tipName, speciesNumber);
            }
        }
        
        int i = 0;
        for (final Tree geneTree : geneTreesInput.get()) {
            int [] localTipNumberMap = new int[geneTree.getLeafNodeCount()];
            for (int j = 0; j < geneTree.getLeafNodeCount(); j++) {
            	final Node geneTreeLeafNode = geneTree.getNode(j);
            	final String geneTreeLeafName = geneTreeLeafNode.getID();
            	final int geneTreeLeafNumber = geneTreeLeafNode.getNr();
            	localTipNumberMap[geneTreeLeafNumber] = tipNumberMap.get(geneTreeLeafName);
            }

            m_taxonMap[i++] = localTipNumberMap;
        }

        nrOfGeneTrees = geneTreesInput.get().size();
        nrOfSpecies = treeInput.get().getLeafNodeCount();
    }

    @Override
    public double proposal() {
        final Tree tree = treeInput.get();
        m_nodes = tree.getNodesAsArray();
        final int nodeCount = tree.getNodeCount();
        // randomly change left/right order
        tree.startEditing(this);  // we change the tree
        reorder(tree.getRoot());
        // collect heights
        final double[] heights = new double[nodeCount];
        final int[] reverseOrder = new int[nodeCount];
        collectHeights(tree.getRoot(), heights, reverseOrder, 0);
        // change height of an internal, non-sampled-ancestor node
        int nodeIndex = Randomizer.nextInt(heights.length);
        Node node = m_nodes[reverseOrder[nodeIndex]];
        while (node.isLeaf() || node.isFake()) {
            nodeIndex = Randomizer.nextInt(heights.length);
            node = m_nodes[reverseOrder[nodeIndex]];
        }
        double maxHeight = calcMaxHeight(reverseOrder, nodeIndex);
        double minHeight = calcMinHeight(node);
        // debugging code to ensure the new maxHeight equals the original one
        /*
        double maxHeight2 = calcMaxHeight2(reverseOrder, nodeIndex);
        if (Math.abs(maxHeight - maxHeight2) > 1e-10) {
            maxHeight = calcMaxHeight(reverseOrder, nodeIndex);
            maxHeight2 = calcMaxHeight2(reverseOrder, nodeIndex);
        	
        }
        */

        // sample a new height compatible with the gene trees
        heights[nodeIndex] = minHeight + (Randomizer.nextDouble() * (maxHeight - minHeight));
        node.setHeight(heights[nodeIndex]);

        // leave sampled ancestors connected to their fake parents
        boolean[] hasParent = new boolean[heights.length];
        for (int i = 0; i < heights.length; i++)
        	hasParent[i] = m_nodes[reverseOrder[i]].isDirectAncestor();

        // reconstruct tree from heights
        final Node root = reconstructTree(heights, reverseOrder, 0, heights.length, hasParent);

        assert checkConsistency(root, new boolean[heights.length]) ;
        //            System.err.println("Inconsistent tree");
        //        }
        root.setParent(null);
        tree.setRoot(root);
        return 0;
    }

    private double calcMinHeight(Node node) {
		if (node.isLeaf()) {
			return node.getHeight();
		}
		
		double minHeight = 0.0;
		for (Node child: node.getChildren()) {
			minHeight = Math.max(minHeight, calcMinHeight(child));
		}
		
		return minHeight;
	}

	private boolean checkConsistency(final Node node, final boolean[] used) {
        if (used[node.getNr()]) {
            // used twice? tha's bad
            return false;
        }
        used[node.getNr()] = true;
        if ( node.isLeaf() ) {
            return true;
        }
        return checkConsistency(node.getLeft(), used) && checkConsistency(node.getRight(), used);
    }

    /**
     * calculate maximum height that node nodeIndex can become restricted
     * by nodes on the left and right
     */
    private double calcMaxHeight(final int[] reverseOrder, final int nodeIndex) {
        // find maximum height between two species. Only upper right part is populated
        final double[][] maxHeight = new double[nrOfSpecies][nrOfSpecies];
        for (int i = 0; i < nrOfSpecies; i++) {
            Arrays.fill(maxHeight[i], Double.POSITIVE_INFINITY);
        }

        // find species on the left of selected node
        final boolean[] isLowerSpecies = new boolean[nrOfSpecies];
        final Node[] nodes = treeInput.get().getNodesAsArray();
        for (int i = 0; i < nodeIndex; i++) {
            final Node node = nodes[reverseOrder[i]];
            if (node.isLeaf()) {
                isLowerSpecies[node.getNr()] = true;
            }
        }
        // find species on the right of selected node
        final boolean[] isUpperSpecies = new boolean[nrOfSpecies];
        for (int i = nodeIndex + 1; i < nodes.length; i++) {
            final Node node = nodes[reverseOrder[i]];
            if (node.isLeaf()) {
                isUpperSpecies[node.getNr()] = true;
            }
        }

        final boolean[] isUsedSpecies = new boolean[nrOfSpecies];
        for (int i = 0; i < nrOfSpecies; i++) {
            isUsedSpecies[i] = isLowerSpecies[i] || isUpperSpecies[i];
        }
        
        // calculate for every species tree the maximum allowable merge point
        for (int i = 0; i < nrOfGeneTrees; i++) {
            final Tree tree = geneTreesInput.get().get(i);
            findMaximaInGeneTree(tree.getRoot(), m_taxonMap[i], maxHeight, isUsedSpecies);
        }

        // find max
        double max = Double.POSITIVE_INFINITY;
        for (int i = 0; i < nrOfSpecies; i++) {
            if (isLowerSpecies[i]) {
                for (int j = 0; j < nrOfSpecies; j++) {
                    if (j != i && isUpperSpecies[j]) {
                        // final int x = Math.min(i, j);
                        // final int y = Math.max(i, j);
                        max = Math.min(max, maxHeight[i][j]);
                        max = Math.min(max, maxHeight[j][i]);
                    }
                }
            }
        }
        return max;
    } // calcMaxHeight


    /**
     * calculate maximum height that node nodeIndex can become restricted
     * by nodes on the left and right
     * 
     * only used for debugging now
     */
    /* private double calcMaxHeight2(final int[] reverseOrder, final int nodeIndex) {
        // find maximum height between two species. Only upper right part is populated
        final double[][] maxHeight = new double[nrOfSpecies][nrOfSpecies];
        for (int i = 0; i < nrOfSpecies; i++) {
            Arrays.fill(maxHeight[i], Double.POSITIVE_INFINITY);
        }

        // calculate for every species tree the maximum allowable merge point
        for (int i = 0; i < nrOfGeneTrees; i++) {
            final GeneTree tree = geneTreesInput.get().get(i);
            findMaximaInGeneTree(tree.getRoot(), new boolean[nrOfSpecies], m_taxonMap[i], maxHeight);
        }

        // find species on the left of selected node
        final boolean[] isLowerSpecies = new boolean[nrOfSpecies];
        final Node[] nodes = treeInput.get().getNodesAsArray();
        for (int i = 0; i < nodeIndex; i++) {
            final Node node = nodes[reverseOrder[i]];
            if (node.isLeaf()) {
                isLowerSpecies[node.getNr()] = true;
            }
        }
        // find species on the right of selected node
        final boolean[] isUpperSpecies = new boolean[nrOfSpecies];
        for (int i = nodeIndex + 1; i < nodes.length; i++) {
            final Node node = nodes[reverseOrder[i]];
            if (node.isLeaf()) {
                isUpperSpecies[node.getNr()] = true;
            }
        }

        // find max
        double max = Double.POSITIVE_INFINITY;
        for (int i = 0; i < nrOfSpecies; i++) {
            if (isLowerSpecies[i]) {
                for (int j = 0; j < nrOfSpecies; j++) {
                    if (j != i && isUpperSpecies[j]) {
                        final int x = Math.min(i, j);
                        final int y = Math.max(i, j);
                        max = Math.min(max, maxHeight[x][y]);
                    }
                }
            }
        }
        return max;
    } */

    /**
     * for every species in the left on the gene tree and for every species in the right
     * cap the maximum join height by the lowest place the two join in the gene tree
     * 
     * only used for debugging now
     */
    /* private void findMaximaInGeneTree(final Node node, final boolean[] taxonSet, final int [] taxonMap, final double[][] maxHeight) {
        if (node.isLeaf()) {
            final int species = taxonMap[node.getNr()];
            taxonSet[species] = true;
        } else {
            final boolean[] isLeftTaxonSet = new boolean[nrOfSpecies];
            findMaximaInGeneTree(node.getLeft(), isLeftTaxonSet, taxonMap, maxHeight);
            final boolean[] isRightTaxonSet = new boolean[nrOfSpecies];
            findMaximaInGeneTree(node.getRight(), isRightTaxonSet, taxonMap, maxHeight);
            for (int i = 0; i < nrOfSpecies; i++) {
                if (isLeftTaxonSet[i]) {
                    for (int j = 0; j < nrOfSpecies; j++) {
                        if (j != i && isRightTaxonSet[j]) {
                            final int x = Math.min(i, j);
                            final int y = Math.max(i, j);
                            maxHeight[x][y] = Math.min(maxHeight[x][y], node.getHeight());
                        }
                    }
                }
            }
            for (int i = 0; i < nrOfSpecies; i++) {
                taxonSet[i] = isLeftTaxonSet[i] | isRightTaxonSet[i];
            }
        }
    } */

    /**
     * for every species in the left on the gene tree and for every species in the right
     * cap the maximum join height by the lowest place the two join in the gene tree
     */
    private void findMaximaInGeneTree(final Node nodeX, final int [] taxonMap, final double[][] maxHeight, boolean[] isUsedSpecies) {
        Tree tree = nodeX.getTree();
        int nrOfNodes = tree.getNodeCount();
        int [][] speciesList = new int[nrOfNodes][nrOfSpecies];
        int [] speciesCount = new int[nrOfNodes];
        for (Node node : tree.listNodesPostOrder(null, null)) { 	
            if (node.isLeaf()) {
                int nodeNr = node.getNr();
                final int species = taxonMap[nodeNr];
                if (isUsedSpecies[species]) {
                    speciesList[nodeNr][0] = species;
                    speciesCount[nodeNr] = 1;
                }
            } else {
                int left = node.getLeft().getNr();
                int right = node.getRight().getNr();
                for (int i = 0; i < speciesCount[left]; i++) {
                    for (int j = 0; j < speciesCount[right]; j++) {
                        if (speciesList[left][i] != speciesList[right][j]) {
                            int sp1 = speciesList[left][i];
                            int sp2 = speciesList[right][j];
                            final int x;
                            final int y;
                            if (sp1 < sp2) {
                                x = sp1; y = sp2;
                            } else {
                                x = sp2; y = sp1;
                            }
                            maxHeight[x][y] = Math.min(maxHeight[x][y], node.getHeight());
                        }
                    }
                }
                int i = 0;
                int j = 0;
                int k = 0;
                int nodeNr = node.getNr();
                int [] spList = speciesList[nodeNr];
                int [] leftList = speciesList[left];
                int [] rightList = speciesList[right];
                while (i < speciesCount[left] && j < speciesCount[right]) {
                    if (leftList[i] == rightList[j]) {
                        spList[k++] = leftList[i];
                        i++;
                        j++;
                    } else if (leftList[i] < rightList[j]) {
                        spList[k++] = leftList[i++];
                    } else {
                        spList[k++] = rightList[j++];
                    }
                }
                if (i == speciesCount[left]) {
                    while (j < speciesCount[right]) {
                        spList[k++] = rightList[j++];
                    }
                } else if (j == speciesCount[right]) {
                    while (i < speciesCount[left]) {
                        spList[k++] = leftList[i++];
                    }
                }
                speciesCount[node.getNr()] = k;
            }
        }
    }

    private Node reconstructTree(final double[] heights, final int[] reverseOrder, final int from, final int to, final boolean[] hasParent) {
        //nodeIndex = maxIndex(heights, 0, heights.length);
        int nodeIndex = -1;
        double max = Double.NEGATIVE_INFINITY;
        for (int j = from; j < to; j++) {
            if (max < heights[j] && !m_nodes[reverseOrder[j]].isLeaf()) {
                max = heights[j];
                nodeIndex = j;
            }
        }
        if (nodeIndex < 0) {
            return null;
        }
        final Node node = m_nodes[reverseOrder[nodeIndex]];
        final boolean keepLeftChild = node.getLeft().isDirectAncestor();
        final boolean keepRightChild = node.getRight().isDirectAncestor();

        if (!keepLeftChild) {
	        //int left = maxIndex(heights, 0, nodeIndex);
	        int left = -1;
	        max = Double.NEGATIVE_INFINITY;
	        for (int j = from; j < nodeIndex; j++) {
	            if (max < heights[j] && !hasParent[j]) {
	                max = heights[j];
	                left = j;
	            }
	        }
	        /* if (left >= reverseOrder.length || left < 0) {
                System.out.println("(reverseOrder[" + left + "] out of bounds) Node number = " + node.getNr());
                System.out.println(node.toNewick());
            } */
	        node.setLeft(m_nodes[reverseOrder[left]]);
	        node.getLeft().setParent(node);
	        if (node.getLeft().isLeaf()) {
	            heights[left] = Double.NEGATIVE_INFINITY;
	        }
	        hasParent[left] = true;
        }

        if (!keepRightChild) {
	        int right = -1;
	        max = Double.NEGATIVE_INFINITY;
	        for (int j = nodeIndex + 1; j < to; j++) {
	            if (max < heights[j] && !hasParent[j]) {
	                max = heights[j];
	                right = j;
	            }
	        }
	        node.setRight(m_nodes[reverseOrder[right]]);
	        node.getRight().setParent(node);
	        if (node.getRight().isLeaf()) {
	            heights[right] = Double.NEGATIVE_INFINITY;
	        }
	        hasParent[right] = true;
        }

        heights[nodeIndex] = Double.NEGATIVE_INFINITY;
        reconstructTree(heights, reverseOrder, from, nodeIndex, hasParent);
        reconstructTree(heights, reverseOrder, nodeIndex, to, hasParent);
        return node;
    }

    // helper for reconstructTree, to find maximum in range
    //    private int maxIndex(final double[] heights, final int from, final int to) {
    //        int maxIndex = -1;
    //        double max = Double.NEGATIVE_INFINITY;
    //        for (int i = from; i < to; i++) {
    //            if (max < heights[i]) {
    //                max = heights[i];
    //                maxIndex = i;
    //            }
    //        }
    //        return maxIndex;
    //    }

    /**
     ** gather height of each node, and the node index associated with the height.*
     **/
    private int collectHeights(final Node node, final double[] heights, final int[] reverseOrder, int current) {
        if (node.isLeaf()) {
            heights[current] = node.getHeight();
            reverseOrder[current] = node.getNr();
            current++;
        } else {
            current = collectHeights(node.getLeft(), heights, reverseOrder, current);
            heights[current] = node.getHeight();
            reverseOrder[current] = node.getNr();
            current++;
            current = collectHeights(node.getRight(), heights, reverseOrder, current);
        }
        return current;
    }

    /**
     * randomly changes left and right children in every internal node *
     */
    private void reorder(final Node node) {
        if (!node.isLeaf()) {
            if (Randomizer.nextBoolean()) {
                final Node tmp = node.getLeft();
                node.setLeft(node.getRight());
                node.setRight(tmp);
            }
            reorder(node.getLeft());
            reorder(node.getRight());
        }
    }
} // class NodeReheight
