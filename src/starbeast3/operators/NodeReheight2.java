package starbeast3.operators;

import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.Operator;
import beast.core.parameter.RealParameter;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.tree.Node;
import beast.util.Randomizer;
import genekernel.GTKOperator;
import starbeast3.SpeciesTree;

import java.util.List;

/**
 * @author Huw Ogilvie
 */

@Description("Tree operator which randomly changes the height of a node, " +
        "then reconstructs the tree from node heights.")
public class NodeReheight2 extends GTKOperator {
	
    public final Input<SpeciesTree> treeInput = new Input<>("tree", "the species tree", Validate.REQUIRED);
    public final Input<TaxonSet> taxonSetInput = new Input<>("taxonset", "taxon set describing species tree taxa and their gene trees", Validate.REQUIRED); // left for compatibility with previous StarBEAST2 versions
    public final Input<Double> windowInput = new Input<>("window", "size of the random walk window", 10.0);
    public final Input<RealParameter> originInput = new Input<RealParameter>("origin", "The time when the process started", (RealParameter) null);

    
    
    private enum RelativePosition {LEFT, RIGHT, BOTH}
    
    private int nextIndex;
    private int nodeCount;
    private int geneTreeCount;
    private int[][] leafNodeMaps;
    private RelativePosition[][] leafPositionArrays;
    private int trueBifurcationCount;
    private Node[] canonicalOrder;
    private int[] canonicalMap;
    private int[] trueBifurcations;
    private double[] nodeHeights;
    private Node[] leftChildren;
    private Node[] rightChildren;
    private Node[] parents;
    private boolean superimposedAncestors;
    private double maxHeight;
    private double window;
    private boolean originSpecified;

    @Override
    public void initAndValidate() {
        final SpeciesTree speciesTree = treeInput.get();
        nodeCount = speciesTree.getNodeCount();
        canonicalOrder = new Node[nodeCount];
        canonicalMap = new int[nodeCount];
        trueBifurcations = new int[nodeCount];
        nodeHeights = new double[nodeCount];
        leftChildren = new Node[nodeCount];
        rightChildren = new Node[nodeCount];
        parents = new Node[nodeCount];
        window = windowInput.get();
        originSpecified = originInput.get() != null;

        
        // Get the gene tree distributions
        // If a gene tree kernel is being used, then initialise the gene trees here to those in the kernel
        super.initAndValidate();
        
        // Assumption: if there is a gene tree kernel, all trees in the kernel have the same taxonset and tip number map
        // Therefore only need to cache these maps once and not every time
        this.calculateMaps();
        
    }
    
    
    
    private void calculateMaps() {
    	
    	
    	geneTreeCount = geneTreeDistributions.size();
        leafNodeMaps = new int[geneTreeCount][];
        leafPositionArrays = new RelativePosition[geneTreeCount][];
        for (int i = 0; i < geneTreeCount; i++) {
            leafNodeMaps[i] = geneTreeDistributions.get(i).getTipNumberMap();
            leafPositionArrays[i] = new RelativePosition[leafNodeMaps[i].length];
        }
    	
    }
    

    /* This proposal improves TREE SLIDE, developed by Joseph Heled. See section 3.4.1 of Heled's 2011 PhD thesis
    "Bayesian Computational Inference of Species Trees and Population Sizes". TREE SLIDE was developed for ultrametric
    binary species trees, this proposal has been made compatible with sampled ancestors by enforcing a minimum height
    and disallowing superimposed sampled ancestor nodes. Also uses a random walk window with reflection in order to
    sample the heights of nodes without maximum height constraints. */
    @Override
    public double proposal() {
    	
    	geneTreeDistributions = this.getTreeDistributions(this);
    	//geneTreeCount = geneTreeDistributions.size();
    	calculateMaps();
    	
    	
    	
        final SpeciesTree tree = treeInput.get();
        final Node originalRoot = tree.getRoot();

        // chooseCanonicalOrder also fills in nodeHeights and trueBifurcations
        // the lastIndex will be the last and right-most node index
        trueBifurcationCount = 0;
        nextIndex = 0;
        chooseCanonicalOrder(originalRoot);

        // no nodes can be changed by this operator
        if (trueBifurcationCount == 0) {
            return Double.NEGATIVE_INFINITY;
        }

        // pick a bifurcation at random and change the height
        final int chosenNode = trueBifurcations[Randomizer.nextInt(trueBifurcationCount)];
        final double originalHeight = nodeHeights[chosenNode];

        // As long as the height is above the tips (or sampled ancestors) immediately either side in the canonical
        // order, the species tree seems buildable from the new times. The exception is fake bifurcation nodes of equal
        // height, which can result in the superimposition of those sampled ancestors.
        final double minHeight = Double.max(nodeHeights[chosenNode - 1], nodeHeights[chosenNode + 1]);

        recalculateMaxHeight(chosenNode);

        // Use reflection to avoid invalid heights. Height returns to original position every 2 * (max - min) units,
        // so modulus is used to avoid unnecessary looping if the difference between window size and the tree scale
        // is extreme.
        final double heightDelta = (window * (Randomizer.nextDouble() - 0.5)) % (2.0 * (maxHeight - minHeight));
        double newHeight = originalHeight + heightDelta;
        while (newHeight < minHeight || newHeight > maxHeight) {
            if (newHeight < minHeight) {
                newHeight = minHeight + minHeight - newHeight;
            }
            if (newHeight > maxHeight) {
                newHeight = maxHeight + maxHeight - newHeight;
            }
        }

        nodeHeights[chosenNode] = newHeight;

        superimposedAncestors = false;
        final int rootIndex = rebuildTree(0, nodeCount - 1);
        parents[rootIndex] = null;
        if (superimposedAncestors) {
            return Double.NEGATIVE_INFINITY;
        }

        // wait until after checking for superimposed ancestors before modifying tree
        canonicalOrder[chosenNode].setHeight(newHeight);

        for (int i = 0; i < nodeCount; i++) {
            canonicalOrder[i].setParent(parents[i]);

            if (i % 2 == 1) { // internal node
                canonicalOrder[i].setLeft(leftChildren[i]);
                canonicalOrder[i].setRight(rightChildren[i]);
            }
        }

        // for some reason if the root is not reset - even if the root node is the same node as before! - the
        // morphological likelihood will be radically wrong (idk why)
        final Node newRoot = canonicalOrder[rootIndex];
        tree.setRoot(newRoot);

        assert checkVisitedCounts(tree);

        return 0.0;
    }

    private void recalculateMaxHeight(final int centerIndex) {
        if (originSpecified) {
            maxHeight = originInput.get().getValue();
        } else {
            maxHeight = Double.POSITIVE_INFINITY;
        }

        for (int i = 0; i < geneTreeCount; i++) {
            final int[] leafNodeMap = leafNodeMaps[i];
            final RelativePosition[] leafPositions = leafPositionArrays[i];
            for (int j = 0; j < leafNodeMap.length; j++) {
                final int speciesNodeNumber = leafNodeMap[j];
                final int speciesIndex = canonicalMap[speciesNodeNumber];
                if (speciesIndex < centerIndex) {
                    leafPositions[j] = RelativePosition.LEFT;
                } else {
                    leafPositions[j] = RelativePosition.RIGHT;
                }
            }

            final Node geneTreeRoot = geneTreeDistributions.get(i).getGeneTree().getRoot();
            recurseMaxHeight(geneTreeRoot, leafPositions);
        }
    }

    private RelativePosition recurseMaxHeight(final Node node, final RelativePosition[] leafPositions) {
        final Node leftChild = node.getLeft();
        final Node rightChild = node.getRight();

        RelativePosition leftDescendantPosition;
        if (leftChild.isLeaf()) {
            leftDescendantPosition = leafPositions[leftChild.getNr()];
        } else {
            leftDescendantPosition = recurseMaxHeight(leftChild, leafPositions);
        }

        RelativePosition rightDescendantPosition;
        if (rightChild.isLeaf()) {
            rightDescendantPosition = leafPositions[rightChild.getNr()];
        } else {
            rightDescendantPosition = recurseMaxHeight(rightChild, leafPositions);
        }

        if (leftDescendantPosition == rightDescendantPosition) {
            return leftDescendantPosition;
        } else {
            // if all descendants of one child are on the left, and all descendants of the other child are on the right
            if (leftDescendantPosition != RelativePosition.BOTH && rightDescendantPosition != RelativePosition.BOTH) {
                maxHeight = Double.min(maxHeight, node.getHeight());
            }
            return RelativePosition.BOTH;
        }
    }

    /* Performs an in-order traversal of the species tree, randomly shuffling left and right nodes, to produce
       a canonical order in the sense of Mau et al 1999. Also identify which nodes are true bifurcations
       (not fake nodes used for sampled ancestors) */
    private double chooseCanonicalOrder(final Node node) {
        Node canonicalLeft;
        Node canonicalRight;

        if (Randomizer.nextBoolean()) {
            canonicalLeft = node.getLeft();
            canonicalRight = node.getRight();
        } else {
            canonicalLeft = node.getRight();
            canonicalRight = node.getLeft();
        }

        double leftChildHeight;
        if (canonicalLeft.isLeaf()) {
            final int leftChildIndex = nextIndex;
            nextIndex++;

            canonicalMap[canonicalLeft.getNr()] = leftChildIndex;
            canonicalOrder[leftChildIndex] = canonicalLeft;

            leftChildHeight = canonicalLeft.getHeight();
            nodeHeights[leftChildIndex] = leftChildHeight;
        } else {
            leftChildHeight = chooseCanonicalOrder(canonicalLeft);
        }

        final int thisIndex = nextIndex;
        nextIndex++;

        canonicalMap[node.getNr()] = thisIndex;
        canonicalOrder[thisIndex] = node;

        final double thisHeight = node.getHeight();
        nodeHeights[thisIndex] = thisHeight;

        double rightChildHeight;
        if (canonicalRight.isLeaf()) {
            final int rightChildIndex = nextIndex;
            nextIndex++;

            canonicalMap[canonicalRight.getNr()] = rightChildIndex;
            canonicalOrder[rightChildIndex] = canonicalRight;

            rightChildHeight = canonicalRight.getHeight();
            nodeHeights[rightChildIndex] = rightChildHeight;
        } else {
            rightChildHeight = chooseCanonicalOrder(canonicalRight);
        }

        if (thisHeight > leftChildHeight && thisHeight > rightChildHeight) {
            trueBifurcations[trueBifurcationCount] = thisIndex;
            trueBifurcationCount++;
        }

        return thisHeight;
    }

    /* from and to are inclusive */
    private int rebuildTree(final int from, final int to) {
        double thisHeight = 0.0;
        int nodeIndex = -1;

        /* Only check internal nodes, which are odd numbered (leaves are even numbered). If there are multiple highest
           internal nodes in the range, they are likely fake bifurcations, and connecting
           them will result in multiple sampled ancestors at the same point in time along the same lineage.
           In this case we reject the move.
           This is similar to the following behaviour of LeafToSampledAncestorJump (see lines 68-70):
           if (getOtherChild(parent, leaf).getHeight() >= leaf.getHeight()) return Double.NEGATIVE_INFINITY; */
        for (int i = from + 1; i < to; i = i + 2) {
            if (nodeHeights[i] > thisHeight) {
                thisHeight = nodeHeights[i];
                nodeIndex = i;
            } else if (nodeHeights[i] == thisHeight) {
                superimposedAncestors = true;
            }
        }

        int leftNodeIndex;
        if (from == nodeIndex - 1) {
            leftNodeIndex = from;
        } else {
            leftNodeIndex = rebuildTree(from, nodeIndex - 1);
        }

        parents[leftNodeIndex] = canonicalOrder[nodeIndex];
        leftChildren[nodeIndex] = canonicalOrder[leftNodeIndex];

        int rightNodeIndex;
        if (nodeIndex + 1 == to) {
            rightNodeIndex = to;
        } else {
            rightNodeIndex = rebuildTree(nodeIndex + 1, to);
        }

        parents[rightNodeIndex] = canonicalOrder[nodeIndex];
        rightChildren[nodeIndex] = canonicalOrder[rightNodeIndex];

        return nodeIndex;
    }

    // for debugging, only called when assertions are enabled
    private boolean checkVisitedCounts(SpeciesTree tree) {
        int[] visitedCounts = new int[nodeCount];
        recurseVisitedCounts(tree.getRoot(), visitedCounts);
        for (int i = 0; i < nodeCount; i++) {
            if (visitedCounts[i] != 1) {
                return false;
            }
        }
        return true;
    }

    // for debugging, only called when assertions are enabled
    private void recurseVisitedCounts(Node node, int[] visitedCounts) {
        visitedCounts[node.getNr()]++;
        final List<Node> children = node.getChildren();
        if (!node.isLeaf()) {
            assert children.size() == 2;
            final Node leftChild = children.get(0);
            final Node rightChild = children.get(1);
            assert leftChild.getParent() == node;
            assert rightChild.getParent() == node;
            assert leftChild.getHeight() <= node.getHeight();
            assert rightChild.getHeight() <= node.getHeight();
            recurseVisitedCounts(leftChild, visitedCounts);
            recurseVisitedCounts(rightChild, visitedCounts);
        }
    }
}
