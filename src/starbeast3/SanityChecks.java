
/**
 * @author Huw A. Ogilvie
 * 
 */


package starbeast3;

import java.util.List;

import beast.evolution.tree.Node;

final class SanityChecks {
    protected static boolean checkTreeSanity(Node node) {
        final double nodeHeight = node.getHeight();
        final List<Node> children = node.getChildren();
        assert children.size() == 2;

        for (Node childNode: children) {
            final double childHeight = childNode.getHeight();
            assert childNode.getParent() == node;

            // direct ancestor branches have zero height
            // so equal height allowed in this case
            assert childHeight <= nodeHeight;

            if (!childNode.isLeaf()) checkTreeSanity(childNode);
        }

        return true;
    }
}