package starbeast3.util;

import java.util.Arrays;

import beast.evolution.tree.Node;

public class TreeTopologySerialiser {
	
	private static int[] dummy = new int[1];

	// Returns a unique identifier of the tree topology
	// The identifier will not describe branch lengths or other meta data, topology only.
	public static synchronized String serialiseNode(Node node) {
		TreeTopologySerialiser.dummy[0] = 1;
		return serialiseNode(node, dummy);
	}
	
    private static synchronized String serialiseNode(Node node, int[] maxNodeInClade) {
        StringBuilder buf = new StringBuilder();

        if (!node.isLeaf()) {

            if (node.getChildCount() <= 2) {
                // Computationally cheap method for special case of <=2 children

                buf.append("(");
                String child1 = TreeTopologySerialiser.serialiseNode(node.getChild(0), maxNodeInClade);
                int child1Index = maxNodeInClade[0];
                if (node.getChildCount() > 1) {
                    String child2 = TreeTopologySerialiser.serialiseNode(node.getChild(1), maxNodeInClade);
                    int child2Index = maxNodeInClade[0];
                    if (child1Index > child2Index) {
                        buf.append(child2);
                        buf.append(",");
                        buf.append(child1);
                    } else {
                        buf.append(child1);
                        buf.append(",");
                        buf.append(child2);
                        maxNodeInClade[0] = child1Index;
                    }
                } else {
                    buf.append(child1);
                }
                buf.append(")");
                if (node.getID() != null) {
                    buf.append(node.getNr()+1);
                }

            } else {
                // General method for >2 children

                String[] childStrings = new String[node.getChildCount()];
                int[] maxNodeNrs = new int[node.getChildCount()];
                Integer[] indices = new Integer[node.getChildCount()];
                for (int i = 0; i < node.getChildCount(); i++) {
                    childStrings[i] = TreeTopologySerialiser.serialiseNode(node.getChild(i), maxNodeInClade);
                    maxNodeNrs[i] = maxNodeInClade[0];
                    indices[i] = i;
                }

                Arrays.sort(indices, (i1, i2) -> {
                    if (maxNodeNrs[i1] < maxNodeNrs[i2])
                        return -1;

                    if (maxNodeNrs[i1] > maxNodeNrs[i2])
                        return 1;

                    return 0;
                });

                maxNodeInClade[0] = maxNodeNrs[maxNodeNrs.length - 1];

                buf.append("(");
                for (int i = 0; i < indices.length; i++) {
                    if (i > 0)
                        buf.append(",");

                    buf.append(childStrings[indices[i]]);
                }

                buf.append(")");

                if (node.getID() != null) {
                    buf.append(node.getNr() + 1);
                }
            }

        } else {
            maxNodeInClade[0] = node.getNr();
            buf.append(node.getNr() + 1);
        }

        return buf.toString();
        
    }


}
