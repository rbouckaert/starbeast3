/*
 * FixedHeightSubtreePruneRegraftOperator.java
 *
 * Copyright (c) 2002-2015 Alexei Drummond, Andrew Rambaut and Marc Suchard
 *
 * This file is part of BEAST.
 * See the NOTICE file distributed with this work for additional
 * information regarding copyright ownership and licensing.
 *
 * BEAST is free software; you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as
 * published by the Free Software Foundation; either version 2
 * of the License, or (at your option) any later version.
 *
 *  BEAST is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with BEAST; if not, write to the
 * Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
 * Boston, MA  02110-1301  USA
 */

package starbeast3.operators;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.evolution.operator.TreeOperator;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.util.Randomizer;

/**
 * Implements the fixed height subtree prune regraft move. Described by Sebastian Hoehna et al
 *
 * @author Andrew Rambaut
 * @version $Id$
 */
@Description("Implements the fixed height subtree prune regraft move. Described by Sebastian Hoehna et al")
public class FixedHeightSubtreePruneRegraftOperator extends TreeOperator {

    private TreeInterface tree;


    /**
     * Constructor
     * @param tree
     * @param weight
     */
    public FixedHeightSubtreePruneRegraftOperator() {
    }
    
    public FixedHeightSubtreePruneRegraftOperator(Tree tree, double weight) {
    	initByName("tree", tree, "weight", weight);
    }
    
    @Override
    public void initAndValidate() {
    	this.tree = treeInput.get();
    	
    }
    
    /**
     * Do a subtree jump move.
     *
     * @return the log-transformed hastings ratio
     */
    @Override
    public double proposal() {

        final Node root = tree.getRoot();

//        double  maxHeight = tree.getNodeHeight(root);

        Node i;
        Node iP = null;
        Node CiP = null;
        Node PiP = null;
        List<Node> destinations = null;

        do {
            // 1. choose a random node avoiding root or child of root
            i = tree.getNode(Randomizer.nextInt(tree.getNodeCount()));

        } while (root == i || i.getParent() == root);

        iP = i.getParent();
        CiP = iP.getLeft() == i ? iP.getRight() : iP.getLeft();
        PiP = iP.getParent();

        // get the height of the parent
        double parentHeight = iP.getHeight();

        // get a list of all edges that intersect this height
        destinations = getIntersectingEdges(tree, parentHeight);

        if (destinations.size() == 0) {
            // if there are no destinations available then reject the move
            return Double.NEGATIVE_INFINITY;
        }

        int r = Randomizer.nextInt(destinations.size());

        // remove the target node and its sibling (shouldn't be there because their parent's height is exactly equal to the target height).
        destinations.remove(i);
        destinations.remove(CiP);

        final Node j = destinations.get(r);
        final Node jP = j.getParent();

        // remove the parent of i by connecting its sibling to its grandparent.
        iP.removeChild( CiP);
        PiP.removeChild( iP);
        PiP.addChild( CiP);

        // remove destination edge j from its parent
        jP.removeChild( j);

        // add destination edge to the parent of i
        iP.addChild( j);

        // and add the parent of i as a child of the former parent of j
        jP.addChild( iP);


        return 0.0;
    }

    /**
     * Gets a list of edges that extend the given height
     * @param tree
     * @param height
     * @return
     */
    private List<Node> getIntersectingEdges(TreeInterface tree, double height) {

        List<Node> intersectingEdges = new ArrayList<Node>();

        for (int i = 0; i < tree.getNodeCount(); i++) {
            final Node node = tree.getNode(i);
            final Node parent = node.getParent();;

            // The original node and its sibling will not be included because their height is exactly equal to the target height
            if (parent != null && node.getHeight() < height && parent.getHeight() > height) {
                intersectingEdges.add(node);
            }
        }
        return intersectingEdges;
    }

    public String getPerformanceSuggestion() {
        return null;
    }


    public String getOperatorName() {
        return this.getClass().getSimpleName() + "(" + tree.getID() + ")";
    }
}
