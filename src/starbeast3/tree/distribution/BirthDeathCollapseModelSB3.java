

package starbeast3.tree.distribution;

import java.util.ArrayList;
import java.util.List;
import java.util.Random;

import beast.core.State;
import beast.core.StateNode;
import beast.core.parameter.RealParameter;
import beast.evolution.speciation.YuleModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
//import stacey.BirthDeathCollapseModel;

public class BirthDeathCollapseModelSB3 { // extends BirthDeathCollapseModel {

	/*
	@Override
	public void sample(State state, Random random) {
		
		
		// Assumes Yule model
		if (relativeDeathRate.get().getValue() > 0) {
			throw new UnsupportedOperationException("Cannot directly sample from a birth-death model. Please set relativeDeathRate to 0.");
		}
		


        if (sampledFlag) return;
        sampledFlag = true;
        
        

        // Cause conditional parameters to be sampled
        sampleConditions(state, random);

        Tree tree = (Tree) treeInput.get();
        RealParameter birthRate = birthDiffRate.get();
        double w = collapseWeight.get().getValue();

        // Simulate tree conditional on new parameters
        List<Node> activeLineages = new ArrayList<>();
        for (Node oldLeaf : tree.getExternalNodes()) {
            Node newLeaf = new Node(oldLeaf.getID());
            newLeaf.setNr(oldLeaf.getNr());
            newLeaf.setHeight(0.0);
            activeLineages.add(newLeaf);
        }

        int nextNr = activeLineages.size();

        double t = 0.0;
        while (activeLineages.size() > 1) {
            int k = activeLineages.size();
            
        	double u = Randomizer.nextDouble();
    		double t_;
    		
    		
    		
    		// Sample from spike (uniformly between 0 and collapseHeight)
    		if (u < w) {
    			t_ = Randomizer.nextDouble() * collapseHeight.get();
    		}
    		
    		// Sample from Yule (exponential distribution)
    		else {
    			double a = birthRate.getValue() * k;
    			t_ = -Math.log(random.nextDouble())/a;
    		}
            


            t += t_;

            Node node1 = activeLineages.get(random.nextInt(k));
            Node node2;
            do {
                node2 = activeLineages.get(random.nextInt(k));
            } while (node2.equals(node1));

            Node newParent = new Node();
            newParent.setNr(nextNr++);
            newParent.setHeight(t);
            newParent.addChild(node1);
            newParent.addChild(node2);

            activeLineages.remove(node1);
            activeLineages.remove(node2);
            activeLineages.add(newParent);
        }

        tree.assignFromWithoutID(new Tree(activeLineages.get(0)));
		
		
		
		
		
		
		
	}
	
	*/
}
