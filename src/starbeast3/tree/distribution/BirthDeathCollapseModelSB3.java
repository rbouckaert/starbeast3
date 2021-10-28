

package starbeast3.tree.distribution;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Random;

import beast.core.State;
import beast.core.StateNode;
import beast.core.parameter.RealParameter;
import beast.evolution.speciation.YuleModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import stacey.BirthDeathCollapseModel;


public class BirthDeathCollapseModelSB3  extends BirthDeathCollapseModel {

	
	
	
    @Override
    public List<String> getConditions() {
    	List<String> conditions = new ArrayList<>();
    	conditions.add(treeInput.get().getID());
    	conditions.add(collapseWeight.get().getID());
    	conditions.add(birthDiffRate.get().getID());
    	conditions.add(relativeDeathRate.get().getID());
    	return conditions;
    }
	
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
        double epsilon = collapseHeight.get();

        // Simulate tree conditional on new parameters
        List<Node> activeLineages = new ArrayList<>();
        for (Node oldLeaf : tree.getExternalNodes()) {
            Node newLeaf = new Node(oldLeaf.getID());
            newLeaf.setNr(oldLeaf.getNr());
            newLeaf.setHeight(0.0);
            activeLineages.add(newLeaf);

        }
        
        
        
        // How many of the internal nodes will be below epsilon? Binomial(n-1, w) distribution. Assuming binary tree
        List<Double> collapseHeights = new ArrayList<>();
        int n = activeLineages.size();
        for (int i = 0; i < n-1; i ++) {
        	if (random.nextDouble() < w) {
        	
        		// Sample a collapse height
        		double h = random.nextDouble() * epsilon;
        		collapseHeights.add(h);
        		
        	}
        	
        }
	
	
        int nextNr = activeLineages.size();
	
		// Create collapse epoch
		Collections.sort(collapseHeights);
		for (int i = 0; i < collapseHeights.size(); i ++) {
			
			
			double t = collapseHeights.get(i);
			int k = activeLineages.size();
			
			// Sample 2 nodes
			Node node1 = activeLineages.get(random.nextInt(k));
            Node node2;
            do {
                node2 = activeLineages.get(random.nextInt(k));
            } while (node2.equals(node1));
            
            
            
            // Join them
            Node newParent = new Node();
            newParent.setNr(nextNr++);
            newParent.setHeight(t);
            newParent.addChild(node1);
            newParent.addChild(node2);

            activeLineages.remove(node1);
            activeLineages.remove(node2);
            activeLineages.add(newParent);
  
			
		}

        

        
		// Sample from Yule beginning at time epsilon
        double t = epsilon;
        while (activeLineages.size() > 1) {
            int k = activeLineages.size();
            
    
    		
  
    		// Sample from Yule (exponential distribution)
			double a = birthRate.getValue() * k;
			double t_ = -Math.log(random.nextDouble())/a;
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
	
	
}
