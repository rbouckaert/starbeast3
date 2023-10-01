package starbeast3.operators;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.util.Randomizer;
import sa.evolution.operators.SampledNodeDateRandomWalker;
import sa.evolution.tree.SamplingDate;
import sa.evolution.tree.TreeWOffset;
import starbeast3.evolution.speciation.GeneTreeForSpeciesTreeDistribution;


@Description("Moves a tip in the species tree, and also ensures the gene leaves remain at the bottom of the tip")
public class SampledNodeDateRandomWalkerSB3 extends SampledNodeDateRandomWalker {
	public final Input<List<GeneTreeForSpeciesTreeDistribution>> geneTreesInput = new Input<>("gene", "list of gene trees that constrain species tree movement", new ArrayList<>());
	
	boolean useNodeNumbers;
    List<String> samplingDateTaxonNames = new ArrayList<>();
    TreeWOffset combinedTree;

    @Override
    public void initAndValidate() {
    	combinedTree = treeWOffsetInput.get();
        if(combinedTree == null) {
        	combinedTree = new TreeWOffset();
        	combinedTree.setInputValue("tree", treeInput.get());
        	combinedTree.initAndValidate();
        }
    	
        windowSize = windowSizeInput.get();
        useGaussian = useGaussianInput.get();

        for (SamplingDate taxon:samplingDatesInput.get()) {
            samplingDateTaxonNames.add(taxon.taxonInput.get().getID());
        }

        // determine taxon set to choose from
        if (m_taxonsetInput.get() != null) {
            useNodeNumbers = false;
            List<String> sTaxaNames = new ArrayList<String>();
            for (String sTaxon : treeInput.get().getTaxaNames()) {
            	Log.warning("Found taxon " + sTaxon);
                sTaxaNames.add(sTaxon);
            }

            List<String> set = m_taxonsetInput.get().asStringList();
            int nNrOfTaxa = set.size();
            taxonIndices = new int[nNrOfTaxa];
            int k = 0;
            for (String sTaxon : set) {
                int iTaxon = sTaxaNames.indexOf(sTaxon);
                if (iTaxon < 0) {
                    throw new IllegalArgumentException("Cannot find taxon " + sTaxon + " in tree " + treeInput.get().getID());
                }
                taxonIndices[k++] = iTaxon;
            }
        } else {
            useNodeNumbers = true;
        }
    }
    

    @Override
    public double proposal() {

        // randomly select a leaf node
        Tree tree = combinedTree.getTree();

        Node node;
        if (useNodeNumbers) {
            int leafNodeCount = tree.getLeafNodeCount();
            int i = Randomizer.nextInt(leafNodeCount);
            node = tree.getNode(i);
        }  else {
            int i = Randomizer.nextInt(taxonIndices.length);
            node = tree.getNode(taxonIndices[i]);
        }

        double value = combinedTree.getHeightOfNode(node.getNr());

        if (value == 0.0) {
            return Double.NEGATIVE_INFINITY;
        }
        double newValue = value;

        boolean drawFromDistribution = samplingDateTaxonNames.contains(node.getID());
        if (drawFromDistribution) {
            SamplingDate taxonSamplingDate = samplingDatesInput.get().get(samplingDateTaxonNames.indexOf(node.getID()));
            double range = taxonSamplingDate.getUpper() - taxonSamplingDate.getLower();
            newValue = taxonSamplingDate.getLower() + Randomizer.nextDouble() * range;
        }  else {
            if (useGaussian) {
                newValue += Randomizer.nextGaussian() * windowSize;
            } else {
                newValue += Randomizer.nextDouble() * 2 * windowSize - windowSize;
            }
        }


        Node fake = null;
        double lower, upper;

        if ((node).isDirectAncestor()) {
            fake = node.getParent();
            lower = combinedTree.getHeightOfNode(getOtherChild(fake, node).getNr());
            if (fake.getParent() != null) {
                upper = combinedTree.getHeightOfNode(fake.getParent().getNr());
            } else upper = Double.POSITIVE_INFINITY;
        } else {
            //lower = Double.NEGATIVE_INFINITY;
            lower = 0.0;
            upper = combinedTree.getHeightOfNode(node.getParent().getNr());
        }

        if (newValue < lower || newValue > upper) {
            return Double.NEGATIVE_INFINITY;
        }

        if (newValue == value) {
            // this saves calculating the posterior
            return Double.NEGATIVE_INFINITY;
        }

        
        
        // Move all gene tree nodes that were mapped to this node
        double dy = newValue - value;
        for (GeneTreeForSpeciesTreeDistribution gene : geneTreesInput.get()) {
        	
        	Tree geneTree = (Tree) gene.getGeneTree();
        	for (Node n : geneTree.getNodesAsArray()) {
        		if (gene.mapGeneNodeToSpeciesNode(n.getNr()) == node) {
        			
        			// Need to move this node
        			Log.warning("Need to adjust node by " + dy);
        			n.setHeight(n.getHeight() + dy);
        			
        		}
        	}
        	
        	
        	// Double check for negative branch lengths
        	for (Node n : geneTree.getNodesAsArray()) {
        		if (n.getLength() < 0) {
        			Log.warning("Unexpected tip date adjustemnt: negative branch length");
        			return Double.NEGATIVE_INFINITY;
        		}
        	}
        	
        }
        
        if (fake != null) {
        	combinedTree.setHeightOfNode(fake.getNr(), newValue);
        }
        combinedTree.setHeightOfNode(node.getNr(), newValue);
        
        return 0.0;
    }
	
	
    /**
     * Gets the list of gene trees which this operator may act on
     */
    public List<Tree> getTrees() {

		List<Tree> trees = new ArrayList<Tree>();
		for (GeneTreeForSpeciesTreeDistribution t : this.geneTreesInput.get() ) {
			trees.add((Tree) t.getGeneTree());
		}
		return trees;
    	
    	
    	
    }
    
    
    @Override
    public List<StateNode> listStateNodes() {
  
        // Pick up all inputs that are stateNodes that are estimated
        final List<StateNode> stateNodes = super.listStateNodes();
    	for (Tree tree : this.getTrees()) {
    		if (!stateNodes.contains(tree)) {
    			stateNodes.add(tree);
    		}
    	}
        	
        return stateNodes;
    }

}
