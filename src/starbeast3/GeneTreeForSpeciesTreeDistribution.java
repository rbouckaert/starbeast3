package starbeast3;



import java.util.Arrays;
import java.util.List;
import java.util.PriorityQueue;
import java.util.Random;

import beast.app.beauti.Beauti;
import beast.core.Description;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.State;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import beast.evolution.speciation.*;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;



@Description("Calculates probability of gene tree conditioned on a species tree (multi-species coalescent)"
		+ "assuming a constant population size for each branch")
public class GeneTreeForSpeciesTreeDistribution extends TreeDistribution {

    final public Input<TreeInterface> speciesTreeInput =
            new Input<>("speciesTree", "species tree containing the associated gene tree", Validate.REQUIRED);

    final public Input<Double> ploidyInput =
            new Input<>("ploidy", "ploidy (copy number) for this gene, typically a whole number or half (default 2 for autosomal_nuclear)", 2.0);

    
    final public Input<SpeciesTreePrior> speciesTreePriorInput =
            new Input<>("speciesTreePrior", "defines population function and its parameters", Validate.REQUIRED);

    // intervals for each of the species tree branches
    private PriorityQueue<Double>[] intervalsInput;
    // count nr of lineages at the bottom of species tree branches
    private int[] nrOfLineages;
    // maps gene tree leaf nodes to species tree leaf nodes. Indexed by node number.
    private int[] nrOfLineageToSpeciesMap;
    private RealParameter popSizesBottom;

    // Ploidy is a constant - cache value of input here
    private double ploidy;
    double [][] alltimes;
    boolean uptodate = false;
    
    public double getPloidy() {return ploidy;}

    public GeneTreeForSpeciesTreeDistribution() {
        treeInput.setRule(Validate.REQUIRED);
    }

	@Override
    public void initAndValidate() {
    	ploidy = ploidyInput.get();

        final Node[] gtNodes = treeInput.get().getNodesAsArray();
        final int gtLineages = treeInput.get().getLeafNodeCount();
        final Node[] sptNodes = speciesTreeInput.get().getNodesAsArray();
        final int speciesCount = speciesTreeInput.get().getNodeCount();


        if (Beauti.isInBeauti()) {
            // we are in BEAUti, so do not initialise
            return;
        }

        // reserve memory for priority queues
        intervalsInput = new PriorityQueue[speciesCount];
        for (int i = 0; i < speciesCount; i++) {
            intervalsInput[i] = new PriorityQueue<>();
        }

        // sanity check lineage nodes are all at height=0
        for (int i = 0; i < gtLineages; i++) {
            if (gtNodes[i].getHeight() != 0) {
                throw new IllegalArgumentException("Cannot deal with taxon " + gtNodes[i].getID() +
                        ", which has non-zero height + " + gtNodes[i].getHeight());
            }
        }
        
        // set up nrOfLineageToSpeciesMap
        nrOfLineageToSpeciesMap = new int[gtLineages];
        Arrays.fill(nrOfLineageToSpeciesMap, -1);
        for (int i = 0; i < gtLineages; i++) {
            final String speciesID = getSpeciesID(gtNodes[i].getID());
            // ??? can this be a startup check? can this happen during run due to tree change?
            if (speciesID == null) {
                throw new IllegalArgumentException("Cannot find species for lineage " + gtNodes[i].getID());
            }
            for (int species = 0; species < speciesCount; species++) {
                if (speciesID.equals(sptNodes[species].getID())) {
                    nrOfLineageToSpeciesMap[i] = species;
                    break;
                }
            }
            if (nrOfLineageToSpeciesMap[i] < 0) {
                throw new IllegalArgumentException("Cannot find species with name " + speciesID + " in species tree");
            }
        }

        // calculate nr of lineages per species
        nrOfLineages = new int[speciesCount];

        final SpeciesTreePrior popInfo = speciesTreePriorInput.get();
        popSizesBottom = popInfo.popSizesBottomInput.get();
        
        alltimes = new double[sptNodes.length][];
        for (int i = 0; i < sptNodes.length; i++) {
        	alltimes[i] = new double[0];
        }
    }

    /**
     * @param lineageID
     * @return species ID to which the lineage ID belongs according to the TaxonSets
     */
    private String getSpeciesID(final String lineageID) {
        final TaxonSet taxonSuperset = speciesTreePriorInput.get().taxonSetInput.get();
        final List<Taxon> taxonSets = taxonSuperset.taxonsetInput.get();
        for (final Taxon taxonSet : taxonSets) {
            final List<Taxon> taxa = ((TaxonSet) taxonSet).taxonsetInput.get();
            for (final Taxon aTaxa : taxa) {
                if (aTaxa.getID().equals(lineageID)) {
                    return taxonSet.getID();
                }
            }
        }
        return null;
    }

    @Override
    public double calculateLogP() {
        logP = 0;
        for (final PriorityQueue<Double> m_interval : intervalsInput) {
            m_interval.clear();
        }

        Arrays.fill(nrOfLineages, 0);

        final TreeInterface stree = speciesTreeInput.get();
        final Node[] speciesNodes = stree.getNodesAsArray();

        traverseLineageTree(speciesNodes, treeInput.get().getRoot());

//        System.err.println(getID());
//		for (int i = 0; i < m_intervals.length; i++) {
//			System.err.println(m_intervals[i]);
//		}

        // if the gene tree does not fit the species tree, logP = -infinity by now
        if (logP == 0) {
            traverseSpeciesTree(stree.getRoot());
        }
        uptodate = true;
        return logP;
    }

    /**
     * calculate contribution to logP for each of the branches of the species tree
     *
     * @param node*
     */
    private void traverseSpeciesTree(final Node node) {
        if (!node.isLeaf()) {
            traverseSpeciesTree(node.getLeft());
            traverseSpeciesTree(node.getRight());
        }
        // calculate contribution of a branch in the species tree to the log probability
        final int nodeIndex = node.getNr();

        // k = nr of intervals, as defined in the *BEAST paper
        final int k = intervalsInput[nodeIndex].size();
        if (alltimes[nodeIndex].length != k + 2) {
        	alltimes[nodeIndex] = new double[k + 2];
        }
        final double[] times = alltimes[nodeIndex];
        times[0] = node.getHeight();
        for (int i = 1; i <= k; i++) {
            times[i] = intervalsInput[nodeIndex].poll();
        }
        if (!node.isRoot()) {
        	// time at top of the branch
            times[k + 1] = node.getParent().getHeight();
        } else {
        	// time at root of gene tree
        	// RRB: why consider node.getHeight(), which is already what times[0] was set to,
        	// which means treeInput.get().getRoot().getHeight() is always > node.getHeight()?
            times[k + 1] = Math.max(node.getHeight(), treeInput.get().getRoot().getHeight());
        }
        // sanity check
        for (int i = 0; i <= k; i++) {
            if (times[i] > times[i + 1]) {
            	Log.warning.println("invalid times");
                calculateLogP();
            }
        }

        final int lineagesBottom = nrOfLineages[nodeIndex];

        calcConstantPopSizeContribution(lineagesBottom, popSizesBottom.getValue(nodeIndex), times, k);
    }

    /* the contribution of a branch in the species tree to
      * the log probability, for constant population function.
      */
    private void calcConstantPopSizeContribution(final int lineagesBottom, final double popSize2,
                                                 final double[] times, final int k) {
        final double popSize = popSize2 * ploidy;
        logP += -k * Math.log(popSize);
        for (int i = 0; i <= k; i++) {
            logP += -((lineagesBottom - i) * (lineagesBottom - i - 1.0) / 2.0) * (times[i + 1] - times[i]) / popSize;
        }
    }



    /**
     * collect intervals for each of the branches of the species tree
     * as defined by the lineage tree.
     *
     * @param speciesNodes
     * @param node
     * @return
     */
    private int traverseLineageTree(final Node[] speciesNodes, final Node node) {
        if (node.isLeaf()) {
            final int species = nrOfLineageToSpeciesMap[node.getNr()];
            nrOfLineages[species]++;
            return species;
        } else {
            int speciesLeft = traverseLineageTree(speciesNodes, node.getLeft());
            int speciesRight = traverseLineageTree(speciesNodes, node.getRight());
            final double height = node.getHeight();

            while (!speciesNodes[speciesLeft].isRoot() && height > speciesNodes[speciesLeft].getParent().getHeight()) {
                speciesLeft = speciesNodes[speciesLeft].getParent().getNr();
                nrOfLineages[speciesLeft]++;
            }
            while (!speciesNodes[speciesRight].isRoot() && height > speciesNodes[speciesRight].getParent().getHeight()) {
                speciesRight = speciesNodes[speciesRight].getParent().getNr();
                nrOfLineages[speciesRight]++;
            }
            // validity check
            if (speciesLeft != speciesRight) {
                // if we got here, it means the gene tree does
                // not fit in the species tree
                logP = Double.NEGATIVE_INFINITY;
            }
            intervalsInput[speciesRight].add(height);
            return speciesRight;
        }
    }

    @Override
    public void restore() {
    	uptodate = false;
    	super.restore();
    }
    

    @Override
    public boolean requiresRecalculation() {
        return true;
    }

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(final State state, final Random random) {
    }

    // number of coalescent events in branch i
	public int getCoalescentCount(int i) {
		if (!uptodate) {
			calculateLogP();			
		}
		return alltimes[i].length - 2;
	}
	
	public double [] getTimes(int i) {
		if (!uptodate) {
			calculateLogP();			
		}
		return alltimes[i];
	}
}
