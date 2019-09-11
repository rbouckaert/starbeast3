package starbeast3;



import java.util.Arrays;
import java.util.List;
import java.util.Map;
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
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;



@Description("Calculates probability of gene tree conditioned on a species tree (multi-species coalescent)"
		+ "assuming a constant population size for each branch")
public class GeneTreeForSpeciesTreeDistribution extends TreeDistribution {

	//TreeInterfaceSB3
    final public Input<SpeciesTree> speciesTreeInput =
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
    boolean logPuptodate = false;
    boolean clockuptodate = false;
    
    
    // Node counts
    private int geneTreeLeafNodeCount;
    private int geneTreeNodeCount;
    private int speciesNodeCount;
    
    
    // The following are matrices associated with each branch of the species tree
    // They are flattened to arrays for optimal java performance
    protected double[] coalescentTimes; // the coalescent event times for this gene tree for all species tree branches
    protected double[] storedCoalescentTimes; // the coalescent event times for this gene tree for all species tree branches
    int coalescentTimesLength; // length of coalescentTimes array
    protected int[] coalescentCounts; // the number of coalescent events in each branch
    protected int[] storedCoalescentCounts; // stored version of coalescentCounts
    final static int DELTA_BLOCK_SIZE = 4;
    private int blocksize = DELTA_BLOCK_SIZE; // size of blocks for storing coalescentTimes, may grow (and shrink) throughout the MCMC
    int maxCoalescentCounts, storedMaxCoalescentCounts; // maximum number of coalescent events in a branch -- blocksize must always be at least as large

    
    
    protected int[] coalescentLineageCounts; // The number of lineages at the tipward end of each branch
    protected int[] storedCoalescentLineageCounts; // The number of lineages at the tipward end of each branch

    
    
    protected double[] speciesOccupancy; // A [1D] matrix describing the length of overlap between each gene tree node's 
    									 // branch with each species tree node's branch
    protected int[] geneNodeSpeciesAssignment;
    protected int[] storedGeneNodeSpeciesAssignment;
    protected double[] storedSpeciesOccupancy;
    
    
    // Pre-calculated lineage counts and node assignments for gene tree leaf nodes
    int[] leafCoalescentLineageCounts;
    int[] leafGeneNodeSpeciesAssignment;
    
    int updateCount = 0;
    boolean stopPopping = false;
    
    
    // Maps gene tree tip numbers to species tree tip number
    private int[] localTipNumberMap;

    
    
    public double getPloidy() {
    	return ploidy;
    }


    public GeneTreeForSpeciesTreeDistribution() {
        treeInput.setRule(Validate.REQUIRED);
    }
    
    

    

	@Override
    public void initAndValidate() {
    	ploidy = ploidyInput.get();
    	logP = 0.0;
    	clockuptodate = false;

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
        

        // Calculate node counts
        geneTreeNodeCount = treeInput.get().getNodeCount();
        geneTreeLeafNodeCount = treeInput.get().getLeafNodeCount();
        speciesNodeCount = speciesTreeInput.get().getNodeCount();
        
        geneNodeSpeciesAssignment = new int[geneTreeNodeCount];
        storedGeneNodeSpeciesAssignment = new int[geneTreeNodeCount];
        
        
      
        // Generate map of species tree tip node names to node numbers
        final Map<String, Integer> tipNumberMap = speciesTreeInput.get().getTipNumberMap();
        localTipNumberMap = new int[treeInput.get().getLeafNodeCount()];
        for (int i = 0; i < treeInput.get().getLeafNodeCount(); i++) {
        	final Node geneTreeLeafNode = treeInput.get().getNode(i);
        	final String geneTreeLeafName = geneTreeLeafNode.getID();
        	final int geneTreeLeafNumber = geneTreeLeafNode.getNr();

        	if (tipNumberMap.containsKey(geneTreeLeafName)) // not in BEAUTi
        	    localTipNumberMap[geneTreeLeafNumber] = tipNumberMap.get(geneTreeLeafName);
        }
        

        // Allocate memory for coalescent counts, lengths, and times
        coalescentLineageCounts = new int[speciesNodeCount];
        storedCoalescentLineageCounts = new int[speciesNodeCount];
        coalescentCounts = new int[speciesNodeCount];
        storedCoalescentCounts = new int[speciesNodeCount];
        coalescentTimesLength = speciesNodeCount * blocksize;
        coalescentTimes = new double[coalescentTimesLength + geneTreeNodeCount];
        storedCoalescentTimes = new double[coalescentTimesLength + geneTreeNodeCount];
        
        
        // Allocate memory for species occupancies
        speciesOccupancy = new double[geneTreeNodeCount * speciesNodeCount];
        storedSpeciesOccupancy = new double[geneTreeNodeCount * speciesNodeCount];
        

        
        leafCoalescentLineageCounts = new int[speciesNodeCount];
        leafGeneNodeSpeciesAssignment = new int[geneTreeNodeCount];
        Arrays.fill(leafGeneNodeSpeciesAssignment, -1);
        
        for (int geneTreeLeafNumber = 0; geneTreeLeafNumber < geneTreeLeafNodeCount; geneTreeLeafNumber++) {
            final Node geneTreeLeafNode = treeInput.get().getNode(geneTreeLeafNumber);
            final int speciesTreeLeafNumber = localTipNumberMap[geneTreeLeafNode.getNr()];
            leafCoalescentLineageCounts[speciesTreeLeafNumber]++;            
            leafGeneNodeSpeciesAssignment[geneTreeLeafNumber] = speciesTreeLeafNumber;
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
    	
    	
        assert SanityChecks.checkTreeSanity(speciesTreeInput.get().getRoot());
        
        
        logP = 0.0;
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

        // If the gene tree does not fit the species tree, logP = -infinity by now
        if (logP == 0) {
            traverseSpeciesTree(stree.getRoot());
        }
        logPuptodate = true;
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
    public void store() {
    	super.store();
    	
    	/*
      
        System.arraycopy(coalescentCounts, 0, storedCoalescentCounts, 0, coalescentCounts.length);
        System.arraycopy(coalescentTimes, 0, storedCoalescentTimes, 0, coalescentTimesLength);
        System.arraycopy(coalescentLineageCounts, 0, storedCoalescentLineageCounts, 0, coalescentLineageCounts.length);

        System.arraycopy(geneNodeSpeciesAssignment, 0, storedGeneNodeSpeciesAssignment, 0, geneNodeSpeciesAssignment.length);
        System.arraycopy(speciesOccupancy, 0, storedSpeciesOccupancy, 0, speciesOccupancy.length);

        storedMaxCoalescentCounts = maxCoalescentCounts;
        */
    }
    
  


    @Override
    public void restore() {
    	clockuptodate = false;
    	super.restore();
    	
    	/*
    	
    	double[] tmpCoalescentTimes = coalescentTimes;
    	int[] tmpCoalescentCounts = coalescentCounts;
    	int[] tmpCoalescentLineageCounts = coalescentLineageCounts;
    	int[] tmpGeneNodeSpeciesAssignment = geneNodeSpeciesAssignment;
    	double[] tmpSpeciesOccupancy = speciesOccupancy;

    	coalescentTimes = storedCoalescentTimes;
    	coalescentCounts = storedCoalescentCounts;
    	coalescentLineageCounts = storedCoalescentLineageCounts;
    	speciesOccupancy = storedSpeciesOccupancy;
    	geneNodeSpeciesAssignment = storedGeneNodeSpeciesAssignment;

    	storedCoalescentTimes = tmpCoalescentTimes;
    	storedCoalescentCounts = tmpCoalescentCounts;
    	storedCoalescentLineageCounts = tmpCoalescentLineageCounts;
    	storedSpeciesOccupancy = tmpSpeciesOccupancy;
    	storedGeneNodeSpeciesAssignment = tmpGeneNodeSpeciesAssignment;

    	maxCoalescentCounts = storedMaxCoalescentCounts;
		 */
    	
    	
    }
    

    @Override
    public boolean requiresRecalculation() {
    	clockuptodate = false;
    	logPuptodate = false;
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


    // number of lineages at the bottom of branch i
	public int getLineageCount(int i) {
		return nrOfLineages[i];
	}

    // number of coalescent events in branch i
	public int getCoalescentCount(int i) {
		if (!logPuptodate) {
			calculateLogP();			
		}
		return alltimes[i].length - 2;
	}
	
	public double [] getTimes(int i) {
		if (!logPuptodate) {
			calculateLogP();			
		}
		return alltimes[i];
	}

	public int getNodeCount() {
		return geneTreeNodeCount;
	}

	
    public double[] getSpeciesOccupancy() {
        if (!clockuptodate) update();
		return speciesOccupancy;
    }


    // Updating required for clock model
	private void update() {
		
    	synchronized (this) {
			if (!clockuptodate) {
				
				updateCount++;

				// shrink memory reservation for coalescent times?
				if (! stopPopping &&  (updateCount & 0x7fff) == 0 && maxCoalescentCounts < blocksize - 4) {
					// ensure stored coalescent times are valid, so that a restore gives proper times
	            	double [] stmp = new double[speciesNodeCount * (blocksize - 4) + geneTreeNodeCount];
	            	for (int i = 0; i < speciesNodeCount; i++) {
	            		System.arraycopy(storedCoalescentTimes, i * blocksize, stmp, i * (blocksize - 4), blocksize - 4);
	            	}
            		System.arraycopy(stmp, 0, storedCoalescentTimes, 0, speciesNodeCount * (blocksize - 4));
	            	
					blocksize -= 4;
	            	coalescentTimesLength = speciesNodeCount * blocksize;
	            	// System.err.print("pop");
				}

	        Arrays.fill(speciesOccupancy, 0);
	        
	        // reset arrays as these values need to be recomputed after any changes to the species or gene tree
	        //Arrays.fill(geneNodeSpeciesAssignment, -1); // -1 means no species assignment for that gene tree node has been made yet
	        System.arraycopy(leafGeneNodeSpeciesAssignment, 0, geneNodeSpeciesAssignment, 0, geneTreeNodeCount);
	
	        
	        // Arrays.fill(coalescentLineageCounts, 0);
	        System.arraycopy(leafCoalescentLineageCounts, 0, coalescentLineageCounts, 0, speciesNodeCount);
	        Arrays.fill(coalescentCounts, 0);
	        
	
	        final TreeInterface geneTree = treeInput.get();
	        for (int geneTreeLeafNumber = 0; geneTreeLeafNumber < geneTreeLeafNodeCount; geneTreeLeafNumber++) {
	            final Node geneTreeLeafNode = geneTree.getNode(geneTreeLeafNumber);
	            final int speciesTreeLeafNumber = localTipNumberMap[geneTreeLeafNode.getNr()];
	            final Node speciesTreeLeafNode = speciesTreeInput.get().getNode(speciesTreeLeafNumber);	
	            final Node firstCoalescenceNode = geneTreeLeafNode.getParent();
	            final int firstCoalescenceNumber = firstCoalescenceNode.getNr();
	            final double lastHeight = 0.0;
	
	            if (!collateCoalescenceEvents(geneTreeLeafNumber, lastHeight,
	            		firstCoalescenceNode, firstCoalescenceNumber, 
	            		speciesTreeLeafNode, speciesTreeLeafNumber)) {
	                // this gene tree IS NOT compatible with the species tree
	            	clockuptodate = true;
	                return;
	            }
	        }
	
	        maxCoalescentCounts = 0;
	        for (int j : coalescentCounts) {
	        	if (j > maxCoalescentCounts) {maxCoalescentCounts = j;}
	        }
            if (maxCoalescentCounts > blocksize) {
            	// grow memory reservation for coalescent times
            	int DELTA_BLOCK_SIZE = 4*((maxCoalescentCounts+3)/4) - blocksize;
            	coalescentTimesLength = speciesNodeCount * (blocksize + DELTA_BLOCK_SIZE);
            	double [] tmp = new double[coalescentTimesLength + geneTreeNodeCount];
            	double [] stmp = new double[coalescentTimesLength + geneTreeNodeCount];
            	for (int i = 0; i < speciesNodeCount; i++) {
            		//System.arraycopy(coalescentTimes, i * blocksize, tmp, i * (blocksize + DELTA_BLOCK_SIZE), blocksize);
            		System.arraycopy(storedCoalescentTimes, i * blocksize, stmp, i * (blocksize + DELTA_BLOCK_SIZE), blocksize);
            	}
            	coalescentTimes = tmp;
            	storedCoalescentTimes = stmp;
            	blocksize += DELTA_BLOCK_SIZE;
            	// System.err.print("blocksize = " + blocksize + " ");
            	
            	// do calculation again, this time with properly sized array
            	// we only get here very occasionally (only when blocksize is updated)
            	update();
            	
            	if (updateCount > 0x7fff) {
            		stopPopping = true;
            	}
            	return;
            }


            clockuptodate = true;
			}
    	}
    	
    	
    	
    }
	
	private boolean coalescentTimesChanged(int i) {
	    	int k = i * blocksize;
	    	for (int j = 0; j < coalescentLineageCounts[i]; j++) {
	    		if (coalescentTimes[k] != storedCoalescentTimes[k]) {
	    			return true;
	    		}
	    		k++;
	    	}
			return false;
	}
	

    // Non-recursively populates speciesOccupancy, 
    private boolean collateCoalescenceEvents(int lastGeneTreeNodeNumber, double lastHeight, Node geneTreeNode, int geneTreeNodeNumber, Node speciesTreeNode, int speciesTreeNodeNumber) {
        while (true) {
            final double geneTreeNodeHeight = geneTreeNode.getHeight();

            // Check if the next coalescence event occurs in an ancestral branch
            while (!speciesTreeNode.isRoot() && geneTreeNodeHeight >= speciesTreeNode.getParent().getHeight()) {
                /* if (geneTreeNode.isDirty() != Tree.IS_CLEAN )
                    speciesBranchIsDirty[speciesTreeNodeNumber] = true;
                */
                final Node speciesTreeParentNode = speciesTreeNode.getParent();
                final double speciesTreeParentHeight = speciesTreeParentNode.getHeight();
                final int speciesTreeParentNodeNumber = speciesTreeParentNode.getNr();

                speciesOccupancy[lastGeneTreeNodeNumber * speciesNodeCount + speciesTreeNodeNumber] = speciesTreeParentHeight - lastHeight;
                coalescentLineageCounts[speciesTreeParentNodeNumber]++;

                speciesTreeNode = speciesTreeParentNode;
                speciesTreeNodeNumber = speciesTreeParentNodeNumber;
                lastHeight = speciesTreeParentHeight;
            }

            // This code executes if the next coalescence event occurs within the current branch
            speciesOccupancy[lastGeneTreeNodeNumber * speciesNodeCount + speciesTreeNodeNumber] = geneTreeNodeHeight - lastHeight;
            final int existingSpeciesAssignment = geneNodeSpeciesAssignment[geneTreeNodeNumber];
            if (existingSpeciesAssignment == -1) {
                geneNodeSpeciesAssignment[geneTreeNodeNumber] = speciesTreeNodeNumber;

                coalescentTimes[speciesTreeNodeNumber * blocksize + coalescentCounts[speciesTreeNodeNumber]++] = geneTreeNodeHeight;

                final Node nextGeneTreeNode = geneTreeNode.getParent();
                if (nextGeneTreeNode == null) {
                    // This is the root of the gene tree and no incompatibilities were detected
                    return true;
                } else {
                    // If this is not the root of the gene tree, check the subsequent (back in time) coalescence event
                    lastGeneTreeNodeNumber = geneTreeNodeNumber;
                    lastHeight = geneTreeNodeHeight;
                    geneTreeNode = nextGeneTreeNode;
                    geneTreeNodeNumber = nextGeneTreeNode.getNr();
                }
            } else if (existingSpeciesAssignment == speciesTreeNodeNumber) {
                return true; // Gene tree OK up to here, but stop evaluating because deeper nodes have already been traversed
            } else {
                return false; // This gene tree IS NOT compatible with the species tree
            }
        }
    }
    
    
    
	public double[] getCoalescentTimes(int nodeI) {
        if (!clockuptodate) update();

        final Node speciesNode = speciesTreeInput.get().getNode(nodeI);
        final Node parentNode = speciesNode.getParent();

        final double speciesEndTime = speciesNode.getHeight();
        final double speciesStartTime = (parentNode == null) ? Double.POSITIVE_INFINITY : parentNode.getHeight();
        final int branchEventCount = coalescentCounts[nodeI];

		final double[] branchCoalescentTimes = new double[branchEventCount + 2];
		branchCoalescentTimes[0] = speciesEndTime;
        branchCoalescentTimes[branchEventCount + 1] = speciesStartTime;

		System.arraycopy(coalescentTimes, nodeI * blocksize, branchCoalescentTimes, 1, branchEventCount);
		Arrays.sort(branchCoalescentTimes);

		return branchCoalescentTimes;
	}
    

	
	
	
}
