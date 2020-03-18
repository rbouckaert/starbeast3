package starbeast3;



import java.util.ArrayList;
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
import beast.evolution.alignment.Taxon;
import beast.evolution.alignment.TaxonSet;
import starbeast3.SpeciesTreePrior;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeDistribution;
import beast.evolution.tree.TreeInterface;
import starbeast3.evolution.speciation.PopulationModel;



@Description("Calculates probability of gene tree conditioned on a species tree (multi-species coalescent)")
public class GeneTreeForSpeciesTreeDistribution extends TreeDistribution {

	//TreeInterfaceSB3
    final public Input<SpeciesTree> speciesTreeInput =
            new Input<>("speciesTree", "species tree containing the associated gene tree", Validate.REQUIRED);

    final public Input<Double> ploidyInput =
            new Input<>("ploidy", "ploidy (copy number) for this gene, typically a whole number or half (default 2 for autosomal_nuclear)", 2.0);
    
    final public Input<Boolean> samplingAlignmentInput =
            new Input<>("sampling", "Set to true if using this class for simulating sequences", false);
    
    final public Input<SpeciesTreePrior> speciesTreePriorInput =
            new Input<>("speciesTreePrior", "defines population function and its parameters");
    
    final public Input<PopulationModel> popModelInput = 
    		new Input<>("populationModel", "Population model used to infer the multispecies coalescent probability for this gene");
    
    
    final public Input<TaxonSet> taxonSetInput = new Input<>("taxonset", "set of taxa mapping lineages to species");
    
    
    private SpeciesTree speciesTree;
    
    
    // intervals for each of the species tree branches
    private PriorityQueue<Double>[] intervalsInput;
    
    // Count nr of lineages at the bottom of species tree branches
    private int[] nrOfLineages;
    private int[] nrOfLineagesStored;
    
    
    // Maps gene tree leaf nodes to species tree leaf nodes. Indexed by node number.
    private int[] nrOfLineageToSpeciesMap;
    
    
   
    
    private RealParameter popSizesBottom;
    private PopulationModel popModel;

    // Ploidy is a constant - cache value of input here
    private double ploidy;
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


    
    
    protected double[] speciesOccupancy; // A [1D] matrix describing the length of overlap between each gene tree node's 
    									 // branch with each species tree node's branch
    protected int[] geneNodeSpeciesAssignment;
    protected int[] storedGeneNodeSpeciesAssignment;
    protected boolean[] speciesTreeNodeGeneNodeAssignment;
    protected boolean[] storedSpeciesTreeNodeGeneNodeAssignment;
    protected double[] storedSpeciesOccupancy;
    protected boolean geneTreeCompatible;
    protected boolean storedGeneTreeCompatible;
    
    
    int[] leafCoalescentLineageCounts;
    
    int updateCount = 0;
    boolean stopPopping = false;
    
    
    // Maps gene tree tip numbers to species tree tip number
    private int[] localTipNumberMap;
    
    
    
    // For storing branch-wise prior probabilities
    private boolean[] speciesBranchIsDirty;
    private double[] perBranchLogP;
    private double[] storedPerBranchLogP;
    
    
    // For sampling a gene tree
    private int currentGeneTreeNodeNumber;
    

    
    
    public TreeInterface getGeneTree() {
    	return treeInput.get();
    }
    
    
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
    	popModel = popModelInput.get();

    	
    	// Coerce to SpeciesTree if it is not already (this enables user to parse TreeParser)
    	//if (speciesTreeInput.get() instanceof SpeciesTree) {
    		//speciesTree = (SpeciesTree) speciesTreeInput.get();
    	//}else {
    		//speciesTree = new SpeciesTree(speciesTreeInput.get().getRoot());
    	//}
    	speciesTree = speciesTreeInput.get();
    	
        final Node[] gtNodes = treeInput.get().getNodesAsArray();
        final int gtLineages = treeInput.get().getLeafNodeCount();
        final Node[] sptNodes = speciesTree.getNodesAsArray();
        final int speciesCount = speciesTree.getNodeCount();
        

        // Calculate node counts
        geneTreeNodeCount = treeInput.get().getNodeCount();
        geneTreeLeafNodeCount = treeInput.get().getLeafNodeCount();
        speciesNodeCount = speciesTree.getNodeCount();


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
        nrOfLineageToSpeciesMap = new int[geneTreeNodeCount];
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
        nrOfLineagesStored = new int[speciesCount];

        if (!samplingAlignmentInput.get()) {
        	final SpeciesTreePrior popInfo = speciesTreePriorInput.get();
        	popSizesBottom = popInfo.popSizesBottomInput.get();
        }
        
        geneNodeSpeciesAssignment = new int[geneTreeNodeCount];
        storedGeneNodeSpeciesAssignment = new int[geneTreeNodeCount];
        
        
        speciesTreeNodeGeneNodeAssignment = new boolean[speciesNodeCount * geneTreeNodeCount];
        storedSpeciesTreeNodeGeneNodeAssignment = new boolean[speciesNodeCount * geneTreeNodeCount];
        
      
        // Generate map of species tree tip node names to node numbers
        final Map<String, Integer> tipNumberMap = speciesTree.getTipNumberMap();
        localTipNumberMap = new int[treeInput.get().getLeafNodeCount()];
        for (int i = 0; i < treeInput.get().getLeafNodeCount(); i++) {
        	final Node geneTreeLeafNode = treeInput.get().getNode(i);
        	final String geneTreeLeafName = geneTreeLeafNode.getID();
        	final int geneTreeLeafNumber = geneTreeLeafNode.getNr();

        	if (tipNumberMap.containsKey(geneTreeLeafName)) // not in BEAUTi
        	    localTipNumberMap[geneTreeLeafNumber] = tipNumberMap.get(geneTreeLeafName);
        }
        
        
        geneTreeCompatible = false;
        storedGeneTreeCompatible = false;
        

        // Allocate memory for coalescent counts, lengths, and times
        coalescentCounts = new int[speciesNodeCount];
        storedCoalescentCounts = new int[speciesNodeCount];
        coalescentTimesLength = speciesNodeCount * blocksize;
        coalescentTimes = new double[coalescentTimesLength + geneTreeNodeCount];
        storedCoalescentTimes = new double[coalescentTimesLength + geneTreeNodeCount];
        
        
        // Allocate memory for species occupancies
        speciesOccupancy = new double[geneTreeNodeCount * speciesNodeCount];
        storedSpeciesOccupancy = new double[geneTreeNodeCount * speciesNodeCount];
        
        
        // Allocate memory for per-species-tree-branch probabilities
        perBranchLogP = new double[speciesNodeCount];
        storedPerBranchLogP = new double[speciesNodeCount];

        speciesBranchIsDirty = new boolean[speciesNodeCount];
        Arrays.fill(speciesBranchIsDirty, true);
        

        
        leafCoalescentLineageCounts = new int[speciesNodeCount];
        
        for (int geneTreeLeafNumber = 0; geneTreeLeafNumber < geneTreeLeafNodeCount; geneTreeLeafNumber++) {
            final Node geneTreeLeafNode = treeInput.get().getNode(geneTreeLeafNumber);
            final int speciesTreeLeafNumber = localTipNumberMap[geneTreeLeafNode.getNr()];
            leafCoalescentLineageCounts[speciesTreeLeafNumber]++;            
        }
        
        
        
        
    }

    /**
     * @param lineageID
     * @return species ID to which the lineage ID belongs according to the TaxonSets
     */
    private String getSpeciesID(final String lineageID) {
        final TaxonSet taxonSuperset = taxonSetInput.get() != null ? taxonSetInput.get() : speciesTreePriorInput.get().taxonSetInput.get();
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
    
    
    // Returns an array of gene tree nodes which map to this species tree node
    public Node[] mapSpeciesNodeToGeneTreeNodes(Node species) {
    	
    	int speciesStartIndex = geneTreeNodeCount * species.getNr();
    	int speciesEndIndex = geneTreeNodeCount * (species.getNr()+1);
    	
    	// Find how many nodes are mapped to this species node
    	int numMapped = 0;
    	for (int i = speciesStartIndex; i < speciesEndIndex; i ++) {
    		if (speciesTreeNodeGeneNodeAssignment[i]) numMapped++;
    	}
    	
    	// Return an array of that length with the gene nodes stored in it
    	Node[] mappedNodes = new Node[numMapped];
    	int j = 0;
    	for (int i = speciesStartIndex; i < speciesEndIndex; i ++) {
    		if (speciesTreeNodeGeneNodeAssignment[i]) {
    			int geneNodeNr = i % geneTreeNodeCount;
    			mappedNodes[j] = treeInput.get().getNode(geneNodeNr);
    			j++;
    		}
    	}

    	return mappedNodes;
    	

        
    }
    

    @Override
    public double calculateLogP() {

    	
        assert SanityChecks.checkTreeSanity(speciesTree.getRoot());
        
        update();

        if (!geneTreeCompatible) {
            logP = Double.NEGATIVE_INFINITY;
            return logP;
        }
        
        
        logP = 0.0;
        
        // Recompute log priors for all branches which are dirty
        final Node[] speciesTreeNodes = speciesTree.getNodesAsArray();
        for (int speciesNodeI = 0; speciesNodeI < speciesNodeCount; speciesNodeI++) {
            Node speciesNode = speciesTreeNodes[speciesNodeI];
            if (isDirtyBranch(speciesNodeI) || popModel.isDirtyBranch(speciesNode)) {
                final int lineagesBottom = nrOfLineages[speciesNodeI];
                final int k = coalescentCounts[speciesNodeI];
                final double[] branchCoalescentTimes = getCoalescentTimes(speciesNodeI);
                
                perBranchLogP[speciesNodeI] = popModel.calculateBranchLogP(lineagesBottom, ploidy, popSizesBottom.getValue(speciesNodeI), branchCoalescentTimes, k); 
            }

            // System.out.println(String.format("%s-%d: %f", getID(), nodeI, logP));
            logP += perBranchLogP[speciesNodeI];
        }

        // System.out.println(String.format("%s-%d: %f", getID(), speciesNodeCount, logP));
        logPuptodate = true;
        return logP;
    
        
    }
    
    
    // Calculates the prior density contribution from this branches lineage history (without using the population size)
    public double calculatePartialLogPBranch(Node speciesNode) {
    	
    	
    	update();
    	
    	final int lineagesBottom = nrOfLineages[speciesNode.getNr()];
        final int k = coalescentCounts[speciesNode.getNr()];
        final double[] branchCoalescentTimes = getCoalescentTimes(speciesNode.getNr());
        
        return popModel.calculatePartialLogPBranch(lineagesBottom, branchCoalescentTimes, k);
    	
    }
    
    
    
    protected boolean isDirtyBranch(int nodeNr) {
		return speciesBranchIsDirty[nodeNr];
	}



	public int[] getTipNumberMap() {
		return localTipNumberMap;
	}


    

    @Override
    public void store() {
    	super.store();
    	
    	
      
        System.arraycopy(coalescentCounts, 0, storedCoalescentCounts, 0, coalescentCounts.length);
        System.arraycopy(coalescentTimes, 0, storedCoalescentTimes, 0, coalescentTimesLength);
        System.arraycopy(nrOfLineages, 0, nrOfLineagesStored, 0, nrOfLineages.length);

        System.arraycopy(geneNodeSpeciesAssignment, 0, storedGeneNodeSpeciesAssignment, 0, geneNodeSpeciesAssignment.length);
        System.arraycopy(speciesTreeNodeGeneNodeAssignment, 0, storedSpeciesTreeNodeGeneNodeAssignment, 0, speciesTreeNodeGeneNodeAssignment.length);
        
        
        
        System.arraycopy(speciesOccupancy, 0, storedSpeciesOccupancy, 0, speciesOccupancy.length);

        System.arraycopy(perBranchLogP, 0, storedPerBranchLogP, 0, perBranchLogP.length);
        
        storedGeneTreeCompatible = geneTreeCompatible;
        storedMaxCoalescentCounts = maxCoalescentCounts;
        
    }
    
  


    @Override
    public void restore() {
    	logPuptodate = false;
    	super.restore();
    	
    	
    	
    	double[] tmpCoalescentTimes = coalescentTimes;
    	int[] tmpCoalescentCounts = coalescentCounts;
    	int[] tmpCoalescentLineageCounts = nrOfLineages;
    	int[] tmpGeneNodeSpeciesAssignment = geneNodeSpeciesAssignment;
    	boolean[] tmpSpeciesTreeNodeGeneNodeAssignment = speciesTreeNodeGeneNodeAssignment;
    	double[] tmpSpeciesOccupancy = speciesOccupancy;
    	double[] tmpPerBranchLogP = perBranchLogP;
    	boolean tmpGeneTreeCompatible = geneTreeCompatible;

    	coalescentTimes = storedCoalescentTimes;
    	coalescentCounts = storedCoalescentCounts;
    	nrOfLineages = nrOfLineagesStored;
    	speciesOccupancy = storedSpeciesOccupancy;
    	geneNodeSpeciesAssignment = storedGeneNodeSpeciesAssignment;
    	speciesTreeNodeGeneNodeAssignment = storedSpeciesTreeNodeGeneNodeAssignment;
    	perBranchLogP = storedPerBranchLogP;
    	geneTreeCompatible = storedGeneTreeCompatible;

    	storedCoalescentTimes = tmpCoalescentTimes;
    	storedCoalescentCounts = tmpCoalescentCounts;
    	nrOfLineagesStored = tmpCoalescentLineageCounts;
    	storedSpeciesOccupancy = tmpSpeciesOccupancy;
    	storedGeneNodeSpeciesAssignment = tmpGeneNodeSpeciesAssignment;
    	storedSpeciesTreeNodeGeneNodeAssignment = tmpSpeciesTreeNodeGeneNodeAssignment;
    	storedPerBranchLogP = tmpPerBranchLogP;
    	storedGeneTreeCompatible = tmpGeneTreeCompatible;

    	maxCoalescentCounts = storedMaxCoalescentCounts;
		
    	
    	
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
        List<String> arguments = new ArrayList<>();
        arguments.add(speciesTreePriorInput.get().getID());
        return arguments;
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
		return getCoalescentTimes(i).length - 2;
	}
	
	public double [] getTimes(int i) {
		if (!logPuptodate) {
			calculateLogP();			
		}
		return getCoalescentTimes(i);
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
		        System.arraycopy(nrOfLineageToSpeciesMap, 0, geneNodeSpeciesAssignment, 0, geneTreeNodeCount);
		        
		        for (int i = 0; i < speciesTreeNodeGeneNodeAssignment.length; i ++) {
		        	speciesTreeNodeGeneNodeAssignment[i] = false;
		        }
		      
		
		        
		        // Arrays.fill(coalescentLineageCounts, 0);
		        System.arraycopy(leafCoalescentLineageCounts, 0, nrOfLineages, 0, speciesNodeCount);
		        Arrays.fill(coalescentCounts, 0);
		        
		        
		        
		        Arrays.fill(speciesBranchIsDirty, false);
		        
		
		        final TreeInterface geneTree = treeInput.get();
		        for (int geneTreeLeafNumber = 0; geneTreeLeafNumber < geneTreeLeafNodeCount; geneTreeLeafNumber++) {
		            final Node geneTreeLeafNode = geneTree.getNode(geneTreeLeafNumber);
		            final int speciesTreeLeafNumber = localTipNumberMap[geneTreeLeafNode.getNr()];
		            final Node speciesTreeLeafNode = speciesTree.getNode(speciesTreeLeafNumber);	
		            final Node firstCoalescenceNode = geneTreeLeafNode.getParent();
		            final int firstCoalescenceNumber = firstCoalescenceNode.getNr();
		            final double lastHeight = 0.0;
		
		            if (!collateCoalescenceEvents(geneTreeLeafNumber, lastHeight,
		            		firstCoalescenceNode, firstCoalescenceNumber, 
		            		speciesTreeLeafNode, speciesTreeLeafNumber)) {
		                // this gene tree IS NOT compatible with the species tree
		            	geneTreeCompatible = false;
		            	clockuptodate = true;
		                return;
		            }
	
		            
		        }
		        
		        
		        // Reverse map the geneToSpecies map into the speciesToGene map
		        for (int geneTreeNodeNr = 0; geneTreeNodeNr < geneTreeNodeCount; geneTreeNodeNr++) {
		        	final int speciesNodeNrMappedTo = geneNodeSpeciesAssignment[geneTreeNodeNr];
		        	final int mapIndex = speciesNodeNrMappedTo * geneTreeNodeCount + geneTreeNodeNr;
		        	speciesTreeNodeGeneNodeAssignment[mapIndex] = true;
		        	
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
	            
	            
	            for (int i = 0; i < speciesNodeCount; i++) {
		        	if (nrOfLineages[i] != nrOfLineagesStored[i] ||
		        		coalescentCounts[i] != storedCoalescentCounts[i]) {
		        		speciesBranchIsDirty[i] = true;
		        	} else {
		        		Node node = speciesTree.getNode(i);
		        		if (node.isDirty() != Tree.IS_CLEAN ||
		        			(!node.isRoot() && node.getParent().isDirty() != Tree.IS_CLEAN) ||
		        			coalescentTimesChanged(i)) speciesBranchIsDirty[i] = true;
		        	}
		        }
	
	            geneTreeCompatible = true;
	            clockuptodate = true;
	            
			}
    	}
    	
    	
    	
    }
	
	private boolean coalescentTimesChanged(int i) {
	    	int k = i * blocksize;
	    	for (int j = 0; j < nrOfLineages[i]; j++) {
	    		if (coalescentTimes[k] != storedCoalescentTimes[k]) {
	    			return true;
	    		}
	    		k++;
	    	}
			return false;
	}
	

    // Iteratively populates speciesOccupancy, coalescentTimes, and geneNodeSpeciesAssignment
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
                nrOfLineages[speciesTreeParentNodeNumber]++;

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

        final Node speciesNode = speciesTree.getNode(nodeI);

        final double speciesEndTime = speciesNode.getHeight();
        
        final double speciesStartTime = speciesNode.isRoot() ? Math.max(speciesNode.getHeight(), treeInput.get().getRoot().getHeight()) 
        													 : speciesNode.getParent().getHeight();
        
        final int branchEventCount = coalescentCounts[nodeI];

		final double[] branchCoalescentTimes = new double[branchEventCount + 2];
		branchCoalescentTimes[0] = speciesEndTime;
        branchCoalescentTimes[branchEventCount + 1] = speciesStartTime;

		System.arraycopy(coalescentTimes, nodeI * blocksize, branchCoalescentTimes, 1, branchEventCount);
		Arrays.sort(branchCoalescentTimes);

		return branchCoalescentTimes;
	}
	
	
	
	@Override
    public void sample(final State state, final Random random) {
    	
    	if (sampledFlag) return;
        sampledFlag = true;
        
        
        // Sample 
        sampleConditions(state, random);
        
        
        // Get branch population size
        SpeciesTreePrior speciesTreePrior = speciesTreePriorInput.get();
        RealParameter popSizes = speciesTreePrior.getPopulationSizes();
        
        currentGeneTreeNodeNumber = treeInput.get().getLeafNodeCount();
        
        List<Node> root = sampleBranchRecursive(speciesTree.getRoot(), random, popSizes);
        final Tree geneTree = (Tree) treeInput.get();
        assert root.size() == 1;
        
        Tree sampledGeneTree = new Tree(root.get(0));
        assert sampledGeneTree.getNodeCount() == treeInput.get().getNodeCount();
        geneTree.assignFromWithoutID(sampledGeneTree);
        
    	
    	
    }
    
    
    private List<Node> sampleBranchRecursive(final Node node, final Random random, final RealParameter popSizes) {
    	
    	
    	List<Node> geneTreeNodesInLineage = new ArrayList<Node>();
    	
    	// Population size of this species branch
    	final double popSizeBranch = popSizes.getValue(node.getNr()) * ploidy;
    	
    	
    	// Sample coalescent events of children first
    	if (!node.isLeaf()) {
    		
    		for (Node c : node.getChildren()) {
    			List<Node> childsTopNodes = sampleBranchRecursive(c, random, popSizes);
    			geneTreeNodesInLineage.addAll(childsTopNodes);
    			
    		}
    		
    	}
    	
    	
    	// Find all the gene tree leaves which map to this species node
    	final Map<String, Integer> tipNumberMap = speciesTree.getTipNumberMap();
        for (int i = 0; i < treeInput.get().getLeafNodeCount(); i++) {
        	final Node geneTreeLeafNode = treeInput.get().getNode(i);
        	final String geneTreeLeafName = geneTreeLeafNode.getID();

        	if (tipNumberMap.get(geneTreeLeafName) == node.getNr()) {
        		geneTreeNodesInLineage.add(geneTreeLeafNode);
        	}
        	    
        } 
    	
    	
    		
        // Simulate coalescence along this branch
        // Stop when either there is only 1 more lineage or there is no time remaining
        int numLineages = geneTreeNodesInLineage.size();
        double time = node.getHeight();
        double parentTime = node.isRoot() ? Double.POSITIVE_INFINITY : node.getParent().getHeight();
        while (numLineages > 1) {
        	
        	
        	// Sample the time to the next coalescence
        	double rateOfCoalescence = numLineages * (numLineages-1) / 2 / popSizeBranch;
        	double timeToCoalescence = -Math.log(random.nextDouble())/rateOfCoalescence;
        	time += timeToCoalescence;
        	if (time > parentTime) break;
        	
        	
        	// Uniformly at random select two nodes to coalesce
        	int leftChildNr = random.nextInt(numLineages);
        	int rightChildNr = random.nextInt(numLineages);
        	while (leftChildNr == rightChildNr) rightChildNr = random.nextInt(numLineages);
        	Node parentNode = new Node();
        	parentNode.setNr(currentGeneTreeNodeNumber);
        	parentNode.setHeight(time);
        	parentNode.addChild(geneTreeNodesInLineage.get(leftChildNr));
        	parentNode.addChild(geneTreeNodesInLineage.get(rightChildNr));
        	currentGeneTreeNodeNumber++;
        	
        	
        	// Remove children from list and add parent
        	geneTreeNodesInLineage.remove(Math.max(leftChildNr, rightChildNr));
        	geneTreeNodesInLineage.remove(Math.min(leftChildNr, rightChildNr));
        	geneTreeNodesInLineage.add(parentNode);
        	numLineages--;
        	
        }
    	
    	
    	
    	return geneTreeNodesInLineage;
    	
	    	
    }
	
	
	
}
