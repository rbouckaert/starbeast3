package genekernel;

import java.util.ArrayList;
import java.util.List;

import beast.core.BEASTInterface;
import beast.core.Input;
import beast.core.Operator;
import beast.core.StateNode;
import beast.evolution.operators.KernelDistribution;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import starbeast3.GeneTreeForSpeciesTreeDistribution;


/**
 * An class for any Operators which act on gene trees, and also has support for Bactrian kernels
 * This abstract class serves as a reminder to probe the kernel for gene trees at the start of each proposal
 * And to ensure that if a kernel is provided, then gene trees are not also provided
 * @author Jordan Douglas
 */
public abstract class GTKOperator extends Operator {

	// Gene trees or a gene tree kernel
	public final Input<List<GeneTreeForSpeciesTreeDistribution>> geneTreesInput = new Input<>("gene", "list of gene trees that constrain species tree movement", new ArrayList<>());
	public final Input<GTKPrior> geneTreeKernelPriorInput = new Input<>("kernel", "the kernel of gene trees");

	public final Input<KernelDistribution> kernelDistributionInput = new Input<>("kernelDistribution", "provides sample distribution for proposals", 
	KernelDistribution.newDefaultKernelDistribution());
    protected KernelDistribution kernelDistribution;
	
	
	protected List<Tree> geneTrees;
	protected List<GeneTreeForSpeciesTreeDistribution> geneTreeDistributions;
	protected GTKPrior geneTreeKernelPrior = null;
    
    
	
	@Override
	public void initAndValidate() {
		
		
		geneTreeKernelPrior = geneTreeKernelPriorInput.get();
		kernelDistribution = kernelDistributionInput.get();
		
		if ((geneTreesInput.get().size() == 0 && geneTreeKernelPrior == null) ||
			(geneTreesInput.get().size() > 0 && geneTreeKernelPrior != null) ) {
					throw new IllegalArgumentException("Must specify either 'gene' or 'kernel' (but not both)");
		}

		
		this.geneTreeDistributions = this.getTreeDistributions(null);
		
		// If using a list of genes, add the tree objects to the list now. If using a kernel, this needs to be done on the fly
		this.geneTrees = new ArrayList<Tree>();
		if (geneTreeKernelPrior == null) {
			for (GeneTreeForSpeciesTreeDistribution t : this.geneTreeDistributions ) {
				this.geneTrees.add((Tree) t.getGeneTree());
			}
		}  
			
		
	}
	
	
	
	/**
	 * Uniformly at random samples a tree from either the gene tree kernel, or the list of gene trees
	 * Method assumes that this tree will be edited
	 * @return a random tree to operate on
	 */
	public Tree sampleTree(Operator operator) {
		this.geneTrees = this.getTrees(operator);
		int index = Randomizer.nextInt(geneTrees.size());
		Tree tree = this.geneTrees.get(index);
		//tree.startEditing(operator);
		

		
		return tree;
	}
	
	
	
	
	/**
	 * Gets the number of gene trees. This number is subject to change if a gene tree kernel is used, otherwise will be constant.
	 * 
	 * Use this function independently at the start of each proposal. Do not rely on the Input elements.
	 * 
	 * @return Current number of gene trees
	 */
	public int getGeneTreeCount(Operator operator) {
		if (geneTreeKernelPrior == null) {
    		return geneTreesInput.get().size();
    	}else {
    		return geneTreeKernelPrior.getGeneTreeDistributions(operator).size();
    	}
	}
	
	
	
    /**
     * Gets the list of GeneTreeForSpeciesTreeDistribution (which point to trees) that this operator may act on
     * If using a gene tree kernel, then this list of trees changes throughout the MCMC
     * Otherwise, the list is constant
     * 
     * Use this function independently at the start of each proposal. Do not rely on the Input elements.
     * 
     * @return the set of GeneTreeForSpeciesTreeDistribution objects that this operator may act on
     */
    public List<GeneTreeForSpeciesTreeDistribution> getTreeDistributions(Operator operator) {
    	if (geneTreeKernelPrior == null) {
    		return geneTreesInput.get();
    	}else {
    		return geneTreeKernelPrior.getGeneTreeDistributions(operator);
    	}
    }
    
    
    /**
     * Gets the list of Trees which this operator may act on
     * If using a gene tree kernel, then this list of trees changes throughout the MCMC
     * Otherwise, the list is constant
     * 
     * Use this function independently at the start of each proposal. Do not rely on the Input elements.
     * 
     * @return the set of Trees which this operator may act on
     */
    public List<Tree> getTrees(Operator operator) {
    	
    	// If using a list of gene trees, the list is constant
    	if (geneTreeKernelPrior == null) {
    		return this.geneTrees;
    	}
    	
    	// If using a gene tree kernel, must generate a new list every time since the list of trees is always changing
    	else {
    		List<Tree> trees = new ArrayList<Tree>();
    		for (GeneTreeForSpeciesTreeDistribution t : this.getTreeDistributions(operator) ) {
    			trees.add((Tree) t.getGeneTree());
    		}
    		return trees;
    	}
    	
    	
    }
    
    
    @Override
    public List<StateNode> listStateNodes() {
  
        // Pick up all inputs that are stateNodes that are estimated
        final List<StateNode> stateNodes = super.listStateNodes();
        
        // Ensure that BEAST2 understands this operator changes the kernel and the trees pointing to the kernel
        if (geneTreeKernelPrior != null) {
        	stateNodes.add(this.geneTreeKernelPrior.getKernel());
        	for (GTKPointerTree tree : this.geneTreeKernelPrior.getPointerTrees()) {
        		stateNodes.add(tree);
        	}
        	
        }
        return stateNodes;
    }

    
    
	@Override
    public void accept() {
        super.accept();
    }
	
	@Override
    public void reject() {
		super.reject();
    }
	
	@Override
	public void reject(final int reason) {
		super.reject(reason);
    }

	
	
	
}



