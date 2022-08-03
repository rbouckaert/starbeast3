package starbeast3.genekernel;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.List;
import java.util.Random;
import java.util.Set;

import beast.base.inference.Distribution;
import beast.base.core.Input;
import beast.base.inference.State;
import beast.base.inference.StateNode;
import beast.base.core.Input.Validate;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.core.Log;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import starbeast3.GeneTreeForSpeciesTreeDistribution;
import starbeast3.SpeciesTree;
import starbeast3.SpeciesTreePrior;
import starbeast3.evolution.speciation.PopulationModel;

public class GTKPrior extends Distribution {

	
	
	final public Input<List<GTKPointerTree>> genesInput = 
			new Input<>("pointer", "Gene trees which point to the kernel trees", new ArrayList<>());
	
	final public Input<GeneTreeKernel> kernelInput = 
			new Input<>("kernel", "The gene tree kernel which contains a set of trees", Input.Validate.REQUIRED);
	
	final public Input<IntegerParameter> geneKernelSizeInput = 
			new Input<>("geneKernelSize", "The parameter governing the gene kernel size", Input.Validate.REQUIRED);
	
	final public Input<IntegerParameter> indicatorInput = 
			new Input<>("indicator", "A parameter which points each observed gene tree to one in the kernel", Input.Validate.REQUIRED);
	
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
    
    final public Input<TaxonSet> taxonSetInput = 
    		new Input<>("taxonset", "set of taxa mapping lineages to species");
    
    
    final public Input<Boolean> pointerPriorInput = 
    		new Input<>("pointerPrior", "If true, only tree pointers have prior density (except for illegal kernal trees which are rejected). If false,"
    				+ "only kernel trees have prior density. ", true);
	
    
	List<GeneTreeForSpeciesTreeDistribution> geneDistributions;
	List<GeneTreeForSpeciesTreeDistribution> geneDistributions_stored;
	
	List<GTKPointerTree> pointers;
	
	GeneTreeKernel kernel;
	IntegerParameter geneKernelSize;
	IntegerParameter indicator;
	boolean needsUpdating;
	boolean pointerPrior;
	
	
	
	// Storing
	GeneTreeForSpeciesTreeDistribution deletedPrior;
	int deletedPriorIndex;
	boolean expanded;
	
	
	@Override
	public void initAndValidate() {
		this.pointers = genesInput.get();
		if (this.pointers.size() == 0) throw new IllegalArgumentException("Please ensure that there is at least 1 gene tree.");
		this.kernel = kernelInput.get();
		this.pointerPrior = pointerPriorInput.get();
		
		this.geneDistributions_stored = new ArrayList<GeneTreeForSpeciesTreeDistribution>();
		
		// Ensure that the kernel has room for at least 1 tree
		this.geneKernelSize = geneKernelSizeInput.get();
		int kernelSize = this.geneKernelSize.getValue();
		if (this.geneKernelSize.getLower() < 1) this.geneKernelSize.setLower(1);
		if (kernelSize < this.geneKernelSize.getLower()) this.geneKernelSize.setValue(this.geneKernelSize.getLower());
		kernelSize = this.geneKernelSize.getValue();
		
		// Initialise the indicator (one pointer per gene tree)
		indicator = indicatorInput.get();
		Integer[] oldValues = indicator.getValues();
		indicator.setDimension(this.pointers.size());
		indicator.setLower(0);
		indicator.setUpper(this.geneKernelSize.getValue()-1);
		
		// Randomly sample the indicator if it was not already set
		if (oldValues.length != indicator.getDimension()) {
			for (int i = 0; i < indicator.getDimension(); i ++) {
				indicator.setValue(i, Randomizer.nextInt(this.geneKernelSize.getValue()));
			}
		}
		
		
		// Ensure that all gene trees have the same taxonset
		Tree geneTree1 = this.pointers.get(0);
		String[] taxaNames1 = geneTree1.getTaxaNames();
		List<String> taxaNamesSorted1 = Arrays.asList(taxaNames1);
		Collections.sort(taxaNamesSorted1);
		
		for (int i = 1; i < this.pointers.size(); i ++) {
			Tree geneTree2 = this.pointers.get(i);
			String[] taxaNames2 = geneTree2.getTaxaNames();
			List<String> taxaNamesSorted2 = Arrays.asList(taxaNames2);
			Collections.sort(taxaNamesSorted2);
			
			// Compare all strings
			if (taxaNamesSorted1.size() != taxaNamesSorted2.size()) {
				throw new IllegalArgumentException("Please ensure that all gene trees share the same taxa. " + geneTree1.getID() + " != " + geneTree2.getID() + ".");
			}
			for (int j = 0; j < taxaNamesSorted1.size(); j ++) {
				String taxon1 = taxaNamesSorted1.get(j);
				String taxon2 = taxaNamesSorted2.get(j);
				if (!taxon1.equals(taxon2)) {
					throw new IllegalArgumentException("Please ensure that all gene trees share the same taxa. " + geneTree1.getID() + " != " + geneTree2.getID() + ".");
				}
			}
		}
		
		
		
		// Create a taxonset which contains all gene tree taxa but without any data
		List<String> taxaNames = this.pointers.get(0).m_taxonset.get().asStringList();
		TaxonSet taxonset = new TaxonSet(Taxon.createTaxonList(taxaNames));
		
		
		// Initialise the kernel trees using the taxonset
		for (int i = 0; i < kernelSize; i ++) {
			GTKGeneTree kernelTree = new GTKGeneTree();
			kernelTree.initByName("taxonset", taxonset);
			kernel.addTree(kernelTree);
		}
		
		
		
		// Point each gene tree to a kernel
		for (int i = 0; i < this.pointers.size(); i ++) {
			this.pointers.get(i).updateTree();
		}
		
		
		
		// Create the GeneTreeForSpeciesTreeDistribution priors
		this.geneDistributions = new ArrayList<GeneTreeForSpeciesTreeDistribution>();
		for (Tree kernelTree : kernel.getTrees()) {
			this.addGeneTreeDistribution(kernelTree);
			
		}
		
		
		this.deletedPrior = null;
		this.expanded = false;
		this.deletedPriorIndex = -1;
		
		needsUpdating = false;
		
	}
	
	/**
	 * @return The list of GeneTreeForSpeciesTreeDistribution priors associated with this kernel
	 */
	public List<GeneTreeForSpeciesTreeDistribution> getGeneTreeDistributions() {
		update();
		return this.geneDistributions;
	}
	
	
	/**
	 * Same as getGeneTreeDistributions() but with this difference that the State can manage
	 * whether to make a copy and register the operator.
	 * <p/>
	 * Only Operators should call this method.
	 * @param operator
	 * @return
	 */
	public List<GeneTreeForSpeciesTreeDistribution> getGeneTreeDistributions(final Operator operator) {
		update();
		this.getKernel().startEditing(operator);
		return this.geneDistributions;
	}

	
	/**
	 * @param index of the tree pointer
	 * @return The GeneTreeForSpeciesTreeDistribution prior at this index
	 */
	public GeneTreeForSpeciesTreeDistribution getGeneTreeDistributions(int index) {
		update();
		return this.geneDistributions.get(index);
	}
	
	
	/**
	 * @return The list of trees pointing to the kernel trees
	 */
	public List<GTKPointerTree> getPointerTrees(){
		return this.pointers;
	}




	
	
	@Override
	public double calculateLogP() {
		
		update();
		logP = 0;
		
		// If the kernel size is not equal to the size it should be, return -Inf
		if (this.getKernel().getDimension() != this.geneKernelSize.getValue() || this.getKernel().getDimension() != this.getGeneTreeDistributions().size()) {
			Log.warning("WError: GTKPrior the size of the kernel should be equal to 'geneKernelSize' and the number of gene prior distributions. " + this.getKernel().getDimension() + " != " + this.geneKernelSize.getValue() + " != " + this.getGeneTreeDistributions().size());
			logP = Double.NEGATIVE_INFINITY;
			return logP;
		}
		
		// Gene tree distribution prior
		int nKernelTrees = this.geneKernelSize.getValue();
		boolean anIllegalTree = false;
		for (int kernelIndex = 0; kernelIndex < this.geneDistributions.size(); kernelIndex ++) {
			
			GeneTreeForSpeciesTreeDistribution prior = this.geneDistributions.get(kernelIndex);
			double priorLogP = prior.calculateLogP();
			
			
			// Always reject an illegal gene trees even if no-one points to it
			if (priorLogP == Double.NEGATIVE_INFINITY) anIllegalTree = true;
			
			
			// Weight each count gene tree prior by how many trees are pointing to it
			if (this.pointerPrior) {
				
				int numPointingTo = 0;
				for (int pointerIndex = 0; pointerIndex < this.indicator.getDimension(); pointerIndex ++) {
					if (this.indicator.getValue(pointerIndex) == kernelIndex) numPointingTo ++;
				}
				
				logP += numPointingTo * priorLogP;
				
			}
			
			// Count all gene tree priors even if they are not being pointed to
			else {
				logP += priorLogP;
			}
			
			
			//if (priorLogP > 0) System.out.println("prior: " + priorLogP + ":" + prior.treeInput.get().getRoot().toNewick());
			//if (Double.isNaN(priorLogP)) throw new ArrayIndexOutOfBoundsException("prior: " + priorLogP + ":" + prior.treeInput.get().getRoot().toNewick());
			
		}
		
		
		if (anIllegalTree) {
			logP = Double.NEGATIVE_INFINITY;
			return logP;
		}
		
		// A uniform prior over the number of trees in the kernel
		// If there are p observed genes and q kernel genes, then the prior is (1/q)^p
		logP += -this.pointers.size() * Math.log(nKernelTrees);
		return logP;
	}
	

	/**
	 * Adds a new gene tree distribution to the end of the list
	 * @param kernelTree
	 */
	public void addGeneTreeDistribution(Tree kernelTree) {
		this.expanded = true;
		this.geneDistributions.add(createGeneTreeDistribution(kernelTree));
	}
	
	/**
	 * @param kernelTree
	 * @return A 'GeneTreeForSpeciesTreeDistribution' object created from 'kernelTree' and the Inputs of this class
	 */
	private GeneTreeForSpeciesTreeDistribution createGeneTreeDistribution(Tree kernelTree) {
		GeneTreeForSpeciesTreeDistribution prior = new GeneTreeForSpeciesTreeDistribution();
		prior.initByName(	"speciesTree", speciesTreeInput.get(), 
							"tree", kernelTree,
							"ploidy", ploidyInput.get(),
							"taxonset", taxonSetInput.get(),
							"sampling", samplingAlignmentInput.get(),
							"speciesTreePrior", speciesTreePriorInput.get(),
							"populationModel", popModelInput.get()
						);
		return prior;
	}
	
	
	/**
	 * Deletes the gene tree distribution at the specified position within the list
	 * @param index
	 */
	public void removeGeneTreeDistribution(int index) {
		this.deletedPrior = this.geneDistributions.get(index);
		this.deletedPriorIndex = index;
		this.geneDistributions.remove(index);
	}
	
	
	
	private boolean needsUpdating() {

		if (needsUpdating) return true;
		if (kernel.isDirtyCalculation()) return true;
		if (geneKernelSize.isDirtyCalculation()) return true;
		if (indicator.isDirtyCalculation()) return true;
		
		return false;
	}
	
	
	private void update() {
		
		
		 synchronized (this) {
			 if (!this.needsUpdating()) return;
			 
			 
			 for (GTKPointerTree pointer: this.pointers) {
				 //if (indicator.isDirtyCalculation()) pointer.updateTree();
				 //else pointer.update();
				 pointer.update();
			 }
			 
			 
			 needsUpdating = false;
		 }
		 
	 }
	
	
	 @Override
	 public void store() {
		 super.store();
		 
		 //this.geneDistributions_stored.clear();
		 for (GeneTreeForSpeciesTreeDistribution prior : this.geneDistributions) {
			 prior.store();
			// geneDistributions_stored.add(prior);
		 }
		 
		 
		 
	 }
	
	
	 @Override
	 public void restore() {
		 super.restore();
		 
		 // If a new tree prior was added, then remove it
		 if (this.expanded) {
			 this.geneDistributions.remove(this.geneDistributions.size()-1);
		 }
		 
		 // If a tree prior was deleted, then add it back
		 else if (this.deletedPrior != null) {
			 this.geneDistributions.add(this.deletedPriorIndex, this.deletedPrior);
		 }
		 
		 this.expanded = false;
		 this.deletedPrior = null;
		 this.deletedPriorIndex = -1;
		 
		 
		 // Restore the priors
		 for (GeneTreeForSpeciesTreeDistribution prior : this.geneDistributions) {
			 prior.restore();
		 }
		 
		 // Ensure that the maximum indicator value is the size of the gene tree kernel
		 this.indicator.setUpper(this.geneKernelSize.getValue()-1);

		 needsUpdating = true;
		 update();
		 
	 }
	 
	 @Override
	 protected void accept() {
		 this.expanded = false;
		 this.deletedPrior = null;
		 this.deletedPriorIndex = -1;
	 }
	
	
    @Override
    public boolean requiresRecalculation() {
    	needsUpdating = true;
    	for (GeneTreeForSpeciesTreeDistribution prior : this.geneDistributions) {
			 prior.requiresRecalculation();
		}
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

	
	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub
		
	}

	public GeneTreeKernel getKernel() {
		return this.kernel;
	}

	
	

}
