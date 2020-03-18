package starbeast3.operators;

import genekernel.GTKGeneTree;
import genekernel.GTKOperator;
import genekernel.GTKPointerTree;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.distribution.PoissonDistribution;
import org.apache.commons.math.distribution.PoissonDistributionImpl;
import org.apache.commons.math.util.MathUtils;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.util.Log;
import beast.evolution.tree.Tree;
import beast.math.distributions.Poisson;
import beast.util.Randomizer;


/**
 * Adds/removes a gene tree to/from the kernel and assigns/unassigns gene pointers accordingly
 * @author Jordan Douglas
 *
 */
public class GTKExpander extends GTKOperator {
	
	final public Input<IntegerParameter> indicatorInput = 
			new Input<>("indicator", "A parameter which points each observed gene tree to one in the kernel", Input.Validate.REQUIRED);
	
	final public Input<IntegerParameter> geneKernelSizeInput = 
			new Input<>("geneKernelSize", "A parameter which controls the size of the gene kernel", Input.Validate.REQUIRED);
	
	final public Input<List<GTKPointerTree>> genesInput = 
			new Input<>("pointer", "Gene trees which point to the kernel trees", new ArrayList<>());
	
	final public Input<Double> poissonScaleInput = 
			new Input<>("poisson", "The number of trees pointers to change is Poisson distributed, with a mean of 'poisson' * 'num pointers' / 'kernel size' ", 1.0);

	
	
	final boolean DEBUG = false;
	List<GTKPointerTree> pointers;
	IntegerParameter indicator;
	IntegerParameter geneKernelSize;
	double poissonScale;
	
	public GTKExpander() {
		geneTreesInput.setRule(Input.Validate.FORBIDDEN);
	}
	
	@Override
	public void initAndValidate() {
		this.indicator = indicatorInput.get();
		this.pointers = genesInput.get();
		this.geneKernelSize = geneKernelSizeInput.get();
		this.poissonScale = poissonScaleInput.get();
		super.initAndValidate();
	}
	
	
	@Override
	public double proposal() {
		
		
		
		geneTreeDistributions = this.getTreeDistributions(this);
		
		double log_q_expand = 0;
		double log_q_contract = 0;
		int originalKernelSize = this.geneTreeKernelPrior.getKernel().getDimension();
		
		// Sanity check
		if (originalKernelSize != this.geneKernelSize.getValue()) {
			return Double.NEGATIVE_INFINITY;
		}
		
		
		// Expanding or contracting? 
		boolean expanding = Randomizer.nextBoolean();
		if (originalKernelSize == this.geneKernelSize.getLower()) expanding = true;
		else if (originalKernelSize == this.geneKernelSize.getUpper()) expanding = false;
		
		//expanding = true;
		
		// Calculations
		int proposedKernelSize = originalKernelSize + (expanding ? 1 : -1);
		double poissonRate = this.poissonScale * this.pointers.size() / (expanding ? originalKernelSize : proposedKernelSize);
		int numPointersReassigned = 0;
		

		if (DEBUG) System.out.println((expanding ? "Expanding" : "Contracting ") + " the gene tree kernel from " + originalKernelSize + " to " + proposedKernelSize);
		
		
		// Add a tree
		if (expanding) {
			
			//if (true) return Double.NEGATIVE_INFINITY;

			int proposedTreeIndex = proposedKernelSize - 1;
			
			// Sample a tree to copy to the new one
			int kernelTreeToCopyIndex = Randomizer.nextInt(originalKernelSize);
			Tree treeToCopy = this.geneTreeKernelPrior.getKernel().getTree(kernelTreeToCopyIndex);
			if (treeToCopy.getRoot().getNodeCount() != 51) {
				int x = 5;
				int y = x;
				
			}
			GTKGeneTree proposedKernelTree = new GTKGeneTree(this.geneTreeKernelPrior.getKernel().getTree(kernelTreeToCopyIndex).getRoot().copy());
					
			
			
			// Add the new tree to the kernel
			this.geneTreeKernelPrior.addGeneTreeDistribution(proposedKernelTree);
			this.geneTreeKernelPrior.getKernel().addTree(proposedKernelTree);
			
			
			// Update the range of the indicator
			this.indicator.setUpper(proposedKernelSize - 1);
			
			
			// Sample pointers and allocate them to the new tree. 
			// The number of pointers sampled follows a Poisson(poissonRate) distribution
			for (GTKPointerTree pointer : this.pointers) {
				
				if (Randomizer.nextExponential(poissonRate) < 1) {
					
					if (DEBUG) Log.warning("Changing pointer from " + this.indicator.getValue(pointer.getTreeIndex()) + " to " + proposedTreeIndex);
					
					this.indicator.setValue(pointer.getTreeIndex(), proposedTreeIndex);
					numPointersReassigned ++;
					
					
					
					
				}
				
			}
			
			
			
			
			
			// Probability of selecting a given tree to add / delete
			log_q_expand += -Math.log(originalKernelSize);
			log_q_contract += -Math.log(proposedKernelSize);
			
			
			// Probability of reallocating pointers due to their tree being deleted (= 'prob assign pointer to gene' ^ 'num pointers reassigned')
			log_q_contract += -numPointersReassigned*Math.log(originalKernelSize);
			
		}

		// Delete a tree
		else {
			
			
			//if (true) return Double.NEGATIVE_INFINITY;
			
			// Sample a tree to delete
			int kernelIndexToDelete = Randomizer.nextInt(originalKernelSize);
			if (DEBUG) Log.warning("Deleting tree " + kernelIndexToDelete);
			
			// Find all genes which are pointing to this
			for (GTKPointerTree pointer : this.pointers) {

				int pointerIndicator = pointer.getIndicatorValue();
				if (pointerIndicator >= kernelIndexToDelete) {
					
					int proposedGTKIndex;
					
					// If the gene is pointing to the deleted tree, randomly reassign it to another tree
					if (pointerIndicator == kernelIndexToDelete) {
						proposedGTKIndex = (kernelIndexToDelete + 1 + Randomizer.nextInt(originalKernelSize - 1)) % originalKernelSize;
						
						// Decrement the indicator if it is greater than the one being deleted
						if (proposedGTKIndex > kernelIndexToDelete) proposedGTKIndex --;
						numPointersReassigned ++;
					}
					
					// If the gene is pointing to an index that is greater than the deleted tree, then decrement its index
					else {
						proposedGTKIndex = pointerIndicator - 1;
					}
					
					this.indicator.setValue(pointer.getTreeIndex(), proposedGTKIndex);
					
					
					if (DEBUG) Log.warning("Changing pointer from " + pointerIndicator + " to " + proposedGTKIndex);
					
				}

				
			}
			
			

			// Delete the tree from the kernel
			this.geneTreeKernelPrior.removeGeneTreeDistribution(kernelIndexToDelete);
			this.geneTreeKernelPrior.getKernel().deleteTree(kernelIndexToDelete);
			
			// Update the range of the indicator
			this.indicator.setUpper(proposedKernelSize - 1);
			
			

			// Probability of selecting a given tree to add / delete
			log_q_expand += -Math.log(proposedKernelSize);
			log_q_contract += -Math.log(originalKernelSize);
			
			
			// Probability of reallocating pointers due to their tree being deleted (= 'prob assign pointer to gene' ^ 'num pointers reassigned')
			log_q_contract += -numPointersReassigned*Math.log(proposedKernelSize);
			
			
			
			
			
			
		}
		

		
		
		// Update gene kernel size parameter
		this.geneKernelSize.setValue(proposedKernelSize);
		
		
		// Probability of sampling this many trees to add, under the Poisson(poissonRate) distribution
		log_q_expand += numPointersReassigned*Math.log(poissonRate) - poissonRate - MathUtils.factorialLog(numPointersReassigned);
		
		
		
		double logHR = expanding ? log_q_contract-log_q_expand : log_q_expand-log_q_contract;
		//if (DEBUG) System.out.println(this.getClass() + " Hastings ratio: " + logHR);
		
		
		// Return hastings ratio
		return logHR;
		
	}
	
	
	@Override
    public void accept() {
		if (DEBUG) System.out.println(this.getClass() + " accepted");
        super.accept();
    }

	@Override
	public void reject(final int reason) {
		if (DEBUG) System.out.println(this.getClass() + " rejected " + reason);
		super.reject(reason);
    }


	

	
	
	
	
}
