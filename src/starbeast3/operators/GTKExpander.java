package starbeast3.operators;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.distribution.PoissonDistribution;
import org.apache.commons.math.distribution.PoissonDistributionImpl;
import org.apache.commons.math.util.MathUtils;

import beast.base.core.Input;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.core.Log;
import beast.base.evolution.tree.Tree;
import beast.base.inference.distribution.Poisson;
import beast.base.util.Randomizer;
import starbeast3.genekernel.GTKGeneTree;
import starbeast3.genekernel.GTKOperator;
import starbeast3.genekernel.GTKPointerTree;


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
	
	final public Input<Double> pInput = 
			new Input<>("p", "The number of trees pointers to change is Binomial distributed, where n = 'num pointers' and p = min('numpointers / kernelsize * p', 'max')", 1.0);

	final public Input<Double> maxInput = 
			new Input<>("max", "The number of trees pointers to change is Binomial distributed, where n = 'num pointers' and p = min('numpointers / kernelsize * p', 'max')", 0.95);

	
	
	boolean DEBUG = false;
	IntegerParameter indicator;
	IntegerParameter geneKernelSize;
	double p;
	double max;
	
	public GTKExpander() {
		geneTreesInput.setRule(Input.Validate.FORBIDDEN);
	}
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
		this.indicator = indicatorInput.get();
		this.geneKernelSize = geneKernelSizeInput.get();
		this.p = pInput.get();
		this.max = maxInput.get();
		
	}
	
	
	@Override
	public double proposal() {
		
		//geneTreeDistributions = this.getTreeDistributions(this);
		
		double log_q_expand = 0;
		double log_q_contract = 0;
		int originalKernelSize = this.geneTreeKernelPrior.getKernel().getDimension();
		
		// Sanity check
		if (originalKernelSize != this.geneKernelSize.getValue()) {
			return Double.NEGATIVE_INFINITY;
		}
		
		
		// Expanding or contracting? 
		boolean expanding = Randomizer.nextBoolean();
		if (originalKernelSize == this.geneKernelSize.getLower() && !expanding) return Double.NEGATIVE_INFINITY;
		else if (originalKernelSize == this.geneKernelSize.getUpper() & expanding) return Double.NEGATIVE_INFINITY;
		
		
		// Calculations
		int proposedKernelSize = originalKernelSize + (expanding ? 1 : -1);
		double binomialP = Math.min(this.max, this.p * this.indicator.getDimension() / (expanding ? originalKernelSize : proposedKernelSize));
		int numPointersReassigned = 0;
		
		//System.out.println(this.indicator.getDimension()  + " / " + (expanding ? originalKernelSize : proposedKernelSize) + " = " + binomialP);
		

		if (DEBUG) System.out.println((expanding ? "Expanding" : "Contracting") + " the gene tree kernel from " + originalKernelSize + " to " + proposedKernelSize);
		
		
		// Add a tree
		if (expanding) {
			
			if (true) return Double.NEGATIVE_INFINITY;

			int proposedTreeIndex = proposedKernelSize - 1;
			
			// Sample a tree to copy to the new one
			int kernelTreeToCopyIndex = Randomizer.nextInt(originalKernelSize);

			GTKGeneTree proposedKernelTree = new GTKGeneTree(this.geneTreeKernelPrior.getKernel().getTree(kernelTreeToCopyIndex).getRoot().copy());
					
			
			
			// Add the new tree to the kernel
			this.geneTreeKernelPrior.addGeneTreeDistribution(proposedKernelTree);
			this.geneTreeKernelPrior.getKernel().addTree(proposedKernelTree);
			
			
			// Update the range of the indicator
			this.indicator.setUpper(proposedKernelSize - 1);
			
			
			// Sample pointers and allocate them to the new tree. 
			// The number of pointers sampled follows a Poisson(poissonRate) distribution
			for (int treeIndex = 0; treeIndex < this.indicator.getDimension(); treeIndex ++) {
				
				if (Randomizer.nextFloat() < binomialP) {
					
					if (DEBUG) Log.warning("Changing pointer from " + this.indicator.getValue(treeIndex) + " to " + proposedTreeIndex);
					
					this.indicator.setValue(treeIndex, proposedTreeIndex);
					numPointersReassigned ++;
					
				}
				
			}
			
			
			// Probability of selecting a given kernel tree to add/delete
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
			for (int treeIndex = 0; treeIndex < this.indicator.getDimension(); treeIndex ++) {

				int pointerIndicator = this.indicator.getValue(treeIndex);
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
					
					this.indicator.setValue(treeIndex, proposedGTKIndex);
					
					
					if (DEBUG) Log.warning("Changing pointer from " + pointerIndicator + " to " + proposedGTKIndex);
					
				}

				
			}

			// Delete the tree from the kernel
			this.geneTreeKernelPrior.removeGeneTreeDistribution(kernelIndexToDelete);
			this.geneTreeKernelPrior.getKernel().deleteTree(kernelIndexToDelete);
			
			// Update the range of the indicator
			this.indicator.setUpper(proposedKernelSize - 1);

			// Probability of selecting a given kernel tree to add/delete
			log_q_expand += -Math.log(proposedKernelSize);
			log_q_contract += -Math.log(originalKernelSize);
			
			
			// Probability of reallocating pointers due to their tree being deleted (= 'prob assign pointer to kernel gene' ^ 'num pointers reassigned')
			 // (= 1 / num genes in shrunk kernel) ^ num pointers reassigned
			log_q_contract += -numPointersReassigned*Math.log(proposedKernelSize);
			
			
		}
		
		
		// Update gene kernel size parameter
		this.geneKernelSize.setValue(proposedKernelSize);
		
		
		// Probability of selecting this set of pointers to reassign (upon expansion) as Bernoulli trials
		log_q_expand += numPointersReassigned*Math.log(binomialP) + (this.indicator.getDimension() - numPointersReassigned)*Math.log(1 - binomialP);
		
		
		
		double logHR = expanding ? log_q_contract-log_q_expand : log_q_expand-log_q_contract;
		if (DEBUG) Log.warning(this.getClass() + " Hastings ratio: " + logHR);
		
		
		// Return hastings ratio
		return -100000000; // Double.NEGATIVE_INFINITY;
		//return logHR;
		
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
