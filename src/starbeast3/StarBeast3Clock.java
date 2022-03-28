

/**
 * Species tree relaxed clock model. The per-branch rates in the species tree are multiplied 
 * by the per-tree rates of each gene tree to get the overall rate of a gene tree branch region
 * 
 * @author Huw A. Ogilvie
 * 
 * Adapted by Jordan Douglas for starbeast3
 */

package starbeast3;


import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.util.Log;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import genekernel.GTKPointerTree;
import genekernel.GTKPrior;
import starbeast3.evolution.branchratemodel.BranchRateModelSB3;
import starbeast3.evolution.branchratemodel.SharedSpeciesClockModel;

public class StarBeast3Clock extends BranchRateModel.Base {
	
    public Input<BranchRateModelSB3> speciesTreeRatesInput = new Input<>("speciesTreeRates", "The real branch rate model of the species tree");
    public Input<SharedSpeciesClockModel> sharedRateModelInput =  new Input<>("sharedRateModel", "Clock model for species tree (instead of speciesTreeRates)");
    
    
    final public Input<GeneTreeForSpeciesTreeDistribution> geneTreeInput = new Input<>("geneTree", "The gene tree this relaxed clock is associated with.", Input.Validate.OPTIONAL);
	final public Input<GTKPrior> geneTreeKernelPriorInput = new Input<>("kernel", "the kernel of gene trees", Input.Validate.OPTIONAL);
	final public Input<GTKPointerTree> geneTreePointerInput = new Input<>("pointer", "the tree which points to the kernel", Input.Validate.OPTIONAL);
  
	
    
    protected int geneNodeCount;
    protected double[] branchRates;
    protected double[] storedBranchRates;
    private boolean needsUpdate;

    GTKPrior kernel;
    GTKPointerTree pointer;
    
    RealParameter meanRate;
    BranchRateModelSB3 speciesTreeRatesX;
    GeneTreeForSpeciesTreeDistribution geneTree;
    
    @Override
    public void initAndValidate() {
        meanRate = meanRateInput.get();
        
        if (speciesTreeRatesInput.get() != null) {
        	speciesTreeRatesX = speciesTreeRatesInput.get();
        }else if (sharedRateModelInput != null) {
        	speciesTreeRatesX = sharedRateModelInput.get().getClockModel();	
        }else {
        	speciesTreeRatesX = null;
        }
        
        
        
        // Must specify either A) a gene tree prior, or B) a gene tree kernel AND the pointer
        if (this.geneTreeInput.get() == null) {
        	if (geneTreeKernelPriorInput.get() == null || this.geneTreePointerInput.get() == null) {
        		Log.warning("StarBeast3Clock: Must specify either A) 'geneTree', or B) 'kernel' AND 'pointer', but not both A and B.");
        		return;
        	}
        } 
        else if (geneTreeKernelPriorInput.get() != null || this.geneTreePointerInput.get() != null) {
        	Log.warning("StarBeast3Clock: Must specify either A) 'geneTree', or B) 'kernel' AND 'pointer', but not both A and B.");
        	return;
        }
        
        this.kernel = geneTreeKernelPriorInput.get();
        this.pointer = geneTreePointerInput.get();
        geneTree = this.getGeneTreePrior();
    
        geneNodeCount = geneTree.getNodeCount();
        branchRates = new double[geneNodeCount];
        storedBranchRates = new double[geneNodeCount];
        needsUpdate = true;
    }
    
    
    
    public GeneTreeForSpeciesTreeDistribution getGeneTreePrior() {
    	
    	// If using a gene kernel, then the gene tree prior will change as the pointer changes
    	if (this.kernel != null) {
    		int indicatorOfPointer = this.pointer.getIndicatorValue();
    		return this.kernel.getGeneTreeDistributions(indicatorOfPointer);
    	}
    	
    	// Otherwise, use the one which was parsed
    	else {
    		return this.geneTreeInput.get();
    	}
    	
    }
    

    @Override
    public boolean requiresRecalculation() {
    	
    	
    	if (speciesTreeRatesInput.get() != null && speciesTreeRatesInput.isDirty()) {
    		needsUpdate = true;
        }else if (sharedRateModelInput != null && sharedRateModelInput.isDirty()) {
        	needsUpdate = true;
        }
        needsUpdate = needsUpdate || meanRateInput.isDirty();
        if (this.kernel == null) needsUpdate = needsUpdate || geneTreeInput.isDirty();
        else needsUpdate = needsUpdate || geneTreeKernelPriorInput.isDirty() || geneTreePointerInput.isDirty();
        
        return needsUpdate;
    }

    @Override
    public void store() {
        System.arraycopy(branchRates, 0, storedBranchRates, 0, branchRates.length);
        super.store();
    }

    @Override
    public void restore() {
        double[] tmpRatesArray = branchRates;
        branchRates = storedBranchRates;
        storedBranchRates = tmpRatesArray;
        super.restore();
    }

    private void update() {
    	geneTree = this.getGeneTreePrior();
        final double geneTreeRate = meanRate.getValue();
        final double[] speciesTreeRates = speciesTreeRatesX.getRatesArray();
        final double[] speciesOccupancy = geneTree.getSpeciesOccupancy();

        final int speciesNodeCount = speciesTreeRates.length;
        for (int i = 0; i < geneNodeCount - 1; i++) {
            double weightedSum = 0.0;
            double branchLength = 0.0;
            for (int j = 0; j < speciesNodeCount; j++) {
                // System.out.println(String.format("%d, %d: %f, %f", i, j, speciesTreeRates[j], speciesOccupancy[i * speciesNodeCount + j]));
                weightedSum += speciesTreeRates[j] * speciesOccupancy[i * speciesNodeCount + j];
                branchLength += speciesOccupancy[i * speciesNodeCount + j];
            }

            branchRates[i] = geneTreeRate * weightedSum / branchLength;
        }
        // set the rate for the root branch of this gene to equal the input mean rate
        branchRates[geneNodeCount - 1] = geneTreeRate;

        needsUpdate = false;
    }

    @Override
    public double getRateForBranch(Node node) {
        if (needsUpdate) {
            update();
        }
        return branchRates[node.getNr()];
    }
    
    
    
}





