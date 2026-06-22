

/**
 * Species tree relaxed clock model. The per-branch rates in the species tree are multiplied 
 * by the per-tree rates of each gene tree to get the overall rate of a gene tree branch region
 * 
 * @author Huw A. Ogilvie
 * 
 * Adapted by Jordan Douglas for starbeast3
 */

package starbeast3.evolution.branchratemodel;



import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beast.base.inference.util.InputUtil;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.evolution.branchratemodel.Base;
import beast.base.spec.type.RealScalar;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import starbeast3.evolution.speciation.GeneTreeForSpeciesTreeDistribution;

public class StarBeast3Clock extends Base {
	
    public Input<BranchRateModelSB3> speciesTreeRatesInput = new Input<>("speciesTreeRates", "The real branch rate model of the species tree");
    public Input<SharedSpeciesClockModel> sharedRateModelInput =  new Input<>("sharedRateModel", "Clock model for species tree (instead of speciesTreeRates)");
    
    
    public Input<Tree> treeInput =  new Input<>("tree", "the gene tree", Input.Validate.OPTIONAL);
    final public Input<GeneTreeForSpeciesTreeDistribution> geneTreeInput = new Input<>("geneTree", "The gene tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
  
	
    
    protected int geneNodeCount;
    protected double[] branchRates;
    protected double[] storedBranchRates;
    protected boolean needsUpdate;

    
    protected RealScalar<PositiveReal> meanRate;
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
        
        
        
        this.geneTree = geneTreeInput.get();
        
        // Get prior from tree?
        if (this.geneTree == null && this.treeInput.get() != null) {
        	
        	Tree tree = this.treeInput.get();
        	
        	// Get prior
    		GeneTreeForSpeciesTreeDistribution prior = null;
    		for (BEASTInterface obj : tree.getOutputs()){
    			if (obj instanceof GeneTreeForSpeciesTreeDistribution){
    				
    				if (prior != null && prior != obj){
    					throw new IllegalArgumentException("Found more than 1 'GeneTreeForSpeciesTreeDistribution' for" + tree.getID());
    				}
    				prior = (GeneTreeForSpeciesTreeDistribution) obj;
    			}
    		
    		}
    		
    		if (prior != null) {
    			geneTree = prior; 
    			geneTreeInput.set(prior);
    		}
        	
        }
        
        

        geneTree = this.getGeneTreePrior();
    
        geneNodeCount = geneTree.getNodeCount();
        branchRates = new double[geneNodeCount];
        storedBranchRates = new double[geneNodeCount];
        needsUpdate = true;
        
        
    }
    
    
    
    public GeneTreeForSpeciesTreeDistribution getGeneTreePrior() {
    	return this.geneTree;
    }
    

    @Override
    public boolean requiresRecalculation() {
    	
    	
    	if (speciesTreeRatesInput.get() != null && InputUtil.isDirty(speciesTreeRatesInput))  {
    		needsUpdate = true;
        }else if (sharedRateModelInput != null && InputUtil.isDirty(sharedRateModelInput)) {
        	needsUpdate = true;
        }
        needsUpdate = needsUpdate || InputUtil.isDirty(meanRateInput);
        needsUpdate = needsUpdate || InputUtil.isDirty(geneTreeInput);
        
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

    protected void update() {
    	geneTree = this.getGeneTreePrior();
        final double geneTreeRate = meanRate.get();
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





