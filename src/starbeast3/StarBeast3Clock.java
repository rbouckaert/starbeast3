

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
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import starbeast3.evolution.branchratemodel.BranchRateModelSB3;

public class StarBeast3Clock extends BranchRateModel.Base {
    public Input<GeneTreeForSpeciesTreeDistribution> geneTreeInput = new Input<>("geneTree", "The gene tree this relaxed clock is associated with.", Input.Validate.REQUIRED);
    public Input<BranchRateModelSB3> speciesTreeRatesInput = new Input<>("speciesTreeRates", "The per-branch rates for the species tree", Input.Validate.REQUIRED);

    private int geneNodeCount;
    private double[] branchRates;
    private double[] storedBranchRates;
    private boolean needsUpdate;

    RealParameter meanRate;
    BranchRateModelSB3 speciesTreeRatesX;
    GeneTreeForSpeciesTreeDistribution geneTree;
    
    @Override
    public void initAndValidate() {
        meanRate = meanRateInput.get();
        speciesTreeRatesX = speciesTreeRatesInput.get();
        geneTree = geneTreeInput.get();
    
        geneNodeCount = geneTree.getNodeCount();
        branchRates = new double[geneNodeCount];
        storedBranchRates = new double[geneNodeCount];
        needsUpdate = true;
    }

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = geneTreeInput.isDirty() || speciesTreeRatesInput.isDirty() || meanRateInput.isDirty();
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
