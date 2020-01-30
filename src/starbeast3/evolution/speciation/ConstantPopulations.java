package starbeast3.evolution.speciation;


import java.text.DecimalFormat;
import java.util.Arrays;

import beast.core.CalculationNode;
import beast.core.Input;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import starbeast3.SpeciesTree;

/**
* @author Huw Ogilvie
 */

public class ConstantPopulations extends CalculationNode implements PopulationModel {
	
	
    final public Input<SpeciesTree> speciesTreeInput =
            new Input<>("speciesTree", "The species tree this population model is assoicated with.", Validate.REQUIRED);
	
    public Input<RealParameter> popSizesInput = new Input<RealParameter>("populationSizes", "Constant per-branch population sizes.", Validate.REQUIRED);

    private SpeciesTree speciesTree;

    private boolean needsUpdate;
    private boolean[] speciesBranchStatus;
    private int speciesNodeCount;

    @Override
    public boolean requiresRecalculation() {
        needsUpdate = true;
        return needsUpdate;
    }

    @Override
    public void initAndValidate() {
        speciesTree = speciesTreeInput.get();
        speciesNodeCount = speciesTree.getNodeCount();
        popSizesInput.get().setDimension(speciesNodeCount);
        speciesBranchStatus = new boolean[speciesNodeCount];
        needsUpdate = true;
    }

    
    /* The contribution of a branch in the species tree to 
     * the log probability, for constant population function.
     */
    @Override
    public double calculateBranchLogP(final int lineagesBottom, final double ploidy, final double popSize2, final double[] times, final int k) {
    	
    	
    	double logPBranch = 0.0;
        final double popSize = popSize2 * ploidy;
        logPBranch += -k * Math.log(popSize);
        for (int i = 0; i <= k; i++) {
        	logPBranch += -((lineagesBottom - i) * (lineagesBottom - i - 1.0) / 2.0) * (times[i + 1] - times[i]) / popSize;
        }
        return logPBranch;
    }
    
    
	@Override
	public double calculatePartialLogPBranch(int lineagesBottom, double[] times, int k) {
		double partialLogP = 0.0;
        for (int i = 0; i <= k; i++) {
        	partialLogP += -((lineagesBottom - i) * (lineagesBottom - i - 1.0) / 2.0) * (times[i + 1] - times[i]);
        }
        return partialLogP;
	}


    @Override
    public void initPopSizes(double popInitial) {
        final RealParameter popSizes = popSizesInput.get();
        final double lower = popSizes.getLower();
        final double upper = popSizes.getUpper();

        if (popSizes.isEstimatedInput.get() && popInitial > lower && popInitial < upper) {
	        for (int i = 0; i < popSizes.getDimension(); i++) {
	            popSizes.setValue(i, popInitial);
	        }
        }
    }

    @Override
    public void serialize(Node speciesTreeNode, StringBuffer buf, DecimalFormat df) {
        final RealParameter popSizes = popSizesInput.get();
        final int speciesTreeNodeNumber = speciesTreeNode.getNr();
        final double branchPopSize = popSizes.getValue(speciesTreeNodeNumber);

        buf.append("dmv={");
        if (df == null) {
            buf.append(branchPopSize);
        } else {
            buf.append(df.format(branchPopSize));
        }
        buf.append("}");
    }

    @Override
    public boolean isDirtyBranch(Node speciesNode) {
        if (needsUpdate) {
            final RealParameter popSizes = popSizesInput.get();

            Arrays.fill(speciesBranchStatus, false);

            for (int nodeI = 0; nodeI < speciesNodeCount; nodeI++)
                speciesBranchStatus[nodeI] = popSizes.isDirty(nodeI);

            needsUpdate = false;
        }

        return speciesBranchStatus[speciesNode.getNr()];
    }




}



