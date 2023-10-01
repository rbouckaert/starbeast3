package starbeast3.evolution.speciation;


import java.text.DecimalFormat;
import java.util.Arrays;

import beast.base.inference.CalculationNode;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import starbeast3.tree.SpeciesTree;
import beast.base.evolution.tree.Node;

/**
 * 
 * @author jdou557
 *
 */
@Description("Gene lineages are assumed to diverge under a speciation process, constrained within species tree")
public class BirthProcess extends CalculationNode implements PopulationModel {
	
	
    final public Input<SpeciesTree> speciesTreeInput =
            new Input<>("speciesTree", "The species tree this population model is assoicated with.", Validate.REQUIRED);
   
	
    public Input<RealParameter> birthRateinput = new Input<RealParameter>("birthRate", "Per-branch birth rate.", Validate.REQUIRED);

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
        birthRateinput.get().setDimension(speciesNodeCount);
        speciesBranchStatus = new boolean[speciesNodeCount];
        needsUpdate = true;
    }

    
    @Override
    public double getRootBranchLogP(double rootHeight, double lambda, int ntaxa) {
    	
    	final double mrh = -lambda * rootHeight;
    	final double treeProb = (ntaxa - 1) * Math.log(lambda);
    	return mrh + treeProb;
    }
    
    
    /* The contribution of a branch in the species tree to 
     * the log probability, for constant population function.
     */
    @Override
    public double calculateBranchLogP(final int lineagesBottom, final double ploidy, final double lambda, final double[] times, final int k) {
    	
    	
    	
    	 double logPBranch = 0.0;
    	 for (int i = 0; i <= k; i++) {
    		 
    		 final double height = times[i];
             final double mrh = -lambda * height;
             double l = mrh;
             //if (node.isRoot()) l += mrh;
             logPBranch += l;
    		 
    	 }
    	
         return logPBranch;

    }
    
    
	@Override
	public double calculatePartialLogPBranch(int lineagesBottom, double[] times, int k) {
		throw new IllegalArgumentException("calculatePartialLogPBranch: not implemented");
	}


    @Override
    public void initPopSizes(double popInitial) {
        final RealParameter birthRates = birthRateinput.get();
        final double lower = birthRates.getLower();
        final double upper = birthRates.getUpper();

        if (birthRates.isEstimatedInput.get() && popInitial > lower && popInitial < upper) {
	        for (int i = 0; i < birthRates.getDimension(); i++) {
	            birthRates.setValue(i, popInitial);
	        }
        }
    }

    @Override
    public void serialize(Node speciesTreeNode, StringBuffer buf, DecimalFormat df) {
        final RealParameter birthRates = birthRateinput.get();
        final int speciesTreeNodeNumber = speciesTreeNode.getNr();
        final double branchPopSize = birthRates.getValue(speciesTreeNodeNumber);

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
            final RealParameter birthRates = birthRateinput.get();

            Arrays.fill(speciesBranchStatus, false);

            for (int nodeI = 0; nodeI < speciesNodeCount; nodeI++)
                speciesBranchStatus[nodeI] = birthRates.isDirty(nodeI);

            needsUpdate = false;
        }

        return speciesBranchStatus[speciesNode.getNr()];
    }




}



