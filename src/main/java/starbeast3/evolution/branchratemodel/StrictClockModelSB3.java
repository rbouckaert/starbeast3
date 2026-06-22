package starbeast3.evolution.branchratemodel;

import beast.base.core.Input;
import beast.base.inference.util.InputUtil;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;
import beast.base.spec.domain.PositiveReal;
import beast.base.spec.evolution.branchratemodel.Base;
import beast.base.spec.type.RealScalar;

public class StrictClockModelSB3 extends Base implements BranchRateModelSB3 {

	
	final public Input<TreeInterface> treeInput = new Input<>("tree", "(Species) tree to apply per-branch rates to.", Input.Validate.REQUIRED);
	
	RealScalar<PositiveReal> muParameter;
	private double mu = 1.0;
	private int numSpecies;
	double[] ratesArray;

    @Override
    public void initAndValidate() {
        muParameter = meanRateInput.get();
        if (muParameter != null) {
            mu = muParameter.get();
        }
        this.numSpecies = treeInput.get().getNodeCount();
        this.ratesArray = new double[this.numSpecies];
        
        System.out.println("Number of species " +  this.numSpecies);
    }

    @Override
    public double getRateForBranch(final Node node) {
        return mu;
    }

    @Override
    public boolean requiresRecalculation() {
    	if (muParameter == null) return false;
    	mu = muParameter.get();
    	return InputUtil.isDirty(meanRateInput);
    }
    
    @Override
    protected void restore() {
    	if (muParameter != null) mu = muParameter.get();
        super.restore();
    }

    @Override
    protected void store() {
    	if (muParameter != null)  mu = muParameter.get();
        super.store();
    }

    @Override
    public double[] getRatesArray() {
    	for (int i = 0; i < this.numSpecies; i ++) {
    		ratesArray[i] = mu;
    	}
    	return ratesArray;
    }
    

}
