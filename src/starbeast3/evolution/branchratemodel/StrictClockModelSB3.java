package starbeast3.evolution.branchratemodel;

import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;

public class StrictClockModelSB3 extends BranchRateModel.Base implements BranchRateModelSB3 {

	
	final public Input<TreeInterface> treeInput = new Input<>("tree", "(Species) tree to apply per-branch rates to.", Input.Validate.REQUIRED);
	
	RealParameter muParameter;
	private double mu = 1.0;
	private int numSpecies;
	double[] ratesArray;

    @Override
    public void initAndValidate() {
        muParameter = meanRateInput.get();
        if (muParameter != null) {
            muParameter.setBounds(Math.max(0.0, muParameter.getLower()), muParameter.getUpper());
            mu = muParameter.getValue();
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
    	mu = muParameter.getValue();
    	return muParameter.isDirtyCalculation();
    }
    
    @Override
    protected void restore() {
    	if (muParameter != null) mu = muParameter.getValue();
        super.restore();
    }

    @Override
    protected void store() {
    	if (muParameter != null)  mu = muParameter.getValue();
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
