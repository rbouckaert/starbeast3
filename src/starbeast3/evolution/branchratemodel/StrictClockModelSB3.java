package starbeast3.evolution.branchratemodel;

import beast.core.parameter.RealParameter;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.tree.Node;

public class StrictClockModelSB3 extends BranchRateModel.Base implements BranchRateModelSB3 {

	
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
        mu = muParameter.getValue();
        return true;
    }
    
    @Override
    protected void restore() {
        mu = muParameter.getValue();
        super.restore();
    }

    @Override
    protected void store() {
        mu = muParameter.getValue();
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
