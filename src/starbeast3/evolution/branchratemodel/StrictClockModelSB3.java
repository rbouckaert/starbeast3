package starbeast3.evolution.branchratemodel;

import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.StateNode;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeInterface;

public class StrictClockModelSB3 extends BranchRateModel.Base implements BranchRateModelSB3 {

	
	final public Input<TreeInterface> treeInput = new Input<>("tree", "(Species) tree to apply per-branch rates to.", Input.Validate.REQUIRED);
	
	Function muParameter;
	private double mu = 1.0;
	private int numSpecies;
	double[] ratesArray;

    @Override
    public void initAndValidate() {
        muParameter = meanRateInput.get();
        if (muParameter != null) {
        	if (muParameter instanceof RealParameter) {
        		((RealParameter)muParameter).setBounds(Math.max(0.0, ((RealParameter)muParameter).getLower()), ((RealParameter)muParameter).getUpper());
        	}
            
            mu = muParameter.getArrayValue();
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
    	mu = muParameter.getArrayValue();
    	if (muParameter instanceof StateNode) {
    		return ((StateNode)muParameter).isDirtyCalculation();
    	}
    	return false;
    	
    }
    
    @Override
    protected void restore() {
    	if (muParameter != null) mu = muParameter.getArrayValue();
        super.restore();
    }

    @Override
    protected void store() {
    	if (muParameter != null)  mu = muParameter.getArrayValue();
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
