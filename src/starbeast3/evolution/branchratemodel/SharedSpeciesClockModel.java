package starbeast3.evolution.branchratemodel;

import beast.base.inference.CalculationNode;
import beast.base.core.Input;

public class SharedSpeciesClockModel extends CalculationNode {
	
	final public Input<BranchRateModelSB3> branchRateModelInput = new Input<>("branchRateModel", "Species tree branch rate model", Input.Validate.REQUIRED);
	

	

	@Override
	public void initAndValidate() {
		// TODO Auto-generated method stub
		
	}
	
	
	public BranchRateModelSB3 getClockModel() {
		return branchRateModelInput.get();
	}






	

}
