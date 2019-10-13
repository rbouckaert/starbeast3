package starbeast3.operators;

import beast.core.Description;

@Description("Operaor that does proposals that count for one step or more steps in the MCMC")
public interface MultiStepOperator {
	
	/** number of steps to be performed by operator **/
	public default int stepCount() {
		return 1;
	}
}
