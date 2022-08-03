package starbeast3.core;

import beast.base.inference.OperatorSchedule;

public class OperatorScheduleRecalculator extends OperatorSchedule {

	
	/**
	 * Recompute all operator weights
	 */
	public void reweight() {
		this.reweightOperators();
	}
	
	
}
