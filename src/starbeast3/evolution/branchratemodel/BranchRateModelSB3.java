
/**
 * @author Huw A. Ogilvie
 * Adapted for starbeast3 by Jordan Douglas
 */

package starbeast3.evolution.branchratemodel;

import beast.core.Input;
import beast.evolution.tree.TreeInterface;

public interface BranchRateModelSB3 {
	
	final public Input<TreeInterface> treeInput = new Input<>("tree", "(Species) tree to apply per-branch rates to.", Input.Validate.REQUIRED);
    double[] getRatesArray();
}
