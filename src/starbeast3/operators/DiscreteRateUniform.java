package starbeast3.operators;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.Input.Validate;
import beast.util.Randomizer;

/**
 *
 * @author Huw A. Ogilvie
 * 
 */

public class DiscreteRateUniform extends AdaptiveOperator {
    final public Input<IntegerParameter> treeRatesInput = new Input<>("treeRates", "The branch rates.", Validate.REQUIRED);

    private int nNodes;
    private int lowerBound;
    private int integerRange;

    @Override
    public void initAndValidate() {
        final IntegerParameter treeRates = treeRatesInput.get();
        nNodes = treeRates.getDimension();
        lowerBound = treeRates.getLower();
        integerRange = 1 + treeRates.getUpper() - lowerBound;

        setLimits(1, nNodes);
        super.initAndValidate();
    }

    // symmetric proposal distribution
    @Override
    public double proposal() {
        final IntegerParameter treeRates = treeRatesInput.get();
        final int[] cycle = chooseK(nNodes);

        for (int i = 0; i < discreteK; i++) {
            final int nodeNumber = cycle[i];
            final int newRate = lowerBound + Randomizer.nextInt(integerRange);
            treeRates.setValue(nodeNumber, newRate);
        }

        return 0.0;
    }
}
