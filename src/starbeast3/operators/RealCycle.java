package starbeast3.operators;

import beast.base.core.Input;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Input.Validate;

/**
 *
 * @author Huw A. Ogilvie
 * 
 */

public class RealCycle extends AdaptiveOperator {
    final public Input<RealParameter> parameterInput = new Input<>("parameter", "The branch rates.", Validate.REQUIRED);

    private int nNodes;

    @Override
    public void initAndValidate() {
        final RealParameter parameter = parameterInput.get();
        nNodes = parameter.getDimension();
        setLimits(2, nNodes);
        super.initAndValidate();
    }

    // symmetric proposal distribution
    @Override
    public double proposal() {
        final RealParameter parameter = parameterInput.get();
        final Double[] parameterArray = parameter.getValues();
        final int[] cycle = chooseK(nNodes);

        final int lastNodeNumber = cycle[discreteK - 1];
        double previousRate = parameterArray[lastNodeNumber];
        for (int i = 0; i < discreteK; i++) {
            final int nodeNumber = cycle[i];
            parameter.setValue(nodeNumber, previousRate);
            previousRate = parameterArray[nodeNumber];
        }

        return 0.0;
    }
}
