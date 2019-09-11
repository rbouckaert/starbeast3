/**
 * @author Huw A. Ogilvie
 * Adapted for starbeast3 by Jordan Douglas
 */


package starbeast3.evolution.branchratemodel;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.BooleanParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.evolution.branchratemodel.BranchRateModel;


@Description("Based on RandomLocalClockModel in BEAST v2.4")
public class RandomLocalClockModelSB3 extends BranchRateModel.Base implements BranchRateModelSB3 {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "(Species) tree to apply per-branch rates to.", Input.Validate.REQUIRED);
    final public Input<Boolean> noCacheInput = new Input<>("noCache", "Always recalculate branch rates.", false);
    final public Input<RealParameter> branchRatesInput = new Input<>("rates", "Continuous per-branch rates.", Input.Validate.REQUIRED);
    final public Input<BooleanParameter> indicatorsInput = new Input<>("indicators", "Indicators associated with nodes in the tree for sampling of individual rate changes among branches.", Input.Validate.REQUIRED);

    private boolean noCache;
    private boolean needsUpdate;
    private int rootNodeNumber;
    private double[] ratesArray;
    private double[] storedRatesArray;
    private double strictRatesTotal;
    private double relaxedRatesTotal;

    @Override
    protected boolean requiresRecalculation() {
        needsUpdate = true;
        return needsUpdate;
    }

    @Override
    public void store() {
        System.arraycopy(ratesArray, 0, storedRatesArray, 0, ratesArray.length);
        super.store();
    }

    @Override
    public void restore() {
        double[] tmpRatesArray = ratesArray;
        ratesArray = storedRatesArray;
        storedRatesArray = tmpRatesArray;
        super.restore();
    }

    @Override
    public void initAndValidate() {
        BooleanParameter indicators = indicatorsInput.get();
        TreeInterface tree = treeInput.get();
        int nodeCount = tree.getNodeCount();

        rootNodeNumber = nodeCount - 1;
        ratesArray = new double[nodeCount];
        storedRatesArray = new double[nodeCount];

        noCache = noCacheInput.get().booleanValue();

        RealParameter rates = branchRatesInput.get();
        rates.setDimension(rootNodeNumber);
        indicators.setDimension(rootNodeNumber);

        if (rates.lowerValueInput.get() == null || rates.lowerValueInput.get() < 0.0) {
            rates.setLower(0.0);
        }
        if (rates.upperValueInput.get() == null || rates.upperValueInput.get() < 0.0) {
            rates.setUpper(Double.MAX_VALUE);
        }
        
        needsUpdate = true;
    }

    /**
     * This is a recursive function that does the work of
     * calculating the unscaled branch rates across the tree
     * taking into account the indicator variables.
     *
     * @param node the node
     * @param rate the rate of the parent node
     */
    private void recurseBranchRates(Node node, double parentHeight, double rate, Boolean[] indicators, Double[] branchRates) {
        final int nodeNumber = node.getNr();
        final double nodeHeight = node.getHeight();

        // not the root, and indicator is "on"
        if (nodeNumber < rootNodeNumber && indicators[nodeNumber]) {
            rate = branchRates[nodeNumber];
        }

        double branchLength = parentHeight - nodeHeight;
        strictRatesTotal += branchLength;
        relaxedRatesTotal += branchLength * rate;

        ratesArray[nodeNumber] = rate;

        if (!node.isLeaf()) {
            recurseBranchRates(node.getLeft(), nodeHeight, rate, indicators, branchRates);
            recurseBranchRates(node.getRight(), nodeHeight, rate, indicators, branchRates);
        }
    }

    private void update() {
        final Boolean[] indicators = indicatorsInput.get().getValues();
        final Double[] rates = branchRatesInput.get().getValues();
        final Node treeRoot = treeInput.get().getRoot();
        final double treeHeight = treeRoot.getHeight();

        double estimatedMean;
        final RealParameter estimatedMeanParameter = meanRateInput.get();
        if (estimatedMeanParameter == null) {
            estimatedMean = 1.0;
        } else {
            estimatedMean = estimatedMeanParameter.getValue();
        }

        ratesArray[rootNodeNumber] = estimatedMean;

        strictRatesTotal = 0.0;
        relaxedRatesTotal = 0.0;
        recurseBranchRates(treeRoot, treeHeight, 1.0, indicators, rates);

        // normalize the weighted average of branch rates to equal the mean rate parameter
        double scaleFactor = estimatedMean * strictRatesTotal / relaxedRatesTotal;
        for (int i = 0; i < rootNodeNumber; i++) {
            ratesArray[i] = ratesArray[i] * scaleFactor; 
        }
    }

    @Override
    public double[] getRatesArray() {
        synchronized (this) {
        	if (needsUpdate || noCache) {
                update();
                needsUpdate = false;
            }
        }

        return ratesArray;
    }

    @Override
    public double getRateForBranch(Node node) {
        synchronized (this) {
        	if (needsUpdate || noCache) {
                update();
                needsUpdate = false;
            }
        }

        assert ratesArray[node.getNr()] > 0.0;
        return ratesArray[node.getNr()];
    }
}
