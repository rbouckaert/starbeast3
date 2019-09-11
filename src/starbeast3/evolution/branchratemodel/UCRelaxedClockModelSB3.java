/**
 * @author Huw A. Ogilvie
 * Adapted for starbeast3 by Jordan Douglas
 */



package starbeast3.evolution.branchratemodel;

import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ExponentialDistribution;
import org.apache.commons.math.distribution.ExponentialDistributionImpl;
import org.apache.commons.math.distribution.NormalDistribution;
import org.apache.commons.math.distribution.NormalDistributionImpl;

import beast.core.Input;
import beast.core.parameter.IntegerParameter;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.Node;
import beast.evolution.tree.TreeInterface;
import beast.evolution.branchratemodel.BranchRateModel;

public class UCRelaxedClockModelSB3 extends BranchRateModel.Base implements BranchRateModelSB3 {
    final public Input<TreeInterface> treeInput = new Input<>("tree", "(Species) tree to apply per-branch rates to.", Input.Validate.REQUIRED);
    final public Input<Integer> nBinsInput = new Input<>("nBins", "Number of discrete branch rate bins (default is equal to the number of estimated branch rates).", -1);
    final public Input<Boolean> estimateRootInput = new Input<>("estimateRoot", "Estimate rate of the root branch.", false);
    final public Input<Boolean> noCacheInput = new Input<>("noCache", "Always recalculate branch rates.", false);
    final public Input<RealParameter> stdevInput = new Input<>("stdev", "Standard deviation of the log-normal distribution for branch rates. If not supplied uses exponential.");
    final public Input<IntegerParameter> branchRatesInput = new Input<>("rates", "Discrete per-branch rates.", Input.Validate.REQUIRED);

    private int nBins;
    private double currentLogNormalStdev;
    private double storedLogNormalStdev;
    private double[] binRates;
    private double[] storedBinRates;
    private double[] ratesArray;
    private double[] storedRatesArray;

    private int nEstimatedRates;
    private int rootNodeNumber;
    private boolean estimateRoot;
    private boolean noCache;
    private boolean needsUpdate;
    private boolean binRatesNeedsUpdate;
    private boolean useLogNormal;

    @Override
    public boolean requiresRecalculation() {
        if (useLogNormal) {
            final double proposedLogNormalStdev = stdevInput.get().getValue();
            if (proposedLogNormalStdev != currentLogNormalStdev) {
                binRatesNeedsUpdate = true;
            } else {
                binRatesNeedsUpdate = false;
            }
        }

        needsUpdate = binRatesNeedsUpdate || branchRatesInput.isDirty() || meanRateInput.isDirty();
        return needsUpdate;
    }

    @Override
    public void store() {
        storedLogNormalStdev = currentLogNormalStdev;
        System.arraycopy(binRates, 0, storedBinRates, 0, binRates.length);
        System.arraycopy(ratesArray, 0, storedRatesArray, 0, ratesArray.length);
        super.store();
    }

    @Override
    public void restore() {
        double tmpLogNormalStdev = currentLogNormalStdev;
        double[] tmpBinRates = binRates;
        double[] tmpRatesArray = ratesArray;

        currentLogNormalStdev = storedLogNormalStdev;
        binRates = storedBinRates;
        ratesArray = storedRatesArray;
        
        storedLogNormalStdev = tmpLogNormalStdev;
        storedBinRates = tmpBinRates;
        storedRatesArray = tmpRatesArray;

        super.restore();
    }

    @Override
    public void initAndValidate() {
        final IntegerParameter branchRates = branchRatesInput.get();
        final TreeInterface speciesTree = treeInput.get();
        final Node[] speciesNodes = speciesTree.getNodesAsArray();
        estimateRoot = estimateRootInput.get().booleanValue();
        noCache = noCacheInput.get().booleanValue();
        rootNodeNumber = speciesTree.getRoot().getNr();
        ratesArray = new double[speciesNodes.length];
        storedRatesArray = new double[speciesNodes.length];

        if (estimateRoot) {
            nEstimatedRates = speciesNodes.length;
        } else {
            nEstimatedRates = speciesNodes.length - 1;
        }

        final int nBinsSupplied = nBinsInput.get().intValue();
        nBins = (nBinsSupplied <= 0) ? nEstimatedRates : nBinsSupplied;

        branchRates.setDimension(nEstimatedRates);
        branchRates.setLower(0);
        branchRates.setUpper(nBins - 1);

        currentLogNormalStdev = -1.0;
        storedLogNormalStdev = -1.0;

        binRates = new double[nBins];
        storedBinRates = new double[nBins];

        if (stdevInput.get() == null) {
            useLogNormal = false;

            final ExponentialDistribution exponentialDistr = new ExponentialDistributionImpl(1.0);
            try {
                for (int i = 0; i < nBins; i++) {
                    binRates[i] = exponentialDistr.inverseCumulativeProbability((i + 0.5) / nBins);
                }
            } catch (MathException e) {
                throw new RuntimeException("Failed to compute inverse cumulative probability!");
            }

            binRatesNeedsUpdate = false;
        } else {
            useLogNormal = true;
            binRatesNeedsUpdate = true;
        }

        needsUpdate = true;
    }

    private void update() {
        if (useLogNormal && (binRatesNeedsUpdate || noCache)) {
            // set the mean in real space to equal 1
            currentLogNormalStdev = stdevInput.get().getValue();
            final double newMean = -(0.5 * currentLogNormalStdev * currentLogNormalStdev);
            final NormalDistribution normalDistr = new NormalDistributionImpl(newMean, currentLogNormalStdev);

            try {
                for (int i = 0; i < nBins; i++) {
                    binRates[i] = Math.exp(normalDistr.inverseCumulativeProbability((i + 0.5) / nBins));
                }
            } catch (MathException e) {
                throw new RuntimeException("Failed to compute inverse cumulative probability!");
            }
        }

        Double estimatedMean;
        final RealParameter estimatedMeanParameter = meanRateInput.get();
        if (estimatedMeanParameter == null) {
            estimatedMean = 1.0;
        } else {
            estimatedMean = estimatedMeanParameter.getValue();
        }

        final Integer[] branchRatePointers = branchRatesInput.get().getValues();
        for (int i = 0; i < nEstimatedRates; i++) {
            int b = branchRatePointers[i];
            ratesArray[i] = binRates[b] * estimatedMean;
        }

        if (!estimateRoot) ratesArray[rootNodeNumber] = estimatedMean;

        /* StringBuffer x = new StringBuffer();
        x.append(treeInput.get().getID());
        for (int i = 0; i < ratesArray.length; i++) {
            x.append(" ");
            x.append(ratesArray[i]);
        }
        System.out.println(x); */
    }

    @Override
    public double[] getRatesArray() {
        if (needsUpdate || noCache) {
            synchronized (this) {
                update();
                needsUpdate = false;
            }
        }

        return ratesArray;
    }

    @Override
    public double getRateForBranch(Node node) {
        if (needsUpdate || noCache) {
            synchronized (this) {
                update();
                needsUpdate = false;
            }
        }

        assert ratesArray[node.getNr()] > 0.0;
        return ratesArray[node.getNr()];
    }
}
