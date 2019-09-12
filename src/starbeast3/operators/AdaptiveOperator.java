package starbeast3.operators;

import beast.core.Input;
import beast.core.Operator;
import beast.util.Randomizer;


/**
 * Common methods for executing the same operator multiple times in one step.
 * Useful for cases where acceptance ratios are way too high.
 *
 * @author Huw A. Ogilvie
 */
public abstract class AdaptiveOperator extends Operator {
    public final Input<Integer> kInput = new Input<>("k", "Number of operations to perform per step.", 1);
    public final Input<Boolean> optimiseInput = new Input<>("optimise", "Adjust 'k' during the MCMC run to improve mixing.", true);

    protected boolean optimise;
    protected int discreteK;
    protected double continuousK;

    protected double lower;
    protected double upper;

    @Override
    public void initAndValidate() {
        setCoercableParameterValue(kInput.get());
        optimise = optimiseInput.get();
    }

    @Override
    public double getCoercableParameterValue() {
        return continuousK;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        continuousK = value;
        if (continuousK < lower) continuousK = lower;
        else if (continuousK > upper) continuousK = upper;
        discreteK = (int) Math.round(continuousK);
    }

    @Override
    public void optimize(final double logAlpha) {
        if (optimise) {
            final double currentK = continuousK;
            final double delta = calcDelta(logAlpha);
            final double logK = delta + Math.log(currentK);
            setCoercableParameterValue(Math.exp(logK));
            // System.out.println(String.format("%f = exp(%f) = exp(cd(%g) + ln(%f)) = exp(%f + ln(%f))", continuousK, logK, logAlpha, currentK, delta, currentK));
        }
    }

    @Override
    public final String getPerformanceSuggestion() {
        final double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        final double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        final double newContinuousK = continuousK * ratio;
        final int newDiscreteK = (int) Math.round(newContinuousK); 

        if ((newDiscreteK != discreteK) && (newDiscreteK >= lower) && (newDiscreteK <= upper) && (prob < 0.10 || prob > 0.40)) {
            return String.format("k = %d... Try setting it to about %d", discreteK, newDiscreteK);
        } else {
            return String.format("k = %d", discreteK);
        }
    }

    protected void setLimits(final int lower, final int upper) {
        this.lower = (double) lower;
        this.upper = (double) upper;
        this.lower -= 0.4999999999;
        this.upper += 0.4999999999;
    }

    // chooses 'k' numbers from between 0 and n - 1
    // i.e., random sampling without replacement
    // this is a Durstenfeld shuffle which terminates after k loops
    protected int[] chooseK(final int n) {
        final int[] sequential = new int[n];

        for (int i = 0; i < n; i++) sequential[i] = i;

        for (int i = 0; i < discreteK; i++) {
            final int j = Randomizer.nextInt(n - i);
            if (j > 0) { // swap [i] with [i + j]
                final int i_temp = sequential[i];
                final int iplusj = i + j;
                sequential[i] = sequential[iplusj];
                sequential[iplusj] = i_temp;
            }
        }

        return sequential;
    }
}
