package starbeast3.math.distributions;


import org.apache.commons.math.MathException;
import org.apache.commons.math.distribution.ContinuousDistribution;
import org.apache.commons.math.distribution.Distribution;
import org.apache.commons.math.distribution.GammaDistributionImpl;

import beast.core.Description;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.math.distributions.ParametricDistribution;



@Description("Inverse Gamma distribution, used as prior.    for x>0  f(x; alpha, beta) = \frac{beta^alpha}{Gamma(alpha)} (1/x)^{alpha + 1}exp(-beta/x) " +
        "If the input x is a multidimensional parameter, each of the dimensions is considered as a " +
        "separate independent component.")
public class InverseGamma extends ParametricDistribution {
    final public Input<RealParameter> alphaInput = new Input<>("alpha", "shape parameter, defaults to 2");
    final public Input<RealParameter> betaInput = new Input<>("beta", "scale parameter, defaults to 2");

    InverseGammaImpl dist = new InverseGammaImpl(2, 2);
    GammaDistributionImpl gamma_dist = new GammaDistributionImpl(2, 1.0/2.0);

    @Override
    public void initAndValidate() {
        refresh();
    }

    /**
     * ensure internal state is up to date *
     */
    void refresh() {
        double alpha;
        double beta;
        if (alphaInput.get() == null) {
            alpha = 2;
        } else {
            alpha = alphaInput.get().getValue();
        }
        if (betaInput.get() == null) {
            beta = 2;
        } else {
            beta = betaInput.get().getValue();
        }
        dist.setAlphaBeta(alpha, beta);
        
        // If X ~ Gamma(shape=a, scale=1/b)
        // Then 1/X ~ InvGamma(shape=a, scale=b)
        gamma_dist.setAlpha(alpha);
        gamma_dist.setBeta(1/beta);
    }

    @Override
    public Distribution getDistribution() {
        refresh();
        return dist;
    }

    class InverseGammaImpl implements ContinuousDistribution {
        double m_fAlpha;
        double m_fBeta;
        // log of the constant beta^alpha/Gamma(alpha)
        double C;

        InverseGammaImpl(double alpha, double beta) {
            setAlphaBeta(alpha, beta);
        }

        void setAlphaBeta(double alpha, double beta) {
            m_fAlpha = alpha;
            m_fBeta = beta;
            C = m_fAlpha * Math.log(m_fBeta) - org.apache.commons.math.special.Gamma.logGamma(m_fAlpha);
        }

        @Override
        public double cumulativeProbability(double x) throws MathException {
            throw new MathException("Not implemented yet");
        }

        @Override
        public double cumulativeProbability(double x0, double x1) throws MathException {
            throw new MathException("Not implemented yet");
        }

        @Override
        public double inverseCumulativeProbability(double p) throws MathException {
            double offset = getOffset();
            double X =  gamma_dist.inverseCumulativeProbability(p);
            return offset + 1/X;
        }

        @Override
        public double density(double x) {
            double logP = logDensity(x);
            return Math.exp(logP);
        }

        @Override
        public double logDensity(double x) {
            double logP = -(m_fAlpha + 1.0) * Math.log(x) - (m_fBeta / x) + C;
            return logP;
        }
    } // class OneOnXImpl


    @Override
    public double getMeanWithoutOffset() {
    	RealParameter alpha = alphaInput.get();
    	RealParameter beta = betaInput.get();
    	return (beta != null ? beta.getArrayValue() : 2.0) / (alpha != null ? (alpha.getArrayValue() - 1.0) : 2.0);
    }
    
} // class InverseGamma
