package starbeast3.operators;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.GammaDistribution;
import org.apache.commons.math3.exception.NotStrictlyPositiveException;

import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.inference.Operator;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.util.Randomizer;
import beast.base.core.Log;
import beast.base.evolution.tree.TreeInterface;
import beast.base.evolution.tree.TreeIntervals;
import beast.base.inference.distribution.Gamma;
import beast.base.inference.distribution.InverseGamma;
import beast.base.inference.distribution.ParametricDistribution;
import starbeast3.evolution.speciation.GeneTreeForSpeciesTreeDistribution;
import starbeast3.genekernel.GTKOperator;
import starbeast3.genekernel.GTKPrior;
import starbeast3.util.*;

/*
gibbs coalescent constant operator gamma prior
                                                coverage Mean ESS Min ESS
TreeHeight                                      95	   4331.22  3578.89
kappa                                           94	   4343.84  3371.71
gammaShape                                      95	   4368.21  3664.47
popSize                                         93	   2809.96  1047.55
CoalescentConstant                              97	   3842.19  1994.09
parameter.hyperInverseGamma-beta-PopSizePrior   91	   1270.48  725.18
HyperPrior.hyperInverseGamma-beta-PopSizePrior  91	   1595.78  1033.85
logP(mrca(root))                                98	   4357.48  3563.65
mrca.age(root)                                  95	   4331.22  3578.89
clockRate                                       0	   4274.31  3262.76
freqParameter.1                                 92	   4338.11  3434.69
freqParameter.2                                 93	   4312.99  2573.32
freqParameter.3                                 94	   4354.94  3273.80
freqParameter.4                                 92	   4337.21  3295.10


gibbs coalescent constant operator inverse gamma prior
                                                coverage Mean ESS Min ESS
kappa                                           86	   4330.76  3398.64
gammaShape                                      7	   4334.10  3685.89
TreeHeight                                      95	   3076.44  2233.29
popSize                                         94	   577.20  331.78
CoalescentConstant                              91	   1620.76  787.30
parameter.hyperInverseGamma-beta-PopSizePrior   94	   545.00  319.01
HyperPrior.hyperInverseGamma-beta-PopSizePrior  91	   556.71  331.29
logP(mrca(root))                                97	   4320.70  3328.88
mrca.age(root)                                  95	   3076.44  2233.29
clockRate                                       0	   3046.64  2174.60
freqParameter.1                                 98	   4332.76  3388.90
freqParameter.2                                 97	   4337.93  3334.29
freqParameter.3                                 96	   4378.30  3462.73
freqParameter.4                                 92	   4348.83  3316.36

 */
@Description("Samples population sizes of constant populations for multi species coalescent model")
public class PopSizeGibbsSampler extends GTKOperator {
	final public Input<RealParameter> popSizesInput = new Input<>("popSizes", "constant population size parameter, one dimension for each branch of the species tree", Validate.REQUIRED);
	final public Input<ParametricDistribution> priorInput = new Input<>("gammaprior", "gamma distributed prior for population sizes", Validate.REQUIRED);
	final public Input<TreeIntervals> treeIntervalsInput = new Input<>("intervals", "tree intervals for use with single tree -- should not be used if gene-attribute is used");
	
	
	RealParameter popSizes;
	Function priorAlpha, priorBeta;
	ParametricDistribution prior;
	TreeIntervals treeIntervals;

	// used to sample a gamma distribution
	MyRandomizer myRandomizer = new MyRandomizer();
	
	@Override
	public void initAndValidate() {
		popSizes = popSizesInput.get();
		prior = priorInput.get();
		treeIntervals = treeIntervalsInput.get();
		priorAlpha = (Function) prior.getInput("alpha").get();
		priorBeta = (Function) prior.getInput("beta").get();
		if (treeIntervals != null) { 
			return;
		}
		
		
		super.initAndValidate();
		if (geneTreeDistributions != null && geneTreeDistributions.size() > 0) {
		
			TreeInterface speciesTree = geneTreeDistributions.get(0).speciesTreeInput.get();
			if (speciesTree.getNodeCount() != popSizes.getDimension()) {
				throw new IllegalArgumentException("The dimension of the population size parameter (" + popSizes.getDimension()+ ") "
						+ "should be equal to the number of branches in the species tree ( " + speciesTree.getNodeCount() + ")");
			}
			
		}else {
			
			Log.warning("PopSizeGibbsSampler: Please provide at least one gene tree");
		}
		
	}
	

	@Override
	public double proposal() {
		
		geneTreeDistributions = this.getTreeDistributions(this);
		
		
		double[] newPopSizes = new double[popSizes.getDimension()];
		double lower = popSizes.getLower();
		double upper = popSizes.getUpper();
		for (int i = 0; i < popSizes.getDimension(); i++) {
			
			try {
				double newValue = samplePopSize(i);
				if (newValue < lower) {
					newValue = lower;
				} else if (newValue > upper) {
					newValue = upper;
				}
			
				newPopSizes[i] = newValue;
				
			}catch(Exception e) {
				return Double.NEGATIVE_INFINITY;
			}
			
			
			
		}
		
		for (int i = 0; i < popSizes.getDimension(); i++) {
			popSizes.setValue(i, newPopSizes[i]);
		}
		return Double.POSITIVE_INFINITY;
		
	}

	/** 
	 * sample population size for branch b 
	 * @param branch
	 * @return
	 */
	private double samplePopSize(int branch) throws Exception {
		
		if (treeIntervals != null) {
			return constantCoalescentSample();
		}
		
		double a = 0; // = sum_j k_{jb}
		for (GeneTreeForSpeciesTreeDistribution prior : geneTreeDistributions) {
			a += prior.getCoalescentCount(branch);
		}
		
		double b = 0; // = sum_j 1/ploidy \sum_i c_jbi(2 choose (n_jb - i))
		for (GeneTreeForSpeciesTreeDistribution prior : geneTreeDistributions) {
			double c = 0;
			double [] times = prior.getTimes(branch);
			double nbj = prior.getLineageCount(branch);
			for (int i = 0; i < times.length - 1; i++) {
				c += (times[i+1] - times[i]) * (nbj - i) * (nbj - i - 1) / 2;
			}
			double ploidy = prior.getPloidy();
			c /= ploidy;
			b += c;
		}
		
		
		double alpha = priorAlpha.getArrayValue() + a;
		double beta = priorBeta.getArrayValue() + b;
		
		GammaDistribution g = new GammaDistribution(myRandomizer, alpha, 1.0/beta, GammaDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
		double newN = 1.0/g.sample();
		return newN;
		
	}
/*
	public static void main(String[] args) {
		MyRandomizer myRandomizer = new MyRandomizer();
		double alpha = 2;
		double beta = 0.0015;
		int n = 1000;
		double [] samples = new double[n];
		for (int i = 0; i < n; i++) {
			GammaDistribution g = new GammaDistribution(myRandomizer, alpha, beta, GammaDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
			samples[i] = g.sample();
		}
		double sum = 0;
		for (double d : samples) {
			sum +=d;
		}
		System.out.println(sum / n);
	}
*/


	private double constantCoalescentSample() {
		double a = 0; // = sum_j k_{jb}
		TreeInterface tree = treeIntervals.treeInput.get();
		a = tree.getInternalNodeCount();
		
		double b = 0; // = sum_j 1/ploidy \sum_i c_jbi(2 choose (n_jb - i))

		double c = 0;
		double [] times = treeIntervals.getCoalescentTimes(null);
		double nbj = tree.getLeafNodeCount(); //prior.getLineageCount(branch);
		c += (times[0]) * (nbj) * (nbj - 1) / 2;
		nbj--;
		for (int i = 0; i < times.length - 1; i++) {
			c += (times[i+1] - times[i]) * (nbj - i) * (nbj - i - 1) / 2;
		}
		double ploidy = 1; // prior.getPloidy();
		c /= ploidy;
		b += c;
		
		
		double alpha = priorAlpha.getArrayValue() + a;
		double beta = priorBeta.getArrayValue() + b;
		
		GammaDistribution g = new GammaDistribution(myRandomizer, alpha, 1.0/beta, GammaDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
		double newN = 1.0/g.sample();
		return newN;
	}	
}
