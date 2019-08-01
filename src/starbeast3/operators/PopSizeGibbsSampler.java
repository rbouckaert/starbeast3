package starbeast3.operators;

import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math3.distribution.GammaDistribution;

import beast.core.Description;
import beast.core.Input;
import beast.core.Operator;
import beast.core.Input.Validate;
import beast.core.parameter.RealParameter;
import beast.evolution.tree.TreeInterface;
import beast.evolution.tree.coalescent.TreeIntervals;
import beast.math.distributions.Gamma;
import beast.math.distributions.InverseGamma;
import starbeast3.GeneTreeForSpeciesTreeDistribution;
import starbeast3.util.*;

@Description("Samples population sizes of constant populations for multi species coalescent model")
public class PopSizeGibbsSampler extends Operator {
	final public Input<RealParameter> popSizesInput = new Input<>("popSizes", "constant population size parameter, one dimension for each branch of the species tree", Validate.REQUIRED);
	final public Input<Gamma> priorInput = new Input<>("gammaprior", "gamma distributed prior for population sizes", Validate.REQUIRED);
	final public Input<List<GeneTreeForSpeciesTreeDistribution>> genePriorsInput = new Input<>("gene", "gene tree for species tree distribution for each of the genes", new ArrayList<>());
	final public Input<TreeIntervals> treeIntervalsInput = new Input<>("intervals", "tree intervals for use with single tree -- should not be used if gene-attribute is used");

	RealParameter popSizes, priorAlpha, priorBeta;
	Gamma prior;
	List<GeneTreeForSpeciesTreeDistribution> genePriors;
	TreeIntervals treeIntervals;

	// used to sample a gamma distribution
	MyRandomizer myRandomizer = new MyRandomizer();
	
	@Override
	public void initAndValidate() {
		popSizes = popSizesInput.get();
		prior = priorInput.get();
		genePriors = genePriorsInput.get();
		treeIntervals = treeIntervalsInput.get();
		priorAlpha = prior.alphaInput.get();
		priorBeta = prior.betaInput.get();
		if (treeIntervals != null) { 
			return;
		}
		
		if (genePriors.size() == 0) {
			throw new IllegalArgumentException("At least one gene needs to be prodided (or intervals specified)");
		}
		
		TreeInterface speciesTree = genePriors.get(0).speciesTreeInput.get();
		if (speciesTree.getNodeCount() != popSizes.getDimension()) {
			throw new IllegalArgumentException("The dimension of the population size parameter (" + popSizes.getDimension()+ ") "
					+ "should be equal to the number of branches in the species tree ( " + speciesTree.getNodeCount() + ")");
		}
	}
	

	@Override
	public double proposal() {
		double lower = popSizes.getLower();
		double upper = popSizes.getUpper();
		for (int i = 0; i < popSizes.getDimension(); i++) {
			double newValue = samplePopSize(i);
			if (newValue < lower) {
				newValue = lower;
			} else if (newValue > upper) {
				newValue = upper;
			}
			popSizes.setValue(i, newValue);
		}
		return Double.POSITIVE_INFINITY;
	}

	/** 
	 * sample population size for branch b 
	 * @param branch
	 * @return
	 */
	private double samplePopSize(int branch) {
		if (treeIntervals != null) {
			return constantCoalescentSample();
		}
		
		double a = 0; // = sum_j k_{jb}
		for (GeneTreeForSpeciesTreeDistribution prior : genePriors) {
			a += prior.getCoalescentCount(branch);
		}
		
		double b = 0; // = sum_j 1/ploidy \sum_i c_jbi(2 choose (n_jb - i))
		for (GeneTreeForSpeciesTreeDistribution prior : genePriors) {
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
		
		
		double alpha = prior.alphaInput.get().getValue() + a;
		double beta = prior.betaInput.get().getValue() + b;
		
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
		
		
		double alpha = prior.alphaInput.get().getValue() + a;
		double beta = prior.betaInput.get().getValue() + b;
		
		GammaDistribution g = new GammaDistribution(myRandomizer, alpha, 1.0/beta, GammaDistribution.DEFAULT_INVERSE_ABSOLUTE_ACCURACY);
		double newN = 1.0/g.sample();
		return newN;
	}	
}
