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
import beast.math.distributions.Gamma;
import starbeast3.GeneTreeForSpeciesTreeDistribution;
import starbeast3.util.*;

@Description("Samples population sizes of constant populations for multi species coalescent model")
public class PopSizeGibbsSampler extends Operator {
	final public Input<RealParameter> popSizesInput = new Input<>("popSizes", "constant population size parameter, one dimension for each branch of the species tree", Validate.REQUIRED);
	final public Input<Gamma> priorInput = new Input<>("gammaprior", "gamma distributed prior for population sizes", Validate.REQUIRED);
	final public Input<List<GeneTreeForSpeciesTreeDistribution>> genePriorsInput = new Input<>("gene", "gene tree for species tree distribution for each of the genes", new ArrayList<>());

	RealParameter popSizes, priorAlpha, priorBeta;
	Gamma prior;
	List<GeneTreeForSpeciesTreeDistribution> genePriors;

	// used to sample a gamma distribution
	MyRandomizer myRandomizer = new MyRandomizer();
	
	@Override
	public void initAndValidate() {
		popSizes = popSizesInput.get();
		prior = priorInput.get();
		genePriors = genePriorsInput.get();
		if (genePriors.size() == 0) {
			throw new IllegalArgumentException("At least one gene needs to be prodided");
		}
		
		TreeInterface speciesTree = genePriors.get(0).speciesTreeInput.get();
		if (speciesTree.getNodeCount() != popSizes.getDimension()) {
			throw new IllegalArgumentException("The dimension of the population size parameter (" + popSizes.getDimension()+ ") "
					+ "should be equal to the number of branches in the species tree ( " + speciesTree.getNodeCount() + ")");
		}
		priorAlpha = prior.alphaInput.get();
		priorBeta = prior.betaInput.get();
	}
	

	@Override
	public double proposal() {
		for (int i = 0; i < popSizes.getDimension(); i++) {
			popSizes.setValue(i, samplePopSize(i));
		}
		return Double.POSITIVE_INFINITY;
	}

	/** 
	 * sample population size for branch b 
	 * @param b
	 * @return
	 */
	private double samplePopSize(int b) {
		double q = 0; // = sum_j k_{jb}
		for (GeneTreeForSpeciesTreeDistribution prior : genePriors) {
			q += prior.getCoalescentCount(b);
		}
		
		double gamma = 0; // = sum_j 1/ploidy \sum_i c_jbi(2 choose (n_jb - i))
		for (GeneTreeForSpeciesTreeDistribution prior : genePriors) {
			double c = 0;
			double [] times = prior.getTimes(b);
			double nbj = times.length - 2;
			for (int i = 0; i < times.length - 1; i++) {
				c += (times[i+1] - times[i]) * (nbj - i) * (nbj - i - 1) / 2;
			}
			double ploidy = prior.getPloidy();
			c /= ploidy;
			gamma += c;
		}
		
		
		double alpha = prior.alphaInput.get().getValue() + q + 1.0;
		double beta = prior.betaInput.get().getValue() + gamma;
		
		GammaDistribution g = new GammaDistribution(myRandomizer, alpha, beta);
		double newN = g.sample();
		return newN;
	}

}
