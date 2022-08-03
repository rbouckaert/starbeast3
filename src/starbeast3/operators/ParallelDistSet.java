package starbeast3.operators;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.BEASTObject;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.tree.Tree;

@Description("A set of parallel tree distributions, which will all be added to the same thread in a ParallelMCMCOperator")
public class ParallelDistSet extends BEASTObject implements Comparable<ParallelDistSet> {
	
	
	final public Input<List<ParallelMCMCTreeOperatorTreeDistribution>> distributionInput = new Input<>("distribution", 
			"Distribution on a tree conditionally independent from all other distributions given the state of the rest of parameter space. ",
			new ArrayList<>());

	
	public ParallelDistSet() {
		
	}
	
	public ParallelDistSet(ParallelDistSet dist1, ParallelDistSet dist2) {
		List<ParallelMCMCTreeOperatorTreeDistribution> dists = new ArrayList<>();
		dists.addAll(dist1.getDists());
		dists.addAll(dist2.getDists());
		distributionInput.setValue(dists, this);
	}

	@Override
	public void initAndValidate() {
		
		if (distributionInput.get().isEmpty()) {
			throw new IllegalArgumentException("Please ensure there is at least 1 distribution");
		}
		
	}

	@Override
	public int compareTo(ParallelDistSet arg0) {
		// TODO Auto-generated method stub
		return 0;
	}
	
	
	
	/**
	 * Return the distributions
	 * @return
	 */
	public List<ParallelMCMCTreeOperatorTreeDistribution> getDists(){
		return distributionInput.get();
	}
	
	
	/**
	 * Number of taxa, summed across all trees
	 * @return
	 */
	public int getTaxonCount() {
		int nTaxa = 0;
		for (ParallelMCMCTreeOperatorTreeDistribution dist : distributionInput.get()) {
			nTaxa += dist.getTaxonCount();
		}
		return nTaxa;
	}
	
	/**
	 * Number of site patterns, summer across all likelihoods
	 * @return
	 */
	public int getNumberPatterns() {
		int nPatterns = 0;
		for (ParallelMCMCTreeOperatorTreeDistribution dist : distributionInput.get()) {
			nPatterns += dist.getNumberPatterns();
		}
		return nPatterns;
	}

	
	/**
	 * Get the trees
	 */
	public List<Tree> getTrees() {
		List<Tree> trees = new ArrayList<>();
		for (ParallelMCMCTreeOperatorTreeDistribution dist : distributionInput.get()) {
			trees.addAll(dist.getTrees());
		}
		return trees;
	}
	
	


}







