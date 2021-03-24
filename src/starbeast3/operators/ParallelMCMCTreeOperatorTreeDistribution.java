package starbeast3.operators;

import beast.core.BEASTObject;
import beast.core.Description;
import beast.core.Distribution;
import beast.core.Param;
import beast.evolution.likelihood.GenericTreeLikelihood;
import beast.evolution.likelihood.TreeLikelihood;
import beast.evolution.tree.Tree;
import starbeast3.GeneTreeForSpeciesTreeDistribution;

@Description("Distribution on a tree conditinionally independent from all other distributions given the state of the rest of parameter space")
public class ParallelMCMCTreeOperatorTreeDistribution extends BEASTObject implements Comparable<ParallelMCMCTreeOperatorTreeDistribution> {
	Tree tree;
	Distribution treelikelihood;
	GeneTreeForSpeciesTreeDistribution geneprior;
	
	public Tree getTree() {return tree;}
	public void setTree(Tree tree) {this.tree = tree;}
	public Distribution getTreelikelihood() {return treelikelihood;}
	public void setTreelikelihood(Distribution treelikelihood) {this.treelikelihood = treelikelihood;}
	public GeneTreeForSpeciesTreeDistribution getGeneprior() {return geneprior;}
	public void setGeneprior(GeneTreeForSpeciesTreeDistribution geneprior) {this.geneprior = geneprior;}

	public ParallelMCMCTreeOperatorTreeDistribution(@Param(name="tree", description="tree for which") Tree tree,
			@Param(name="treelikelihood", description="treelikelihood part of the distribution") Distribution treelikelihood,
			@Param(name="geneprior", description="prior on the gene tree") GeneTreeForSpeciesTreeDistribution geneprior) {
		this.tree = tree;
		this.treelikelihood = treelikelihood;
		this.geneprior = geneprior;
	}
	
	@Override
	public void initAndValidate() {
	}
	
	
	/**
	 * More site patterns = greater
	 */
	@Override
	public int compareTo(ParallelMCMCTreeOperatorTreeDistribution other) {
		int patterns1 = this.getNumberPatterns();
		int patterns2 = other.getNumberPatterns();
		if (patterns1 < patterns2) return -1;
		if (patterns1 > patterns2) return  1;
		return 0;
	}
	
	/**
	 * Number of site patterns
	 * @return
	 */
	public int getNumberPatterns() {
		if (this.getTreelikelihood() instanceof GenericTreeLikelihood) {
			return ((GenericTreeLikelihood) this.getTreelikelihood()).dataInput.get().getPatternCount();
		}
		return 0;
	}
	
}
