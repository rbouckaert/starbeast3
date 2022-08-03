package starbeast3.operators;

import java.util.ArrayList;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import starbeast3.GeneTreeForSpeciesTreeDistribution;

@Description("Distribution on a tree conditinionally independent from all other distributions given the state of the rest of parameter space")
public class ParallelMCMCTreeOperatorTreeDistribution extends ParallelDistSet {
	
	
	
	
	final public Input<Tree> treeInput = new Input<>("tree", "the tree", Validate.REQUIRED);
	final public Input<GenericTreeLikelihood> treelikelihoodInput = new Input<>("treelikelihood", "treelikelihood part of the distribution", Validate.REQUIRED);
	final public Input<GeneTreeForSpeciesTreeDistribution> genepriorInput = new Input<>("geneprior", "prior on the gene tree", Validate.REQUIRED);
	final public Input<List<RealParameter>> includeInput = new Input<>("include", "optional additional parameters to include", new ArrayList<>());
	
	
	Tree tree;
	//Distribution treelikelihood;
	GenericTreeLikelihood treelikelihood;
	GeneTreeForSpeciesTreeDistribution geneprior;
	List<RealParameter> includeExtraParameters;
	
	
	
	 
	public Tree getTree() {return tree;}
	public void setTree(Tree tree) {this.tree = tree;}
	public GenericTreeLikelihood getTreelikelihood() {return treelikelihood;}
	public void setTreelikelihood(GenericTreeLikelihood treelikelihood) {this.treelikelihood = treelikelihood;}
	public GeneTreeForSpeciesTreeDistribution getGeneprior() {return geneprior;}
	public void setGeneprior(GeneTreeForSpeciesTreeDistribution geneprior) {this.geneprior = geneprior;}
	public List<RealParameter> getInclude() {return includeExtraParameters;}
	public void setInclude(List<RealParameter> parameters) {this.includeExtraParameters = parameters;}

	
	public ParallelMCMCTreeOperatorTreeDistribution() {
		
	}
	/*
	public ParallelMCMCTreeOperatorTreeDistribution(@Param(name="tree", description="tree for which") Tree tree,
			@Param(name="treelikelihood", description="treelikelihood part of the distribution") GenericTreeLikelihood treelikelihood,
			@Param(name="geneprior", description="prior on the gene tree") GeneTreeForSpeciesTreeDistribution geneprior){
		this.tree = tree;
		this.treelikelihood = treelikelihood;
		this.geneprior = geneprior;
		this.includeExtraParameters = new ArrayList<>();
	}
	
	
	public ParallelMCMCTreeOperatorTreeDistribution(@Param(name="tree", description="tree for which") Tree tree,
			@Param(name="treelikelihood", description="treelikelihood part of the distribution") GenericTreeLikelihood treelikelihood,
			@Param(name="geneprior", description="prior on the gene tree") GeneTreeForSpeciesTreeDistribution geneprior,
			@Param(name="include", description="additional real parameter to include") List<RealParameter> include) {
		this.tree = tree;
		this.treelikelihood = treelikelihood;
		this.geneprior = geneprior;
		this.includeExtraParameters = include;
	}
	*/

	@Override
	public void initAndValidate() {
		this.tree = treeInput.get();
		this.treelikelihood = treelikelihoodInput.get();
		this.geneprior = genepriorInput.get();
		this.includeExtraParameters = includeInput.get();
	}
	
	
	/**
	 * More site patterns = greater
	 */
	@Override
	public int compareTo(ParallelDistSet other) {
		int patterns1 = this.getNumberPatterns();
		int patterns2 = other.getNumberPatterns();
		if (patterns1 < patterns2) return -1;
		if (patterns1 > patterns2) return  1;
		return 0;
	}
	
	
	/**
	 * Return the actual distribution as a list
	 * @return
	 */
	@Override
	public List<ParallelMCMCTreeOperatorTreeDistribution> getDists(){
		List<ParallelMCMCTreeOperatorTreeDistribution> dists = new ArrayList<>();
		dists.add(this);
		return dists;
	}
	
	
	
	/**
	 * Get the tree, as well as the tree of the likelihood and prior
	 */
	@Override
	public List<Tree> getTrees() {
		List<Tree> trees = new ArrayList<>();
		trees.add(this.tree);
		TreeInterface priorTree = this.getGeneprior().treeInput.get();
		TreeInterface likelihoodTree = this.treelikelihood.treeInput.get();
		if (!trees.contains(priorTree) && priorTree instanceof Tree) trees.add((Tree)priorTree);
		if (!trees.contains(likelihoodTree) && likelihoodTree instanceof Tree) trees.add((Tree)likelihoodTree);
		return trees;
	}
	
	
	/**
	 * Number of taxa
	 * @return
	 */
	@Override
	public int getTaxonCount() {
		return this.tree.getTaxaNames().length;
	}
	
	/**
	 * Number of site patterns
	 * @return
	 */
	@Override
	public int getNumberPatterns() {
		return this.getTreelikelihood().dataInput.get().getPatternCount();
	}
	
	
	/**
	 * Call before beginning parallel MCMC
	 */
	public void startThreading() {
		
	}
	
	
	/**
	 * Call after finishing parallel MCMC
	 */
	public void stopThreading() {
		
	}

	
}
