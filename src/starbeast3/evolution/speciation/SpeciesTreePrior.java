package starbeast3.evolution.speciation;


import java.util.ArrayList;
import java.util.Arrays;
import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.State;
import beast.base.inference.parameter.RealParameter;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.TreeDistribution;
import beast.base.inference.distribution.Gamma;
import beast.base.inference.distribution.InverseGamma;



@Description("Species tree prior for *BEAST analysis")
public class SpeciesTreePrior extends TreeDistribution {
    //public Input<Tree> m_speciesTree = new Input<>("speciesTree", "species tree containing the associated gene tree", Validate.REQUIRED);

    protected enum TreePopSizeFunction {constant, linear, linear_with_constant_root}

    public final Input<TreePopSizeFunction> popFunctionInput = new Input<>("popFunction", "Population function. " +
            "This can be " + Arrays.toString(TreePopSizeFunction.values()) + " (default 'constant')", TreePopSizeFunction.constant, TreePopSizeFunction.values());

    public final Input<RealParameter> popSizesBottomInput = new Input<>("bottomPopSize", "population size parameter for populations at the bottom of a branch. " +
            "For linear population function, this is the same at the top of the branch.", Validate.REQUIRED);
    public final Input<RealParameter> popSizesTopInput = new Input<>("topPopSize", "population size parameter at the top of a branch. " +
            "Ignored for constant population function, but required for linear population function.");

    public final Input<RealParameter> gammaParameterInput = new Input<>("gammaParameter", "scale parameter of the gamma distribution over population sizes. "
    		+ "This makes this parameter half the expected population size on all branches for constant population function, "
    		+ "but a quarter of the expected population size for tip branches only for linear population functions.", Validate.REQUIRED);

    final public Input<PopulationModel> popModelInput = 
    		new Input<>("populationModel", "Population model used to infer the multispecies coalescent probability", Validate.REQUIRED);
    
    final public Input<TreeDistribution> treePriorInput =
    		new Input<>("treePrior", "Prior distribution behind the species tree");
    
    
    
    
//	public Input<RealParameter> m_rootHeightParameter = new Input<>("rootBranchHeight","height of the node above the root, representing the root branch", Validate.REQUIRED);
    /**
     * m_taxonSet is used by GeneTreeForSpeciesTreeDistribution *
     */
    final public Input<TaxonSet> taxonSetInput = new Input<>("taxonset", "set of taxa mapping lineages to species", Validate.REQUIRED);


    private TreePopSizeFunction popFunction;
    private RealParameter popSizesBottom;
    private RealParameter popSizesTop;

    private InverseGamma gamma2Prior;
    private Gamma gamma4Prior;

    @Override
    public void initAndValidate() {
        popFunction = popFunctionInput.get();
        popSizesBottom = popSizesBottomInput.get();
        popSizesTop = popSizesTopInput.get();

        // set up sizes of population functions
        final int speciesCount = treeInput.get().getLeafNodeCount();
        final int nodeCount = treeInput.get().getNodeCount();
        switch (popFunction) {
            case constant:
                popSizesBottom.setDimension(nodeCount);
                break;
            case linear:
                if (popSizesTop == null) {
                    throw new IllegalArgumentException("topPopSize must be specified");
                }
                popSizesBottom.setDimension(speciesCount);
                popSizesTop.setDimension(nodeCount);
                break;
            case linear_with_constant_root:
                if (popSizesTop == null) {
                    throw new IllegalArgumentException("topPopSize must be specified");
                }
                popSizesBottom.setDimension(speciesCount);
                popSizesTop.setDimension(nodeCount - 1);
                break;
        }

        // bottom prior = Gamma(2,Psi)
        gamma2Prior = new InverseGamma();
        gamma2Prior.betaInput.setValue(gammaParameterInput.get(), gamma2Prior);

        // top prior = Gamma(4,Psi)
        gamma4Prior = new Gamma();
        final RealParameter parameter = new RealParameter(new Double[]{4.0});
        gamma4Prior.alphaInput.setValue(parameter, gamma4Prior);
        gamma4Prior.betaInput.setValue(gammaParameterInput.get(), gamma4Prior);

        if (popFunction != TreePopSizeFunction.constant && gamma4Prior == null) {
            throw new IllegalArgumentException("Top prior must be specified when population function is not constant");
        }
        // make sure the m_taxonSet is a set of taxonsets
// HACK to make Beauti initialise: skip the check here
//		for (Taxon taxon : m_taxonSet.get().m_taxonset.get()) {
//			if (!(taxon instanceof TaxonSet)) {
//				throw new IllegalArgumentException("taxonset should be sets of taxa only, not individual taxons");
//			}
//		}
    }
    
    
    public RealParameter getPopulationSizes() {
    	return popSizesBottomInput.get();
    }
    

    @Override
    public double calculateLogP() {
        logP = 0;
        // make sure the root branch length is positive
//		if (m_rootHeightParameter.get().getValue() < m_speciesTree.get().getRoot().getHeight()) {
//			logP = Double.NEGATIVE_INFINITY;
//			return logP;
//		}

        final Node[] speciesNodes = treeInput.get().getNodesAsArray();
        try {
            switch (popFunction) {
                case constant:
                    // constant pop size function
                    //logP += gamma2Prior.calcLogP(popSizesBottom);
//			for (int i = 0; i < speciesNodes.length; i++) {
//				double popSize = m_fPopSizesBottom.getValue(i);
//				logP += m_bottomPrior.logDensity(popSize); 
//			}
                    break;
                case linear:
                    // linear pop size function
//			int speciesCount = m_tree.get().getLeafNodeCount();
//			m_fPopSizesBottom.setDimension(speciesCount);
//			logP += m_gamma4Prior.calcLogP(m_fPopSizesBottom);
//			int nodeCount = m_tree.get().getNodeCount();
//			m_fPopSizesTop.setDimension(nodeCount-1);
//			logP += m_gamma2Prior.calcLogP(m_fPopSizesTop);

                    for (int i = 0; i < speciesNodes.length; i++) {
                        final Node node = speciesNodes[i];
                        final double popSizeBottom;
                        if (node.isLeaf()) {
                            // Gamma(4, psi) prior
                            popSizeBottom = popSizesBottom.getValue(i);
                            logP += gamma4Prior.logDensity(popSizeBottom);
                        }
                        final double popSizeTop = popSizesTop.getValue(i);
                        logP += gamma2Prior.logDensity(popSizeTop);
                    }
                    break;
                case linear_with_constant_root:
//			logP += m_gamma4Prior.calcLogP(m_fPopSizesBottom);
//			logP += m_gamma2Prior.calcLogP(m_fPopSizesTop);
//			int rootNr = m_tree.get().getRoot().getNr();
//			double popSize = m_fPopSizesTop.getValue(rootNr);
//			logP -= m_gamma2Prior.logDensity(popSize); 

                    for (int i = 0; i < speciesNodes.length; i++) {
                        final Node node = speciesNodes[i];
                        if (node.isLeaf()) {
                            final double popSizeBottom = popSizesBottom.getValue(i);
                            logP += gamma4Prior.logDensity(popSizeBottom);
                        }
                        if (!node.isRoot()) {
                            if (i < speciesNodes.length - 1) {
                                final double popSizeTop = popSizesTop.getArrayValue(i);
                                logP += gamma2Prior.logDensity(popSizeTop);
                            } else {
                                final int nodeIndex = treeInput.get().getRoot().getNr();
                                final double popSizeTop = popSizesTop.getArrayValue(nodeIndex);
                                logP += gamma2Prior.logDensity(popSizeTop);
                            }
                        }
                    }
                    break;
            }
        } catch (Exception e) {
            // exceptions can be thrown by the gamma priors
            e.printStackTrace();
            return Double.NEGATIVE_INFINITY;
        }
        return logP;
    }

    @Override
    protected boolean requiresRecalculation() {
        return true;
    }

    @Override
    public List<String> getArguments() {
    	List<String> arguments = new ArrayList<>();
        return arguments;
    }

    @Override
    public List<String> getConditions() {
        List<String> arguments = new ArrayList<>();
        arguments.add(treeInput.get().getID());
        if (popSizesBottomInput.get() != null) arguments.add(popSizesBottomInput.get().getID());
        if (popSizesTopInput.get() != null) arguments.add(popSizesTopInput.get().getID());
        return arguments;
    }

    @Override
    public void sample(final State state, final Random random) {
    	
    	if (sampledFlag) return;
        sampledFlag = true;
        
        
        // Sample the species tree
        sampleConditions(state, random);
        
        
    }
}
