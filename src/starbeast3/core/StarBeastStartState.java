package starbeast3.core;

import static java.lang.Math.abs;
import static java.lang.Math.max;
import static java.lang.Math.min;

import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.LinkedHashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.MathException;

import beast.base.core.BEASTInterface;
import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.inference.parameter.RealParameter;
import beastlabs.evolution.tree.ConstrainedClusterTree;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.alignment.Taxon;
import beast.base.evolution.alignment.TaxonSet;
import beast.base.evolution.distance.Distance;
import beast.base.evolution.distance.JukesCantorDistance;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.speciation.CalibratedYuleModel;
import beast.base.evolution.speciation.CalibrationPoint;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.coalescent.RandomTree;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeParser;
import beast.base.evolution.tree.TreeUtils;
import beast.base.evolution.tree.coalescent.ConstantPopulation;
import beast.base.evolution.tree.MRCAPrior;
import beast.base.evolution.tree.ClusterTree;
import starbeast3.evolution.branchratemodel.BranchRateModelSB3;
import starbeast3.evolution.branchratemodel.SharedSpeciesClockModel;
import starbeast3.evolution.branchratemodel.UCRelaxedClockModelSB3;
import starbeast3.evolution.speciation.GeneTreeForSpeciesTreeDistribution;
import starbeast3.evolution.speciation.SpeciesTreePrior;
import starbeast3.genekernel.GTKPrior;
import starbeast3.tree.SpeciesTree;

/**
* @author Joseph Heled
 */

@Description("Set a starting point for a *BEAST analysis from gene alignment data.")
@Citation(value="Douglas, Jordan, Cinthy L. Jimenez-Silva, and Remco Bouckaert. StarBeast3: Adaptive Parallelised Bayesian Inference under the Multispecies Coalescent. Systematic Biology (2022).", DOI="10.1093/sysbio/syac010")
public class StarBeastStartState extends Tree implements StateNodeInitialiser {

    static enum Method {
        POINT("point-estimate"),
        FIXED("fixed"),
        ALL_RANDOM("random");

        Method(final String name) {
            this.ename = name;
        }

        @Override
		public String toString() {
            return ename;
        }

        private final String ename;
    }
    final public Input<Method> initMethod = new Input<>("method", "Initialise either with a totally random " +
            "state or a point estimate based on alignments data (default point-estimate)",
            Method.POINT, Method.values());

    final public Input<Tree> speciesTreeInput = new Input<>("speciesTree", "The species tree to initialize");
    final public Input<TreeParser> fixedInput = new Input<>("fixed", "Optionally provide a newick of the species tree to use at start state");
    
    
    final public Input<List<Tree>> genesInput = new Input<>("gene", "Gene trees to initialize", new ArrayList<>());
    //,
    //        Validate.REQUIRED);

    final public Input<CalibratedYuleModel> calibratedYule = new Input<>("calibratedYule",
            "The species tree (with calibrations) to initialize", Validate.XOR, speciesTreeInput);

    final public Input<RealParameter> popMean = new Input<>("popMean",
            "Population mean hyper prior to initialse");
    
    
    final public Input<List<RealParameter>> originInput = new Input<>("origin",
            "Species tree origin height", new ArrayList<>());
    
    final public Input<List<RealParameter>> originBranchInput = new Input<>("originLength",
            "Species tree origin branch length", new ArrayList<>());

    final public Input<RealParameter> birthRate = new Input<>("birthRate",
            "Tree prior birth rate to initialize");

    final public Input<SpeciesTreePrior> speciesTreePriorInput =
            new Input<>("speciesTreePrior", "Population size parameters to initialise");

    final public Input<Function> muInput = new Input<>("baseRate",
            "Main clock rate used to scale trees (default 1).");
    
    final public Input<Function> rootHeightInput = new Input<>("rootHeight",
            "Initial root height (default at substitution length)");
    
    
    final public Input<BranchRateModelSB3> speciesTreeRatesInput =
            new Input<>("speciesTreeRates", "Clock model for species tree");
    
    final public Input<SharedSpeciesClockModel> sharedRateModelInput =
            new Input<>("sharedRateModel", "Clock model for species tree (instead of speciesTreeRates)");
    
    
    // A gene tree kernel
    final public Input<GTKPrior> geneKernelPriorInput =
            new Input<>("kernel", "Clock model for species tree", Input.Validate.OPTIONAL);
    
    
    List<Tree> genes;
    private boolean hasCalibrations;


    @Override
    public void initAndValidate() {
        // what does this do and is it dangerous to call it or not to call it at the start or at the end??????
        super.initAndValidate();
        hasCalibrations = calibratedYule.get() != null;
        genes = genesInput.get();
        
        // Fix the intial tree?
        if (fixedInput.get() != null) {
        	initMethod.set(Method.FIXED);
        }
        if (initMethod.get() == Method.FIXED) {
        	if (fixedInput.get() == null) {
        		throw new IllegalArgumentException("Please provide a starting 'newick' if using the fixed method");
        	}
        }

        
        // Get clock model
        BranchRateModelSB3 s = null;
        if (speciesTreeRatesInput.get() != null) {
        	s = speciesTreeRatesInput.get();
        }else if (sharedRateModelInput.get() != null) {
        	s = sharedRateModelInput.get().getClockModel();	
        }
    	if (s != null && s instanceof UCRelaxedClockModelSB3) {
    		UCRelaxedClockModelSB3 s2 = (UCRelaxedClockModelSB3) s;
    		RealParameter p = s2.realRatesInput.get();
    		if (p != null) {
    			rates = p;
    			lowerRate = 0.1; // p .getLower();
    		}
    	}
    	
    	
        
    }
    
    
    

    @Override
    public void initStateNodes() {

    	
    	// Find calibrations
        final Set<BEASTInterface> treeOutputs = speciesTreeInput.get().getOutputs();
        List<MRCAPrior> calibrations = new ArrayList<>();
        for (final Object plugin : treeOutputs ) {
            if( plugin instanceof MRCAPrior ) {
            	MRCAPrior o = (MRCAPrior) plugin;
            	if (!calibrations.contains(o)) {
            		calibrations.add(o);
            	}
            }
        }
        
        
        
        // If a gene tree kernel is being used, then initialise the gene trees here to those in the kernel
        if (geneKernelPriorInput.get() != null) {
        	GTKPrior kernel = geneKernelPriorInput.get();
        	
        	// There should not be any genes added to the gene list
        	if (genes.size() > 0) {
        		throw new IllegalArgumentException("Make sure you do not provide any gene trees when using a gene tree kernel");
        	}
        	
        	// Set the genes tree list as the kernel gene tree list
        	genes = kernel.kernelInput.get().getCoercedTrees();
        	
        }
        
        // Calibrations
        if( hasCalibrations ) {
            if( calibrations.size() > 0 ) {
                throw new IllegalArgumentException("Not implemented: mix of calibrated yule and MRCA priors: " +
                        "place all priors in the calibrated Yule");
            }
            try {
				initWithCalibrations();
			} catch (MathException e) {
				throw new IllegalArgumentException(e);
			}
        } else {
            if( calibrations.size() > 0 )  {
                initWithMRCACalibrations(calibrations);
                //return;
            }

            final Method method = initMethod.get();

            switch( method ) {
                case POINT:
                    fullInit(calibrations);
                    break;
                case ALL_RANDOM:
                    randomInit();
                    break;
                case FIXED:
                    fixedInit();
                    break;
            }
        }
        
        
        
        
        // Ensure that all gene tree tips are the same height as their species
        for (Tree gtree : genes) {
        	this.resetGeneTreeTipHeights((SpeciesTree) speciesTreeInput.get(), gtree);
        }
        

        if (rates != null) {
        	// rates.setLower(lowerRate);
        }
        
        
        // Set origin branch length to slightly above tallest tree
        if (!originInput.get().isEmpty() || !originBranchInput.get().isEmpty()) {
        	
        	double maxGeneTreeHeight = 0;
        	for (Tree gene : genes) {
        		maxGeneTreeHeight = Math.max(gene.getRoot().getHeight(), maxGeneTreeHeight);
        	}
        	double speciesHeight = speciesTreeInput.get().getRoot().getHeight();
        	
        	for (RealParameter originParam : originInput.get()) {
        		originParam.setValue(2 * (maxGeneTreeHeight - speciesHeight) + speciesHeight);
        	}
        	
        	for (RealParameter originLengthParam : originBranchInput.get()) {
        		originLengthParam.setValue(2 * (maxGeneTreeHeight - speciesHeight));
        	}
        	
        }
        
    }
    
    

	RealParameter rates;
    double lowerRate;

    private double[] firstMeetings(final Tree gtree, final Map<String, Integer> tipName2Species, final int speciesCount) {
        final Node[] nodes = gtree.listNodesPostOrder(null, null);
        @SuppressWarnings("unchecked")
		final Set<Integer>[] tipsSpecies = new Set[nodes.length];
        for(int k = 0; k < tipsSpecies.length; ++k) {
            tipsSpecies[k] = new HashSet<>();
        }
        // d[i,j] = minimum height of node which has tips belonging to species i and j
        // d is is upper triangular
        final double[] dmin = new double[(speciesCount*(speciesCount-1))/2];
        Arrays.fill(dmin, Double.MAX_VALUE);

        for (final Node n : nodes) {
            if (n.isLeaf()) {
                tipsSpecies[n.getNr()].add(tipName2Species.get(n.getID()));
            } else {
                assert n.getChildCount() == 2;
                @SuppressWarnings("unchecked")
				final Set<Integer>[] sps = new Set[2];
                sps[0] = tipsSpecies[n.getChild(0).getNr()];
                sps[1] = tipsSpecies[n.getChild(1).getNr()];
                final Set<Integer> u = new HashSet<>(sps[0]);
                u.retainAll(sps[1]);
                sps[0].removeAll(u);
                sps[1].removeAll(u);

                for (final Integer s1 : sps[0]) {
                    for (final Integer s2 : sps[1]) {
                        final int i = getDMindex(speciesCount, s1, s2);
                        dmin[i] = min(dmin[i], n.getHeight());
                    }
                }
                u.addAll(sps[0]);
                u.addAll(sps[1]);
                tipsSpecies[n.getNr()] = u;
            }
        }
        return dmin;
    }

    private int getDMindex(final int speciesCount, final int s1, final int s2) {
        final int mij = min(s1,s2);
        return (mij*(2*speciesCount-1 - mij))/2 + (abs(s1-s2)-1);
    }


    private void fullInit(List<MRCAPrior> mrcapriors) {
    	
        // Build gene trees from  alignments
    	if (geneKernelPriorInput.get() != null) {
    		throw new IllegalArgumentException("Point-estimates are currently not supported when using a gene tree kernel. Please use method='random'.");
    	}
    	

        final Function muInput = this.muInput.get();
        final double mu =  (muInput != null )  ? muInput.getArrayValue() : 1;

        final Tree stree = speciesTreeInput.get();
        final TaxonSet species = stree.m_taxonset.get();
        final List<String> speciesNames = species.asStringList();
        final int speciesCount = speciesNames.size();

        final List<Tree> geneTrees = genes;

        //final List<Alignment> alignments = genes.get();
        //final List<Tree> geneTrees = new ArrayList<>(alignments.size());
        double maxNsites = 0;
        //for( final Alignment alignment : alignments)  {
        
        
        
        for (final Tree gtree : geneTrees) {
        	
            final Alignment alignment = gtree.m_taxonset.get().alignmentInput.get();
            final ClusterTree ctree = new ClusterTree();
            ctree.initByName("initial", gtree, "clusterType", "upgma", "taxa", alignment);
            gtree.scale(1 / mu);

            maxNsites = max(maxNsites, alignment.getSiteCount());
        }
        
        
        // Build tip map
        final Map<String, Integer> geneTips2Species = new HashMap<>();
        final List<Taxon> taxonSets = species.taxonsetInput.get();
        for(int k = 0; k < speciesNames.size(); ++k) {
            final Taxon nx = taxonSets.get(k);
            final List<Taxon> taxa = ((TaxonSet) nx).taxonsetInput.get();
            for(final Taxon n : taxa ) {
              geneTips2Species.put(n.getID(), k);
            }
        }

        
        final double[] dg = new double[(speciesCount*(speciesCount-1))/2];

        final double[][] genesDmins = new double[geneTrees.size()][];

        for( int ng = 0; ng < geneTrees.size(); ++ng ) {
            final Tree g = geneTrees.get(ng);
            final double[] dmin = firstMeetings(g, geneTips2Species, speciesCount);
            genesDmins[ng] = dmin;

            for(int i = 0; i < dmin.length; ++i) {
                if (dmin[i] == Double.MAX_VALUE) {
                	// this happens when a gene tree has no taxa for some species-tree taxon.
                	// TODO: ensure that if this happens, there will always be an "infinite"
                	// distance between species-taxon 0 and the species-taxon with missing lineages,
                	// so i < speciesCount - 1.
                	// What if lineages for species-taxon 0 are missing? Then all entries will be 'infinite'.
                	String id = (i < speciesCount - 1? stree.getExternalNodes().get(i+1).getID() : "unknown taxon");
                	if (i == 0) {
                		// test that all entries are 'infinite', which implies taxon 0 has lineages missing 
                		boolean b = true;
                		for (int k = 1; b && k < speciesCount - 1; k++) {
                			b = (dmin[k] == Double.MAX_VALUE);
                		}
                		if (b) {
                			// if all entries have 'infinite' distances, it is probably the first taxon that is at fault
                			id = stree.getExternalNodes().get(0).getID();
                		}
                	}
                	// throw new RuntimeException("Gene tree " + g.getID() + " has no lineages for species taxon " + id + " ");
                } else {
                    dg[i] += dmin[i];
                }
            }
        }

        for(int i = 0; i < dg.length; ++i) {
            double d = dg[i] / geneTrees.size();
            if( d == 0 ) {
               d = (0.5/maxNsites) * (1/mu);
            } else {
                // heights to distances
                d *= 2;
            }
            dg[i] = d;
        }

        if (mrcapriors.size() == 0) {
        	// no MRCA calibrations
		    final ClusterTree ctree = new ClusterTree();
		    final Distance distance = new Distance() {
		        @Override
		        public double pairwiseDistance(final int s1, final int s2) {
		            final int i = getDMindex(speciesCount, s1,s2);
		            return dg[i];
		        }
		    };
		    ctree.initByName("initial", stree, "taxonset", species,"clusterType", "upgma", "distance", distance);
        } else {
		    final ConstrainedClusterTree ctree = new ConstrainedClusterTree();
		    final Distance distance = new Distance() {
		        @Override
		        public double pairwiseDistance(final int s1, final int s2) {
		            final int i = getDMindex(speciesCount, s1,s2);
		            return dg[i];
		        }
		    };
		    ctree.initByName("initial", stree, "taxonset", species,"clusterType", "upgma", "distance", distance, "constraint", mrcapriors);
        	
        }

        // Set height?
        if (rootHeightInput.get() != null) {
        	double rootHeight = rootHeightInput.get().getArrayValue();
        	if (rootHeight > 0) {
	        	double scale = rootHeight / stree.getRoot().getHeight();
	        	stree.scale(scale);
	        	System.out.println("Scaling species tree to height " + rootHeight);
        	}
        }
        
        

       // System.out.println(ctree.getRoot().toNewick());
        
        final Map<String, Integer> sptips2SpeciesIndex = new HashMap<>();
        for(int i = 0; i < speciesNames.size(); ++i) {
            sptips2SpeciesIndex.put(speciesNames.get(i), i);
        }
        final double[] spmin = firstMeetings(stree, sptips2SpeciesIndex, speciesCount);

        for( int ng = 0; ng < geneTrees.size(); ++ng ) {
            final double[] dmin = genesDmins[ng];
            boolean compatible = true;
            for(int i = 0; i < spmin.length; ++i) {
                if( dmin[i] <= spmin[i] ) {
                    compatible = false;
                    break;
                }
            }
            if( ! compatible ) {
                final Tree gtree = geneTrees.get(ng);
                final TaxonSet gtreeTaxa = gtree.m_taxonset.get();
                final Alignment alignment = gtreeTaxa.alignmentInput.get();
                final List<String> taxaNames = alignment.getTaxaNames();
                final int taxonCount =  taxaNames.size();
                // speedup
                final Map<Integer,Integer> g2s = new HashMap<>();
                for(int i = 0; i < taxonCount; ++i) {
                    g2s.put(i, geneTips2Species.get(taxaNames.get(i)));
                }

                final JukesCantorDistance jc = new JukesCantorDistance();
                jc.setPatterns(alignment);
                final Distance gdistance = new Distance() {
                    @Override
                    public double pairwiseDistance(final int t1, final int t2) {
                        final int s1 = g2s.get(t1);
                        final int s2 = g2s.get(t2);
                        double d = jc.pairwiseDistance(t1,t2)/mu;
                        if( s1 != s2 ) {
                            final int i = getDMindex(speciesCount, s1,s2);
                            final double minDist = 2 * spmin[i];
                            if( d <= minDist ) {
                                d = minDist * 1.001;
                            }
                        }
                        return d;
                    }
                };
                final ClusterTree gtreec = new ClusterTree();
                gtreec.initByName("initial", gtree, "taxonset", gtreeTaxa,
                        "clusterType", "upgma", "distance", gdistance);
            }
        }

        {
            final RealParameter lambda = birthRate.get();
            if( lambda != null ) {
                final double rh = stree.getRoot().getHeight();
                double l = 0;
                for(int i = 2; i < speciesCount+1; ++i) {
                    l += 1./i;
                }
                setParameterValue(lambda, (1 / rh) * l);
            }

            double totBranches = 0;
            final Node[] streeNodeas = stree.getNodesAsArray();
            for( final Node n : streeNodeas ) {
                if( ! n.isRoot() ) {
                    totBranches += n.getLength();
                }
            }
            totBranches /= 2* (streeNodeas.length - 1);
            final RealParameter popm = popMean.get();
            if( popm != null ) {
            	setParameterValue(popm, totBranches);
            }
            final SpeciesTreePrior speciesTreePrior = speciesTreePriorInput.get();
            if( speciesTreePrior != null ) {
                final RealParameter popb = speciesTreePrior.popSizesBottomInput.get();
                if( popb != null ) {
                    for(int i = 0; i < popb.getDimension(); ++i) {
                    	setParameterValue(popb, i, 2*totBranches);
                    }
                }
                final RealParameter popt = speciesTreePrior.popSizesTopInput.get();
                if( popt != null ) {
                    for(int i = 0; i < popt.getDimension(); ++i) {
                    	setParameterValue(popt, i, totBranches);
                    }
                }
            }
        }
    }
    
    /** set parameter value taking bounds in account: if out of bounds, use closest boundary value instead **/
    private void setParameterValue(RealParameter p, double value) {
    	setParameterValue(p, 0, value);
    }
    
    private void setParameterValue(RealParameter p, int index, double value) {
    	if (value < p.getLower()) {
    		value = p.getLower();
    	}
    	if (value > p.getUpper()) {
    		value = p.getUpper();
    	}
    	p.setValue(index, value);
    }

    private void randomInitGeneTrees(double speciesTreeHeight) {
      final List<Tree> geneTrees = genes;
        for (final Tree gtree : geneTrees) {
            gtree.makeCaterpillar(speciesTreeHeight, speciesTreeHeight/gtree.getInternalNodeCount(), true);
        }
    }
    
    
    /**
     * Fix the starting species tree
     */
    private void fixedInit() {
    	
    	//Log.warning("before: " +  speciesTreeInput.get().getRoot().toNewick());
    	TreeParser tp = fixedInput.get();
    	tp.initStateNodes();
    	//Log.warning("after: " +  speciesTreeInput.get().getRoot().toNewick());
    	
    	// Fit the gene trees into the species tree
    	
    	for (Tree gtree : genes) {

	    	// Find GeneTreeForSpeciesTreeDistribution for this gene tree
	        final Set<BEASTInterface> treeOutputs = gtree.getOutputs();
	        GeneTreeForSpeciesTreeDistribution prior = null;
	        for (final Object plugin : treeOutputs ) {
	            if( plugin instanceof GeneTreeForSpeciesTreeDistribution ) {
	            	prior = (GeneTreeForSpeciesTreeDistribution) plugin;
	            	break;
	            }
	        }
	        
	        if (prior == null) {
	        	continue;
	        }
	        
	        if (prior.speciesTreeInput.get() != speciesTreeInput.get()) {
	        	continue;
	        }
	        
	        prior.requiresRecalculation();
	    	
	        
	        
	    	// Ensure that all gene tree nodes are above the species they are mapped to
	        double ddt = speciesTreeInput.get().getRoot().getHeight() * 1e-4;
	        double dt = ddt;
	        double rootHeight = speciesTreeInput.get().getRoot().getHeight();
	        
	        
	       
	        
	        

	        // Place all coalescent events slightly above the root initially in case the following fails
	        for (Node geneNode: gtree.getInternalNodes()) {
	        	
	        
	        	// Set its height
	        	//Log.warning("Setting node heigt from " + geneNode.getHeight() + " to " + (speciesNode.getHeight() + dt) );
	        	geneNode.setHeight(rootHeight + dt);
	        	dt += ddt;
	        }
	        
	        
	        
	        placeGeneTreeWithinSpeciesTree(speciesTreeInput.get().getRoot(), prior, gtree, gtree.getLeafNodeCount());
	       // if (newRoot != null) {
	        	
	        	//Tree newTree = new Tree(newRoot);
	        	//if (gtree.m_taxonset.get() != null) {
	        		//newTree.m_taxonset.setValue(gtree.m_taxonset.get(), gtree);
	        	//}
	        	//Log.warning(gtree.getRoot().toNewick());
		        //gtree.assignFromWithoutID(newTree);
		        prior.requiresRecalculation();
		        
		        
	        //}

	
    	}
    	
    	
    }
    
    
    private int placeGeneTreeWithinSpeciesTree(Node speciesNode, GeneTreeForSpeciesTreeDistribution prior, Tree gTree, int internalNodeNr) {
    	
    	
    	// Create a clade within the species leaf
    	if (speciesNode.isLeaf()) {
    		
    		
    		// Find all of the gene leaves that belong to this species tree leaf
    		Set<String> leaves = prior.getLineagesInSpeciesLeaf(speciesNode.getID());
        	if (leaves.isEmpty()) {
        		Log.warning("Unexpected: cannot find gene leaves for " + speciesNode.getID());
        		return -1;
        	}
        	List<Node> leafNodes = new ArrayList<>();
        	for (Node leaf : gTree.getExternalNodes()) {
        		if (leaves.contains(leaf.getID()) )leafNodes.add(leaf);
        	}
        	if (leafNodes.isEmpty()) {
        		Log.warning("Unexpected: cannot find gene nodes for " + speciesNode.getID());
        		return -1;
        	}
        	
        	
        	// Rearrange the gene tree into a caterpillar
        	double ddt = speciesNode.getLength() / leaves.size();
        	double dt = speciesNode.getHeight() + ddt;
        	Node mrca = null;
        	for (Node leaf : leafNodes) {
        		
        		leaf.setHeight(speciesNode.getHeight());
        		if (mrca == null) {
        			mrca = leaf;
        		}else {
        			
        			Node newMrca = gTree.getNode(internalNodeNr);
        			newMrca.setHeight(dt);
        			newMrca.removeAllChildren(true);
        			newMrca.addChild(mrca);
        			newMrca.addChild(leaf);
        			mrca = newMrca;
        			dt += ddt;
        			
        			internalNodeNr++;
        		}
        	}
        	
        	return internalNodeNr;
        	
    	}
    	
    	
    	// Species tree internal node: take the two child nodes and coalesce them within the common ancestor
    	int internalNodeNrLeft = placeGeneTreeWithinSpeciesTree(speciesNode.getChild(0), prior, gTree, internalNodeNr)-1;
    	if (internalNodeNrLeft > 0) {
    		internalNodeNr = internalNodeNrLeft+1;
    	}
    	int internalNodeNrRight = placeGeneTreeWithinSpeciesTree(speciesNode.getChild(1), prior, gTree, internalNodeNr)-1;
    	if (internalNodeNrRight > 0) {
    		internalNodeNr = internalNodeNrRight+1;
    	}
    	
    	Node mrca = gTree.getNode(internalNodeNr);
    	internalNodeNr++;
    	
    	double height = speciesNode.isRoot() ? speciesNode.getHeight()*1.1 : (speciesNode.getHeight() + speciesNode.getLength()*0.5);
    	mrca.setHeight(height);
    	mrca.removeAllChildren(true);
    	mrca.addChild(gTree.getNode(internalNodeNrLeft));
    	mrca.addChild(gTree.getNode(internalNodeNrRight));
    	return internalNodeNr;
    	
    }
    
    
    

    private void randomInit() {
        double lam = 1;
        final RealParameter lambda = birthRate.get();
        if( lambda != null ) {
            lam = lambda.getArrayValue();
        }
        final Tree stree = speciesTreeInput.get();
        final TaxonSet species = stree.m_taxonset.get();
        final int speciesCount = species.asStringList().size();
        double s = 0;
        for(int k = 2; k <= speciesCount; ++k) {
            s += 1.0/k;
        }
        final double rootHeight = (1/lam) * s;
        stree.scale(rootHeight/stree.getRoot().getHeight());
        randomInitGeneTrees(rootHeight);
//        final List<Tree> geneTrees = genes.get();
//        for (final Tree gtree : geneTrees) {
//            gtree.makeCaterpillar(rootHeight, rootHeight/gtree.getInternalNodeCount(), true);
//        }
    }

    private void initWithCalibrations() throws MathException {
        final CalibratedYuleModel cYule = calibratedYule.get();
        final Tree spTree = (Tree) cYule.treeInput.get();

        final List<CalibrationPoint> cals = cYule.calibrationsInput.get();

        final CalibratedYuleModel cym = new CalibratedYuleModel();
        
        cym.getOutputs().addAll(cYule.getOutputs());

        for( final CalibrationPoint cal : cals ) {
          cym.setInputValue("calibrations", cal);
        }
        cym.setInputValue("tree", spTree);
        cym.setInputValue("type", CalibratedYuleModel.Type.NONE);
        cym.initAndValidate();

        final Tree t = cym.compatibleInitialTree();
        assert spTree.getLeafNodeCount() == t.getLeafNodeCount();

        spTree.assignFromWithoutID(t);

//        final CalibratedYuleInitialTree ct = new CalibratedYuleInitialTree();
//        ct.initByName("initial", spTree, "calibrations", cYule.calibrationsInput.get());
//        ct.initStateNodes();
        final double rootHeight = spTree.getRoot().getHeight();
        randomInitGeneTrees(rootHeight);

        cYule.initAndValidate();
    }

    private void initWithMRCACalibrations(List<MRCAPrior> calibrations) {
        final Tree spTree = speciesTreeInput.get();
        final RandomTree rnd = new RandomTree();
        rnd.setInputValue("taxonset", spTree.getTaxonset());

        for( final MRCAPrior cal : calibrations ) {
          rnd.setInputValue("constraint", cal);
        }
        ConstantPopulation pf = new ConstantPopulation();
        pf.setInputValue("popSize", new RealParameter("1.0"));

        rnd.setInputValue("populationModel", pf);
        rnd.initAndValidate();
        spTree.assignFromWithoutID((Tree)rnd);

        final double rootHeight = spTree.getRoot().getHeight();
        randomInitGeneTrees(rootHeight);
    }

    @Override
    public void getInitialisedStateNodes(List<StateNode> stateNodes) {
        if( hasCalibrations ) {
            stateNodes.add((Tree) calibratedYule.get().treeInput.get());
        } else {
          stateNodes.add(speciesTreeInput.get());
        }

        for( final Tree g : genes ) {
            stateNodes.add(g);
        }

        final RealParameter popm = popMean.get();
        if( popm != null ) {
            stateNodes.add(popm);
        }
        final RealParameter brate = birthRate.get();
        if( brate != null ) {
            stateNodes.add(brate) ;
        }

        final SpeciesTreePrior speciesTreePrior = speciesTreePriorInput.get();
        if( speciesTreePrior != null ) {
            final RealParameter popb = speciesTreePrior.popSizesBottomInput.get();
            if( popb != null ) {
                stateNodes.add(popb) ;
            }
            final RealParameter popt = speciesTreePrior.popSizesTopInput.get();
            if( popt != null ) {
                stateNodes.add(popt);
            }
        }
    }
    
    
    /**
     * Ensure that all gene tree leaves are the same as their species leaf
     * Then ensure that all internal nodes are above their children
     * @param speciesTree
     * @param gtree
     */
    private void resetGeneTreeTipHeights(SpeciesTree speciesTree, Tree gtree) {
    	
    	
    	// Find GeneTreeForSpeciesTreeDistribution for this gene tree
        final Set<BEASTInterface> treeOutputs = gtree.getOutputs();
        GeneTreeForSpeciesTreeDistribution prior = null;
        for (final Object plugin : treeOutputs ) {
            if( plugin instanceof GeneTreeForSpeciesTreeDistribution ) {
            	prior = (GeneTreeForSpeciesTreeDistribution) plugin;
            	break;
            }
        }
    	
        if (prior == null) {
        	Log.warning("Cannot reset tip dates for " + gtree.getID() + " because it does not have a tree prior");
        	return;
        }
        
        if (prior.speciesTreeInput.get() != speciesTree) {
        	Log.warning("Cannot reset tip dates for " + gtree.getID() + " because its tree prior does not have the matching species tree " + prior.speciesTreeInput.get().getID());
        	return;
        }
    	
    	// Set leaf heights
        for (Node geneLeaf: gtree.getExternalNodes()) {
        	
        	// Find mapped node
        	Node speciesLeaf = prior.mapGeneNodeToSpeciesNode(geneLeaf.getNr());
        	if (!speciesLeaf.isLeaf()) {
        		throw new IllegalArgumentException("Unexpected error: the species tree node mapped to gene leaf " + geneLeaf.getID() + " is not a leaf");
        	}
        	
        	// Set its height
        	geneLeaf.setHeight(speciesLeaf.getHeight());
        }
        
        
       // Recursively set internal node heights above their leaves
        shiftInternalNodes(gtree.getRoot());
        	
        
	}
    
    /**
     * Ensure that all internal nodes are above both children
     * @param node
     */
    private void shiftInternalNodes(Node node) {
    	
    	if (node.isLeaf()) return;
    	
    	
    	// Repeat for children
    	double maxChildHeight = 0;
    	for (Node child : node.getChildren()) {
    		shiftInternalNodes(child);
    		maxChildHeight = Math.max(maxChildHeight, child.getHeight());
    	}
    	
    	// Set this height above the oldest child
    	if (node.getHeight() <= maxChildHeight) {
    		double newHeight = maxChildHeight*1.1 + 1e-8;
    		node.setHeight(newHeight);
    	}
    	
    	
    	
    }
    

}
