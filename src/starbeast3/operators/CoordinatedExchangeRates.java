package starbeast3.operators;

import java.util.ArrayList;
import java.util.List;


import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.parameter.RealParameter;
import beast.base.core.Log;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.distribution.LogNormalDistributionModel;
import beast.base.inference.distribution.ParametricDistribution;
import beast.base.util.Randomizer;
import orc.consoperators.ConsOperatorUtils;
import orc.consoperators.InConstantDistanceOperator;
import orc.consoperators.PiecewiseLinearDistribution;
import starbeast3.GeneTreeForSpeciesTreeDistribution;
import starbeast3.genekernel.GTKPrior;


/**
 * Coordinated exchange rate operator
 * Applies narrow exchange rate on the species tree and changes gene tree topologies where a conflict is created
 * @author jdou557
 *
 */
public class CoordinatedExchangeRates extends InConstantDistanceOperator {
	
	

	public final  Input<ParametricDistribution> distributionInput = new Input<>("distr", "Rate distribution. Required if using the quantile parameterisation.");
	
	
	// Starbeast (optional)
	final public Input<RealParameter> popSizeInput = new Input<>("popSizes", "the constant population sizes associated with nodes in the tree.");
    final public Input<List<GeneTreeForSpeciesTreeDistribution>> geneTreeDistributionsInput = new Input<>("gene", "gene tree for species tree distribution for each of the genes", new ArrayList<>());
    final public Input<GTKPrior> geneTreeKernelPriorInput = new Input<>("geneTreeKernel", "the kernel of gene trees", Input.Validate.XOR, geneTreeDistributionsInput);
    
    private RealParameter popsizes;
    private List<GeneTreeForSpeciesTreeDistribution> geneTreeDistributions;
    private boolean proposeNewPopulationSizes;
    GTKPrior geneTreeKernelPrior;
	
    // Time proposal
    protected double tdp;
    
    // Rate proposals
    protected double rap, rbp, rcp, rdp;
    
    // Quantile proposals
    protected double qap, qbp, qcp, qdp;
    
    
    // Clock model
    protected enum ClockMode {
        quantiles,
        rates,
        strict
    }
    
    
    protected enum TreeNeighbours {
    	left,
    	right,
    }
    TreeNeighbours neighbourOfProposal;
    double[] neighbourSamplingProbabilitiesBefore;
    double[] neighbourSamplingProbabilitiesAfter;
    
    
    // Inputs
    private ClockMode clockMode;
    private RealParameter rates;
    private double twindowSize = 0;
    protected ParametricDistribution rateDistribution;
    private Tree tree;
    protected KernelDistribution proposalKernel;
    
    
    
    public CoordinatedExchangeRates() {
    	twindowSizeInput.setRule(Validate.OPTIONAL);
    	clockModelInput.setRule(Validate.FORBIDDEN);
    	
    	
    	// CNER operator: do not change rates
    	if (this.getClass().equals(CoordinatedExchangeRates.class)) {
    		rateInput.setRule(Validate.OPTIONAL);
    		quantileInput.setRule(Validate.OPTIONAL);
    	}
    	

    }
    
    
	@Override
	public void initAndValidate() {
		
		
		// Real rates or quantile parameterisation?
		if (rateInput.get() != null) {
			rates = rateInput.get();
			clockMode = ClockMode.rates;
		}else if (quantileInput.get() != null) { 
			rates = quantileInput.get();
			clockMode = ClockMode.quantiles;
			rateDistribution = distributionInput.get();
			if (rateDistribution == null) {
				throw new IllegalArgumentException("Please specify the rate distribution 'distr' when using the quantile parameterisation");
			}
		}else {
			rates = null;
			clockMode = ClockMode.strict;
		}
		
		
		// Window size, tree, and sub-operator (use this operator if unspecified)
		twindowSize = twindowSizeInput.get() == null ? 0 : twindowSizeInput.get();
		tree =  treeInput.get();
		proposalKernel = proposalKernelInput.get();
		
		
		// Starbeast. Any gene trees? Are population sizes being proposed?
		geneTreeDistributions = geneTreeDistributionsInput.get();
        proposeNewPopulationSizes = popSizeInput.get() != null && geneTreeDistributions != null && geneTreeDistributions.size() > 0;
        
        
        // Function is out of order due to mathematical errors in the proposal
        proposeNewPopulationSizes = false;
        
        if (proposeNewPopulationSizes) {
        	popsizes = popSizeInput.get();
        }
        
        
        // Get the gene tree distributions
        // If a gene tree kernel is being used, then initialise the gene trees here to those in the kernel
        geneTreeKernelPrior = geneTreeKernelPriorInput.get();
        geneTreeDistributions = geneTreeDistributionsInput.get();
        if (geneTreeKernelPrior != null) {
        	
        	// There should not be any genes added to the gene list
        	if (geneTreeDistributions.size() > 0) {
        		throw new IllegalArgumentException("Make sure you do not provide any gene trees when using a gene tree kernel");
        	}
        	geneTreeDistributions = geneTreeKernelPrior.getGeneTreeDistributions();
        }
        
		
		
	}

	@Override
	public double proposal() {
		
		
		// Get nodes which operator may apply to
        final List<Node> applicableNodesBeforeOperation = getApplicableNodes();
        if (applicableNodesBeforeOperation.size() == 0) return Double.NEGATIVE_INFINITY;
        
        
        
        // Sample a node E uniformly at random
        int Eindex = Randomizer.nextInt(applicableNodesBeforeOperation.size());
        Node E = applicableNodesBeforeOperation.get(Eindex);
        

		// Access to the child nodes of E: C and D. D must not be a leaf.
        int nodeDNum = Randomizer.nextInt(2);
        Node D = E.getChild(nodeDNum);
        Node C = E.getChild(1 - nodeDNum);
        if (D.getHeight() < C.getHeight()) {
        	nodeDNum = 1 - nodeDNum;
            D = E.getChild(nodeDNum);
            C = E.getChild(1 - nodeDNum);
        }
        
        
        // Tree with dated tips
        if (D.isLeaf()) {
            return Double.NEGATIVE_INFINITY;
        }
		
        
        
        // Tree guiding. Compute the scores of the trees obtained by:
        //	1) moving the left child of D
        //  2) moving the right child of D
        //  3) not changing anything
        // And sample one of these three from the scores
        neighbourOfProposal = this.sampleNeighbour(C, D, E);
        
        
        Node A = null;
        Node B = null;
        switch (neighbourOfProposal) {
        
	        case left: {
	        	
	        	// Perform narrow exchange on the 1st child of D
	        	A = D.getChild(1);
	        	B = D.getChild(0);
	        	break;
	        }
	        
	        case right: {
	        	
	        	// Perform narrow exchange on the 2nd child of D
	        	A = D.getChild(0);
	        	B = D.getChild(1);
	        	break;
	        }
	        	
        
        }
        
        
        
        // Get original node times
        double ta = A.getHeight(); // Fixed
        double tb = B.getHeight(); // Fixed
        double tc = C.getHeight(); // Fixed
        double td = D.getHeight(); // Free
        double te = E.getHeight(); // Fixed
        
        // Make the proposal
        double logJD = 0;
        switch (clockMode) {
        
	        case rates: {
	        	
	            // Get node rates 
	            double ra = rates.getValue(A.getNr()); // Free
	            double rb = rates.getValue(B.getNr()); // Free
	            double rc = rates.getValue(C.getNr()); // Free
	            double rd = rates.getValue(D.getNr()); // Free
	            
	            logJD = this.proposalRates(twindowSize, ta, tb, tc, td, te, ra, rb, rc, rd);
	            
	            // Ensure that the constraints have not been broken
	            if (!this.validateProposalRates(ta, tb, tc, td, te)) return Double.NEGATIVE_INFINITY;
	            break;
	        	
	        }
	        
	        case quantiles: {
	        	
	        	// Get node quantiles 
	            double qa = rates.getValue(A.getNr()); // Free
	            double qb = rates.getValue(B.getNr()); // Free
	            double qc = rates.getValue(C.getNr()); // Free
	            double qd = rates.getValue(D.getNr()); // Free
	            
	            logJD = this.proposalQuantiles(twindowSize, ta, tb, tc, td, te, qa, qb, qc, qd);
	            
	            // Ensure that the constraints have not been broken
	            if (!this.validateProposalQuantiles(ta, tb, tc, td, te)) return Double.NEGATIVE_INFINITY;
	            break;
	        }
	        
	        case strict:{
	        	this.tdp = td;
	        	break;
	        }
        
        }
        
        // Irreversible proposal
        if (logJD == Double.NEGATIVE_INFINITY) return Double.NEGATIVE_INFINITY;
        
        
        
        
        
        // Regrafting gene trees? Get the gene nodes which need to move
        List<List<Node>> nodesToRegraft_original = null;
        if (geneTreeDistributions != null && geneTreeDistributions.size() > 0) {
        	nodesToRegraft_original = getGeneNodesToRegraft(A, B, C, D, E);
        }
        
        
        // Proposing new population sizes? Calculate the contribution to log prior of the current system
        double partialLogPriorPopulation_A_original = 0;
        double partialLogPriorPopulation_C_original = 0;
        double partialLogPriorPopulation_D_original = 0;
        if (proposeNewPopulationSizes) {
        	
        	
        	
        	// The A, C, and D branches will be affected 
        	for (int g = 0; g < geneTreeDistributions.size(); g++) {
        		partialLogPriorPopulation_A_original += geneTreeDistributions.get(g).calculatePartialLogPBranch(A);
        		partialLogPriorPopulation_C_original += geneTreeDistributions.get(g).calculatePartialLogPBranch(C);
        		partialLogPriorPopulation_D_original += geneTreeDistributions.get(g).calculatePartialLogPBranch(D);
        	}
        	
        }
        
        

        // Rearrange the tree
        exchangeNodes(A, C, D, E);
        
        
        
        // Regraft gene trees and calculate Hastings ratio
        double logHastingsRatioGeneTrees = 0;
        if (geneTreeDistributions != null && geneTreeDistributions.size() > 0) {
        	logHastingsRatioGeneTrees = regraftGeneTrees(nodesToRegraft_original, A, D);
        }
        
        
        // Propose new population sizes
        double logHastingsRatioPopulationSize = 0;
        if (proposeNewPopulationSizes) {
        	
        	
        	// Recompute population sizes such that their tree partial contribution to the prior remains constant
        	double partialLogPriorPopulation_A_proposal = 0;
            double partialLogPriorPopulation_C_proposal = 0;
            double partialLogPriorPopulation_D_proposal = 0;
            
            
            // The A, C, and D branches have been affected by the proposal
        	for (int g = 0; g < geneTreeDistributions.size(); g++) {
        		partialLogPriorPopulation_A_proposal += geneTreeDistributions.get(g).calculatePartialLogPBranch(A);
        		partialLogPriorPopulation_C_proposal += geneTreeDistributions.get(g).calculatePartialLogPBranch(C);
        		partialLogPriorPopulation_D_proposal += geneTreeDistributions.get(g).calculatePartialLogPBranch(D);
        	}
            
            
        	Node X;
        	double partial_original, partial_proposal, scale, N, Np;
        	for (int i = 0; i < 3; i ++) {
        		
        		// The A branch has become longer and may have more coalescent events. 
        		if (i == 0) {
        			X = A;
        			partial_original = partialLogPriorPopulation_A_original;
        			partial_proposal = partialLogPriorPopulation_A_proposal;
        		}
        	
        		// The C branch has become shorter and may have fewer coalescent events. 
        		else if (i == 1) {
        			X = C;
        			partial_original = partialLogPriorPopulation_C_original;
        			partial_proposal = partialLogPriorPopulation_C_proposal;
        		}
        		
        		// The D branch has the same length (unless subjected to a random walk) but different coalescent events. 
        		else {
        			X = D;
        			partial_original = partialLogPriorPopulation_D_original;
        			partial_proposal = partialLogPriorPopulation_D_proposal;
        		}
        		
        		
        		// Scale the population size by the ratio between the old and new partial prior contribution
        		if (partial_original == 0 || partial_proposal == 0) scale = 1;
        		else scale = partial_proposal / partial_original;
        		
        		
        		N = popsizes.getArrayValue(X.getNr());
        		Np = N * scale;
        		popsizes.setValue(X.getNr(), Np);
        		logHastingsRatioPopulationSize += Math.log(scale);
        		
        		
        	}
            
            
        }
        
        
      
        
        
        // Set the new times + rates/quantiles
        D.setHeight(this.tdp);
        
        switch (clockMode) {
	        case rates: {
	        	
	            // Set the new rates
	        	rates.setValue(A.getNr(), this.rap);
	        	rates.setValue(B.getNr(), this.rbp);
	        	rates.setValue(C.getNr(), this.rcp);
	        	rates.setValue(D.getNr(), this.rdp);
	        	break;
	        }
	        case quantiles: {
	        	
	        	// Set the new quantiles
	        	rates.setValue(A.getNr(), this.qap);
	        	rates.setValue(B.getNr(), this.qbp);
	        	rates.setValue(C.getNr(), this.qcp);
	        	rates.setValue(D.getNr(), this.qdp);
	        	break;
	        	
	        }
	        
	        case strict:{
	        	break;
	        }
    
        }
        

        
        
        // Hastings ratio
        double proposalForward = 0;
        double proposalBackward = 0;
        final List<Node> applicableNodesAfterOperation = getApplicableNodes();
    	proposalForward = -Math.log(applicableNodesBeforeOperation.size());
    	proposalBackward = -Math.log(applicableNodesAfterOperation.size());
        double logHastingsRatioExchange = proposalBackward - proposalForward;
		
        
        

        
		//return Double.POSITIVE_INFINITY;
        return logHastingsRatioExchange + logJD + logHastingsRatioGeneTrees + logHastingsRatioPopulationSize;
        
	}
	
	
	


	private TreeNeighbours sampleNeighbour(Node C, Node D, Node E) {
		if (Randomizer.nextBoolean()) return TreeNeighbours.left;
		return TreeNeighbours.right;
	}
    
    
    
    


    // p1 becomes the parent of c2 
    // p2 becomes the parent of c1
	private void exchangeNodes(Node c1, Node c2, Node p1, Node p2) {
        replace(p1, c1, c2);
        replace(p2, c2, c1);
    }

	
	
	// Returns a list of all nodes which this operator may apply to
    private List<Node> getApplicableNodes(){
    	
    	
    	// Get all internal/root nodes which have grandchildren
    	List<Node> nodes = new ArrayList<Node>();

    	
    	
        
        // Iterate through all internal nodes
        for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
        	
        	Node E = tree.getNode(i);
        	
        	// Ensure that at least 1 child is not a leaf
        	if (E.getChildCount() == 2 && (!E.getChild(0).isLeaf() || !E.getChild(1).isLeaf())) {
        		nodes.add(E);
        	}
        	
        }
    	
    	return nodes;
    	
    }
    
    
    
    // Sample a random walk step size
    protected double getRandomWalkStepSize(double windowSize) {
    	
    	
    	// Sample a random walk step size
    	if (proposalKernel != null) return proposalKernel.getRandomDelta(1, windowSize);
    	
    	// No random walk
    	return 0;
    	
    	
    }
    
	
	
	
	/**
	 * Makes proposals: td', ra', rb', rc', and rd' and returns the hastings ratio
	 * Assumes real rate parameterisation
	 * The proposal variables tdp, rap, rbp, rcp, and rdp contain the new values
	 * @param rWindowSize - proposal width size
	 * @param ta - original time of node A
	 * @param tb - original time of node B
	 * @param tc - original time of node C
	 * @param td - original time of node D
	 * @param te - original time of node E
	 * @param ra - original rate of node A
	 * @param rb - original rate of node B
	 * @param rc - original rate of node C
	 * @param rd - original rate of node D
	 * @return logJD - the natural logarithm of the determinant of the Jacobian
	 */
	protected double proposalRates(double rWindowSize, double ta, double tb, double tc, double td, double te, 
												  double ra, double rb, double rc, double rd) {
		
		
		double r_td = this.getRandomWalkStepSize(rWindowSize); 
		
		
		// Null proposal: keep all rates constant. Keep tD constant unless a random walk is applied.
		this.tdp = td + r_td;
		this.rap = ra;
		this.rbp = rb;
		this.rcp = rc;
		this.rdp = rd;
		return 0;
		
	}
	
	/**
	 * Makes proposals: td', qa', qb', qc', and qd' and returns the hastings ratio
	 * Assumes quantile parameterisation
	 * The proposal variables tdp, qap, qbp, qcp, and qdp contain the new values
	 * @param rWindowSize - proposal width size
	 * @param ta - original time of node A
	 * @param tb - original time of node B
	 * @param tc - original time of node C
	 * @param td - original time of node D
	 * @param te - original time of node E
	 * @param qa - original quantile of node A
	 * @param qb - original quantile of node B
	 * @param qc - original quantile of node C
	 * @param qd - original quantile of node D
	 * @return logJD - the natural logarithm of the determinant of the Jacobian
	 */
	protected double proposalQuantiles(double rWindowSize, double ta, double tb, double tc, double td, double te, 
			  double qa, double qb, double qc, double qd) {
	
	
		double logJD = 0;
		try {
			
			
			// Convert quantiles to rates using the i-cdf
			double ra = rateDistribution.inverseCumulativeProbability(qa);
			double rb = rateDistribution.inverseCumulativeProbability(qb);
			double rc = rateDistribution.inverseCumulativeProbability(qc);
			double rd = rateDistribution.inverseCumulativeProbability(qd);
			
			
			// Propose new rates + times
			logJD = this.proposalRates(rWindowSize, ta,  tb,  tc, td, te, ra, rb, rc, rd);
			if (logJD == Double.NEGATIVE_INFINITY) return logJD;
			
			
			// Convert proposed rates into to proposed quantiles using the cdf
			this.qap = rateDistribution.cumulativeProbability(this.rap);
			this.qbp = rateDistribution.cumulativeProbability(this.rbp);
			this.qcp = rateDistribution.cumulativeProbability(this.rcp);
			this.qdp = rateDistribution.cumulativeProbability(this.rdp);
			
			
			// Contribution of the cdf to the hastings ratio
			if (this.qap != qa) logJD += getQuantileHastingsRatioContribution(this.rap, qa, this.qap);
			if (this.qbp != qb) logJD += getQuantileHastingsRatioContribution(this.rbp, qb, this.qbp);
			if (this.qcp != qc) logJD += getQuantileHastingsRatioContribution(this.rcp, qc, this.qcp);
			if (this.qdp != qd) logJD += getQuantileHastingsRatioContribution(this.rdp, qd, this.qdp);

		
		} catch (Exception e) {
			return Double.NEGATIVE_INFINITY;
		}
		
		return logJD;
		
		
	}
	
	
	
	// Get the contribution of a q -> r -> r' -> q' transformation to the hastings ratio
	protected double getQuantileHastingsRatioContribution(double rNew, double qOld, double qNew) {
		
		
		double logHR = 0;
		if (rateDistribution instanceof LogNormalDistributionModel) {
			logHR = ConsOperatorUtils.getHRForLN(rNew, qOld, rateDistribution);
        }

        else if (rateDistribution instanceof PiecewiseLinearDistribution) {
        	logHR = ConsOperatorUtils.getHRForPieceWise(rNew, qOld, qNew, rateDistribution);
        }

        else {
        	logHR = ConsOperatorUtils.getHRUseNumericApproximation(rNew, qOld, rateDistribution);
        }
		
		return logHR;
		
		
	}
	

	
    @Override
    public void optimize(double logAlpha) {
    	 double delta = calcDelta(logAlpha);
         delta += Math.log(twindowSize);
         twindowSize = Math.exp(delta);
    }
	
	
	
	// Validate proposal for rates
	protected boolean validateProposalRates(double ta, double tb, double tc, double td, double te) {
		
		if (this.tdp > te || this.tdp < tc || this.tdp < tb || this.tdp < ta) return false;
		if (this.rap <= 0) return false;
		if (this.rbp <= 0) return false;
		if (this.rcp <= 0) return false;
		if (this.rdp <= 0) return false;
		return true;
		
	}
	
	
	// Validate proposal for quantiles
	protected boolean validateProposalQuantiles(double ta, double tb, double tc, double td, double te) {
		
		if (this.tdp > te || this.tdp < tc || this.tdp < tb || this.tdp < ta) return false;
		if (this.qap <= 0 || this.qap >= 1) return false;
		if (this.qbp <= 0 || this.qbp >= 1) return false;
		if (this.qcp <= 0 || this.qcp >= 1) return false;
		if (this.qdp <= 0 || this.qdp >= 1) return false;
		
		
		// Ensure that proposed quantiles are not associated with rates which go outside the rate boundaries
		if (rateDistribution instanceof PiecewiseLinearDistribution) {
            PiecewiseLinearDistribution piecewise = (PiecewiseLinearDistribution) rateDistribution;
            
            try {
            
	            double rmin = piecewise.getRangeMin();
	            double rmax = piecewise.getRangeMax();
	            if (this.rap < rmin || this.rap > rmax) return false;
	            if (this.rbp < rmin || this.rbp > rmax) return false;
	            if (this.rcp < rmin || this.rcp > rmax) return false;
	            if (this.rdp < rmin || this.rdp > rmax) return false;
	            
            } catch (Exception e) {
            	e.printStackTrace();
            	return false;
            }
        }
		
		
		return true;
		
	}
	
	


	// Returns all gene tree branches which must be regrafted
	private List<List<Node>> getGeneNodesToRegraft(Node A, Node B, Node C, Node D, Node E) {
		
		
		List<List<Node>> nodesToRegraft = new ArrayList<List<Node>>();
		
	 	// Find all gene tree nodes that were in species D and need moving
		for (int g = 0; g < geneTreeDistributions.size(); g++) {
			
			
			List<Node> nodes = new ArrayList<Node>();
			
			// Find gene nodes from this gene tree which map to D
			GeneTreeForSpeciesTreeDistribution geneTree = geneTreeDistributions.get(g);
			Node[] geneNodeMap_D = geneTree.mapSpeciesNodeToGeneTreeNodes(D);
			Node[] geneNodeMap_A = geneTree.mapSpeciesNodeToGeneTreeNodes(A);
			
			// Get the gene nodes which have exactly 1 child which only have descendents in B
			for (int i = 0; i < geneNodeMap_D.length; i++) {
				Node geneNode = geneNodeMap_D[i];
				if (geneNodeNeedsToMove(geneNode, geneNodeMap_A)) nodes.add(geneNode);
			}
			
			nodesToRegraft.add(nodes);
			
			
		}
		
		return nodesToRegraft;
		
	}


	// Checks whether the gene node needs to be regrafted by checking that 
	// exactly one child of 'geneNode' has 1 or more descendent in A (and the other has no descendents in A)
	private boolean geneNodeNeedsToMove(Node geneNode, Node[] geneNodeMap_A) {
		
		
		boolean leftChildHasDescendentsInA = false;
		boolean rightChildHasDescendentsInA = false;
		for (int i = 0; i < 2; i ++) {
			
			Node child = geneNode.getChild(0);
			
			if (i == 0) leftChildHasDescendentsInA = geneSubtreeMapsToSpecies(child, geneNodeMap_A);
			else rightChildHasDescendentsInA = geneSubtreeMapsToSpecies(child, geneNodeMap_A);
			
		}
		
		
		
		// Return the XOR
		return (leftChildHasDescendentsInA || rightChildHasDescendentsInA) && !(leftChildHasDescendentsInA && rightChildHasDescendentsInA);
	}

	
	// Check whether any of the nodes in this subtree are in the node list
	private boolean geneSubtreeMapsToSpecies(Node geneSubtree, Node[] geneNodeMap) {
		
		
		
		// Check if this node is in the list
		for (int i = 0; i < geneNodeMap.length; i ++) {
			if (geneNodeMap[i].getNr() == geneSubtree.getNr()) return true;
		}
		
		
		// Recurse to children
		for (int c = 0; c < geneSubtree.getChildCount(); c++) {
			if (geneSubtreeMapsToSpecies(geneSubtree.getChild(c), geneNodeMap)) return true;
		}
		
		
		return false;
		
	}
	
	
	// Sample locations and regraft all of the gene nodes in nodesToRegraft_original
	// And calculate the hastings ratio
	// This is called after the exchange has occurred, so node A is now adjacent to node E
    private double regraftGeneTrees(List<List<Node>> nodesToRegraft, Node A, Node D) {
    	
    	double logHR = 0;
    	for (int g = 0; g < nodesToRegraft.size(); g++) {
    		
    		List<Node> nodesToMove_g = nodesToRegraft.get(g);
    		GeneTreeForSpeciesTreeDistribution geneTree = geneTreeDistributions.get(g);
    		
    		Node[] geneNodeMap_A = geneTree.mapSpeciesNodeToGeneTreeNodes(A);
    		Node[] geneNodeMap_D = geneTree.mapSpeciesNodeToGeneTreeNodes(D);
    		List<Node> destinationBranches_forward;
    		Node candidate, destination;
    		for (int i = 0; i < nodesToMove_g.size(); i ++) {
    			
    			Node nodeToMove = nodesToMove_g.get(i);
    			
    			
    			
    			
    			// Find the places the branch can move to (forward)
    			destinationBranches_forward = new ArrayList<Node>();
    			for (int j = 0; j < geneNodeMap_D.length; j ++) {
    				candidate = geneNodeMap_D[j];
    				if (candidate.getHeight() < nodeToMove.getHeight() && candidate.getParent().getHeight() > nodeToMove.getHeight()) {
    					destinationBranches_forward.add(candidate);
    				}
    			}
    			int numberOfDestinations_forward = destinationBranches_forward.size();
    			
    			
    			// Find the [number of] places the branch can move to (reverse)
    			int numberOfDestinations_reverse = 0;
    			for (int j = 0; j < geneNodeMap_A.length; j ++) {
    				candidate = geneNodeMap_A[j];
    				if (candidate.getHeight() < nodeToMove.getHeight() && candidate.getParent().getHeight() > nodeToMove.getHeight()) {
    					numberOfDestinations_reverse ++;
    				}
    			}
    			
    			
    			// Hastings ratio
    			logHR += Math.log(numberOfDestinations_forward) + Math.log(numberOfDestinations_reverse);
    			
    			
    			// Sample a branch to move to
    			destination = destinationBranches_forward.get(Randomizer.nextInt(numberOfDestinations_forward));
    			
    			
    			// Regraft the branch to its randomly sampled destination
    			exchangeNodes(nodeToMove, nodeToMove.getParent(), destination, destination.getParent());
    			
    			
    		}
    		
    		
    	}
    	
		return logHR;
	}

	
	
	
	
}












