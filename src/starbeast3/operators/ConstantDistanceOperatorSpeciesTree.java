package starbeast3.operators;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.UCRelaxedClockModel;
import beast.evolution.operators.KernelDistribution;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.evolution.tree.TreeInterface;
import beast.math.distributions.PiecewiseLinearDistribution;
import beast.util.Randomizer;
import starbeast3.GeneTreeForSpeciesTreeDistribution;
import starbeast3.StarBeast3Clock;
import starbeast3.evolution.branchratemodel.BranchRateModelSB3;
import starbeast3.evolution.branchratemodel.UCRelaxedClockModelSB3;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.MathException;

@Description("For internal nodes: propose a new node time")
public class ConstantDistanceOperatorSpeciesTree extends TreeOperator {
	final public Input<Double> twindowSizeInput = new Input<>("twindowSize", "the size of the window when proposing new node time", Input.Validate.REQUIRED);
    //public final Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            //"A model describing the rates on the branches of the beast.tree.");
   // final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.");
    //final public Input<RealParameter> quantilesInput = new Input<>("quantiles", "The real branch rates of the species tree, parameterised as quantiles", Input.Validate.XOR, rateInput);
    final public Input<RealParameter> popSizeInput = new Input<>("popsizes", "the constant population sizes associated with nodes in the tree.");
    final public Input<List<GeneTreeForSpeciesTreeDistribution>> geneTreeDistributionsInput = new Input<>("gene", "gene tree for species tree distribution for each of the genes", new ArrayList<>());
    final public Input<Boolean> proportionalToBranchLengthInput = new Input<>("proportionalToBranchLength", "Set proposal step sizes proportional to branch length (true) or a constant (false)", false);
    final public Input<UCRelaxedClockModelSB3> clockModelInput = new Input<>("clock", "the relaxed clock model associated with species tree brancg rates.", Input.Validate.REQUIRED);
    final public Input<KernelDistribution> proposalKernelInput = new Input<>("kernelDistribution", "Proposal kernel for a random walk on the internal node height.");
	
    
    
    UCRelaxedClockModelSB3 clockModel;
    
    private double twindowSize;
    private RealParameter rates;
    private RealParameter quantiles;
    private RealParameter popsizes;
    private List<GeneTreeForSpeciesTreeDistribution> geneTreeDistributions;
    private boolean proposeNewPopulationSizes;
    
    private Node[] geneNodeMap_x;
    private Node[] geneNodeMap_L;
    private Node[] geneNodeMap_R;
    
    
    // Quantiles
    PiecewiseLinearDistribution piecewise = null;
    
    // Proposal kernel
    private KernelDistribution kernel;
   

    //protected BranchRateModel.Base branchRateModel;
    //protected UCRelaxedClockModel branchRateModel;
    //JacobianMatrixDeterminant JD = new JacobianMatrixDeterminant();

    @Override
    public void initAndValidate() {
        twindowSize = twindowSizeInput.get();
        clockModel = clockModelInput.get();
        
        
        // Ensure that either rates or quantiles are used and not categories
        if (clockModel.getRateMode() != UCRelaxedClockModelSB3.Mode.rates && clockModel.getRateMode() != UCRelaxedClockModelSB3.Mode.quantiles) {
        	throw new IllegalArgumentException("Clock model must parameterise rates as real numbers or quantiles! This operator does not work with categories.");
        }
        
        
        // Get rate parameter
        if (clockModel.getRateMode() == UCRelaxedClockModelSB3.Mode.rates) {
        	rates = clockModel.realRatesInput.get();
        }else {
        	quantiles = clockModel.quantilesInput.get();
        }
        
        
        // Are population sizes being proposed?
        proposeNewPopulationSizes = popSizeInput.get() != null;
        if (proposeNewPopulationSizes) popsizes = popSizeInput.get();
       
        
        // Get the gene tree distributions
        geneTreeDistributions = geneTreeDistributionsInput.get();
        
        kernel = proposalKernelInput.get();
        
        
    }

    @Override
    public double proposal() {
    	
    	
    	
        final Tree tree = treeInput.get(this);
        int nodeCount = tree.getNodeCount(); //return the number of nodes in the tree
        int branchCount = nodeCount - 1; //the number of branches of the tree

        //the chosen node to work on
        Node node;

        // Original node times
        double t_x, t_L, t_R;
        
        // Original rates
        double r_x, r_R, r_L;

        //the proposed node time
        double t_x_;

       // Step 1: randomly select an internal, denoted by node x. x can be a root node as long as there are gene trees
       final int firstNonLeafNr = tree.getLeafNodeCount();
       final int lastNodeNr = geneTreeDistributions.size() == 0 ? nodeCount - 1 : nodeCount;
       final int nodeNr = firstNonLeafNr + Randomizer.nextInt(lastNodeNr - firstNonLeafNr);
       node = tree.getNode(nodeNr);
       
       
       //Step 2: Access to the child nodes of this node
       // Left child
       Node leftNode = node.getChild(0);//get the left child of this node
       t_L = leftNode.getHeight();//node time of son

       int leftNr = leftNode.getNr();// node number of son
       if (leftNr == branchCount) {
           leftNr = leftNode.getTree().getRoot().getNr();
        }
       
       
       // Right child
       Node rightNode = node.getChild(1);//get the right child of this node
       t_R = rightNode.getHeight();//node time of right child

       int rightNr = rightNode.getNr(); // node time of right child
       if (rightNr == branchCount) {
            rightNr = rightNode.getTree().getRoot().getNr();
       }
       
      
       // Original node times
       t_x = node.getHeight();
       
       
       

       // Rates
       switch(clockModel.getRateMode()) {
       
	       case rates: {
	    	   r_x = rates.getValues()[nodeNr];
	    	   r_L = rates.getValues()[leftNr]; // Rate of branch above left child
	    	   r_R = rates.getValues()[rightNr]; // Rate of branch above right child
	    	   break;
	       }
	       
	       case quantiles: {
	    	   piecewise = clockModel.getPiecewiseQuantileApproximation();
	    	   try {
					r_x = piecewise.inverseCumulativeProbability(quantiles.getValues()[nodeNr]);
					r_L = piecewise.inverseCumulativeProbability(quantiles.getValues()[leftNr]);
					r_R = piecewise.inverseCumulativeProbability(quantiles.getValues()[rightNr]);
				} catch (MathException e) {
					e.printStackTrace();
					return Double.NEGATIVE_INFINITY;
				}
	    	   break;
	       }
	       
	       default: {
	    	   return Double.NEGATIVE_INFINITY;
	       }
       
       }


       // Compute lower and upper bounds
       final double lower = Math.max(t_L, t_R);
       double upper = 0 ;
       if (node.isRoot()) {
    	   
    	   // If this is the root node then the upper limit is the maximum gene tree height
    	   for (int i = 0; i < geneTreeDistributions.size(); i ++) {
    		   upper = Math.max(upper, geneTreeDistributions.get(i).getGeneTree().getRoot().getHeight());
    	   }
    	   
       }else {
    	   upper = node.getParent().getHeight();
       }
       
       
       
       
       // Step3-4: propose a new node time for this node
       double alpha;
       if (kernel != null) alpha = kernel.getRandomDelta(twindowSize);
       else alpha = Randomizer.uniform(-twindowSize, twindowSize);
       if (proportionalToBranchLengthInput.get()) {
    	   
    	   // Proposal size is proportional to branch length
    	   double beta = alpha * (upper - lower);
    	   t_x_ = t_x + beta;
       }else {
    	   
    	   // Constant proposal width
    	   t_x_ = t_x + alpha;
       }
       
       
       
       // Reject the proposal if exceeds the boundary
       if (t_x_<= lower || t_x_ >= upper) {
            return Double.NEGATIVE_INFINITY;
        }
       


       // Step5: propose new rates - r_x, r_L, r_R
       double r_x_ = r_x * (upper - t_x) / (upper - t_x_);
       double r_L_ = r_L * (t_x - t_L) / (t_x_ - t_L);
       double r_R_ = r_R * (t_x - t_R) / (t_x_ - t_R);

       // Set the proposed new rates (and calculate the Jacobian contribution from quantiles if applicable)
       double logJD_quantiles = 0;
       switch(clockModel.getRateMode()) {
       
	       case rates: {
	           rates.setValue(nodeNr, r_x_);
	           rates.setValue(leftNr, r_L_);
	           rates.setValue(rightNr, r_R_);
	           logJD_quantiles = 0;
	    	   break;
	       }
	       
	       case quantiles: {
	    	   try {
	    		   
   
					// Ensure that proposed rates are within the piecewise approximation's range
					double rmin = piecewise.getRangeMin();
					double rmax = piecewise.getRangeMax();
					if (r_x_ <= rmin || r_x_ >= rmax) return Double.NEGATIVE_INFINITY;
					if (r_L_ <= rmin || r_L_ >= rmax) return Double.NEGATIVE_INFINITY;
					if (r_R_ <= rmin || r_R_ >= rmax) return Double.NEGATIVE_INFINITY;
					
					
					//System.out.println(quantiles.getValues()[nodeNr] + "," + r_x + "," + r_x_ + "," + r_L_ + "," + r_R_);
					
					
					// Calculate quantiles from rates
					double q_x_ = piecewise.cumulativeProbability(r_x_);
					double q_L_ = piecewise.cumulativeProbability(r_L_);
					double q_R_ = piecewise.cumulativeProbability(r_R_);
					   
					// Jacobian contribution from the icdf derivative
					double q_x = quantiles.getValues()[nodeNr];
					double q_L = quantiles.getValues()[leftNr];
					double q_R = quantiles.getValues()[rightNr];
					logJD_quantiles += Math.log(piecewise.getDerivativeAtQuantile(q_x));
					logJD_quantiles += Math.log(piecewise.getDerivativeAtQuantile(q_L));
					logJD_quantiles += Math.log(piecewise.getDerivativeAtQuantile(q_R));
					   
					// Jacobian contribution from the cdf derivative
					// double dqx = piecewise.getDerivativeAtQuantile(q_x_);
					//double drx = piecewise.getDerivativeAtQuantileInverse(r_x_, q_x_);
					logJD_quantiles += Math.log(piecewise.getDerivativeAtQuantileInverse(r_x_, q_x_));
					logJD_quantiles += Math.log(piecewise.getDerivativeAtQuantileInverse(r_L_, q_L_));
					logJD_quantiles += Math.log(piecewise.getDerivativeAtQuantileInverse(r_R_, q_R_));
					   
					// Set new quantiles
					quantiles.setValue(nodeNr, q_x_);
					quantiles.setValue(leftNr, q_L_);
					quantiles.setValue(rightNr, q_R_);
		    	   
				} catch (MathException e) {
					e.printStackTrace();
					return Double.NEGATIVE_INFINITY;
				}
	    	   break;
	    	   
	       }
	       
	       default: {
	    	   
	       }
	   
       }

       

       
       
       // Step6: propose new population sizes
       if (proposeNewPopulationSizes) {
    	   
	       double N_x = popsizes.getValues()[nodeNr];
	       double N_L = popsizes.getValues()[leftNr];
	       double N_R = popsizes.getValues()[rightNr];
	       
	
	       double N_x_ = N_x * (upper - t_x_) / (upper - t_x);
	       double N_L_ = N_L * (t_x_ - t_L) / (t_x - t_L);
	       double N_R_ = N_R * (t_x_ - t_R) / (t_x - t_R);
	       
	       popsizes.setValue(nodeNr, N_x_);
	       popsizes.setValue(leftNr, N_L_);
	       popsizes.setValue(rightNr, N_R_);
       
       }
       
       
       // Iterate through gene trees
       Node geneTreeNode;
       double t_g;
       double t_g_;
       
       
       // Count the number of nodes mapped to each branch (for computing Green ratio). Only count the nodes which change heights (ie. not the ultrametric leaves)
       int numNodesMappedX = 0;
       int numNodesMappedL = 0;
       int numNodesMappedR = 0;


       
       
       for (int i = 0; i < geneTreeDistributions.size(); i ++) {
    	   
    	   
    	   // Get the nodes in this gene tree which map to species tree branch x, L, and R 
    	   geneNodeMap_x = geneTreeDistributions.get(i).mapSpeciesNodeToGeneTreeNodes(node);
    	   geneNodeMap_L = geneTreeDistributions.get(i).mapSpeciesNodeToGeneTreeNodes(leftNode);
    	   geneNodeMap_R = geneTreeDistributions.get(i).mapSpeciesNodeToGeneTreeNodes(rightNode);
    	   

    	   
    	   
    	   /* -------------------------------
 	   	  	  ---------- Proposals ----------
 	          -------------------------------  */
    	   
    	   // Propose new time for each gene tree node which mapped to branch above x
    	   for (int j = 0; j < geneNodeMap_x.length; j ++) {
    		   geneTreeNode = geneNodeMap_x[j];
    		   t_g = geneTreeNode.getHeight();
    		   t_g_ = upper - (r_x / r_x_) * (upper - t_g);
    		   geneTreeNode.setHeight(t_g_);
    		   if(t_g != t_g_) numNodesMappedX ++;
    	   }
    	   
    	   
    	   //  Propose new time for each gene tree node which mapped to branch above L
    	   for (int j = 0; j < geneNodeMap_L.length; j ++) {
    		   geneTreeNode = geneNodeMap_L[j];
    		   t_g = geneTreeNode.getHeight();
    		   t_g_ = t_L + (r_L / r_L_) * (t_g - t_L);
    		   geneTreeNode.setHeight(t_g_);
    		   if(t_g != t_g_) numNodesMappedL ++;
    	   }
    	   
    	   
    	   //  Propose new time for each gene tree node which mapped to branch above R
    	   for (int j = 0; j < geneNodeMap_R.length; j ++) {
    		   geneTreeNode = geneNodeMap_R[j];
    		   t_g = geneTreeNode.getHeight();
    		   t_g_ = t_R + (r_R / r_R_) * (t_g - t_R);
    		   geneTreeNode.setHeight(t_g_);
    		   if(t_g != t_g_) numNodesMappedR ++;
    	   }
    	   
    	   

       }
       
       
       // Change the node time of the species tree AFTER changing the gene tree nodes 
       // or the species to gene tree mapper will no longer apply
       node.setHeight(t_x_);
       

       
       if (!proposeNewPopulationSizes) {
    	   numNodesMappedX -= 1;
    	   numNodesMappedL -= 1;
    	   numNodesMappedR -= 1;
       }
       
       
       // Calculate Green ratio
       double logJD = 0;
       logJD += numNodesMappedX * (Math.log(r_x) - Math.log(r_x_));
       logJD += numNodesMappedL * (Math.log(r_L) - Math.log(r_L_));
       logJD += numNodesMappedR * (Math.log(r_R) - Math.log(r_R_));
       logJD += logJD_quantiles;


 
       return logJD;
       

       
    }
    
    

    
    /*
    Tuning the parameter: twindowsize represents the range of Uniform distribution
     */
    @Override
    public double getCoercableParameterValue() {
        return twindowSize;
    }


    @Override
    public void setCoercableParameterValue(double value) {
        twindowSize = value;
    }

    /**
     * called after every invocation of this operator to see whether
     * a parameter can be optimised for better acceptance hence faster
     * mixing
     *
     * @param logAlpha difference in posterior between previous state & proposed state + hasting ratio
     */

    @Override
    public void optimize(double logAlpha) {
        // must be overridden by operator implementation to have an effect
        double delta = calcDelta(logAlpha);

        delta += Math.log(twindowSize);
        twindowSize = Math.exp(delta);
    }

    @Override
    public final String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = twindowSize * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else if (prob > 0.40) {
            return "Try setting window size to about " + formatter.format(newWindowSize);
        } else return "";
    }


}

