package starbeast3.operators;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;
import beast.evolution.branchratemodel.BranchRateModel;
import beast.evolution.branchratemodel.UCRelaxedClockModel;
import beast.evolution.operators.TreeOperator;
import beast.evolution.tree.Node;
import beast.evolution.tree.Tree;
import beast.util.Randomizer;
import starbeast3.GeneTreeForSpeciesTreeDistribution;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

@Description("For internal nodes: propose a new node time")
public class ConstantDistanceOperatorSpeciesTree extends TreeOperator {
	final public Input<Double> twindowSizeInput = new Input<>("twindowSize", "the size of the window when proposing new node time", Input.Validate.REQUIRED);
    //public final Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            //"A model describing the rates on the branches of the beast.tree.");
    final public Input<RealParameter> rateInput = new Input<>("rates", "the rates associated with nodes in the tree for sampling of individual rates among branches.", Input.Validate.REQUIRED);
    final public Input<RealParameter> popSizeInput = new Input<>("popsizes", "the constant population sizes associated with nodes in the tree.", Input.Validate.REQUIRED);
    final public Input<CompoundDistribution> coalescentModelInput = new Input<>("coalescentModel", "The coalescent model containing gene tree priors.", Input.Validate.REQUIRED);
    
    
    private double twindowSize;
    private RealParameter rates;
    private RealParameter popsizes;
    private List<GeneTreeForSpeciesTreeDistribution> geneTreeDistributions;
    
    
    private Node[] geneNodeMap_x;
    private Node[] geneNodeMap_L;
    private Node[] geneNodeMap_R;

    //protected BranchRateModel.Base branchRateModel;
    //protected UCRelaxedClockModel branchRateModel;
    //JacobianMatrixDeterminant JD = new JacobianMatrixDeterminant();

    @Override
    public void initAndValidate() {
        twindowSize = twindowSizeInput.get();
        //branchRateModel = branchRateModelInput.get();
        rates = rateInput.get();
        popsizes = popSizeInput.get();
        
        
        // Get the gene tree distributions
        geneTreeDistributions = new ArrayList<GeneTreeForSpeciesTreeDistribution>();
        List<Distribution> distributions = coalescentModelInput.get().pDistributions.get();
        for (int i = 0; i < distributions.size(); i ++) {
        	if (distributions.get(i) instanceof GeneTreeForSpeciesTreeDistribution) {
        		geneTreeDistributions.add((GeneTreeForSpeciesTreeDistribution) distributions.get(i));
        	}
        }
        
        
    }

    @Override
    public double proposal() {
    	
    	
        final Tree tree = treeInput.get(this);
        int nodeCount = tree.getNodeCount(); //return the number of nodes in the tree
        int branchCount = nodeCount - 1; //the number of branches of the tree

        //the chosen node to work on
        Node node;

        
        double t_P;
        
        //the original node times
        double t_x;
        double t_L;
        double t_R;
        //the original distances
        double d_L;
        double d_R;
        //the original rates
        double r_R;
        double r_L;

        //the proposed node time
        double t_x_;

       // /Step 1: randomly select an internal node, denoted by node x
       do {
            final int nodeNr = nodeCount / 2 + 1 + Randomizer.nextInt(nodeCount / 2);
            node = tree.getNode(nodeNr);
       } while (node.isRoot() || node.isLeaf());

       // the number of this node
        int nodeNr = node.getNr();
        if (nodeNr == branchCount) {
            nodeNr = node.getTree().getRoot().getNr();
        }

       //rate and time for this node
       t_x = node.getHeight();
       double r_x = rates.getValues()[nodeNr];
       t_P = node.getParent().getHeight();
       

       //Step 2: Access to the child nodes of this node
       // Left child
       Node leftNode = node.getChild(0);//get the left child of this node
       t_L = leftNode.getHeight();//node time of son

       int leftNr = leftNode.getNr();// node number of son
       if (leftNr == branchCount) {
           leftNr = leftNode.getTree().getRoot().getNr();
        }

       r_L = rates.getValues()[leftNr]; // rate of branch above left child
       d_L = r_L * (t_x - t_L); // distance of branch above left child


       // Right child
       Node rightNode = node.getChild(1);//get the right child of this node
       t_R = rightNode.getHeight();//node time of right child

       int rightNr = rightNode.getNr(); // node time of right child
       if (rightNr == branchCount) {
            rightNr = rightNode.getTree().getRoot().getNr();
       }

       r_R = rates.getValues()[rightNr];// rate of branch above right child
       d_R = r_R * (t_x - t_R);// distance of branch above right child



       //Step3-4: to propose a new node time for this node
       double alpha = Randomizer.uniform(-twindowSize, twindowSize);
       t_x_ = t_x + alpha;
       
       //double beta = alpha * tree.getRoot().getHeight();
       //t_x_ = t_x + beta;

       //reject the proposal if exceeds the boundary
       double upper = node.getParent().getHeight();
       double lower = Math.max(t_L, t_R);

       if (t_x_<= lower || t_x_ >= upper) {
            return Double.NEGATIVE_INFINITY;
        }
        node.setHeight(t_x_);


       //Step5: propose the new rates
       //there are three rates in total
       //r_x, r_L, r_R
       double r_x_ = r_x * (upper - t_x) / (upper - t_x_);
       double r_L_ = d_L / (t_x_ - t_L);
       double r_R_ = d_R / (t_x_ - t_R);

       // set the proposed new rates
       rates.setValue(nodeNr, r_x_);
       rates.setValue(leftNr, r_L_);
       rates.setValue(rightNr, r_R_);
       
       
       // Step6: propose new population sizes
       double N_x = popsizes.getValues()[nodeNr];
       double N_L = popsizes.getValues()[leftNr];
       double N_R = popsizes.getValues()[rightNr];
       
       double N_x_ = N_x * (upper - t_x_) / (upper - t_x);
       double N_L_ = N_L * (t_x_ - t_L) / (t_x - t_L);
       double N_R_ = N_R * (t_x_ - t_R) / (t_x - t_R);
       
       popsizes.setValue(nodeNr, N_x_);
       popsizes.setValue(leftNr, N_L_);
       popsizes.setValue(rightNr, N_R_);
       
       
       
       
       
       
       // Iterate through gene trees
       Node geneTreeNode;
       double t_g;
       double t_g_;
       int numNodesMappedX = 0;
       int numNodesMappedL = 0;
       int numNodesMappedR = 0;
       for (int i = 0; i < geneTreeDistributions.size(); i ++) {
    	   
    	   
    	   // Get the nodes in this gene tree which map to species tree branch x, L, and R
    	   geneNodeMap_x = geneTreeDistributions.get(i).mapSpeciesNodeToGeneTreeNodes(node);
    	   geneNodeMap_L = geneTreeDistributions.get(i).mapSpeciesNodeToGeneTreeNodes(leftNode);
    	   geneNodeMap_R = geneTreeDistributions.get(i).mapSpeciesNodeToGeneTreeNodes(rightNode);
    	   
    	   
    	   // Count the number of nodes mapped to each branch (for computing Green ratio)
    	   numNodesMappedX += geneNodeMap_x.length;
    	   numNodesMappedL += geneNodeMap_L.length;
    	   numNodesMappedR += geneNodeMap_R.length;
    	   
    	   
    	   
    	   // Propose new time for each gene tree node which mapped to branch above x
    	   for (int j = 0; j < geneNodeMap_x.length; j ++) {
    		   geneTreeNode = geneNodeMap_x[j];
    		   t_g = geneTreeNode.getHeight();
    		   t_g_ = upper - (r_x / r_x_) * (upper - t_g);
    		   
    		   
    		   double d_g = r_x * (t_P - t_g);
    		   double d_g_ = r_x_ * (t_P - t_g_);
    		   
    		   if (Math.abs(d_g - d_g_) > 0.000001) {
    			   System.out.println("X distance has changed from " + d_g + " to " + d_g);    		   
    		   }
    		   
    		   geneTreeNode.setHeight(t_g_);
    	   }
    	   
    	   
    	   //  Propose new time for each gene tree node which mapped to branch above L
    	   for (int j = 0; j < geneNodeMap_L.length; j ++) {
    		   geneTreeNode = geneNodeMap_L[j];
    		   t_g = geneTreeNode.getHeight();
    		   t_g_ = t_L + (r_L / r_L_) * (t_g - t_L);
    		   
    		   
    		   double d_g, d_g_;
    		   double d_g_p = geneTreeNode.getParent().getHeight();
    		   if (false && d_g_p > t_x) {
    			   d_g = 0;
    			   d_g_ = 0;
    			   
    		   }else {
    			   d_g = r_L * (t_g - t_L);
    			   d_g_ = r_L_ * (t_g_ - t_L);
    		   }

    		   
    		   if (Math.abs(d_g - d_g_) > 0.000001) {
    			   System.out.println("L distance has changed from " + d_g + " to " + d_g);    		   
    		   }
    		   
    		   geneTreeNode.setHeight(t_g_);
    	   }
    	   
    	   
    	   //  Propose new time for each gene tree node which mapped to branch above R
    	   for (int j = 0; j < geneNodeMap_R.length; j ++) {
    		   geneTreeNode = geneNodeMap_R[j];
    		   t_g = geneTreeNode.getHeight();
    		   t_g_ = t_R + (r_R / r_R_) * (t_g - t_R);
    		   
    		   
    		   
    		   double d_g, d_g_;
    		   double d_g_p = geneTreeNode.getParent().getHeight();
    		   if (false && d_g_p > t_x) {
    			   d_g = 0;
    			   d_g_ = 0;
    			   
    		   }else {
    			   d_g = r_R * (t_g - t_R);
    			   d_g_ = r_R_ * (t_g_ - t_R);
    		   }

    		   
    		   if (Math.abs(d_g - d_g_) > 0.000001) {
    			   System.out.println("R distance has changed from " + d_g + " to " + d_g);    		   
    		   }
    		   
    		   
    		   geneTreeNode.setHeight(t_g_);
    	   }
    	   

       }


       
       // Calculate Green ratio
       double logJD = 0;
       logJD += numNodesMappedX * (Math.log(r_x) - Math.log(r_x_));
       logJD += numNodesMappedL * (Math.log(r_L) - Math.log(r_L_));
       logJD += numNodesMappedR * (Math.log(r_R) - Math.log(r_R_));
       
       //logJD += Math.log(tree.getRoot().getHeight());
       
       //return 0;
       return logJD;
    }
    
    
    //private double getGeneticDistance(Node node) {
    	
    	
    	//return 0;
    	
    //}
    
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

