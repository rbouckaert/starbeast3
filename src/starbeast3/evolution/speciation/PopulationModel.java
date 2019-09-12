package starbeast3.evolution.speciation;

import java.text.DecimalFormat;
import beast.evolution.tree.Node;

public interface PopulationModel {
	
	
    // Calculate the truncated coalescent probability for a single species tree branch and gene
    double calculateBranchLogP(final int lineagesBottom, final double ploidy, final double popSize2, final double[] times, final int k);

    // Sets model-compatible default population sizes
    // To successfully begin a run, this must be called from a StateNodeInitializer
    void initPopSizes(final double initialPopSizes);

    // Per-branch population size information which will be added to a Newick string.
    void serialize(Node speciesTreeNode, StringBuffer buf, DecimalFormat df);

    // Checks if the population function for a given branch has been changed 
    boolean isDirtyBranch(Node speciesTreeNode);
    
}
