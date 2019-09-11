package starbeast3;

import java.util.List;
import java.util.Random;

import beast.core.Description;
import beast.core.Distribution;
import beast.core.Input;
import beast.core.parameter.RealParameter;
import beast.core.util.CompoundDistribution;
import beast.core.State;

/**
* @author Remco Bouckaert
* @author Joseph Heled
* @author Huw Ogilvie
 */

@Description("Calculates probability of gene trees conditioned on a species tree (the multi-species coalescent).")
public class MultispeciesCoalescent extends CompoundDistribution {
    final public Input<RealParameter> populationShapeInput = new Input<>("populationShape", "Shape of the inverse gamma prior distribution on population sizes.");
    final public Input<RealParameter> populationMeanInput = new Input<>("populationMean", "Mean of the inverse gamma prior distribution on population sizes.");

    private RealParameter invGammaShape;
    private RealParameter invGammaMean;

    private int nGeneTrees;
    private int speciesNodeCount;
    private double[] perGenePloidy;

    private double alpha;
    private double beta;

    private double storedAlpha;
    private double storedBeta;

    private int[] allLineageCounts;
    private int[] allEventCounts;
    private double[][] allCoalescentTimes;
    private double[] perBranchLogP;

    private int[] storedLineageCounts;
    private int[] storedEventCounts;
    private double[][] storedCoalescentTimes;
    private double[] storedPerBranchLogP;

    private boolean dontCalculate;

    @Override
    public void store() {
        super.store();
        if (dontCalculate) return;

        storedAlpha = alpha;
        storedBeta = beta;

        System.arraycopy(allLineageCounts, 0, storedLineageCounts, 0, allLineageCounts.length);
        System.arraycopy(allEventCounts, 0, storedEventCounts, 0, allEventCounts.length);
        System.arraycopy(allCoalescentTimes, 0, storedCoalescentTimes, 0, allCoalescentTimes.length);
        System.arraycopy(perBranchLogP, 0, storedPerBranchLogP, 0, perBranchLogP.length);
    }

    @Override
    public void restore() {
        super.restore();
        if (dontCalculate) return;

        double tmpAlpha = alpha;
        double tmpBeta = beta;
        int[] tmpLineageCounts = allLineageCounts;
        int[] tmpEventCounts = allEventCounts;
        double[][] tmpCoalescentTimes = allCoalescentTimes;
        double[] tmpPerBranchLogP = perBranchLogP;

        alpha = storedAlpha;
        beta = storedBeta;
        allLineageCounts = storedLineageCounts;
        allEventCounts = storedEventCounts;
        allCoalescentTimes = storedCoalescentTimes;
        perBranchLogP = storedPerBranchLogP;

        storedAlpha = tmpAlpha;
        storedBeta = tmpBeta;
        storedLineageCounts = tmpLineageCounts;
        storedEventCounts = tmpEventCounts;
        storedCoalescentTimes = tmpCoalescentTimes;
        storedPerBranchLogP = tmpPerBranchLogP;
    }

    @Override
    public void initAndValidate() {
        super.initAndValidate();

        if (populationShapeInput.get() == null ^ populationMeanInput.get() == null) {
            throw new IllegalArgumentException("Either specify both population size prior parameters for analytical integration,"
                    + "or neither for MCMC integration of population sizes.");
        } else if (populationShapeInput.get() == null) {
            dontCalculate = true;
            return;
        }

        final List<Distribution> geneTrees = pDistributions.get();

        dontCalculate = false;
        checkHyperparameters(true);
        nGeneTrees = geneTrees.size();
        perGenePloidy = new double[nGeneTrees];
        speciesNodeCount = -1;
        for (int geneI = 0; geneI < nGeneTrees; geneI++) {
            final Distribution pDist = geneTrees.get(geneI);
            if (pDist instanceof GeneTreeForSpeciesTreeDistribution) {
                final GeneTreeForSpeciesTreeDistribution gt = (GeneTreeForSpeciesTreeDistribution) pDist;
                perGenePloidy[geneI] = gt.getPloidy();
                if (speciesNodeCount == -1)
                    speciesNodeCount = gt.speciesTreeInput.get().getNodeCount();
            } else { // check that all input distributions are gene trees
                throw new IllegalArgumentException("Input distributions must all be of class GeneTree.");
            }
        }

        if (speciesNodeCount != -1) { // not BEAUTi
            allLineageCounts = new int[speciesNodeCount*nGeneTrees];
            allEventCounts = new int[speciesNodeCount*nGeneTrees];
            allCoalescentTimes = new double[speciesNodeCount*nGeneTrees][];
            perBranchLogP = new double[speciesNodeCount];
    
            storedLineageCounts = new int[speciesNodeCount*nGeneTrees];
            storedEventCounts = new int[speciesNodeCount*nGeneTrees];
            storedCoalescentTimes = new double[speciesNodeCount*nGeneTrees][];
            storedPerBranchLogP = new double[speciesNodeCount];
        }
    }

    private boolean checkHyperparameters(final boolean force) {
        invGammaShape = populationShapeInput.get();
        invGammaMean = populationMeanInput.get();
        final double currentAlpha = invGammaShape.getValue();
        final double currentBeta = invGammaMean.getValue() * (alpha - 1.0);

        if (force || currentAlpha != alpha || currentBeta != beta) {
            alpha = currentAlpha;
            beta = currentBeta;
            return true;
        }

        return false;
    }

    @Override
	public double calculateLogP() {
        super.calculateLogP();
        // System.out.println(tmpLogP + " -> " + logP);
        if (dontCalculate || Double.isInfinite(logP) || Double.isNaN(logP)) return logP;

        // need to recompute all branches if the parameters of the prior distribution have changed
        final boolean updatedPrior = checkHyperparameters(false);

        final int[] branchLineageCounts = new int[nGeneTrees];
        final int[] branchEventCounts = new int[nGeneTrees];
        final double[][] branchCoalescentTimes = new double[nGeneTrees][];

        final List<Distribution> pDists = pDistributions.get();
        final GeneTreeForSpeciesTreeDistribution[] geneTrees = new GeneTreeForSpeciesTreeDistribution[pDists.size()];

        int nodeGeneI = 0;
        for (int nodeI = 0; nodeI < speciesNodeCount; nodeI++) {
            boolean dirtyBranch = false;
            for (int geneI = 0; geneI < nGeneTrees; geneI++) {
                if (nodeI == 0) geneTrees[geneI] = (GeneTreeForSpeciesTreeDistribution) pDists.get(geneI);
                final GeneTreeForSpeciesTreeDistribution geneTree = geneTrees[geneI];

               // TODO if (geneTree.isDirtyBranch(nodeI)) {
                    dirtyBranch = true;

                    final double[] geneBranchCoalescentTimes = geneTree.getCoalescentTimes(nodeI);
                    final int geneBranchLineageCount = geneTree.coalescentLineageCounts[nodeI];
                    final int geneBranchEventCount = geneTree.coalescentCounts[nodeI];

                    allLineageCounts[nodeGeneI] = geneBranchLineageCount;
                    allEventCounts[nodeGeneI] = geneBranchEventCount;
                    allCoalescentTimes[nodeGeneI] = geneBranchCoalescentTimes;
               // }

                branchLineageCounts[geneI] = allLineageCounts[nodeGeneI];
                branchEventCounts[geneI] = allEventCounts[nodeGeneI];
                branchCoalescentTimes[geneI] = allCoalescentTimes[nodeGeneI];

                nodeGeneI++;
            }

            if (updatedPrior || dirtyBranch)
                perBranchLogP[nodeI] = analyticalLogP(alpha, beta, perGenePloidy, branchCoalescentTimes, branchLineageCounts, branchEventCounts);

            logP += perBranchLogP[nodeI];
        }

        return logP;
    }

    static private double analyticalLogP(double alpha, double beta, double[] perGenePloidy, double[][] branchCoalescentTimes, int[] branchLineageCounts, int[] branchEventCounts) {
        final int nGenes = perGenePloidy.length;

        int branchQ = 0;
        double branchLogR = 0.0;
        double branchGamma = 0.0;

        for (int j = 0; j < nGenes; j++) {
            final int geneN = branchLineageCounts[j];
            final double[] geneCoalescentTimes = branchCoalescentTimes[j];
            final int geneK = branchEventCounts[j];
            final double genePloidy = perGenePloidy[j]; 
            branchLogR -= geneK * Math.log(genePloidy);
            branchQ += geneK;

            double partialGamma = 0.0;
            for (int i = 0; i < geneK; i++) {
                partialGamma += (geneCoalescentTimes[i + 1] - geneCoalescentTimes[i]) * (geneN - i) * (geneN - (i + 1.0)) / 2.0;
            }
            
            if (geneN - geneK > 1) {
                partialGamma += (geneCoalescentTimes[geneK + 1] - geneCoalescentTimes[geneK]) * (geneN - geneK) * (geneN - (geneK + 1.0)) / 2.0;
            }

            branchGamma += partialGamma / genePloidy;
        }

        double logGammaRatio = 0.0;
        for (int i = 0; i < branchQ; i++) {
            logGammaRatio += Math.log(alpha + i);
        }

        final double logP = branchLogR + (alpha * Math.log(beta)) - ((alpha + branchQ) * Math.log(beta + branchGamma)) + logGammaRatio;

        return logP;
    }

    /*@Override
    public double getCurrentLogP() {
        return calculateLogP();
    }*/

    @Override
    public List<String> getArguments() {
        return null;
    }

    @Override
    public List<String> getConditions() {
        return null;
    }

    @Override
    public void sample(State state, Random random) {
    }
}
