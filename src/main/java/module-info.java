open module starbeast3 {
	requires beast.pkgmgmt;
    requires beast.base;
    requires static javafx.controls;
    requires static beast.fx;
    requires static orc;
    requires static org.apache.commons.statistics.distribution;
    requires static commons.math3;
    

    exports starbeast3.app.beauti;
    exports starbeast3.core;
    exports starbeast3.evolution.branchratemodel;
    exports starbeast3.evolution.speciation;
    exports starbeast3.evolution.substitutionmodel;
    exports starbeast3.math.distributions;
    exports starbeast3.operators;
    exports starbeast3.simulation;
    exports starbeast3.tree;


    provides beast.base.core.BEASTInterface with
        
        starbeast3.app.beauti.StarBeastAlignmentProvider3,
        starbeast3.app.beauti.StarBeast3MorphModelAlignmentProvider,
        starbeast3.simulation.DirectSimulator,
        starbeast3.core.OperatorScheduleRecalculator,
        starbeast3.core.ParallelMCMC,
        starbeast3.evolution.branchratemodel.SharedSpeciesClockModel,
        starbeast3.evolution.branchratemodel.StrictClockModelSB3,
        starbeast3.evolution.branchratemodel.UCRelaxedClockModelSB3,
        starbeast3.evolution.speciation.ConstantPopulations,
        starbeast3.evolution.substitutionmodel.LewisMK,
        starbeast3.tree.StarBeast3TaxonSet,
        starbeast3.evolution.speciation.GeneTreeForSpeciesTreeDistribution,
        starbeast3.math.distributions.MRCAPriorSB3,
        starbeast3.operators.AVMNCubeOperator,
        starbeast3.operators.ConstantDistanceOperatorSpeciesTree,
        starbeast3.operators.CoordinatedExchangeRates,
        starbeast3.operators.CoordinatedExponential,
        starbeast3.operators.CoordinatedUniform,
        starbeast3.operators.SampledNodeDateRandomWalkerSB3,
        starbeast3.operators.NEROperator_dAE_dBE_dCE,
        starbeast3.operators.NodeReheight2,
        starbeast3.operators.ParallelDistSet,
        starbeast3.operators.ParallelMCMCRealParameterOperator,
        starbeast3.operators.ParallelMCMCTreeOperator,
        starbeast3.operators.ParallelMCMCTreeOperatorTreeDistribution,
        starbeast3.operators.PopSizeGibbsSampler,
        starbeast3.tree.SpeciesTree,
        starbeast3.core.SpeciesTreeLogger,
        starbeast3.evolution.speciation.SpeciesTreePrior,
        starbeast3.evolution.branchratemodel.StarBeast3Clock,
        starbeast3.core.StarBeastStartState;
    
    
    provides beastfx.app.inputeditor.InputEditor with
    	starbeast3.app.beauti.MRCAPriorInputEditorSB3,
    	starbeast3.app.beauti.StarBeast3ClockInputEditor,
    	starbeast3.app.beauti.StarBeastTipDatesInputEditor;
 
    provides beastfx.app.inputeditor.BeautiAlignmentProvider with
    	starbeast3.app.beauti.StarBeastAlignmentProvider3,
    	starbeast3.app.beauti.StarBeast3MorphModelAlignmentProvider;
}
