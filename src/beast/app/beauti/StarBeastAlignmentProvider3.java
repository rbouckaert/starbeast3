package beast.app.beauti;

import java.io.File;
import java.util.List;

import beast.core.BEASTInterface;
import beast.evolution.operators.DeltaExchangeOperator;

public class StarBeastAlignmentProvider3 extends BeautiAlignmentProvider {

	@Override
	public List<BEASTInterface> getAlignments(BeautiDoc doc, File[] files) {
		final List<BEASTInterface> newAlignments = super.getAlignments(doc, files);
		final int alignmentCount = newAlignments.size();

		doc.autoSetClockRate = true;
		doc.beauti.autoSetClockRate.setSelected(true);
		doc.autoUpdateFixMeanSubstRate = false;

		System.out.println(String.format("N_ALIGNMENTS = %d", doc.alignments.size()));
		
		
		
		
		// initialize delta exchange operator in order to increase weight to something more sensible
		DeltaExchangeOperator operator = (DeltaExchangeOperator) doc.pluginmap.get("FixMeanMutationRatesOperator");
		if (operator == null) {
			operator = new DeltaExchangeOperator();
			try {
				operator.setID("FixMeanMutationRatesOperator");
				operator.initByName("weight", (double) alignmentCount, "delta", 0.75);
			} catch (Throwable e1) {
				// ignore initAndValidate exception
			}
			doc.addPlugin(operator);
		} else {
			final double updatedWeight = doc.alignments.size() + alignmentCount;
			operator.setInputValue("weight", updatedWeight);
		}
		
		
		

		return newAlignments;
	}
}
