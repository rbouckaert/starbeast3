package starbeast3.app.beauti;

import java.io.File;
import java.util.List;

import beast.base.core.BEASTInterface;
import beast.base.inference.operator.kernel.BactrianDeltaExchangeOperator;
import beastfx.app.inputeditor.BeautiAlignmentProvider;
import beastfx.app.inputeditor.BeautiDoc;

public class StarBeastAlignmentProvider3 extends BeautiAlignmentProvider {

	
	
	
	@Override
	public List<BEASTInterface> getAlignments(BeautiDoc doc, File[] files) {
		if (doc.getPartitions("Partitions").size() == 0) {
			// initial partition to be added
			// so set up default mode
			doc.autoSetClockRate = false;
			doc.beauti.get_autoSetClockRate().setSelected(false);
			doc.autoUpdateFixMeanSubstRate = false;
			doc.beauti.get_autoUpdateFixMeanSubstRate().setSelected(false);
			// initialValueOfAutoSet = doc.beauti.get_autoSetClockRate().isSelected();
		}

		final List<BEASTInterface> newAlignments = super.getAlignments(doc, files);
		final int alignmentCount = newAlignments.size();

		System.out.println(String.format("N_ALIGNMENTS = %d", doc.alignments.size()));
		
		
		// initialize delta exchange operator in order to increase weight to something more sensible
		BactrianDeltaExchangeOperator operator = (BactrianDeltaExchangeOperator) doc.pluginmap.get("FixMeanMutationRatesOperator");
		if (operator == null) {
			operator = new BactrianDeltaExchangeOperator();
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
