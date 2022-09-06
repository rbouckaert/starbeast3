package starbeast3.app.beauti;

import java.io.File;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import beast.base.core.BEASTInterface;
import beast.base.core.Log;
import beast.base.inference.StateNode;
import beast.base.inference.operator.kernel.BactrianDeltaExchangeOperator;
import beastfx.app.inputeditor.BeautiAlignmentProvider;
import beastfx.app.inputeditor.BeautiDoc;
import starbeast3.evolution.branchratemodel.StarBeast3Clock;

public class StarBeastAlignmentProvider3 extends BeautiAlignmentProvider {

	
	
	
	@Override
	public List<BEASTInterface> getAlignments(BeautiDoc doc, File[] files) {
		if (doc.getPartitions("Partitions").size() == 0) {
			// initial partition to be added
			// so set up default mode
			doc.autoSetClockRate = false;
			doc.beauti.get_autoSetClockRate().setSelected(false);
			doc.autoUpdateFixMeanSubstRate = false;
			initialValueOfAutoSet = doc.beauti.get_autoSetClockRate().isSelected();
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
	
	
	
    /**
     * Set clock rate estimate=true initially
     */
    private static StarBeast3Clock firstBeautiClock = null;
    private static Map<StarBeast3Clock, Integer> beautiClocksChecked = new HashMap<>();
    private static Boolean initialValueOfAutoSet = null;

    public static void estimateGeneClockRates(BeautiDoc doc) {
    	if (initialValueOfAutoSet == null) {
        	// Loaded in an old file. get autoSetClockRate value
    		initialValueOfAutoSet = doc.beauti.get_autoSetClockRate().isSelected();
    	}
    	
    	// Loaded in an old file. Do not overwrite
    	if (!initialValueOfAutoSet && !doc.beauti.get_autoSetClockRate().isSelected() && firstBeautiClock == null) {
    		return;
    	}
    	
    	// Reset
    	if (doc.beauti.get_autoSetClockRate().isSelected()) {
    		firstBeautiClock = null;
    		initialValueOfAutoSet = true;
    		beautiClocksChecked.clear();
    		return;
    	}

    	for (String str : doc.pluginmap.keySet()) {
    		
    		// Find the StarBeast3Clocks
    		BEASTInterface obj = doc.pluginmap.get(str);
    		if (obj instanceof StarBeast3Clock) {
    			
    			StarBeast3Clock clock = (StarBeast3Clock)obj;
    			
    			/*
    			// Do not revisit this clock or the users wishes will be overridden
				if (beautiClocksChecked.containsKey(clock) && beautiClocksChecked.get(clock) > 0) {
					Log.warning("StarBeast3Clock Skipping " + clock.getID());
					continue;
				}
    			*/
    			
    			if (clock.meanRateInput.get() != null && clock.meanRateInput.get() instanceof StateNode) {
					
    				StateNode rp = (StateNode)clock.meanRateInput.get();
    				rp.isEstimatedInput.set(true);
    				// The first clock is not estimated
        			if (firstBeautiClock == null || firstBeautiClock == clock) {
        				firstBeautiClock = clock;
        				rp.isEstimatedInput.set(false);
        			}else {
    					rp.isEstimatedInput.set(true);
    					Log.warning("StarBeast3Clock estimating " + rp.getID());
        			}
        			
        			int val = 1;
        			if (beautiClocksChecked.containsKey(clock)) {
        				val += beautiClocksChecked.get(clock);
        			}
        			beautiClocksChecked.put(clock, val);
				}
    			
    		}
    		
    	}
    
    	
    }

}
