package starbeast3.app.beauti;

import beast.base.core.BEASTInterface;
import beast.base.core.Input;
import beastfx.app.inputeditor.ListInputEditor;
import starbeast3.evolution.branchratemodel.StarBeast3Clock;

public class StarBeast3ClockInputEditor extends ListInputEditor {

	@Override
	public Class<?> baseType() {
		return StarBeast3Clock.class;
	}
	
	@Override
	public void init(Input<?> input, BEASTInterface beastObject, int itemNr, ExpandOption isExpandOption,
			boolean addButtons) {
		// using default list input editor suppresses addition of drop down box
		super.init(input, beastObject, itemNr, isExpandOption, addButtons);
	}

}
