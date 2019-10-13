package starbeast3.operators;

import java.util.*;

import beast.core.*;

@Description("Operator schedule that takes in account MultiStepOperators so logging "
		+ "is done at the correct MCMC step")
public class MultiStepOperatorSchedule extends OperatorSchedule {
	final public Input<List<Logger>> loggerInput = new Input<>("logger", "loggers used to infer log-frequencies. Defaults to 1000 if no loggers are specified", new ArrayList<>()); 
	
	
	// TODO: automatically weighting of species tree operators vs others?
	
	int [] logEvery;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
		// collect log frequencies (removes duplicates and those that can divide others)
		List<Integer>  every = new ArrayList<>();
		for (Logger logger : loggerInput.get()) {
			int freq = logger.everyInput.get();
			every.add(freq);
		}
		if (every.size() == 0) {
			every.add(1000);
		}
		for (int i = every.size() - 1; i >= 0; i--) {
			boolean isDivider = false;
			for (int j = 0; j < every.size(); j++) {
				if (j != i && (every.get(j) % every.get(i) == 0)) {
					isDivider = true;
					break;
				}
			}
			if (isDivider) {
				every.remove(i);
			}
		}
		
		// convert List to array		
		logEvery = new int[every.size()];
		int k = 0;
		for (Integer i : every) {
			logEvery[k] = i;
		}
	}
	
    public Operator selectOperator(long sampleNr) {
    	Operator operator = super.selectOperator();
    	
    	if (operator instanceof MultiStepOperator) {
    		int steps = ((MultiStepOperator) operator).stepCount();
    		boolean tooManySteps = false;
    		for (int i : logEvery) {
    			if (sampleNr % i + steps >= i) {
    				tooManySteps = true;
    				break;
    			}
    		}
    		if (tooManySteps) {
    			return selectOperator(sampleNr);
    		}
    	}
    	
    	
        return operator;
    }
    
    @Override
    public Operator selectOperator() {
    	throw new RuntimeException("MultiStepOperatorSchedule should be used with " + beast.core.MCMCsb3.class.getName() + " not with beast.core.MCMC");
    }
}
