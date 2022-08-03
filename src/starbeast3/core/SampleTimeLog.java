package starbeast3.core;

import java.io.PrintStream;
import java.time.LocalDateTime;
import java.time.temporal.ChronoUnit;

import beast.base.inference.CalculationNode;
import beast.base.core.Function;
import beast.base.core.Input;
import beast.base.core.Loggable;
import beast.base.inference.Operator;
import beast.base.inference.OperatorSchedule;
import beast.base.core.Log;
import starbeast3.operators.MultiStepOperator;

public class SampleTimeLog extends CalculationNode implements Loggable, Function {
	
	
	
	final public Input<OperatorSchedule> operatorScheduleInput = new Input<>("schedule", "The oeprator schedule, used for determining number of proposals", Input.Validate.REQUIRED);

  	protected LocalDateTime dateTimeLast_l;
  	protected LocalDateTime dateTimeFirst_l;
  	protected double nproposalsLast;
  
  
	@Override
	public void initAndValidate() {
		LocalDateTime now = LocalDateTime.now();
		this.dateTimeLast_l = now;
		this.dateTimeFirst_l = now;
		nproposalsLast = 0;
	}
		
	
	
	@Override
	public void init(PrintStream out) {
		out.print("MsmplPerMin_i" + "\t" + "MsmplPerMin_c" + "\t");
	}
	
	@Override
	public void close(PrintStream out) {
		
	}

	@Override
	public void log(long sample, PrintStream out) {
		
		
		// Total runtime
		LocalDateTime dateTimeLast = this.dateTimeLast_l;
		LocalDateTime dateTimeFirst = this.dateTimeFirst_l;

		// Current time
		LocalDateTime dateTimeNow = LocalDateTime.now();
		double i_time = ChronoUnit.MILLIS.between(dateTimeLast, dateTimeNow);
		i_time = i_time / 60000.0;
		double c_time = ChronoUnit.MILLIS.between(dateTimeFirst, dateTimeNow);
		c_time = c_time / 60000.0;
				
		
		
		
		// Total n proposals
		double nproposals = 0;
		for (Operator operator : operatorScheduleInput.get().operators ) {
			
			
			
			if (operator instanceof MultiStepOperator) {
				MultiStepOperator multi = (MultiStepOperator) operator;
				nproposals += multi.getNrProposals();
				//Log.warning(operator.getID() + " " + multi.getNrProposals());
			}else {
				nproposals += operator.get_m_nNrAccepted() + operator.get_m_nNrRejected();
				//Log.warning(operator.getID() + " " +(operator.get_m_nNrAccepted() + operator.get_m_nNrRejected()));
			}
		}
		nproposals = nproposals / 1e6;
		
		// Sampling rate
		double smplPrHr_i = 1.0*(nproposals-nproposalsLast) / i_time;
		double smplPrHr_c = 1.0*nproposals / c_time;
		if (sample == 0) {
			smplPrHr_i=0;
			smplPrHr_c=0;
		}
		
		
		//Log.warning(nproposals + " / " + time);
		
		// Log
		out.print(smplPrHr_i + "\t" + smplPrHr_c + "\t");
		
		// Update time
		this.dateTimeLast_l = dateTimeNow;
		nproposalsLast = nproposals;
		
		
		
		
	}


	@Override
	public int getDimension() {
		return 2;
	}


	@Override
	public double getArrayValue(int dim) {
		// TODO Auto-generated method stub
		return 0;
	}

}
