package starbeast3.evolution.substitutionmodel;


import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.datatype.Binary;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.StandardData;
import beast.base.evolution.tree.Node;
import beast.base.evolution.substitutionmodel.*;

import java.util.Arrays;

/**
 * @author Alexandra Gavryushkina
 * @author Joseph Heled
 */
@Description("LewisMK subtitution model: equal rates, equal frequencies")
@Citation(value = "Lewis, P. O. (2001). \n" +
        "A likelihood approach to estimating phylogeny from discrete morphological character data. \n" +
        "Systematic biology, 50(6), 913-925.",
        year = 2001, firstAuthorSurname = "Lewis", DOI="10.1080/106351501753462876")
public class LewisMK extends SubstitutionModel.Base  {

    public Input<Integer> nrOfStatesInput = new Input<Integer>("stateNumber", "the number of character states");

    public Input<DataType> dataTypeInput = new Input<DataType>("datatype", "datatype, used to determine the number of states", Validate.XOR, nrOfStatesInput);

    //public Input<Boolean> isMKvInput = new Input<>("isMkv", "whether to use MKv model or just MK", false);
    boolean hasFreqs;
    private boolean updateFreqs;

    public LewisMK() {
        // this is added to avoid a parsing error inherited from superclass because frequencies are not provided.
        frequenciesInput.setRule(Validate.OPTIONAL);
        try {
            // this call will be made twice when constructed from XML
            // but this ensures that the object is validly constructed for testing purposes.
        	// RRB: for testing you should call initAndValidate() independently
        	// At this stage, the inputs are not valid, hence initAndValidate fails.
            // initAndValidate();
        } catch (Exception e) {
            e.printStackTrace();
            throw new RuntimeException("initAndValidate() call failed when constructing LewisMK()");
        }
    }

    double totalSubRate;
    double[] frequencies;
    EigenDecomposition eigenDecomposition;

    private void setFrequencies() {
        final Frequencies frequencies1 = frequenciesInput.get();
        if( frequencies1 != null ) {
            if( frequencies1.getFreqs().length != nrOfStates ) {
                throw new RuntimeException("number of stationary frequencies does not match number of states.");
            }
            System.arraycopy(frequencies1.getFreqs(), 0, frequencies, 0, nrOfStates);
            totalSubRate = 1;
            for(int k = 0; k < nrOfStates; ++k) {
               totalSubRate -= frequencies[k]*frequencies[k];
            }
            hasFreqs = true;
        }  else {
            Arrays.fill(frequencies, (1.0 / nrOfStates));
            hasFreqs = false;
        }
        updateFreqs = false;
    }

    @Override
    public void initAndValidate() {
    	if (nrOfStatesInput.get() != null) {
    		nrOfStates = nrOfStatesInput.get();
    	} else {
    		nrOfStates = dataTypeInput.get().getStateCount();
    	}
    	
    	if (nrOfStates <= 1) {    		
    		// this goes here so that nrOfStates/(nrOfStates-1) in getTransitionProbabilities()    		
    		// does not throw a division by zero exception
    		throw new IllegalArgumentException("The number of states should be at least 2 but is" + nrOfStates + ".\n"
    				+ "This may be due to a site in the alignment having only 1 state, which can be fixed by "
    				+ "removing the site from the alignment.");
    	}

//        double[] eval = new double[]{0.0, -1.3333333333333333, -1.3333333333333333, -1.3333333333333333};
//        double[] evec = new double[]{1.0, 2.0, 0.0, 0.5, 1.0, -2.0, 0.5, 0.0, 1.0, 2.0, 0.0, -0.5, 1.0, -2.0, -0.5, 0.0};
//        double[] ivec = new double[]{0.25, 0.25, 0.25, 0.25, 0.125, -0.125, 0.125, -0.125, 0.0, 1.0, 0.0, -1.0, 1.0, 0.0, -1.0, 0.0};
//
//        eigenDecomposition = new EigenDecomposition(evec, ivec, eval);
        frequencies = new double[nrOfStates];
        setFrequencies();
    }

    @Override
    public double[] getFrequencies() {
        return frequencies;
    }

    @Override
    public void getTransitionProbabilities(Node node, double fStartTime, double fEndTime, double fRate, double[] matrix) {
        if( updateFreqs ) {
           setFrequencies();
        }

        if( hasFreqs ) {
            final double e1 = Math.exp(-(fStartTime - fEndTime) * fRate/totalSubRate);
            final double e2 = 1 - e1;

            for( int i = 0; i < nrOfStates; ++i ) {
                final int r = i * nrOfStates;
                for( int j = 0; j < nrOfStates; ++j ) {
                    matrix[r + j] = frequencies[j] * e2;
                }
                matrix[r + i] += e1;
            }
        } else {
            double fDelta = (nrOfStates / (nrOfStates - 1.0)) * (fStartTime - fEndTime);
            double fPStay = (1.0 + (nrOfStates - 1) * Math.exp(-fDelta * fRate)) / nrOfStates;
            double fPMove = (1.0 - Math.exp(-fDelta * fRate)) / nrOfStates;
            // fill the matrix with move probabilities
            Arrays.fill(matrix, fPMove);
            // fill the diagonal
            for (int i = 0; i < nrOfStates; i++) {
                matrix[i * (nrOfStates + 1)] = fPStay;
            }
        }
    }

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {
        return eigenDecomposition;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        if (dataType instanceof StandardData || dataType instanceof Binary) {
            return true;
        }
        return false;
        //throw new Exception("Can only handle StandardData and binary data");
    }

    protected boolean requiresRecalculation() {
        if( ! hasFreqs ) {
            return false;
        }
        // we only get here if something is dirty
        updateFreqs = true;
        return true;
    }

    @Override
    protected void restore() {
        updateFreqs = true;
        super.restore();
    }
}
