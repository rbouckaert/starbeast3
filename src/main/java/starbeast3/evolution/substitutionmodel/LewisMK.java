package starbeast3.evolution.substitutionmodel;


import beast.base.core.Citation;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.datatype.Binary;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.StandardData;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.tree.Node;
import beast.base.spec.evolution.substitutionmodel.Base;
import beast.base.spec.evolution.substitutionmodel.Frequencies;

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
public class LewisMK extends Base {

    public Input<Integer> nrOfStatesInput = new Input<Integer>("stateNumber", "the number of character states");

    public Input<DataType> dataTypeInput = new Input<DataType>("datatype", "datatype, used to determine the number of states", Validate.XOR, nrOfStatesInput);

    boolean hasFreqs;
    private boolean updateFreqs;

    public LewisMK() {
        // this is added to avoid a parsing error inherited from superclass because frequencies are not provided.
        frequenciesInput.setRule(Validate.OPTIONAL);
    }

    double totalSubRate;
    double[] freqValues;
    EigenDecomposition eigenDecomposition;

    private void setFrequencies() {
        final Frequencies frequencies1 = frequenciesInput.get();
        if( frequencies1 != null ) {
            if( frequencies1.getFreqs().length != nrOfStates ) {
                throw new RuntimeException("number of stationary frequencies does not match number of states.");
            }
            System.arraycopy(frequencies1.getFreqs(), 0, freqValues, 0, nrOfStates);
            totalSubRate = 1;
            for(int k = 0; k < nrOfStates; ++k) {
               totalSubRate -= freqValues[k]*freqValues[k];
            }
            hasFreqs = true;
        }  else {
            Arrays.fill(freqValues, (1.0 / nrOfStates));
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
    		throw new IllegalArgumentException("The number of states should be at least 2 but is" + nrOfStates + ".\n"
    				+ "This may be due to a site in the alignment having only 1 state, which can be fixed by "
    				+ "removing the site from the alignment.");
    	}

        freqValues = new double[nrOfStates];
        setFrequencies();
    }

    @Override
    public double[] getFrequencies() {
        return freqValues;
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
                    matrix[r + j] = freqValues[j] * e2;
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
    public double[] getRateMatrix(Node node) {
        if (updateFreqs) {
            setFrequencies();
        }

        double[] matrix = new double[nrOfStates * nrOfStates];
        if (hasFreqs) {
            for (int i = 0; i < nrOfStates; i++) {
                double rowSum = 0;
                for (int j = 0; j < nrOfStates; j++) {
                    if (i != j) {
                        matrix[i * nrOfStates + j] = freqValues[j] / totalSubRate;
                        rowSum += matrix[i * nrOfStates + j];
                    }
                }
                matrix[i * nrOfStates + i] = -rowSum;
            }
        } else {
            double offDiag = 1.0 / (nrOfStates - 1);
            for (int i = 0; i < nrOfStates; i++) {
                for (int j = 0; j < nrOfStates; j++) {
                    if (i == j) {
                        matrix[i * nrOfStates + j] = -1.0;
                    } else {
                        matrix[i * nrOfStates + j] = offDiag;
                    }
                }
            }
        }
        return matrix;
    }

    @Override
    public EigenDecomposition getEigenDecomposition(Node node) {
        return eigenDecomposition;
    }

    @Override
    public boolean canHandleDataType(DataType dataType) {
        return dataType instanceof StandardData || dataType instanceof Binary;
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

