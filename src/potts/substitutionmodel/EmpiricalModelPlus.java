package potts.substitutionmodel;

import java.lang.reflect.InvocationTargetException;

import beast.base.core.Input;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.EmpiricalSubstitutionModel;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.substitutionmodel.WAG;
import beast.base.evolution.tree.Node;
import potts.datatype.AminoacidPlus;

public class EmpiricalModelPlus extends GeneralSubstitutionModel {

	final public Input<EmpiricalSubstitutionModel> substModelInput = new Input<>("substmodel", "empirical subsitution model for aminoacids", new WAG());
	
	private EmpiricalSubstitutionModel model;
	
	/** rate mutating into or out of a GAP **/
	final static double GAP_RATE = 1.0;
	
	@Override
	public void initAndValidate() {
		model = substModelInput.get();
		nrOfStates = 21;
        rateMatrix = new double[nrOfStates][nrOfStates];
        relativeRates = new double[21 * 20];
        storedRelativeRates = new double[21 * 20];
		
	    try {
			eigenSystem = createEigenSystem();
		} catch (SecurityException | ClassNotFoundException | InstantiationException | IllegalAccessException | IllegalArgumentException
				| InvocationTargetException e) {
			throw new IllegalArgumentException(e.getMessage());
		}
	}
	
    @Override
    public void setupRelativeRates() {
        double[] empiricalRates = model.getEmpericalRateValues();
        int states = 21;
        for (int i = 0; i < states-1; i++) {
        	System.arraycopy(empiricalRates, i * (nrOfStates-2), relativeRates, i * (nrOfStates-1), (nrOfStates-2));
        	empiricalRates[i * (nrOfStates-1) + nrOfStates - 1] = GAP_RATE;
        }
    }

    @Override
    public void setupRateMatrix() {
        double[] freqs = getFrequencies();
        for (int i = 0; i < nrOfStates; i++) {
            rateMatrix[i][i] = 0;
            for (int j = 0; j < i; j++) {
                rateMatrix[i][j] = relativeRates[i * (nrOfStates - 1) + j];
            }
            for (int j = i + 1; j < nrOfStates; j++) {
                rateMatrix[i][j] = relativeRates[i * (nrOfStates - 1) + j - 1];
            }
        }
        // bring in frequencies
        for (int i = 0; i < nrOfStates; i++) {
            for (int j = i + 1; j < nrOfStates; j++) {
                rateMatrix[i][j] *= freqs[j];
                rateMatrix[j][i] *= freqs[i];
            }
        }
        // set up diagonal
        for (int i = 0; i < nrOfStates; i++) {
            double sum = 0.0;
            for (int j = 0; j < nrOfStates; j++) {
                if (i != j)
                    sum += rateMatrix[i][j];
            }
            rateMatrix[i][i] = -sum;
        }
        // normalise rate matrix to one expected substitution per unit time
        double subst = 0.0;
        for (int i = 0; i < nrOfStates; i++)
            subst += -rateMatrix[i][i] * freqs[i];

        for (int i = 0; i < nrOfStates; i++) {
            for (int j = 0; j < nrOfStates; j++) {
                rateMatrix[i][j] = rateMatrix[i][j] / subst;
            }
        }
    }

    @Override
    public double[] getRateMatrix(Node node) {
        double[][] matrix = model.getEmpiricalRates();
        int states = 21;
        double[] rates = new double[states * states];
        for (int i = 0; i < states-1; i++) {
            for (int j = i + 1; j < states-1; j++) {
                rates[i * states + j] = matrix[i][j];
                rates[j * states + i] = matrix[i][j];
            }
        }
        for (int j = 0; j < states-1; j++) {
            rates[20 * states + j] = GAP_RATE;
            rates[j * states + 20] = GAP_RATE;
        }
        // determine diagonal
        for (int i = 0; i < states; i++) {
            double sum = 0;
            for (int j = i + 1; j < states; j++) {
                sum += rates[i * states + j];
            }
            rates[i * states + i] = -sum;
        }
        return rates;
    }
	
	@Override
	public EigenDecomposition getEigenDecomposition(Node node) {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public boolean canHandleDataType(DataType dataType) {
		return dataType instanceof AminoacidPlus;
	}

}
