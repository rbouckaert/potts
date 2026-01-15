package potts.substitutionmodel;

import java.lang.reflect.InvocationTargetException;
import java.util.Arrays;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.substitutionmodel.EigenDecomposition;
import beast.base.evolution.substitutionmodel.EmpiricalSubstitutionModel;
import beast.base.evolution.substitutionmodel.GeneralSubstitutionModel;
import beast.base.evolution.substitutionmodel.WAG;
import beast.base.evolution.tree.Node;
import potts.datatype.AminoacidPlus;


@Description("Model for explicitly representing a gap state. "
		+ "A constant rate in and out of the gap state is assumed. "
		+ "For the other states, an empirical substitution model for amino acids is assumed.")
public class EmpiricalModelPlus extends GeneralSubstitutionModel {

	final public Input<EmpiricalSubstitutionModel> substModelInput = new Input<>("substModel", "empirical subsitution model for aminoacids", new WAG());
	
	private EmpiricalSubstitutionModel model;
	
	/** rate mutating into or out of a GAP **/
	final static double GAP_RATE = 1.0;
	double GAP_PROPORTION = 1.0/21.0;
	
	
	public EmpiricalModelPlus() {
		ratesInput.setRule(Validate.OPTIONAL);
		frequenciesInput.setRule(Validate.OPTIONAL);
	}
	
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
        if (empiricalRates.length != (nrOfStates-2)*(nrOfStates-1) || relativeRates.length != (nrOfStates-1)*nrOfStates) {
        	int h = 3;
        	h++;
        }
        Arrays.fill(relativeRates, 0);
        for (int i = 0; i < nrOfStates-1; i++) {
        	System.arraycopy(empiricalRates, i * (nrOfStates-2), relativeRates, i * (nrOfStates-1), nrOfStates-2);
        	relativeRates[i * (nrOfStates-1) + nrOfStates - 2] = GAP_RATE;
        }
        
        int offset = (nrOfStates-1)*(nrOfStates-1);
        for (int i = 0; i < nrOfStates-1; i++) {
        	relativeRates[offset++] = GAP_RATE;
        }
        
    }

    
    @Override
    public double[] getFrequencies() {
    	double [] empiricalFreqs = model.getEmpiricalFrequencies();
    	double sum = 1.0 + GAP_PROPORTION;
    	double [] freqs = new double[nrOfStates];
    	for (int i = 0; i < nrOfStates-1; i++) {
    		freqs[i] = empiricalFreqs[i] / sum;
    	}
    	freqs[nrOfStates-1] = GAP_PROPORTION / sum;
    	return freqs;
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
	
//	@Override
//	public EigenDecomposition getEigenDecomposition(Node node) {
//		return eigenDecomposition;
//	}

	@Override
	public boolean canHandleDataType(DataType dataType) {
		return dataType instanceof AminoacidPlus;
	}

}
