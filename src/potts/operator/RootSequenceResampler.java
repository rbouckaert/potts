package potts.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.util.Randomizer;
import potts.datatype.AminoacidPlus;
import potts.dca.DCA;
import potts.dca.DCASequenceSimulator;
import potts.likelihood.PottsSequenceLikelihood;
import potts.likelihood.TreeLikelihoodWithRootStates;

@Description("Gibbs sampler for root sequence in TreeLikelihoodWithRootStates")
public class RootSequenceResampler extends Operator {

	final public Input<PottsSequenceLikelihood> pottsSeqLikelihoodInput = new Input<>("pottsSeqLikelihood", "DCA json file previously trained on an alignment",
			Validate.REQUIRED);

	final public Input<IntegerParameter> sequenceInput = new Input<>("sequence",
			"specifies sequence states to calculated likelihood for", Validate.REQUIRED);

	final public Input<TreeLikelihoodWithRootStates> likelihoodInput = new Input<>("likelihood", 
			"TreeLikelihood With Root States to provide root partials", Validate.REQUIRED);
	
	final public Input<Integer> stepCountInput = new Input<>("stepCount", 
			"average number of times a site is resampled", 20);

	DCA dca;
	IntegerParameter sequence;
	TreeLikelihoodWithRootStates likelihood;
	int totalStepCount, stateCount;
	Alignment data;
	double temperaturFactor;
	
	@Override
	public void initAndValidate() {
		PottsSequenceLikelihood psl = pottsSeqLikelihoodInput.get();
		dca = psl.dcaInput.get();
		temperaturFactor = 1.0 / psl.temperatureInput.get();
		
		sequence = sequenceInput.get();
		likelihood = likelihoodInput.get();
		data = likelihood.dataInput.get();
		totalStepCount = stepCountInput.get() * sequence.getDimension();

		boolean seperateGapState = data.getDataType() instanceof AminoacidPlus;
		stateCount = dca.getStateCount() - (seperateGapState ? 0 : 1);
	}

	@Override
	public double proposal() {
		Integer [] values = sequence.getValues();
		int [] seq = new int[values.length];
		for (int i = 0; i < seq.length; i++) {
			seq[i] = values[i];
		}
		
		double [] rootPartials = likelihood.getRootPartials();
		double [] probs = new double[stateCount];
		
		
		double [] freqs = likelihood.getSubstitutionModel().getFrequencies();

totalStepCount = values.length;

		double logHR = 0;
		for (int i = 0; i < totalStepCount; i++) {
			int site = Randomizer.nextInt(values.length);
			site = i;
			int oldState = seq[site];
			for (int newState = 0; newState < stateCount; newState++) {
				probs[newState] = computeDeltaHamiltonian(dca, seq, site, oldState, newState, rootPartials, freqs);
			}
			
			// find max
			double max = probs[0];
			int iMax = 0;
			for (int k = 0; k < stateCount; k++) {
				if (probs[k] > max) {
					iMax = k;
					max = Math.max(probs[k], max);
				}
			}
			
			// to real space
			for (int k = 0; k < stateCount; k++) {
				probs[k] = Math.exp(probs[k] - max);
			}
			double sum = 0;
			for (double d : probs) {
				sum += d;
			}
			try {
				int newState = (sum > 0) ? Randomizer.randomChoicePDF(probs) : Randomizer.nextInt(freqs.length);
				//int newState = Randomizer.randomChoicePDF(freqs);
				//int newState = Randomizer.nextInt(freqs.length);
				//int newState = iMax;
				seq[site] = newState;
				logHR += Math.log(probs[oldState]) - Math.log(probs[newState]);
			} catch (Error e) {
				int newState = Randomizer.nextInt(freqs.length);
				seq[site] = newState;
			}
		}
		
		
		for (int i = 0; i < values.length; i++) {
			sequence.setValue(i, seq[i]);
		}
		
		return logHR;
	}

	
    private double computeDeltaHamiltonian(DCA dca, int[] seq, int i, int oldState, int newState, 
    		double [] rootPartials, double [] freqs) {
        double delta =  (temperaturFactor > 0) ? 
        		temperaturFactor * DCASequenceSimulator.computeDeltaHamiltonian(dca, seq, i, oldState, newState)
        		: 0;
        
        int patternIndexOffset = data.getPatternIndex(i) * stateCount;

        // Change in transition probability
        delta += Math.log(rootPartials[patternIndexOffset + newState]);// * freqs[newState]);// - Math.log(rootPartials[patternIndexOffset + oldState]);
        		

        return delta;
    }    

}
