package potts.operator;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.inference.Operator;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.util.Randomizer;
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
		stateCount = dca.getStateCount() - 1;
		
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

//totalStepCount = 1;

		for (int i = 0; i < totalStepCount; i++) {
			int site = Randomizer.nextInt(values.length);
			int oldState = seq[site];
			for (int newState = 0; newState < stateCount; newState++) {
				probs[newState] = computeDeltaHamiltonian(dca, seq, site, oldState, newState, rootPartials, freqs);
			}
			
			// find max
			double max = probs[0];
			for (double d : probs) {
				max = Math.max(d, max);
			}
			
			// to real space
			for (int k = 0; k < stateCount; k++) {
				probs[k] = Math.exp(probs[k] - max);
			}
			
			int newState = Randomizer.randomChoicePDF(probs);
			//int newState = Randomizer.randomChoicePDF(freqs);
			//int newState = Randomizer.nextInt(freqs.length);
            seq[site] = newState;
		}
		
		
		for (int i = 0; i < values.length; i++) {
			sequence.setValue(i, seq[i]);
		}
		
		return 0;
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
