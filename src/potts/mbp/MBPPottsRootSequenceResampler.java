package potts.mbp;


import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.util.Randomizer;
import potts.dca.DCA;
import potts.dca.DCASequenceSimulator;
import potts.likelihood.PottsSequenceLikelihood;

@Description("Sampler for root sequence in TreeLikelihoodWithRootStates with middle base paired prior")
public class MBPPottsRootSequenceResampler extends MBPRootSequenceResampler {

	final public Input<PottsSequenceLikelihood> pottsLikelihood1Input = new Input<>("PottsSequenceLikelihood1", 
			"PottsSequenceLikelihood over Root States to provide root partials 1", Validate.REQUIRED);
	final public Input<PottsSequenceLikelihood> pottsLikelihood2Input = new Input<>("PottsSequenceLikelihood2", 
			"PottsSequenceLikelihood over Root States to provide root partials 2", Validate.REQUIRED);
	
	
	final public Input<Integer> stepCountInput = new Input<>("stepCount", 
			"average number of times a site is resampled", 1);

	PottsSequenceLikelihood pottsLikelihood1, pottsLikelihood2;
	DCA dca1, dca2;
	int totalStepCount;
	double temperaturFactor;
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();

		pottsLikelihood1 = pottsLikelihood1Input.get();
		pottsLikelihood2 = pottsLikelihood2Input.get();
		
		dca1 = pottsLikelihood1.dcaInput.get();
		dca2 = pottsLikelihood2.dcaInput.get();
		temperaturFactor = 1.0 / pottsLikelihood1.temperatureInput.get();
		totalStepCount = stepCountInput.get() * sequence1.getDimension();
	}

	
	private boolean first = true;
	
	@Override
	public double proposal() {
		if (first) {
			// skip first proposal in order to guarantee the 
			// first logged rootsequences are the true value
			// in a WCSS
			first = false;
			return Double.POSITIVE_INFINITY;
		}

		Integer [] values = sequence1.getValues();
		int n = values.length;
		int [] seq1 = new int[values.length];
		for (int i = 0; i < seq1.length; i++) {
			seq1[i] = values[i];
		}
		
		values = sequence2.getValues();
		if (values.length != n) {
			throw new IllegalArgumentException("sequence1 and sequence2 have different lengths");
		}
		int [] seq2 = new int[values.length];
		for (int i = 0; i < seq2.length; i++) {
			seq2[i] = values[i];
		}

		
		double [] rootPartials1 = likelihood1.getRootPartials();
		double [] rootPartials2 = likelihood2.getRootPartials();
		double [] probs = new double[stateCount * stateCount];
		
		
		double [] freqs1 = likelihood1.getSubstitutionModel().getFrequencies();
		double [] freqs2 = likelihood2.getSubstitutionModel().getFrequencies();

		//totalStepCount = values.length;

		double logHR = 0;
		for (int i = 0; i < totalStepCount; i++) {
			int site = Randomizer.nextInt(n);
			site = i;
			int oldState1 = seq1[site];
			int oldState2 = seq2[site];
			for (int newState1 = 0; newState1 < stateCount; newState1++) {
				double p1 = computeDeltaHamiltonian(dca1, seq1, site, oldState1, newState1, rootPartials1, freqs1, data1);
				for (int newState2 = 0; newState2 < stateCount; newState2++) {
					probs[newState1 * stateCount + newState2] = 
							  p1 
							+ computeDeltaHamiltonian(dca2, seq2, site, oldState2, newState2, rootPartials2, freqs2, data2) 
							+ mbpPrior.logP(newState1, newState2);
				}
			}
			
			// find max
			double max = probs[0];
			for (int k = 0; k < probs.length; k++) {
				if (probs[k] > max) {
					max = Math.max(probs[k], max);
				}
			}
			
			// to real space
			for (int k = 0; k < probs.length; k++) {
				probs[k] = Math.exp(probs[k] - max);
			}
			double sum = 0;
			for (double d : probs) {
				sum += d;
			}
			try {
				int newState = (sum > 0) ? Randomizer.randomChoicePDF(probs) : Randomizer.nextInt(freqs1.length);
				//int newState = Randomizer.randomChoicePDF(freqs);
				//int newState = Randomizer.nextInt(freqs.length);
				//int newState = iMax;
				seq1[site] = newState / stateCount;
				seq2[site] = newState % stateCount;
				logHR += Math.log(probs[oldState1 * stateCount + oldState2]) - Math.log(probs[newState]);
			} catch (Error e) {
				int newState = Randomizer.nextInt(freqs1.length);
				seq1[site] = newState;
			}
		}
		
		
		for (int i = 0; i < n; i++) {
			sequence1.setValue(i, seq1[i]);
			sequence2.setValue(i, seq2[i]);
		}
		
		return logHR;
	}

	
    private double computeDeltaHamiltonian(DCA dca, int[] seq, int i, int oldState, int newState, 
    		double [] rootPartials, double [] freqs, Alignment data1) {
        double delta =  (temperaturFactor > 0) ? 
        		temperaturFactor * DCASequenceSimulator.computeDeltaHamiltonian(dca, seq, i, oldState, newState)
        		: 0;
        
        int patternIndexOffset = data1.getPatternIndex(i) * stateCount;

        // Change in transition probability
        delta += Math.log(rootPartials[patternIndexOffset + newState]);// * freqs[newState]);// - Math.log(rootPartials[patternIndexOffset + oldState]);
        		

        return delta;
    }    

}
