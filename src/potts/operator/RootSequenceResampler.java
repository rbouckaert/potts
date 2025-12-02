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
import potts.likelihood.TreeLikelihoodWithRootStates;

@Description("Gibbs sampler for root sequence in TreeLikelihoodWithRootStates")
public class RootSequenceResampler extends Operator {

	final public Input<DCA> dcaInput = new Input<>("dca", "DCA json file previously trained on an alignment",
			Validate.REQUIRED);

	final public Input<IntegerParameter> sequenceInput = new Input<>("sequence",
			"specifies sequence states to calculated likelihood for", Validate.REQUIRED);

	final public Input<TreeLikelihoodWithRootStates> likelihoodInput = new Input<>("likelihood", 
			"TreeLikelihood With Root States to provide root partials", Validate.REQUIRED);
	
	final public Input<Integer> stepCountInput = new Input<>("stepCount", 
			"number of times a site is resampled", 10000);

	DCA dca;
	IntegerParameter sequence;
	TreeLikelihoodWithRootStates likelihood;
	int stepCount, stateCount;
	Alignment data;
	
	@Override
	public void initAndValidate() {
		dca = dcaInput.get();
		sequence = sequenceInput.get();
		likelihood = likelihoodInput.get();
		data = likelihood.dataInput.get();
		stepCount = stepCountInput.get();
		stateCount = sequence.getUpper();
	}

	@Override
	public double proposal() {
		Integer [] values = sequence.getValues();
		int [] seq = new int[values.length];
		for (int i = 0; i < seq.length; i++) {
			seq[i] = values[i];
		}
		
		double [] rootPartials = likelihood.getRootPartials();
		
		for (int i = 0; i < stepCount; i++) {
			int site = Randomizer.nextInt(values.length);
			int oldState = seq[site];
			int newState = Randomizer.nextInt(stateCount);
			if (oldState != newState) {
                double deltaH = computeDeltaHamiltonian(dca, seq, site, oldState, newState, rootPartials); 
                
                // Metropolis Criterion
                if (deltaH >= 0 || Randomizer.nextDouble() < Math.exp(deltaH)) {
                    seq[site] = newState; // Accept mutation
                }
			}
		}
		
		
		for (int i = 0; i < values.length; i++) {
			sequence.setValue(i, seq[i]);
		}
		
		return Double.POSITIVE_INFINITY;
	}

	
    private double computeDeltaHamiltonian(DCA dca, int[] seq, int i, int oldState, int newState, 
    		double [] rootPartials) {
        double delta = DCASequenceSimulator.computeDeltaHamiltonian(dca, seq, i, oldState, newState);
        
        int patternIndexOffset = data.getPatternIndex(i) * stateCount;

        // Change in transition probability
        delta += Math.log(rootPartials[patternIndexOffset + newState]) - Math.log(rootPartials[patternIndexOffset + oldState]);
        		

        return delta;
    }    

}
