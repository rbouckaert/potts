package potts.mbp;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.inference.Distribution;
import beast.base.inference.Operator;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.util.Randomizer;
import potts.likelihood.TreeLikelihoodWithRootStates;

@Description("Gibbs sampler for root sequence in TreeLikelihoodWithRootStates")
public class RootSequenceResampler extends Operator {

	final public Input<IntegerParameter> sequence1Input = new Input<>("sequence1",
			"specifies sequence states to calculated likelihood for", Validate.REQUIRED);
	final public Input<IntegerParameter> sequence2Input = new Input<>("sequence2",
			"specifies sequence states to calculated likelihood for", Validate.REQUIRED);

	final public Input<TreeLikelihoodWithRootStates> likelihood1Input = new Input<>("likelihood1", 
			"TreeLikelihood With Root States to provide root partials 1", Validate.REQUIRED);
	final public Input<TreeLikelihoodWithRootStates> likelihood2Input = new Input<>("likelihood2", 
			"TreeLikelihood With Root States to provide root partials 2", Validate.REQUIRED);
	final public Input<MBPPrior> mbpPriorInput = new Input<>("mbpPrior", 
			"Middle base pair prior", Validate.REQUIRED);

	final public Input<State> stateInput = new Input<>("state", "elements of the state space", Validate.REQUIRED);

    final public Input<Distribution> posteriorInput =
            new Input<>("posterior", "probability distribution to sample over (e.g. a posterior)",
                    Input.Validate.REQUIRED);

	IntegerParameter sequence1, sequence2;
	TreeLikelihoodWithRootStates likelihood1, likelihood2;
	MBPPrior mbpPrior;
	Distribution posterior;
	
	int stateCount = 20;
	State state;
	
	@Override
	public void initAndValidate() {
		state = stateInput.get();
		posterior = posteriorInput.get();
		
		sequence1 = sequence1Input.get();
		sequence2 = sequence2Input.get();
		
		likelihood1 = likelihood1Input.get();
		likelihood2 = likelihood2Input.get();
		mbpPrior = mbpPriorInput.get();
	}

	@Override
	public double proposal() {
		int n = sequence1.getValues().length;
		
		double [][] logProbs = new double[n][stateCount*stateCount];
		

		double logHR = 0;

		for (int i = 0; i < 20; i++) {
			for (int j = 0; j < 20; j++) {
				state.store(-1);
	            state.storeCalculationNodes();
	            state.checkCalculationNodesDirtiness();
				for (int k = 0; k < n; k++) {
					sequence1.setValue(k, i);
					sequence2.setValue(k, j);
				}
				posterior.calculateLogP();
				
				double [] rootPartials1 = likelihood1.getRootPartials();
				double [] rootPartials2 = likelihood2.getRootPartials();

				
				for (int k = 0; k < n; k++) {
					logProbs[k][i*stateCount+j] = rootPartials1[k] + rootPartials2[k] 
							+ mbpPrior.logP(i,j);
				}
				
                state.restore();
                state.restoreCalculationNodes();
                state.setEverythingDirty(false);
			}
		}
		
		
		int [] seq = new int[n];
		
		for (int i = 0; i < n; i++) {
			double [] probs = logProbs[i];
			
			// find max
			double max = probs[0];
			int iMax = 0;
			for (int k = 0; k < probs.length; k++) {
				if (probs[k] > max) {
					iMax = k;
					max = Math.max(probs[k], max);
				}
			}
			
			// to real space
			for (int k = 0; k < probs.length; k++) {
				probs[k] = Math.exp(probs[k] - max);
			}
			int newState = Randomizer.randomChoicePDF(probs);
			seq[i] = newState;
		}
		
		
		for (int i = 0; i < n; i++) {
			sequence1.setValue(i, seq[i] / stateCount);
			sequence2.setValue(i, seq[i] % stateCount);
		}
		
		return logHR;
	}
	
}
