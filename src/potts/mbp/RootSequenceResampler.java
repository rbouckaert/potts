package potts.mbp;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.inference.Operator;
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

	IntegerParameter sequence1, sequence2;
	TreeLikelihoodWithRootStates likelihood1, likelihood2;
	Alignment data1, data2;
	MBPPrior mbpPrior;
	
	int stateCount = 20;
	
	@Override
	public void initAndValidate() {
		
		sequence1 = sequence1Input.get();
		sequence2 = sequence2Input.get();
		
		likelihood1 = likelihood1Input.get();
		likelihood2 = likelihood2Input.get();
		
		data1 = likelihood1.dataInput.get();
		data2 = likelihood2.dataInput.get();
		
		mbpPrior = mbpPriorInput.get();
	}

	@Override
	public double proposal() {
		boolean reverse = mbpPrior.reverseInput.get();
		
		int n = sequence1.getValues().length;
		
		double [][] logProbs = new double[n][stateCount*stateCount];
		

		double logHR = 0;

		double [] rootPartials1 = likelihood1.getRootPartials();
		double [] rootPartials2 = likelihood2.getRootPartials();

		for (int i = 0; i < 20; i++) {

	        for (int j = 0; j < 20; j++) {
				
		        for (int k = 0; k < n; k++) {
			        int patternIndexOffset1 = data1.getPatternIndex(k) * stateCount;
			        int patternIndexOffset2 = data2.getPatternIndex(reverse ? n - k - 1: k) * stateCount;
					double logP = 
							Math.log(rootPartials1[patternIndexOffset1 + i]) +
							Math.log(rootPartials2[patternIndexOffset2 + j]) +
							mbpPrior.logP(i, j);
					if (Double.isNaN(logP)) {
						logP = Double.NEGATIVE_INFINITY;
					}
					logProbs[k][i*stateCount + j] = logP;
				}
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
			if (!Double.isNaN(probs[0])) {
				int newState = Randomizer.randomChoicePDF(probs);
				seq[i] = newState;
			}
		}
		
		
		for (int i = 0; i < n; i++) {
			sequence1.setValue(i, seq[i] / stateCount);
			sequence2.setValue(i, seq[i] % stateCount);
		}
		
		logHR = Double.POSITIVE_INFINITY;
		return logHR;
	}
	
}
