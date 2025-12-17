package potts.operator;

import java.text.DecimalFormat;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Operator;
import beast.base.inference.operator.kernel.KernelDistribution;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.util.Randomizer;
import potts.dca.DCA;
import potts.dca.DCASequenceSimulator;
import potts.likelihood.PottsSequenceLikelihood;
import potts.likelihood.TreeLikelihoodWithRootStates;

@Description("Gibbs sampler for root sequence in TreeLikelihoodWithRootStates")
public class RootHeightAndSequenceResampler extends Operator {

	final public Input<PottsSequenceLikelihood> pottsSeqLikelihoodInput = new Input<>("pottsSeqLikelihood", "DCA json file previously trained on an alignment",
			Validate.REQUIRED);

	final public Input<IntegerParameter> sequenceInput = new Input<>("sequence",
			"specifies sequence states to calculated likelihood for", Validate.REQUIRED);

	final public Input<TreeLikelihoodWithRootStates> likelihoodInput = new Input<>("likelihood", 
			"TreeLikelihood With Root States to provide root partials", Validate.REQUIRED);
	
	final public Input<Integer> stepCountInput = new Input<>("stepCount", 
			"average number of times a site is resampled", 20);

	
    public final Input<KernelDistribution> kernelDistributionInput = new Input<>("kernelDistribution", "provides sample distribution for proposals", 
    		KernelDistribution.newDefaultKernelDistribution());

    final public Input<Boolean> optimiseInput = new Input<>("optimise", "flag to indicate that the scale factor is automatically changed in order to achieve a good acceptance rate (default true)", true);

    
    final public Input<Double> scaleUpperLimit = new Input<>("upper", "Upper Limit of scale factor", 1.0 - 1e-8);
    final public Input<Double> scaleLowerLimit = new Input<>("lower", "Lower limit of scale factor", 1e-8);
    public final Input<Double> scaleFactorInput = new Input<>("scaleFactor", "scaling factor: range from 0 to 1. Close to zero is very large jumps, close to 1.0 is very small jumps.", 0.75);

    private double scaleFactor;
    private double upper, lower;

	DCA dca;
	IntegerParameter sequence;
	TreeLikelihoodWithRootStates likelihood;
	int totalStepCount, stateCount;
	Alignment data;
	double temperaturFactor;
    Tree tree;
    protected KernelDistribution kernelDistribution;

	@Override
	public void initAndValidate() {
        scaleFactor = scaleFactorInput.get();
        upper = scaleUpperLimit.get();
        lower = scaleLowerLimit.get();

        PottsSequenceLikelihood psl = pottsSeqLikelihoodInput.get();
		dca = psl.dcaInput.get();
		temperaturFactor = 1.0 / psl.temperatureInput.get();
		
		sequence = sequenceInput.get();
		likelihood = likelihoodInput.get();
		tree = (Tree) likelihood.treeInput.get();
		data = likelihood.dataInput.get();
		totalStepCount = stepCountInput.get() * sequence.getDimension();
		stateCount = dca.getStateCount() - 1;

		kernelDistribution = kernelDistributionInput.get();

	}

	@Override
	public double proposal() {

        final Node root = tree.getRoot();
        final double scale = getScaler(root.getNr(), root.getHeight());
        final double newHeight = root.getHeight() * scale;

        if (newHeight < Math.max(root.getLeft().getHeight(), root.getRight().getHeight())) {
            return Double.NEGATIVE_INFINITY;
        }
        root.setHeight(newHeight);
           		
        // update partials in likelihood
        likelihood.calculateLogP();
		
		
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
//				if (oldState == newState) {
//					probs[newState] = 0;
//				} else {
					probs[newState] = computeDeltaHamiltonian(dca, seq, site, oldState, newState, rootPartials, freqs);
//				}
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

        return Math.log(scale);
	}

	
    private double computeDeltaHamiltonian(DCA dca, int[] seq, int site, int oldState, int newState, 
    		double [] rootPartials, double [] freqs) {
        double delta =  (temperaturFactor > 0) ? 
        		temperaturFactor * DCASequenceSimulator.computeDeltaHamiltonian(dca, seq, site, oldState, newState)
        		: 0;
        
        int patternIndexOffset = data.getPatternIndex(site) * stateCount;

        // Change in transition probability
        delta += Math.log(rootPartials[patternIndexOffset + newState]);// * freqs[newState]);// - Math.log(rootPartials[patternIndexOffset + oldState]);
        		

        return delta;
    }    

    
    @Override
    public double getCoercableParameterValue() {
        return scaleFactor;
    }

    @Override
    public void setCoercableParameterValue(final double value) {
        scaleFactor = Math.max(Math.min(value, upper), lower);
    }

    protected double getScaler(int i, double value) {
        return kernelDistribution.getScaler(i, value, getCoercableParameterValue());
    }

    @Override
    public void optimize(double logAlpha) {
        // must be overridden by operator implementation to have an effect
    	if (optimiseInput.get()) {
	        double delta = calcDelta(logAlpha);
	        double scaleFactor = getCoercableParameterValue();
	        delta += Math.log(scaleFactor);
	        scaleFactor = Math.exp(delta);
	        setCoercableParameterValue(scaleFactor);
    	}
    }
    
    @Override
    public double getTargetAcceptanceProbability() {
    	return 0.3;
    }
    


    @Override
    public String getPerformanceSuggestion() {
        double prob = m_nNrAccepted / (m_nNrAccepted + m_nNrRejected + 0.0);
        double targetProb = getTargetAcceptanceProbability();

        double ratio = prob / targetProb;
        if (ratio > 2.0) ratio = 2.0;
        if (ratio < 0.5) ratio = 0.5;

        // new scale factor
        double newWindowSize = getCoercableParameterValue() * ratio;

        DecimalFormat formatter = new DecimalFormat("#.###");
        if (prob < 0.10 || prob > 0.40) {
            return "Try setting scale factor to about " + formatter.format(newWindowSize);
        } else return "";
    }

}
