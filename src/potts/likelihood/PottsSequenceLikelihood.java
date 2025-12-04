package potts.likelihood;

import java.util.List;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import potts.dca.DCA;

@Description("Contribution of Potts model for one sequence")
public class PottsSequenceLikelihood extends Distribution {

	final public Input<DCA> dcaInput = new Input<>("dca", "DCA json file previously trained on an alignment",
			Validate.REQUIRED);

	final public Input<IntegerParameter> sequenceInput = new Input<>("sequence",
			"specifies sequence states to calculated likelihood for", Validate.REQUIRED);

    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);

	final public Input<Double> temperatureInput = new Input<>("temperature", "temperature balancing the effect of Potts model Pt and substitution model Ps. "
			+ "Mutations are chosen proportional to Pt^1/t * Ps (thus the Potss model weighted by 1/t). "
			+ "If negative, Potts model contributino will be ignored", 1.0);

	/** direct coupling analysis **/
	DCA dca;
	IntegerParameter sequence;
	int siteCount, stateCount;
	double temperatureFactor;

	@Override
	public void initAndValidate() {
		sequence = sequenceInput.get();
		siteCount = sequence.getDimension();
		stateCount = dataInput.get().getMaxStateCount();

		dca = dcaInput.get();

		if (dca.getSiteCount() != siteCount) {
			throw new IllegalArgumentException("Site count of DCA (" + dca.getSiteCount() + ") and alignment ("
					+ siteCount + ") should match");
		}
		if (dca.getStateCount() != stateCount + 1) {
			throw new IllegalArgumentException("State count of DCA (" + dca.getStateCount() + ") and alignment ("
					+ stateCount + ") should match");
		}

		temperatureFactor = 1.0 / temperatureInput.get();
		super.initAndValidate();
	}

	@Override
	public double calculateLogP() {
		logP = 0;
		
		if (temperatureFactor > 0) {
			Integer[] seq = sequence.getValues();
			double[][] h = dca.getH();
			double[][][][] J = dca.getJ();
			for (int i = 0; i < siteCount; i++) {
				int a = seq[i];
				logP += h[i][a];
				for (int j = i + 1; j < siteCount; j++) {
					int b = seq[j];
					logP += J[i][j][a][b];
				}
			}
	
			logP *= temperatureFactor;
		}
		
		return logP;
	}

	@Override
	public List<String> getArguments() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public List<String> getConditions() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void sample(State state, Random random) {
		// TODO Auto-generated method stub

	}

}
