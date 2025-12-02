package potts.likelihood;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Random;

import org.json.JSONException;
import org.json.JSONObject;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;
import beastfx.app.inputeditor.BeautiDoc;
import potts.dca.DCA;

@Description("Contribution of Potts model for one sequence")
public class PottsSequenceLikelihood extends Distribution {

	final public Input<File> dcaInput = new Input<>("dca", "DCA json file previously trained on an alignment",
			Validate.REQUIRED);

	final public Input<IntegerParameter> sequenceInput = new Input<>("sequence",
			"specifies sequence states to calculated likelihood for", Validate.REQUIRED);

    final public Input<Alignment> dataInput = new Input<>("data", "sequence data for the beast.tree", Validate.REQUIRED);

	/** direct coupling analysis **/
	DCA dca;
	IntegerParameter sequence;
	int siteCount, stateCount;

	@Override
	public void initAndValidate() {
		sequence = sequenceInput.get();
		siteCount = sequence.getDimension();
		stateCount = dataInput.get().getMaxStateCount();

		try {
			dca = new DCA();
			String json = BeautiDoc.load(dcaInput.get());
			dca.fromJSON(new JSONObject(json));

			if (dca.getSiteCount() != siteCount) {
				throw new IllegalArgumentException("Site count of DCA (" + dca.getSiteCount() + ") and alignment ("
						+ siteCount + ") should match");
			}
			if (dca.getStateCount() != stateCount + 1) {
				throw new IllegalArgumentException("State count of DCA (" + dca.getStateCount() + ") and alignment ("
						+ stateCount + ") should match");
			}
		} catch (IOException | JSONException e) {
			e.printStackTrace();
		}

		super.initAndValidate();
	}

	@Override
	public double calculateLogP() {
		logP = 0;

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
