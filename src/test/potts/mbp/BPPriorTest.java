package test.potts.mbp;

import java.util.Arrays;

import beast.base.evolution.substitutionmodel.WAG;
import beast.base.inference.parameter.IntegerParameter;
import potts.mbp.BPPrior;
import potts.mbp.MBPPrior;

public class BPPriorTest {

	public static void main(String[] args) {
		MBPPrior prior = new BPPrior();
		MBPPrior prior2 = new MBPPrior();
		IntegerParameter sequence1 = new IntegerParameter();
		sequence1.initByName("lower", 0, "upper", 19, "value", "0");
		IntegerParameter sequence2 = new IntegerParameter();
		sequence2.initByName("lower", 0, "upper", 19, "value", "0");
		WAG wag = new WAG();
		wag.initAndValidate();
		prior.initByName("sequence1", sequence1, "sequence2", sequence2, "substModel", wag);

		for (int i = 0; i < 20; i++) {
			System.out.println(Arrays.toString(prior.match[i]));
		}
		System.out.println();
		prior2.initByName("sequence1", sequence1, "sequence2", sequence2, "substModel", wag);
		for (int i = 0; i < 20; i++) {
			System.out.println(Arrays.toString(prior2.match[i]));
		}

	}
}
