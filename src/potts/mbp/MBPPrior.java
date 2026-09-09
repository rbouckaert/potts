package potts.mbp;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Random;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.inference.Distribution;
import beast.base.inference.State;
import beast.base.inference.parameter.IntegerParameter;

@Description("Middle base pair prior: takes 2 amino acide sequences as input and for every "
		+ "site contribute p for paired and 1-p for unpaired")
public class MBPPrior extends Distribution {
	final public Input<IntegerParameter> sequence1Input = new Input<>("sequence1", 
			  "specifies site specific sequence", Validate.REQUIRED);
	final public Input<IntegerParameter> sequence2Input = new Input<>("sequence2", 
			  "specifies site specific sequence to pair agains", Validate.REQUIRED);
	final public Input<Boolean> reverseInput = new Input<>("reverse", "whether to reverse the second sequence (default) or keep it in order", true);
	final public Input<Double> pInput = new Input<>("p", "probability of basepair matching for a single site", 0.99);

	final public Input<SubstitutionModel> substModelInput =
            new Input<>("substModel", "substitution model along branches in the beast.tree", null, Validate.REQUIRED);
	
	
	// the two sequences to score
	IntegerParameter sequence1, sequence2;
	SubstitutionModel substModel;
	
	// per site contribution for matching and not-matching MBP sites
	double logPmatch, logPmiss;

	// 20x20 matrix to score how well two amino acid characters have a middle base pair in common
	double [][] match;
	
	double [] logFreqs;
	
	boolean reverse;

	@Override
	public void initAndValidate() {
		super.initAndValidate();
		
		reverse = reverseInput.get();
		substModel = substModelInput.get();
		
		sequence1 = sequence1Input.get();
		sequence2 = sequence2Input.get();
		
		if (sequence1.getUpper() < 19 || sequence2.getUpper() < 19) {
			throw new IllegalArgumentException("Expected amino acid sequences");
		}
		
		if (sequence1.getDimension() != sequence2.getDimension()) {
			throw new IllegalArgumentException("Sequences have different length: cannot basepair these sequences");
		}
		
		logPmatch = Math.log(pInput.get());
		logPmiss = Math.log(1.0 - pInput.get());
		
		
		// middle basepair map from amino acid to nucleotide
		Map<String, String> mapAA2N = new HashMap<>();
		
	    // T middle
	    mapAA2N.put("F","T");
	    mapAA2N.put("L","T");
		mapAA2N.put("I","T");
		mapAA2N.put("M","T");
		mapAA2N.put("V","T");

		// C middle (S is handled separately)
	    mapAA2N.put("P","C");
		mapAA2N.put("T","C");
		mapAA2N.put("A","C");
	    
	    // A middle
	    mapAA2N.put("Y","A");
		mapAA2N.put("H","A");
		mapAA2N.put("Q","A");
		mapAA2N.put("N","A");
		mapAA2N.put("K","A");
		mapAA2N.put("D","A");
		mapAA2N.put("E","A");

	    // G middle (S is handled separately)
	    mapAA2N.put("C","G");
		mapAA2N.put("W","G");
		mapAA2N.put("R","G");
		mapAA2N.put("G","G");

		
		Aminoacid aminoacid = new Aminoacid();
		match = new double[20][20];
		for (int i = 0; i < 20; i++) {
			String c1 = aminoacid.getCharacter(i);
			for (int j = 0; j < 20; j++) {
				String c2 = aminoacid.getCharacter(j);
				if (c1.equals("S") || c2.equals("S")) {
					if (!c1.equals("S")) {
						String tmp = c1;
						c1 = c2;
						c2 = tmp;
					}
					if (!c2.equals("S")) {
						if (c2.equals("C")) {
							match[i][j] = 1.0/3.0;
						}
						if (c2.equals("G")) {
							match[i][j] = 2.0/3.0;
						}
					} else {
						match[i][j] = 2.0 * 1.0/3.0 * 2.0/3.0;
					}
				} else {
					String mbp1 = mapAA2N.get(c1);
					String mbp2 = mapAA2N.get(c2);
					if ((mbp1.equals("A") && mbp2.equals("T")) ||
						(mbp1.equals("T") && mbp2.equals("A")) ||
						(mbp1.equals("C") && mbp2.equals("G")) ||
						(mbp1.equals("G") && mbp2.equals("C"))) {
						match[i][j] = 1.0;
					}
				}
				
			}
		}
		
		logFreqs = new double[20];
		double [] freqs = substModel.getFrequencies();
		for (int i = 0; i < freqs.length; i++) {
			logFreqs[i] = Math.log(freqs[i]);
		}
	}
	
	
	@Override
	public double calculateLogP() {
		logP = 0;
		
		double [] freqs = substModel.getFrequencies();
		for (int i = 0; i < freqs.length; i++) {
			logFreqs[i] = Math.log(freqs[i]);
		}
		
		int n = sequence1.getDimension();
		for (int i = 0; i < n; i++) {
			int site1 = sequence1.getValue(i);
			int site2 = sequence2.getValue(reverse ? n - i - 1: i);
			
			// frequencies contribution
			if (site1 >= 0 & site1 <= 20) {
				logP += freqs[site1];
			}
			if (site2 >= 0 & site2 <= 20) {
				logP += freqs[site2];
			}
			
			// base pair matching contribution
			if (site1 >=0 && site1 <= 20 && site2 >= 0 && site2 <= 20) {
				logP += match[site1][site2] * logPmatch + (1 - match[site1][site2]) * logPmiss;
			}
		}
        return logP;
	}
	
	public double logP(int site1, int site2) {
		double logP = 0;
		// frequencies contribution
		if (site1 >= 0 & site1 <= 20) {
			logP += logFreqs[site1];
		}
		if (site2 >= 0 & site2 <= 20) {
			logP += logFreqs[site2];
		}
		
		// base pair matching contribution
		if (site1 >=0 && site1 <= 20 && site2 >= 0 && site2 <= 20) {
			logP += match[site1][site2] * logPmatch + (1 - match[site1][site2]) * logPmiss;
		}
		return logP;
	}
	
	
	@Override
	public List<String> getArguments() {
		List<String> args = new ArrayList<>();
		args.add(sequence1Input.get().getID());
		args.add(sequence2Input.get().getID());
		return args;
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
