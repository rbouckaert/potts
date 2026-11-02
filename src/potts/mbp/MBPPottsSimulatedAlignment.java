package potts.mbp;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.List;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.StateNode;
import beast.base.inference.StateNodeInitialiser;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.util.Randomizer;
import beastfx.app.seqgen.SimulatedAlignment;
import potts.dca.DCA;
import potts.dca.DCASequenceSimulator;
import potts.likelihood.PottsSequenceLikelihood;
import potts.likelihood.TreeLikelihoodWithRootStates;

@Description("Simulator of MBP linked root sequences and their alignments")
public class MBPPottsSimulatedAlignment extends SimulatedAlignment implements StateNodeInitialiser {
	final public Input<MBPPrior> mbpPriorInput = new Input<>("mbpPrior", 
			"Middle base pair prior", Validate.REQUIRED);
	final public Input<Alignment> otherInput = new Input<>("other", "the second alignment to be simulated", Validate.REQUIRED);

	final public Input<Tree> tree2Input = new Input<>("tree2", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

	final public Input<File> rootSeqFileInput = new Input<>("rootSeqFile", "file to store simulated root sequences in tsv format");

	final public Input<PottsSequenceLikelihood> pottsLikelihood1Input = new Input<>("PottsSequenceLikelihood1", 
			"PottsSequenceLikelihood over Root States to provide root partials 1", Validate.REQUIRED);
	final public Input<PottsSequenceLikelihood> pottsLikelihood2Input = new Input<>("PottsSequenceLikelihood2", 
			"PottsSequenceLikelihood over Root States to provide root partials 2", Validate.REQUIRED);
	final public Input<TreeLikelihoodWithRootStates> likelihood1Input = new Input<>("likelihood1", 
			"TreeLikelihood With Root States to provide root partials 1", Validate.REQUIRED);
	final public Input<TreeLikelihoodWithRootStates> likelihood2Input = new Input<>("likelihood2", 
			"TreeLikelihood With Root States to provide root partials 2", Validate.REQUIRED);

	private IntegerParameter sequence1, sequence2;
	private boolean firstSeq, reverse; 
	
	@Override
	public void initAndValidate() {
		
		MBPPrior mbpPrior;

		mbpPrior = mbpPriorInput.get();
		double [] p = new double[20*20];
		for (int i = 0; i < 20; i++) {
			for (int j = 0; j < 20; j++) {
				p[i*20+j] = mbpPrior.logP(i,j);
			}
		}
		double max = p[0];
		for (int i = 0; i < p.length; i++) {
			max = Math.max(max, p[i]);
		}
		for (int i = 0; i < p.length; i++) {
			p[i] = Math.exp(p[i] - max);
		}
		
		sequence1 = mbpPrior.sequence1Input.get();
		sequence2 = mbpPrior.sequence2Input.get();
		int n = sequence1.getDimension();
		reverse = mbpPrior.reverseInput.get();
		for (int i = 0; i < n; i++) {
			int k = Randomizer.randomChoicePDF(p);
			int state1 = k / 20;
			int state2 = k % 20;
			sequence1.setValue(i, state1);
			sequence2.setValue(reverse ? n - i - 1 : i, state2);
		}

		resample(pottsLikelihood1Input.get(),
				pottsLikelihood2Input.get(),
				likelihood1Input.get(),
				likelihood2Input.get(),
				sequence1,
				sequence2,
				mbpPrior
				);
		
		if (rootSeqFileInput.get() != null) {
			try {
				PrintStream out = new PrintStream(rootSeqFileInput.get());
				for (int i = 0; i < sequence1.getDimension(); i++) {
					out.print(sequence1.getValue(i));
					out.print("\t");
				}
				for (int i = 0; i < sequence2.getDimension(); i++) {
					out.print(sequence2.getValue(i));
					out.print("\t");
				}
				out.println();
				out.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		// simulate the other alignment first
		Tree tree = m_treeInput.get();
		Tree tree2 = tree2Input.get();
		m_treeInput.set(tree2);
		firstSeq = false;
		super.initAndValidate();
		
		Alignment second = otherInput.get();
		List<Sequence> seqs2 = second.sequenceInput.get();
		seqs2.clear();
		seqs2.addAll(this.sequences);
		second.initAndValidate();
		
		// now simulate this alignment
		firstSeq = true;
		m_treeInput.set(tree);
		super.initAndValidate();
		Alignment first = m_data.get();
		List<Sequence> seqs1 = first.sequenceInput.get();
		seqs1.clear();
		seqs1.addAll(this.sequences);
		first.initAndValidate();
		
		
		likelihood1Input.get().initAndValidate();
		likelihood2Input.get().initAndValidate();

	}
	
    private void resample(PottsSequenceLikelihood pottsSequenceLikelihood1,
			PottsSequenceLikelihood pottsSequenceLikelihood2, 
			TreeLikelihoodWithRootStates likelihood1,
			TreeLikelihoodWithRootStates likelihood2, 
			IntegerParameter sequence1, 
			IntegerParameter sequence2,
			MBPPrior mbpPrior) {
    	int totalStepCount = sequence1.getDimension();
    	int stateCount = 20;
    	double temperaturFactor = 1.0;
    	
    	DCA dca1 = pottsSequenceLikelihood1.dcaInput.get();
    	DCA dca2 = pottsSequenceLikelihood2.dcaInput.get();

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
				double p1 = computeDeltaHamiltonian(dca1, seq1, site, oldState1, newState1, freqs1, temperaturFactor);
				for (int newState2 = 0; newState2 < stateCount; newState2++) {
					probs[newState1 * stateCount + newState2] = 
							  p1 
							+ computeDeltaHamiltonian(dca2, seq2, site, oldState2, newState2, freqs2, temperaturFactor) 
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
		
	}

	
    private double computeDeltaHamiltonian(DCA dca, int[] seq, int i, int oldState, int newState, 
    		double [] freqs, double temperaturFactor) {
        double delta =  (temperaturFactor > 0) ? 
        		temperaturFactor * DCASequenceSimulator.computeDeltaHamiltonian(dca, seq, i, oldState, newState)
        		: 0;

        return delta;
    }    
    
    
	public void simulate() {
        Node root = m_tree.getRoot();


        double[] categoryProbs = m_siteModel.getCategoryProportions(root);
        int[] category = new int[m_sequenceLength];
        for (int i = 0; i < m_sequenceLength; i++) {
            category[i] = Randomizer.randomChoicePDF(categoryProbs);
        }

        int[] seq = new int[m_sequenceLength];
        for (int i = 0; i < m_sequenceLength; i++) {
            seq[i] = firstSeq ? sequence1.getValue(i) : sequence2.getValue(i);
        }

        traverse(root, seq, category);

    } // simulate

	@Override
	public void initStateNodes() {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void getInitialisedStateNodes(List<StateNode> stateNodes) {
		// TODO Auto-generated method stub
		
	}

}
