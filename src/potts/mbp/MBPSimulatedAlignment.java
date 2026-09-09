package potts.mbp;

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

@Description("Simulator of MBP linked root sequences and their alignments")
public class MBPSimulatedAlignment extends SimulatedAlignment implements StateNodeInitialiser {
	final public Input<MBPPrior> mbpPriorInput = new Input<>("mbpPrior", 
			"Middle base pair prior", Validate.REQUIRED);
	final public Input<Alignment> otherInput = new Input<>("other", "the second alignment to be simulated", Validate.REQUIRED);

	final public Input<Tree> tree2Input = new Input<>("tree2", "phylogenetic beast.tree with sequence data in the leafs", Validate.REQUIRED);

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
            seq[i] = firstSeq ? sequence1.getValue(i) : sequence2.getValue(reverse ? m_sequenceLength-i-1: i);
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
