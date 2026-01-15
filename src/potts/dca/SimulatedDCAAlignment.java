package potts.dca;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.util.Randomizer;
import potts.datatype.AminoacidPlus;

@Description("Simulate alignment on a tree using the Potts model for root sequence")
public class SimulatedDCAAlignment extends Alignment {
	final public Input<DCA> dcaInput = new Input<>("dca", "DCA previously trained on an alignment", Validate.REQUIRED);
	final public Input<Tree> treeInput = new Input<>("tree", "tree file to simulate alignment for", Validate.REQUIRED);

	final public Input<Integer> rootStepCountInput = new Input<>("rootStepCount", "number of times each site for the root sequence is resampled (at initialisation only)", 1000);
	final public Input<Integer> stepCountInput = new Input<>("stepCount", "number of times each site is resampled when sampling sequences conditioned on parent or parents and children", 100);
	final public Input<Double> temperatureInput = new Input<>("temperature", "temperature balancing the effect of Potts model Pt and substitution model Ps. "
			+ "Mutations are chosen proportional to Pt^1/t * Ps (thus 1/t", 1.0);
    final public Input<SiteModel.Base> siteModelInput = new Input<>("siteModel", "site model for leafs in the beast.tree", Validate.REQUIRED);
    final public Input<BranchRateModel.Base> branchRateModelInput = new Input<>("branchRateModel",
            "A model describing the rates on the branches of the beast.tree.");
    final public Input<Long> localSeedInput = new Input<>(
            "seed",
            "Optional local random seed for simulating this alignment. If not set, global seed is used.",
            Input.Validate.OPTIONAL
    );

	int stepCount;
	int rootStepCount;
	double temperatureFactor;
	double [][] probabilities;
	int sequenceLength, categoryCount, stateCount;
	Tree tree;
	
    /**
     * site model used for generating samples *
     */
    protected SiteModel.Base siteModel;
    /**
     * branch rate model used for generating samples *
     */
    protected BranchRateModel branchRateModel;

	@Override
	public void initAndValidate() {
		siteModel = siteModelInput.get();
		categoryCount = siteModel.getCategoryCount();
		stateCount = siteModel.getSubstitutionModel().getStateCount();
		probabilities = new double[categoryCount][(stateCount+1) * stateCount];
		branchRateModel = branchRateModelInput.get();
		if (branchRateModel == null) {
			branchRateModel = new StrictClockModel();
			((BEASTInterface)branchRateModel).initAndValidate();
		}
		
        if (userDataTypeInput.get() != null) {
            m_dataType = userDataTypeInput.get();
        } else {
            initDataType();
        }
				
        Long customSeed = localSeedInput.get();
        long originalSeed = Randomizer.getSeed();
        long seedToUse = customSeed != null ? customSeed : originalSeed;

        if (customSeed != null) {
            Log.info.println();
            Log.info.println("Random number seed for alignment simulation: " + customSeed);
            Log.info.println();
        }
		
		try {
            Randomizer.setSeed(seedToUse);
			simulate();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
        } finally {
            Randomizer.setSeed(originalSeed);
        }
		dataTypeInput.setValue(m_dataType.toString(), this);
		super.initAndValidate();
	}

	public void simulate() throws Exception {

		stepCount = stepCountInput.get();
		rootStepCount = rootStepCountInput.get();
		temperatureFactor = 1.0 / temperatureInput.get();
		tree =  treeInput.get();
		
		DCA dca = dcaInput.get();
		sequenceLength = dca.getSiteCount();
		

		if (m_dataType == null) {
			m_dataType = dca.stateCount == 21 ? new Aminoacid() : new Nucleotide();
		}
		SubstitutionModel substModel = siteModel.getSubstitutionModel();
		if (dca.stateCount -1 != substModel.getStateCount() && !(m_dataType instanceof AminoacidPlus)) {
			throw new IllegalArgumentException("DCA state count should be 1 + model state count");
		} else if (dca.stateCount != substModel.getStateCount() && m_dataType instanceof AminoacidPlus) {
			throw new IllegalArgumentException("DCA state count should equal model state count");
		}
		
        Node root = tree.getRoot();
        double[] categoryProbs = siteModel.getCategoryProportions(root);
        int[] category = new int[sequenceLength];
        for (int i = 0; i < sequenceLength; i++) {
            category[i] = Randomizer.randomChoicePDF(categoryProbs);
        }

        int [] seq = new int[sequenceLength];
		DCASequenceSimulator.sampleRootSequence(seq, dca, rootStepCount, dca.stateCount);
		
		try {
			String fileName = "rootseq" + Randomizer.getSeed() + ".txt";
			PrintStream out = new PrintStream(new File(fileName));
			for (int i = 0; i < seq.length; i++) {
				out.print(seq[i]);
				out.print("\t");
			}
			out.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
        traverse(root, seq, category);
	}
	
    /**
     * recursively walk through the tree top down, and add sequence to alignment whenever
     * a leave node is reached.
     *
     * @param node           reference to the current node, for which we visit all children
     * @param parentSequence randomly generated sequence of the parent node
     * @param category       array of categories for each of the sites
     * @param alignment
     */
    void traverse(Node node, int[] parentSequence, int[] category) {
        for (int childIndex = 0; childIndex < 2; childIndex++) {
            Node child = (childIndex == 0 ? node.getLeft() : node.getRight());
            for (int i = 0; i < categoryCount; i++) {
                getTransitionProbabilities(child, i, probabilities[i]);
            }

            int[] seq = new int[sequenceLength];
            double[] cProb = new double[stateCount];
            for (int i = 0; i < sequenceLength; i++) {
                System.arraycopy(probabilities[category[i]], parentSequence[i] * stateCount, cProb, 0, stateCount);
                seq[i] = Randomizer.randomChoicePDF(cProb);
            }

            if (child.isLeaf()) {
                sequenceInput.setValue(intArray2Sequence(seq, child), this);
            } else {
                traverse(child, seq, category);
            }
        }
    } // traverse

    /**
     * Convert integer representation of sequence into a Sequence
     *
     * @param seq  integer representation of the sequence
     * @param node used to determine taxon for sequence
     * @return Sequence
     */
    Sequence intArray2Sequence(int[] seq, Node node) {
        String seqString = m_dataType.encodingToString(seq);
        
        return new Sequence(node.getID(), seqString);
    } // intArray2Sequence

    /**
     * get transition probability matrix for particular rate category *
     */
    void getTransitionProbabilities(Node node, int rateCategory, double[] probs) {

        Node parent = node.getParent();
        double branchRate = (branchRateModel == null ? 1.0 : branchRateModel.getRateForBranch(node));
        branchRate *= siteModel.getRateForCategory(rateCategory, node);
        siteModel.getSubstitutionModel().getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), branchRate, probs);
        System.arraycopy(siteModel.getSubstitutionModel().getFrequencies(),
        		0, probabilities[rateCategory], stateCount * stateCount, stateCount);
    } // getTransitionProbabilities
	
}
