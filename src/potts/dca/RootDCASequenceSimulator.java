package potts.dca;

import java.io.File;
import java.util.List;

import org.json.JSONObject;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;
import beast.base.evolution.datatype.Nucleotide;
import beast.base.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.substitutionmodel.WAG;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Runnable;
import beast.base.parser.NexusParser;
import beast.base.util.Randomizer;
import beastfx.app.inputeditor.BeautiDoc;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;

@Description("Simulate alignment on a tree using the Potts model for root and tip sequences and WAG model in between")
public class RootDCASequenceSimulator extends Runnable {
	final public Input<File> dcaInput = new Input<>("dca", "DCA json file previously trained on an alignment", Validate.REQUIRED);
	final public Input<File> treeFileInput = new Input<>("treeFile", "NEXUS tree file to simulate alignment for -- if there are "
			+ "multiple trees, output files will be numbered."
			+ "Branch lengths are used as evolutionary distance (i.e. clock rate = 1 and strict clock assumed)", Validate.REQUIRED);
	final public Input<OutFile> outputInput = new Input<>("out", "file to save fasta file", Validate.REQUIRED);

	final public Input<Integer> rootStepCountInput = new Input<>("rootStepCount", "number of times each site for the root sequence is resampled (at initialisation only)", 1000);
	final public Input<Integer> stepCountInput = new Input<>("stepCount", "number of times each site is resampled when sampling sequences conditioned on parent or parents and children", 100);
	final public Input<Integer> resampleCountInput = new Input<>("resampleCount", "number of times each all sites in the alignment are resampled", 1000);

	final public Input<Integer> repeatsInput = new Input<>("repeats", "number of times the sequence is sampled for a taxon."
			+ "The total sequence length equals repeast times sequence lenght in DCA", 1);

	final public Input<Long> seedInput = new Input<>("seed", "random number seed, if not specified, take default seed (time dependent)");

	final public Input<Double> temperatureInput = new Input<>("temperature", "temperature balancing the effect of Potts model Pt and substitution model Ps. "
			+ "Mutations are chosen proportional to Pt^1/t * Ps (thus the Potss model weighted by 1/t). "
			+ "If negative, Potts model is ignored.", 1.0);

	int stepCount;
	int rootStepCount;
	int resampleCount;

	int dataStateCount;
	
	double temperatureFactor;

    final public static String BAR = "|---------|---------|---------|---------|---------|---------|---------|---------|";
	
	
	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		long start = System.currentTimeMillis();
		
		if (seedInput.get() != null ) {
			Randomizer.setSeed(seedInput.get());
		}
		temperatureFactor = 1.0 / temperatureInput.get();

		stepCount = stepCountInput.get();
		rootStepCount = rootStepCountInput.get();
		resampleCount = resampleCountInput.get();
		
		DCA dca = new DCA();
		String json = BeautiDoc.load(dcaInput.get());
		dca.fromJSON(new JSONObject(json));
		
		dataStateCount = dca.stateCount - 1;
		
		NexusParser parser = new NexusParser();
		parser.parseFile(treeFileInput.get());
		List<Tree> trees = parser.trees;

		
		DataType dataType = dataStateCount == 20 ? new Aminoacid() : new Nucleotide();
		SubstitutionModel model = dataStateCount == 20 ? new WAG() : new JukesCantor();
		((BEASTInterface)model).initAndValidate();
		
		if (trees.size() > 1) {
			for (int i = 0; i < trees.size(); i++) {
				Tree tree =  trees.get(i);
				int [][] alignment = simulateAlignment(dca, tree, model);
				
				String path = outputInput.get().getPath();
				if (path.lastIndexOf('.') > 0) {
					int k = path.lastIndexOf('.');
					path = path.substring(0, k) + (i<100?"0":"") + (i<10?"0":"") + i + path.substring(k+1);
				}
				
				DCASequenceSimulator.toFasta(alignment, path, tree, dataType, false);
				for (int k = 1; k < repeatsInput.get(); k++) {
					alignment = simulateAlignment(dca, tree, model);
					DCASequenceSimulator.toFasta(alignment, path, tree, dataType, true);
				}
			}
		} else {
			Tree tree =  trees.get(0);
			int [][] alignment = simulateAlignment(dca, tree, model);
			DCASequenceSimulator.toFasta(alignment,  outputInput.get().getPath(), tree, dataType, false);
			for (int k = 1; k < repeatsInput.get(); k++) {
				alignment = simulateAlignment(dca, tree, model);
				DCASequenceSimulator.toFasta(alignment,  outputInput.get().getPath(), tree, dataType, true);
			}
		}
		
		long end = System.currentTimeMillis();
		Log.warning("Done in " + (end-start) + " ms");

	}
	
	private int[][] simulateAlignment(DCA dca, Tree tree, SubstitutionModel model) {
		// initialise transition probability matrices
		double [][] matrices = new double[tree.getNodeCount()][dca.stateCount * dca.stateCount];
		Node [] nodes = tree.getNodesAsArray();
		for (int i = 0; i < tree.getNodeCount() - 1; i++) {
			model.getTransitionProbabilities(nodes[i], nodes[i].getParent().getHeight(), nodes[i].getHeight(), 1.0, matrices[i]);
		}
		
		// reserve memory for sequences at all internal and external nodes and randomly initialise
		int [][] alignment = new int[tree.getNodeCount()][dca.siteCount];
		for (int i = 0; i < alignment.length; i++) {
			int [] seq = alignment[i];
			for (int j = 0; j < seq.length; j++) {
				seq[j] = Randomizer.nextInt(dataStateCount);
			}
		}
		
		// sample root sequence
		if (temperatureFactor > 0) {
			DCASequenceSimulator.sampleRootSequence(alignment[tree.getRoot().getNr()], dca, rootStepCount, dataStateCount);
		} else {
			int siteCount = alignment[tree.getRoot().getNr()].length;
			double [] freqs = model.getFrequencies();
			for (int i = 0; i < siteCount; i++) {
				alignment[tree.getRoot().getNr()][i] = Randomizer.randomChoicePDF(freqs);
			}
		}
		

		// sample internal nodes & leaf nodes
		traverseDown(alignment, dca, tree.getRoot(), matrices);
		
		// resample all sequences, so dependencies are taken in account
		int reported = 0;
        int lines = Math.max(1, resampleCount / 80);

        Log.warning.println(BAR);
        for (int i = 0; i < resampleCount; i++) {
			reSample(alignment, dca, tree.getRoot(), matrices);
            if (i > 0 && i % lines == 0 && reported < 81) {
				while (10000 * reported < 810000 * (i + 1)/ resampleCount) {
					if (reported % 10 != 0) {
						Log.warning.print("*");
					} else {
						Log.warning.print("|");
					}
	                reported++;
        	    }
            }
		}
        Log.warning.println();
		
		return alignment;
	}

	private void reSample(int[][] alignment, DCA dca, Node node, double [][] matrices) {
		if (node.isRoot()) {
			resampleRoot(alignment, dca, node, matrices);
		} else {
			if (node.isLeaf()) {
				resampleLeaf(alignment, dca, node, matrices);
			} else {
				resampleInternalNode(alignment, dca, node, matrices);
			}
		}
		if (!node.isLeaf()) {
			for (Node child : node.getChildren()) {
				reSample(alignment, dca, child, matrices);
			}
		}
		if (!node.isRoot() && ! node.isLeaf()) {
			resampleInternalNode(alignment, dca, node, matrices);
		}
	}

	private void resampleInternalNode(int[][] alignment, DCA dca, Node node, double[][] matrices) {
		int [] seq = alignment[node.getNr()];
		double [] matrix = matrices[node.getNr()];
		int [] parentseq = alignment[node.getParent().getNr()];
		int [] childseq1 = alignment[node.getLeft().getNr()];
		double [] matrix1 = matrices[node.getLeft().getNr()];
		int [] childseq2 = alignment[node.getRight().getNr()];
		double [] matrix2 = matrices[node.getRight().getNr()];
        double [] probs = new double[dataStateCount];

		
		for (int step = 0; step < stepCount; step++) {
            // Try to mutate every site once (Standard sweep)
            for (int i = 0; i < dca.siteCount; i++) {
                int partenState = parentseq[i];
                int child1State = childseq1[i];
                int child2State = childseq2[i];
                for (int k = 0; k < dataStateCount; k++) {
                	probs[k] = matrix[partenState * dataStateCount + k] *
                			matrix1[k * dataStateCount + child1State] *
                			matrix2[k * dataStateCount + child2State];
                }
                int newState = Randomizer.randomChoicePDF(probs);
                seq[i] = newState;
            }
        }
	}

	private void resampleRoot(int[][] alignment, DCA dca, Node root, double[][] matrices) {
		if (!root.isRoot()) {
			throw new IllegalArgumentException("Expected root node");
		}
		int [] seq = alignment[root.getNr()];
		int [] childseq1 = alignment[root.getLeft().getNr()];
		double [] matrix1 = matrices[root.getLeft().getNr()];
		int [] childseq2 = alignment[root.getRight().getNr()];
		double [] matrix2 = matrices[root.getRight().getNr()];
		for (int step = 0; step < stepCount; step++) {
            // Try to mutate every site once (Standard sweep)
            for (int i = 0; i < dca.siteCount; i++) {
                int oldState = seq[i];
                int newState = Randomizer.nextInt(dca.stateCount);

                if (oldState == newState) continue;

                // Calculate change in Hamiltonian (Score)
                // Delta H = H(new) - H(old)
                // We accept if Delta H > 0 (more probable) or with prob exp(Delta H)

                double deltaH = computeDeltaHamiltonian(dca, seq, i, oldState, newState,
                		childseq1[i], matrix1, childseq2[i], matrix2);

                // Metropolis Criterion
                if (deltaH >= 0 || Randomizer.nextDouble() < Math.exp(deltaH)) {
                    seq[i] = newState; // Accept mutation
                }
            }
        }		
	}
	
	private void resampleLeaf(int[][] alignment, DCA dca, Node node, double[][] matrices) {
		int [] seq = alignment[node.getNr()];
		int [] parentseq = alignment[node.getParent().getNr()];
		double [] matrix = matrices[node.getNr()];
		for (int step = 0; step < stepCount; step++) {
            // Try to mutate every site once (Standard sweep)
            for (int i = 0; i < dca.siteCount; i++) {
                int oldState = seq[i];
                int newState = Randomizer.nextInt(dataStateCount);

                if (oldState == newState) continue;

                // Calculate change in Hamiltonian (Score)
                // Delta H = H(new) - H(old)
                // We accept if Delta H > 0 (more probable) or with prob exp(Delta H)

                double deltaH = computeDeltaHamiltonian(dca, seq, i, oldState, newState,
                		parentseq[i], matrix);

                // Metropolis Criterion
                if (deltaH >= 0 || Randomizer.nextDouble() < Math.exp(deltaH)) {
                    seq[i] = newState; // Accept mutation
                }
            }
        }
	}

	private void traverseDown(int[][] alignment, DCA dca, Node node, double [][] matrices) {
		if (!node.isRoot()) {
			int [] seq = alignment[node.getNr()];
			double [] matrix = matrices[node.getNr()];
			int [] parentseq = alignment[node.getParent().getNr()];
            double [] probs = new double[dataStateCount];
			//for (int step = 0; step < stepCount; step++) {
	            // Try to mutate every site once (Standard sweep)
            for (int i = 0; i < dca.siteCount; i++) {
                int parentState = parentseq[i];
                System.arraycopy(matrix, parentState * dataStateCount, probs, 0, dataStateCount);
                
                int newState = Randomizer.randomChoicePDF(probs);
                seq[i] = newState; // Accept mutation
            }
	        //}
		}
			
		if (!node.isLeaf()) {
			for (Node child : node.getChildren()) {
				traverseDown(alignment, dca, child, matrices);
			}
		} else {
 			// resampleLeaf(alignment, dca, node, matrices);
		}
	}
	



    /**
     * Efficiently calculates the change in energy for a single mutation on a single sequence
     * taking parent state in account
     * deltaH = (h_new - h_old) + Sum_neighbors(J_new_neighbor - J_old_neighbor)
     */
    private double computeDeltaHamiltonian(DCA dca, int[] seq, int i, int oldState, int newState, int parentState, double [] matrix) {
        double delta =
        		temperatureFactor > 0 ?
        		temperatureFactor * DCASequenceSimulator.computeDeltaHamiltonian(dca, seq, i, oldState, newState) :
        		0;

        // Change in transition probability
        delta += Math.log(matrix[parentState * dataStateCount + newState]) - Math.log(matrix[parentState * dataStateCount + oldState]);

        return delta;
    }


    /**
     * Efficiently calculates the change in energy for a single mutation on a single sequence
     * taking parent and child states in account
     * deltaH = (h_new - h_old) + Sum_neighbors(J_new_neighbor - J_old_neighbor)
     */
    private double computeDeltaHamiltonian(DCA dca, int[] seq, int i, int oldState, int newState,
    		int childState1, double [] matrix1,
    		int childState2, double [] matrix2
    		) {
        double delta =
        		temperatureFactor > 0 ?
        		temperatureFactor * DCASequenceSimulator.computeDeltaHamiltonian(dca, seq, i, oldState, newState):
        		0;

        // Change in transition probability
        delta +=
        		+Math.log(matrix1[newState * dataStateCount + childState1]) - Math.log(matrix1[oldState * dataStateCount + childState1])
        		+Math.log(matrix2[newState * dataStateCount + childState2]) - Math.log(matrix2[oldState * dataStateCount + childState2]);
        		

        return delta;
    }

	public static void main(String[] args) throws Exception {
		new Application(new RootDCASequenceSimulator(), "RootDCASequenceSimulator", args);
	}
}