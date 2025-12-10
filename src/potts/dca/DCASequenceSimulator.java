package potts.dca;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.PrintStream;
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

@Description("Simulate alignment on a tree using the Potts model and WAG model")
public class DCASequenceSimulator extends Runnable {
	final public Input<File> dcaInput = new Input<>("dca", "DCA json file previously trained on an alignment", Validate.REQUIRED);
	final public Input<File> treeFileInput = new Input<>("treeFile", "NEXUS tree file to simulate alignment for -- if there are "
			+ "multiple trees, output files will be numbered."
			+ "Branch lengths are used as evolutionary distance (i.e. clock rate = 1 and strict clock assumed)", Validate.REQUIRED);
	final public Input<OutFile> outputInput = new Input<>("out", "file to save fasta file", Validate.REQUIRED);

	final public Input<Integer> rootStepCountInput = new Input<>("rootStepCount", "number of times each site for the root sequence is resampled (at initialisation only)", 1000);
	final public Input<Integer> stepCountInput = new Input<>("stepCount", "number of times each site is resampled when sampling sequences conditioned on parent or parents and children", 100);
	final public Input<Integer> resampleCountInput = new Input<>("resampleCount", "number of times each all sites in the alignment are resampled", 100);
	
	int stepCount;
	int rootStepCount;
	int resampleCount;

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		long start = System.currentTimeMillis();

		stepCount = stepCountInput.get();
		rootStepCount = rootStepCountInput.get();
		resampleCount = resampleCountInput.get();
		
		DCA dca = new DCA();
		String json = BeautiDoc.load(dcaInput.get());
		dca.fromJSON(new JSONObject(json));
		
		NexusParser parser = new NexusParser();
		parser.parseFile(treeFileInput.get());
		List<Tree> trees = parser.trees;

		
		DataType dataType = dca.stateCount == 21 ? new Aminoacid() : new Nucleotide();
		SubstitutionModel model = dca.stateCount == 21 ? new WAG() : new JukesCantor();
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
				
				toFasta(alignment, path, tree, dataType);
			}
		} else {
			Tree tree =  trees.get(0);
			int [][] alignment = simulateAlignment(dca, tree, model);
			toFasta(alignment, outputInput.get().getPath(), tree, dataType);
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
				seq[j] = Randomizer.nextInt(dca.stateCount-1);
			}
		}
		
		// sample root sequence
		sampleRootSequence(alignment[tree.getRoot().getNr()], dca, rootStepCount, dca.stateCount);
		traverseDown(alignment, dca, tree.getRoot(), matrices);
		for (int i = 0; i < resampleCount; i++) {
			reSample(alignment, dca, tree.getRoot(), matrices);
		}
		
		return alignment;
	}

	private void reSample(int[][] alignment, DCA dca, Node node, double [][] matrices) {
		if (!node.isRoot()) {
			if (node.isLeaf()) {
				traverseDown(alignment, dca, node, matrices);
			} else {
				int [] seq = alignment[node.getNr()];
				double [] matrix = matrices[node.getNr()];
				int [] parentseq = alignment[node.getParent().getNr()];
				int [] childseq1 = alignment[node.getLeft().getNr()];
				double [] matrix1 = matrices[node.getLeft().getNr()];
				int [] childseq2 = alignment[node.getRight().getNr()];
				double [] matrix2 = matrices[node.getRight().getNr()];
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
		                		parentseq[i], matrix,
		                		childseq1[i], matrix1,
		                		childseq2[i], matrix2);
		                
		                // Metropolis Criterion
		                if (deltaH >= 0 || Randomizer.nextDouble() < Math.exp(deltaH)) {
		                    seq[i] = newState; // Accept mutation
		                }
		            }
		        }
			}
		}
		if (!node.isLeaf()) {
			for (Node child : node.getChildren()) {
				reSample(alignment, dca, child, matrices);
			}
		}
	}

	private void traverseDown(int[][] alignment, DCA dca, Node node, double [][] matrices) {
		if (!node.isRoot()) {
			int [] seq = alignment[node.getNr()];
			double [] matrix = matrices[node.getNr()];
			int [] parentseq = alignment[node.getParent().getNr()];
			for (int step = 0; step < stepCount; step++) {
	            // Try to mutate every site once (Standard sweep)
	            for (int i = 0; i < dca.siteCount; i++) {
	                int oldState = seq[i];
	                int newState = Randomizer.nextInt(dca.stateCount);
	                
	                if (oldState == newState) continue;
	
	                // Calculate change in Hamiltonian (Score)
	                // Delta H = H(new) - H(old)
	                // We accept if Delta H > 0 (more probable) or with prob exp(Delta H)
	                
	                double deltaH = computeDeltaHamiltonian(dca, seq, i, oldState, newState, parentseq[i], matrix);
	                
	                // Metropolis Criterion
	                if (deltaH >= 0 || Randomizer.nextDouble() < Math.exp(deltaH)) {
	                    seq[i] = newState; // Accept mutation
	                }
	            }
	        }
		}
			
		if (!node.isLeaf()) {
			for (Node child : node.getChildren()) {
				traverseDown(alignment, dca, child, matrices);
			}
		}
	}
	
	
	static void sampleRootSequence(int[] seq, DCA dca, int rootStepCount, int stateCount) {
		for (int step = 0; step < rootStepCount; step++) {
            // Try to mutate every site once (Standard sweep)
            for (int i = 0; i < dca.siteCount; i++) {
                int oldState = seq[i];
                int newState = Randomizer.nextInt(stateCount);
                
                if (oldState == newState) continue;

                // Calculate change in Hamiltonian (Score)
                // Delta H = H(new) - H(old)
                // We accept if Delta H > 0 (more probable) or with prob exp(Delta H)
                
                double deltaH = computeDeltaHamiltonian(dca, seq, i, oldState, newState);
                
                // Metropolis Criterion
                if (deltaH >= 0 || Randomizer.nextDouble() < Math.exp(deltaH)) {
                    seq[i] = newState; // Accept mutation
                }
            }
        }
	}
	
	
    /**
     * Efficiently calculates the change in energy for a single mutation on a single sequence 
	 * No parents or children.
     * deltaH = (h_new - h_old) + Sum_neighbors(J_new_neighbor - J_old_neighbor)
     */
    public static double computeDeltaHamiltonian(DCA dca, int[] seq, int i, int oldState, int newState) {
        // 1. Change in Field Energy
        double delta = dca.h[i][newState] - dca.h[i][oldState];

        // 2. Change in Coupling Energy with all other sites j
        for (int j = 0; j < dca.siteCount; j++) {
            if (i == j) continue;
            
            int neighborState = seq[j];
            
            // Access J respecting i < j storage convention
            if (i < j) {
                delta += dca.J[i][j][newState][neighborState] - dca.J[i][j][oldState][neighborState];
            } else {
                // if i > j, we look up J[j][i][neighbour][target]
                delta += dca.J[j][i][neighborState][newState] - dca.J[j][i][neighborState][oldState];
            }
        }
        return delta;
    }

    
    /**
     * Efficiently calculates the change in energy for a single mutation on a single sequence
     * taking parent state in account
     * deltaH = (h_new - h_old) + Sum_neighbors(J_new_neighbor - J_old_neighbor)
     */
    private double computeDeltaHamiltonian(DCA dca, int[] seq, int i, int oldState, int newState, int parentState, double [] matrix) {
        double delta = computeDeltaHamiltonian(dca, seq, i, oldState, newState); 

        // Change in transition probability
        delta += Math.log(matrix[parentState * (dca.stateCount-1) + newState]) - Math.log(matrix[parentState * (dca.stateCount-1) + oldState]);

        return delta;
    }    

    
    /**
     * Efficiently calculates the change in energy for a single mutation on a single sequence
     * taking parent and child states in account
     * deltaH = (h_new - h_old) + Sum_neighbors(J_new_neighbor - J_old_neighbor)
     */
    private double computeDeltaHamiltonian(DCA dca, int[] seq, int i, int oldState, int newState, 
    		int parentState, double [] matrix,
    		int childState1, double [] matrix1,
    		int childState2, double [] matrix2
    		) {
        double delta = computeDeltaHamiltonian(dca, seq, i, oldState, newState); 

        // Change in transition probability
        delta += Math.log(matrix[parentState * (dca.stateCount-1) + newState]) - Math.log(matrix[parentState * (dca.stateCount-1) + oldState])
        		+Math.log(matrix1[newState * (dca.stateCount-1) + childState1]) - Math.log(matrix1[oldState * (dca.stateCount-1) + childState1])
        		+Math.log(matrix2[newState * (dca.stateCount-1) + childState2]) - Math.log(matrix2[oldState * (dca.stateCount-1) + childState2]);
        		

        return delta;
    }    
    
	static public void toFasta(int[][] alignment, String path, Tree tree, DataType dataType) throws FileNotFoundException {
		toFasta(alignment, path, tree, dataType, false);
	}
	
	static public void toFasta(int[][] alignment, String path, Tree tree, DataType dataType, boolean append) throws FileNotFoundException {
		Log.warning("Output written to " + path);
		PrintStream out = new PrintStream(new FileOutputStream(path, append));
		Node [] nodes = tree.getNodesAsArray();
		for (int i = 0; i < tree.getLeafNodeCount(); i++) {
			out.print(">");
			out.println(nodes[i].getID());
			int [] seq = alignment[i];
			for (int j = 0; j < seq.length; j++) {
				out.print(dataType.getCharacter(seq[j]));
			}
			out.println();
		}
		out.close();
	}

	public static void main(String[] args) throws Exception {
		new Application(new DCASequenceSimulator(), "DCASequenceSimulator", args);
	}
}
