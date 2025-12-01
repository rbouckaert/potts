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

@Description("Simulate alignment on a tree using the Potts model and WAG model")
public class DCASequenceSimulator2 extends Runnable {
	final public Input<File> dcaInput = new Input<>("dca", "DCA json file previously trained on an alignment", Validate.REQUIRED);
	final public Input<File> treeFileInput = new Input<>("treeFile", "NEXUS tree file to simulate alignment for -- if there are "
			+ "multiple trees, output files will be numbered."
			+ "Branch lengths are used as evolutionary distance (i.e. clock rate = 1 and strict clock assumed)", Validate.REQUIRED);
	final public Input<OutFile> outputInput = new Input<>("out", "file to save fasta file", Validate.REQUIRED);

	final public Input<Integer> rootStepCountInput = new Input<>("rootStepCount", "number of times each site for the root sequence is resampled (at initialisation only)", 1000);
	final public Input<Integer> stepCountInput = new Input<>("stepCount", "number of times each site is resampled when sampling sequences conditioned on parent or parents and children", 100);
	final public Input<Double> temperatureInput = new Input<>("temperature", "temperature balancing the effect of Potts model Pt and substitution model Ps. "
			+ "Mutations are chosen proportional to Pt^1/t * Ps (thus 1/t", 1.0);
	
	int stepCount;
	int rootStepCount;
	double temperatureFactor;

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		long start = System.currentTimeMillis();

		stepCount = stepCountInput.get();
		rootStepCount = rootStepCountInput.get();
		temperatureFactor = 1.0 / temperatureInput.get();
		
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
				
				DCASequenceSimulator.toFasta(alignment, path, tree, dataType);
			}
		} else {
			Tree tree =  trees.get(0);
			int [][] alignment = simulateAlignment(dca, tree, model);
			DCASequenceSimulator.toFasta(alignment, outputInput.get().getPath(), tree, dataType);
		}
		
		long end = System.currentTimeMillis();
		Log.warning("Done in " + (end-start) + " ms");

	}
	
	private int[][] simulateAlignment(DCA dca, Tree tree, SubstitutionModel model) {
		double [] rateMatrix = model.getRateMatrix(null);
		
		// reserve memory for sequences at all internal and external nodes and randomly initialise
		int [][] alignment = new int[tree.getNodeCount()][dca.siteCount];
		for (int i = 0; i < alignment.length; i++) {
			int [] seq = alignment[i];
			for (int j = 0; j < seq.length; j++) {
				seq[j] = Randomizer.nextInt(dca.stateCount-1);
			}
		}
		
		// sample root sequence
		DCASequenceSimulator.sampleRootSequence(alignment[tree.getRoot().getNr()], dca, rootStepCount);
		traverseDown(alignment, dca, tree.getRoot(), rateMatrix);
		
		return alignment;
	}

	private void traverseDown(int[][] alignment, DCA dca, Node node, double [] rateMatrix) {
		if (!node.isRoot()) {
			int [] parentseq = alignment[node.getParent().getNr()];
			int [] seq = alignment[node.getNr()];
			System.arraycopy(parentseq, 0, seq, 0, seq.length);
			int mutationCount = (int)(dca.siteCount * node.getLength()) + 1;
			
			
			double [] mutationProb = new double[dca.siteCount * dca.stateCount];
			// introduce mutationCount mutations
			for (int step = 0; step < mutationCount; step++) {
				// calc probability of every possible mutation
	            for (int i = 0; i < dca.siteCount; i++) {
	                int oldState = seq[i];
	                for (int newState = 0; newState < dca.stateCount; newState++) {
	                	if (newState != oldState) {
	                		mutationProb[i * dca.stateCount + newState] = 
	                			Math.log(rateMatrix[oldState * (dca.stateCount-1) + newState]) + 
	                			temperatureFactor * DCASequenceSimulator.computeDeltaHamiltonian(dca, seq, i, oldState, newState);
	                	} else {
	                		mutationProb[i * dca.stateCount + newState] = Double.NEGATIVE_INFINITY;
	                	}
	                };
	            }
	            
	            // find max prob
	            double max = mutationProb[1];
	            for (double d : mutationProb) {
	            	max = Math.max(d, max);
	            }
	            
	            // transform from log space
	            for (int i = 0; i < mutationProb.length; i++) {
	            	mutationProb[i] = Math.exp(mutationProb[i] - max);
	            }
	            
	            // pick a mutation proportional to mutationProb
	            int r = Randomizer.randomChoicePDF(mutationProb);
	            int site = r / dca.stateCount;
	            int newState = r % dca.stateCount;
	            seq[site] = newState;
	        }
		}
	}
	
	
	public static void main(String[] args) throws Exception {
		new Application(new DCASequenceSimulator2(), "DCASequenceSimulator2", args);
	}
}
