package potts.tools;

import java.io.File;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.evolution.alignment.Alignment;
import beast.base.inference.Runnable;
import beastfx.app.tools.Application;
import potts.dca.AnnealedImportanceSampling;
import potts.dca.DCAFromFile;
import potts.dca.TrainDCA;

@Description("Score an alignment for a given DCA by estimating normalisation constant first (by annealed importance sampling)")
public class ScoreAlignment extends Runnable {
	final public Input<File> dcaFileInput = new Input<>("dca", "json file name containing DCA data", new File("[[none]]"));
	final public Input<File> inputInput = new Input<>("data", "sequence file (nexus or fasta) containing sequence alignment", new File("[[none]]"));

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		DCAFromFile dca = new DCAFromFile();
		dca.initByName("fileName", dcaFileInput.get());

		double [][] h = dca.getH();
		double [][][][] J = dca.getJ();
		
		AnnealedImportanceSampling ais = new AnnealedImportanceSampling(h, J, dca.getStateCount());
		
        // Run Estimation
        // 10,000 steps, average over 10 runs
        double logZ = ais.estimateLogZ(10000, 10);

		Alignment data = TrainDCA.getAlignment(inputInput.get());
		int stateCount = data.getMaxStateCount();
		int siteCount = data.getSiteCount();
		int seqCount = data.getTaxonCount();

        // convert BEAST alignment to int matrix
        int [][] matrix = new int[seqCount][siteCount];
        for(int i = 0; i < siteCount; i++) {
        	int [] pattern = data.getPattern(data.getPatternIndex(i));
        	for (int m = 0; m < seqCount; m++) {
        		matrix[m][i] = pattern[m] >= 0 && pattern[m] < stateCount ? pattern[m] : stateCount;
            }
        }

        double score = 0;
        double logProb = 0;

        for (int k = 0; k < seqCount; k++) {
	        int [] seq = matrix[k];
	        for(int i=0; i<siteCount; i++) {
	        	score += h[i][seq[i]];
	        }
	        for(int i=0; i<siteCount; i++)
	           for(int j=i+1; j<siteCount; j++)
	              score += J[i][j][seq[i]][seq[j]];
	        
	        logProb += score - logZ;
        }

        System.out.println("logZ = " + logZ);
        System.out.println("logProb = " + logProb);
	}

	public static void main(String[] args) throws Exception {
		new Application(new ScoreAlignment(), "ScoreAlignment", args);
	}

}