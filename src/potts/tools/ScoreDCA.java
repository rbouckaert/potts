package potts.tools;

import java.io.File;
import java.util.Arrays;

import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.Runnable;
import beastfx.app.tools.Application;
import potts.dca.DCAFromFile;

public class ScoreDCA extends Runnable {
	final public Input<File> dcaFileInput = new Input<>("fileName", "json file name containing DCA data", new File("[[none]]"));

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		
		DCAFromFile dca = new DCAFromFile();
		dca.initByName("fileName", dcaFileInput.get());
		int len = dca.getSiteCount();
		
        Log.info("\nInferred Coupling Scores (Frobenius Norm):");
        // We expect high score between 0 and 1
        for(int i=0; i<len; i++) {
            for(int j=i+1; j<len; j++) {
                double score = dca.getContactScore(i, j);
                System.out.printf("Site %d - %d : %.4f\n", (i+1), (j+1), score);
            }
        }
        
        Log.info("\nInferred Coupling Scores (Frobenius Norm) PCA corrected:");
        double [][] scores = dca.getContactScores();
        for (int i = 0; i < scores.length; i++) {
        	System.out.println(Arrays.toString(scores[i]));
        }

	}

	public static void main(String[] args) throws Exception {
		new Application(new ScoreDCA(), "Score DCA", args);

	}

}
