package potts.dca;


import java.io.File;
import java.io.PrintStream;
import java.util.List;
import java.util.Set;

import beast.base.core.BEASTInterface;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.inference.Runnable;
import beast.pkgmgmt.BEASTClassLoader;
import beastfx.app.inputeditor.AlignmentImporter;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import beastfx.app.util.Utils;

@Description("Train direct-coupling analysis (DCA) and save Potts parameters into file")
public class TrainDCA extends Runnable {
	final public Input<File> inputInput = new Input<>("in", "sequence file (nexus or fasta) containing sequence alignment", Validate.REQUIRED);
	final public Input<OutFile> outputInput = new Input<>("out", "file to save DCA parameters", Validate.REQUIRED);
	final public Input<Double> learningRateInput = new Input<>("learningRate", "Step size for gradient ascent (e.g. 0.01)", 0.01);
	final public Input<Double> regularisationInput = new Input<>("regularisation", "L2 penalty weight (e.g. 0.001)", 0.001);
	final public Input<Integer> stepsCountInput = new Input<>("epochs", "Number of gradient updates to perform", 10000);
	final public Input<Integer> mcStepsCountInput = new Input<>("mcSteps", "Number of Metropolis sweeps per epoch", 10);

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {
		// load alignment
        Log.info("Loading file " + inputInput.get().getPath());
		Alignment data = getAlignment(inputInput.get());
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
        
        // train using Boltzmann DCA
        BoltzmannDCA dca = new BoltzmannDCA(matrix, stateCount+1, learningRateInput.get(), regularisationInput.get());
        dca.train(stepsCountInput.get(), mcStepsCountInput.get()); 

        
        // output results
        Log.info("Output written to " + outputInput.get().getPath());
        PrintStream out = new PrintStream(outputInput.get().getPath());
        out.print(dca.toJSON());
        out.close();
        Log.warning("Done");
	}

	
	private static Alignment getAlignment(File file) {
        Set<String> importerClasses = Utils.loadService(AlignmentImporter.class);        
        for (String _class: importerClasses) {
        	try {
				AlignmentImporter importer = (AlignmentImporter) BEASTClassLoader.forName(_class).newInstance();
				if (importer.canHandleFile(file)) {
					List<BEASTInterface> os = importer.loadFile(file);
					Alignment a = (Alignment) os.get(0);
					return a;
				}
        	} catch (InstantiationException | IllegalAccessException | ClassNotFoundException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
        }

		throw new IllegalArgumentException("Could not find an importer that can handle the file " + file.getPath());
	}

	public static void main(String[] args) throws Exception {
		new Application(new TrainDCA(), "TrainDCA", args);
	}

}
