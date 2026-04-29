package potts.tools;

import java.awt.Color;
import java.io.File;
import java.io.PrintStream;
import java.text.DecimalFormat;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.inference.Runnable;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import potts.dca.TrainDCA;

@Description("Create heatmap from sequence alignment based on mutual information between sites in the alignment")
public class MutualInformation extends Runnable {
	final public Input<File> inputInput = new Input<>("in", "sequence file (nexus or fasta) containing sequence alignment", Validate.REQUIRED);
	final public Input<OutFile> svgInput = new Input<>("svg", "svg file name to print heatmap to, ignored if not specified", new OutFile("[[none]]"));

	@Override
	public void initAndValidate() {
	}

    /**
     * Calculates the Mutual Information between two discrete arrays.
     * 
     * @param x The first array of discrete values.
     * @param y The second array of discrete values.
     * @param n The size of the state space (number of possible discrete values).
     * @return The mutual information in bits.
     */
    public static double calculateMI(int[] x, int[] y, int n) {
        if (x.length != y.length) {
            throw new IllegalArgumentException("Arrays must have the same length.");
        }

        int totalSamples = x.length;

        // 1. Calculate joint and marginal counts
        // jointCounts[i][j] stores how many times (x=i, y=j) occurs
        long[][] jointCounts = new long[n][n];
        long[] countX = new long[n];
        long[] countY = new long[n];

        for (int i = 0; i < totalSamples; i++) {
            int valX = x[i];
            int valY = y[i];
            
            jointCounts[valX][valY]++;
            countX[valX]++;
            countY[valY]++;
        }

        // 2. Compute Mutual Information
        double mi = 0.0;
        double log2 = Math.log(2.0);

        for (int i = 0; i < n; i++) {
            if (countX[i] == 0) continue; // Skip if x never takes value i
            
            for (int j = 0; j < n; j++) {
                if (jointCounts[i][j] > 0) {
                    // P(x, y)
                    double pXY = (double) jointCounts[i][j] / totalSamples;
                    // P(x)
                    double pX = (double) countX[i] / totalSamples;
                    // P(y)
                    double pY = (double) countY[j] / totalSamples;

                    // MI sum: P(x,y) * log2( P(x,y) / (P(x)*P(y)) )
                    mi += pXY * (Math.log(pXY / (pX * pY)) / log2);
                }
            }
        }

        return mi;
    }

	@Override
	public void run() throws Exception {
        Log.info("Loading file " + inputInput.get().getPath());
		Alignment data = TrainDCA.getAlignment(inputInput.get());
		int stateCount = data.getMaxStateCount();
		int siteCount = data.getSiteCount();
		int seqCount = data.getTaxonCount();
        
        // convert BEAST alignment to int matrix
        int [][] matrix = new int[siteCount][seqCount];
        for(int i = 0; i < siteCount; i++) {
        	int [] pattern = data.getPattern(data.getPatternIndex(i));
        	for (int m = 0; m < seqCount; m++) {
        		matrix[i][m] = pattern[m] >= 0 && pattern[m] < stateCount ? pattern[m] : stateCount; 
            }
        }
        
        double [][] scores = new double[siteCount][siteCount];
        for(int i = 0; i < siteCount; i++) {
        	int [] seq1 =  matrix[i];
        	for (int j = 0; j < siteCount; j++) {
            	int [] seq2 =  matrix[j];
            	scores[i][j] = calculateMI(seq1, seq2, stateCount + 1);
        	}
        }

        if (!svgInput.get().getName().equals("[[none]]")) {
        	Log.info("Writing to file " + svgInput.get().getPath());
        	PrintStream out = new PrintStream(svgInput.get());
        	out.print("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
	        out.print("<svg viewBox=\"0 0 " + (100+scores.length * 10) + " " + (50+scores.length * 10) + "\" "
	        		+ "xmlns=\"http://www.w3.org/2000/svg\">\n");
	        out.println("<g transform=\"translate(50,50)\">\n");
	        double min = scores[0][1];
	        double max = min;
	        for (int i = 0; i < scores.length; i++) {
	            for (int j = i + 1; j < scores.length; j++) {
	                double score = scores[i][j];
	            	min = Math.min(min,  score);
	            	max = Math.max(max,  score);
	            }
	        }
	        for (int i = 0; i < scores.length; i++) {
	            for (int j = i + 1; j < scores.length; j++) {
	                double score = (scores[i][j] - min) / (max - min);
	            	int c = (int) ((score - min) * 255/ (max - min));
	            	c = Color.HSBtoRGB(0.1f+(float)(0.9*score), 0.7f, 0.9f);
	            	String s = Integer.toHexString(c).substring(2);
	            	out.println(" <rect "
							+ "style=\"fill:#" + s+ ";\" "
	            			+ "width=\"10\" "
	            			+ "height=\"10\" "
	            			+ "x=\"" + i * 10 + "\" "
	            			+ "y=\"" + j * 10 + "\" />");
	            	out.println(" <rect "
							+ "style=\"fill:#" + s+ ";\" "
	            			+ "width=\"10\" "
	            			+ "height=\"10\" "
	            			+ "x=\"" + j * 10 + "\" "
	            			+ "y=\"" + i * 10 + "\" />");
	            }
	        }
	        out.println("</g>\n");

            for (int j = 0; j < scores.length; j++) {
            	int c = Color.HSBtoRGB(0.1f+((float)(0.9*j) / scores.length), 0.7f, 0.9f);
            	String s = Integer.toHexString(c).substring(2);
            	out.println(" <rect "
						+ "style=\"fill:#" + s+ ";\" "
            			+ "width=\"10\" "
            			+ "height=\"10\" "
            			+ "x=\"" + (50+j * 10) + "\" "
            			+ "y=\"0\" />");
            }
            DecimalFormat df = new DecimalFormat("#.###");
	        out.println("<text x=\"10\" y=\"10\">" + df.format(min) + "</text>");
	        out.println("<text x=\""+(scores.length*10+52)+"\" y=\"10\">" + df.format(max) + "</text>");
	        
	        for (int i = 0; i < scores.length; i += 10) {
            	out.println("<rect x=\"" + (50+i * 10) + "\" y=\"40\" width=\"1\" height=\"10\"/>");	        	
            	out.println("<text x=\"" + (50+i * 10) + "\" y=\"40\">" + i + "</text>");	        	
            	out.println("<rect y=\"" + (50+i * 10) + "\" x=\"40\" width=\"10\" height=\"1\"/>");	        	
            	out.println("<text y=\"" + (50+i * 10) + "\" x=\"30\">" + i + "</text>");	        	
	        }
	        out.print("</svg>\n");
	        out.close();
        }
        
        Log.warning("Done");
	}

	public static void main(String[] args) throws Exception {
		new Application(new MutualInformation(), "MutualInformation", args);

	}

}