package potts.tools;

import java.awt.Color;
import java.io.File;
import java.io.PrintStream;
import java.text.DecimalFormat;
import java.util.Arrays;

import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.inference.Runnable;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import potts.dca.DCAFromFile;

public class ScoreDCA extends Runnable {
	final public Input<File> dcaFileInput = new Input<>("fileName", "json file name containing DCA data", new File("[[none]]"));
	final public Input<OutFile> outputInput = new Input<>("out", "text file name to print to, stdout if not specified", new OutFile("[[none]]"));
	final public Input<OutFile> svgInput = new Input<>("svg", "svg file name to print heatmap to, ignored if not specified", new OutFile("[[none]]"));

	@Override
	public void initAndValidate() {
	}

	@Override
	public void run() throws Exception {

		DCAFromFile dca = new DCAFromFile();
		dca.initByName("fileName", dcaFileInput.get());
		int len = dca.getSiteCount();

        PrintStream out = System.out;
        if (!outputInput.get().getName().equals("[[none]]")) {
        	Log.info("Writing to file " + outputInput.get().getPath());
        	out = new PrintStream(outputInput.get());
        }
        Log.info("\nInferred Coupling Scores (Frobenius Norm):");
        // We expect high score between 0 and 1
        for(int i=0; i<len; i++) {
            for(int j=i+1; j<len; j++) {
                double score = dca.getContactScore(i, j);
                out.printf("Site %d - %d : %.4f\n", (i+1), (j+1), score);
            }
        }
        
        Log.info("\nInferred Coupling Scores (Frobenius Norm) PCA corrected:");
        
        double [][] scores = dca.getContactScores();
        for (int i = 0; i < scores.length; i++) {
        	out.println(Arrays.toString(scores[i]));
        }

        if (!outputInput.get().getName().equals("[[none]]")) {
        	out.close();
        }

        if (!svgInput.get().getName().equals("[[none]]")) {
        	Log.info("Writing to file " + svgInput.get().getPath());
        	out = new PrintStream(svgInput.get());
        	out.print("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
	        out.print("<svg viewBox=\"0 0 " + (100+scores.length * 10) + " " + (50+scores.length * 10) + "\" "
	        		+ "xmlns=\"http://www.w3.org/2000/svg\">\n");
	        out.println("<g transform=\"translate(50,50)\">\n");
	        double min = dca.getContactScore(0, 1);
	        double max = min;
	        for (int i = 0; i < scores.length; i++) {
	            for (int j = i + 1; j < scores.length; j++) {
	                double score = dca.getContactScore(i, j);
	            	min = Math.min(min,  score);
	            	max = Math.max(max,  score);
	            }
	        }
	        for (int i = 0; i < scores.length; i++) {
	            for (int j = i + 1; j < scores.length; j++) {
	                double score = (dca.getContactScore(i, j) - min) / (max - min);
	            	int c = (int) ((score - min) * 255/ (max - min));
	            	//c = Color.HSBtoRGB((float)((scores[i][j] - min)/ (max - min)), 0.7f, 0.9f);
	            	c = Color.HSBtoRGB(0.1f+(float)(0.9*score), 0.7f, 0.9f);
	            	String s = Integer.toHexString(c).substring(2);
	            	// String s = (c < 16 ? "0" : "") + Integer.toHexString(c);
	            	out.println(" <rect "
//	            			+ "style=\"fill:#00" + s + s+ ";\" "
							+ "style=\"fill:#" + s+ ";\" "
	            			+ "width=\"10\" "
	            			+ "height=\"10\" "
	            			+ "x=\"" + i * 10 + "\" "
	            			+ "y=\"" + j * 10 + "\" />");
	            	out.println(" <rect "
//	            			+ "style=\"fill:#00" + s + s+ ";\" "
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
            	// String s = (c < 16 ? "0" : "") + Integer.toHexString(c);
            	out.println(" <rect "
//            			+ "style=\"fill:#00" + s + s+ ";\" "
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
		new Application(new ScoreDCA(), "Score DCA", args);

	}

}
