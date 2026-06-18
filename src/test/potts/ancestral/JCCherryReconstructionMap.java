package test.potts.ancestral;

import java.awt.Color;
import java.io.PrintStream;
import java.text.DecimalFormat;

import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.JukesCantor;
import beast.base.evolution.substitutionmodel.WAG;
import beast.base.evolution.tree.Tree;
import beast.base.inference.Runnable;
import beastfx.app.tools.Application;
import beastfx.app.util.OutFile;
import test.beast.BEASTTestCase;

@Description("Create heatmap from reconstruction at the root of a cherry tree with the same tip states")
public class JCCherryReconstructionMap extends Runnable {
	final public Input<OutFile> svgInput = new Input<>("svg", "svg file name to print heatmap to, ignored if not specified", new OutFile("[[none]]"));
	final public Input<Double> detlaHInput = new Input<>("detlaH", "change in tree height at each step", 0.02);
	final public Input<Integer> nrOfPStepsInput = new Input<>("nrOfPSteps", "number of steps in increasing P (deltap = 1/nrOfPSteps)", 100);
	final public Input<Integer> nrOfHStepsInput = new Input<>("nrOfHSteps", "number of steps in increasing tree height (tree range = [0,..,nrOfHSteps*detlaH])", 100);
	final public Input<Boolean> useAminoAcidInput = new Input<>("useAminoAcid", "flag to use Amino Acids with WAG instead of Nucleotides with Jukes Cantor", false);

	final public Input<Double> rootHeightInput = new Input<>("rootHeight", "ignored if < 0, otherwise fix height of the root at this height", -1.0);

	@Override
	public void initAndValidate() {
	}

    int N = 100; // nr of steps in x direction
    double deltah = 0.02;
    int M = 100; // nr of steps in y direction
    double deltap = 1.0/(M-1);

    
    public double[][] calcRootStateProbs() throws Exception {

        String seq1Probs = new String("1.,0.0,0.0,0.0;");
        String seq2Probs = new String("1.,0.0,0.0,0.0;");
        Sequence seq1 = new Sequence();
		seq1.initByName("taxon","A","value",seq1Probs,"uncertain",true);
		Sequence seq2 = new Sequence();
		seq2.initByName("taxon","B","value",seq2Probs,"uncertain",true);

		
		Sequence seqFixed = new Sequence();
		seqFixed.initByName("taxon","B","value","1.,0.0,0.0,0.0;","uncertain",true);
		seqFixed.initByName("taxon","B","value","0.5,0.161616,0.161616,0.161616;","uncertain",true);
		seqFixed.initByName("taxon","B","value","0.25,0.25,0.25,0.25;","uncertain",true);
		seqFixed.initByName("taxon","B","value","0.0,0.333333,0.333333,0.333333;","uncertain",true);
		
        Alignment data0 = new Alignment();
        data0.initByName("sequence", seq1, "sequence", seq2, "dataType", "nucleotide");
        
        
        Tree tree = BEASTTestCase.getTree(data0, "(A:1.0,B:1.0)");
        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);

        TreeLikelihood likelihood0 = new TreeLikelihood();
        likelihood0.initByName("useTipLikelihoods", true, "data", data0, "tree", tree, "siteModel", siteModel, "scaling", TreeLikelihood.Scaling.none);
        double logP = 0;
        logP = likelihood0.calculateLogP();

        
        
        double [][] scores = new double[M][N];
        
        double h = 1e-10;
        for (int x = 0; x < N; x++) {
        	tree.getRoot().setHeight(h);
        	//tree.getRoot().getLeft().setHeight(h/2);
        	double p = 0;
        	for (int y = 0; y < M; y++) {
        		
                String seqProbs = new String(p+","+(1-p)/3.0+","+(1-p)/3.0+","+(1-p)/3.0+";");
                seq1 = new Sequence();
        		seq1.initByName("taxon","A","value",seqProbs,"uncertain",true);
        		seq2 = new Sequence();
        		seq2.initByName("taxon","B","value",seqProbs,"uncertain",true);
                        
                Alignment data = new Alignment();
                data.initByName("sequence", seq1, "sequence", seq2, "dataType", "nucleotide");
                TreeLikelihood likelihood = new TreeLikelihood();
                likelihood.initByName("useTipLikelihoods", true, "data", data, "tree", tree, "siteModel", siteModel, "scaling", TreeLikelihood.Scaling.none);

                logP = likelihood.calculateLogP();
                double [] partials = likelihood.getRootPartials();
                
                double sum = 0;
                for (double d : partials) {
                	sum += d;
                }
                
                //System.err.println(p + ":" + partials[0]/sum + ":" + Arrays.toString(partials) + " " + sum);
                System.err.println(tree.getRoot().toNewick());
            	scores[y][x] = partials[0]/sum;
            	
            	p += deltap;
        	}
        	h += deltah;
        }
        
        for (int x = 0; x < N; x++) {
        	double p = 0;
        	for (int y = 0; y < M; y++) {
        		//if (Math.abs(scores[y][x]-p) < deltap/2) {
            	if (scores[y][x]>p) {
        			scores[y][x] = - scores[y][x];
        		}
            	p += deltap;
        	}
        }
        
        
        return scores;
    }
     
    
    public double[][] calcRootStateProbs2() throws Exception {
        N = nrOfPStepsInput.get();
        deltap = 1.0/(N-1);
        M = N;
        deltah = deltap;

        String seq1Probs = new String("1.,0.0,0.0,0.0;");
        String seq2Probs = new String("1.,0.0,0.0,0.0;");
        Sequence seq1 = new Sequence();
		seq1.initByName("taxon","A","value",seq1Probs,"uncertain",true);
		Sequence seq2 = new Sequence();
		seq2.initByName("taxon","B","value",seq2Probs,"uncertain",true);
		
        Alignment data0 = new Alignment();
        data0.initByName("sequence", seq1, "sequence", seq2, "dataType", "nucleotide");
        
        Tree tree = BEASTTestCase.getTree(data0, "(A:1.0,B:1.0)");
        JukesCantor JC = new JukesCantor();
        JC.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", JC);

        TreeLikelihood likelihood0 = new TreeLikelihood();
        likelihood0.initByName("useTipLikelihoods", true, "data", data0, "tree", tree, "siteModel", siteModel, "scaling", TreeLikelihood.Scaling.none);
        double logP = 0;
        logP = likelihood0.calculateLogP();

    	tree.getRoot().setHeight(rootHeightInput.get());
        
        
        double [][] scores = new double[M][N];
        
        double q = 0;
        for (int x = 0; x < N; x++) {
        	//tree.getRoot().getLeft().setHeight(h/2);
        	double p = 0;
        	for (int y = 0; y < M; y++) {
        		
                String seqProbs = new String(p+","+(1-p)/3.0+","+(1-p)/3.0+","+(1-p)/3.0+";");
                seq1 = new Sequence();
        		seq1.initByName("taxon","A","value",seqProbs,"uncertain",true);
                String seqProbs2 = new String((1-q)/3.0+","+q+","+(1-q)/3.0+","+(1-q)/3.0+";");
        		seq2 = new Sequence();
        		seq2.initByName("taxon","B","value",seqProbs2,"uncertain",true);
                        
                Alignment data = new Alignment();
                data.initByName("sequence", seq1, "sequence", seq2, "dataType", "nucleotide");
                TreeLikelihood likelihood = new TreeLikelihood();
                likelihood.initByName("useTipLikelihoods", true, "data", data, "tree", tree, "siteModel", siteModel, "scaling", TreeLikelihood.Scaling.none);

                logP = likelihood.calculateLogP();
                double [] partials = likelihood.getRootPartials();
                
                double sum = 0;
                for (double d : partials) {
                	sum += d;
                }
                
                //System.err.println(p + ":" + partials[0]/sum + ":" + Arrays.toString(partials) + " " + sum);
                System.err.println(tree.getRoot().toNewick());
            	scores[y][x] = partials[0]/sum;
            	
            	p += deltap;
        	}
        	q += deltah;
        }
        
        for (int x = 0; x < N; x++) {
        	double p = 0;
        	for (int y = 0; y < M; y++) {
        		//if (Math.abs(scores[y][x]-p) < deltap/2) {
            	if (scores[y][x]>p) {
        			scores[y][x] = - scores[y][x];
        		}
            	p += deltap;
        	}
        }
        
        
        return scores;
    }

    public double[][] calcRootStateProbsAA() throws Exception {
        N = nrOfPStepsInput.get();
        deltap = 1.0/(N-1);
        M = nrOfHStepsInput.get();
        deltah = detlaHInput.get();

        String seq1Probs = new String("1.,0.0,0.0,0.0;");
        String seq2Probs = new String("1.,0.0,0.0,0.0;");
        Sequence seq1 = new Sequence();
		seq1.initByName("taxon","A","value",seq1Probs,"uncertain",true);
		Sequence seq2 = new Sequence();
		seq2.initByName("taxon","B","value",seq2Probs,"uncertain",true);

		
//		Sequence seqFixed = new Sequence();
//		seqFixed.initByName("taxon","B","value","1.,0.0,0.0,0.0;","uncertain",true);
//		seqFixed.initByName("taxon","B","value","0.5,0.161616,0.161616,0.161616;","uncertain",true);
//		seqFixed.initByName("taxon","B","value","0.25,0.25,0.25,0.25;","uncertain",true);
//		seqFixed.initByName("taxon","B","value","0.0,0.333333,0.333333,0.333333;","uncertain",true);
		
        Alignment data0 = new Alignment();
        data0.initByName("sequence", seq1, "sequence", seq2, "dataType", "aminoacid");
        
        Tree tree = BEASTTestCase.getTree(data0, "(A:1.0,B:1.0)");

        WAG WAG = new WAG();
        WAG.initAndValidate();

        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", WAG);
        
        
        double [][] scores = new double[M][N];
        
        double h = 1e-10;
        for (int x = 0; x < N; x++) {
        	tree.getRoot().setHeight(h);
        	tree.getRoot().getLeft().setHeight(h/2);
        	double p = 0;
        	for (int y = 0; y < M; y++) {
        		
                String seqProbs = p+"";
                String q = "," + (1-p)/19.0;
                for (int a = 0; a < 19; a++) {
                	seqProbs += q;
                }
                seqProbs += ";";

//                String seqProbs = p*2/3+","+p/3;
//                String q = "," + (1-p)/18.0;
//                for (int a = 0; a < 18; a++) {
//                	seqProbs += q;
//                }
//                seqProbs += ";";
                
                seq1 = new Sequence();
        		seq1.initByName("taxon","A","value",seqProbs,"uncertain",true);
        		seq2 = new Sequence();
        		seq2.initByName("taxon","B","value",seqProbs,"uncertain",true);
                        
                Alignment data = new Alignment();
                data.initByName("sequence", seq1, "sequence", seq2, "dataType", "aminoacid");
                TreeLikelihood likelihood = new TreeLikelihood();
                likelihood.initByName("useTipLikelihoods", true, "data", data, "tree", tree, "siteModel", siteModel, "scaling", TreeLikelihood.Scaling.none);

                double logP = likelihood.calculateLogP();
                double [] partials = likelihood.getRootPartials();
                
                double sum = 0;
                for (double d : partials) {
                	sum += d;
                }
                
                //System.err.println(p + ":" + partials[0]/sum + ":" + Arrays.toString(partials) + " " + sum);
                System.err.println(tree.getRoot().toNewick());
            	scores[y][x] = partials[0]/sum;
            	
            	p += deltap;
        	}
        	h += deltah;
        }
        
        for (int x = 0; x < N; x++) {
        	double p = 0;
        	for (int y = 0; y < M; y++) {
        		//if (Math.abs(scores[y][x]-p) < deltap/2) {
            	if (scores[y][x]>p) {
        			scores[y][x] = - scores[y][x];
        		}
            	p += deltap;
        	}
        }
        
        
        return scores;
    }
        
    @Override
	public void run() throws Exception {
        double [][] scores = useAminoAcidInput.get() ? calcRootStateProbsAA() : 
        	rootHeightInput.get() < 0 ? calcRootStateProbs() : calcRootStateProbs2();
                
        if (!svgInput.get().getName().equals("[[none]]")) {
        	Log.info("Writing to file " + svgInput.get().getPath());
        	PrintStream out = new PrintStream(svgInput.get() + (rootHeightInput.get() > 0?rootHeightInput.get() + ".svg" : ""));
        	out.print("<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?>\n");
	        out.print("<svg viewBox=\"0 0 " + (100+scores.length * 10) + " " + (50+scores.length * 10) + "\" "
	        		+ "xmlns=\"http://www.w3.org/2000/svg\">\n");
	        out.println("<g transform=\"translate(50,50)\">\n");
	        double min = 0;//Math.abs(scores[0][0]);
	        double max = 1;
	        for (int i = 0; i < scores.length; i++) {
	            for (int j = 0; j < scores[0].length; j++) {
	                double score = Math.abs(scores[i][j]);
	            	min = Math.min(min,  score);
	            	max = Math.max(max,  score);
	            }
	        }
	        for (int i = 0; i < scores.length; i++) {
	            for (int j = 0; j < scores[0].length; j++) {
	                double score = (Math.abs(scores[i][j]) - min) / (max - min);
	            	int c = (int) ((score - min) * 255/ (max - min));
	            	c = Color.HSBtoRGB(0.1f+(float)(0.9*score), 0.7f, 0.9f);
	            	String s = Integer.toHexString(c).substring(2);
	            	out.println(" <rect "
							+ "style=\"fill:#" + s+ ";\" "
	            			+ "width=\"10\" "
	            			+ "height=\"10\" "
	            			+ "x=\"" + i * 10 + "\" "
	            			+ "y=\"" + j * 10 + "\" />");

	            }
	        }
	        
	        for (int i = 0; i < scores.length; i++) {
	            for (int j = 0; j < scores[0].length; j++) {
	                double score = (Math.abs(scores[i][j]) - min) / (max - min);
	            	int c = (int) ((score - min) * 255/ (max - min));
	            	c = Color.HSBtoRGB(0.1f+(float)(0.9*score), 0.7f, 0.9f);
	            	String s = Integer.toHexString(c).substring(2);
	            	if (scores[i][j] < 0) {
		            	out.println(" <rect "
								+ "style=\"fill:#00F;\" "
		            			+ "width=\"4\" "
		            			+ "height=\"4\" "
		            			+ "x=\"" + (i * 10 +4)+ "\" "
		            			+ "y=\"" + (j * 10 +4)+ "\" />");
		            	if (rootHeightInput.get() > 0) {
			            	out.println(" <rect "
									+ "style=\"fill:#F00;\" "
			            			+ "width=\"4\" "
			            			+ "height=\"4\" "
			            			+ "x=\"" + (j * 10 +2)+ "\" "
			            			+ "y=\"" + (i * 10 +2)+ "\" />");
		            		
		            	}
	            		
	            	}
	            }
	        }
	        out.println("</g>\n");

	        // legend
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
            
            // min/max labals for legend
            DecimalFormat df = new DecimalFormat("#.###");
            DecimalFormat df1 = new DecimalFormat("#.#");
	        out.println("<text x=\"10\" y=\"10\">" + df.format(min) + "</text>");
	        out.println("<text x=\""+(scores.length*10+52)+"\" y=\"10\">" + df.format(max) + "</text>");
	        
	        // ticks
	        for (int i = 0; i <= scores.length; i += 10) {
            	out.println("<rect x=\"" + (50+i * 10) + "\" y=\"40\" width=\"1\" height=\"10\"/>");	        	
            	out.println("<text x=\"" + (50+i * 10) + "\" y=\"40\">" + df1.format(i * deltap) + "</text>");
	        }
	        if (rootHeightInput.get() > 0) {
	        	df = df1;
	        }
	        for (int i = 0; i <= scores[0].length; i += 10) {
            	out.println("<rect y=\"" + (50+i * 10) + "\" x=\"40\" width=\"10\" height=\"1\"/>");	        	
            	out.println("<text y=\"" + (50+i * 10) + "\" x=\"30\">" + df.format(i * deltah) + "</text>");	        	
	        }
	        out.print("</svg>\n");
	        out.close();
        }
        
        Log.warning("Done");
	}

	public static void main(String[] args) throws Exception {
		new Application(new JCCherryReconstructionMap(), "MutualInformation", args);

	}

}