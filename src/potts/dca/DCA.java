package potts.dca;

import java.util.Arrays;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

import beast.base.core.BEASTObject;
import beast.base.core.Description;

@Description("Base class for representing a direct-coupling analysis (DCA)")
public class DCA extends BEASTObject {
    protected int stateCount;    // Number of states (e.g., 21 for Amino Acids)
	protected int siteCount;     // Length of sequence
    protected int sequenceCount; // Number of sequences in MSA

    
	@Override
	public void initAndValidate() {
	}

    public int getStateCount() {
		return stateCount;
	}

	public int getSiteCount() {
		return siteCount;
	}

	public int getSequenceCount() {
		return sequenceCount;
	}

    // --- Parameters to Learn ---
    // h[i][a] : Field for site i, state a
    protected double [][] h; 
    // J[i][j][a][b] : Coupling for site i (state a) and site j (state b)
    // Note: We only utilise entries where i < j to save redundant math
    protected double [][][][] J;
    
    
    public double [][] getH() {return h;}
    public double [][][][] getJ() {return J;}

    public String toJSON() {
    	StringBuilder out = new StringBuilder();
        out.append("{" + "\n");
        out.append("\"stateCount\":\"" + stateCount +"\"," + "\n");
        out.append("\"siteCount\":\"" + siteCount +"\"," + "\n");
        out.append("\"seqCount\":\"" + sequenceCount +"\"," + "\n");

        out.append("\"h\":[" + "\n");
        for (int i = 0; i < siteCount; i++) {
        	out.append("\t" + Arrays.toString(h[i]));
        	if (i < siteCount-1) {
        		out.append(",");
        	}
        	out.append("\n");
        }
        out.append("]," + "\n");

        // J[i][j][a][b] : Coupling for site i (state a) and site j (state b)
        out.append("\"J\":[" + "\n");
        for (int i = 0; i < siteCount; i++) {
            out.append("\t[" + "\n");
            for (int j = i + 1; j < siteCount; j++) {
                out.append("\t\t[" + "\n");
            	for (int a = 0; a < stateCount; a++) {
            		out.append("\t\t\t" + Arrays.toString(J[i][j][a]));
                	out.append((a < stateCount-1) ? "," : "" + "\n");            		
            	}
            	out.append("\t\t]" + ((j < siteCount-1) ? "," : "") + "\n");
        	}
        	out.append("\t]" + ((i < siteCount-1) ? "," : "") + "\n");
        }
        out.append("]" + "\n");
        
        out.append("}" + "\n");
        return out.toString();
    }
    
    public void fromJSON(JSONObject o) throws JSONException {
    	stateCount = o.getInt("stateCount");
    	siteCount = o.getInt("siteCount");
    	sequenceCount = o.getInt("seqCount");
    	
        // Initialise Parameters (starts at 0.0)
        this.h = new double[siteCount][stateCount];
        this.J = new double[siteCount][siteCount][stateCount][stateCount];

    	JSONArray hArray = o.getJSONArray("h");
        for (int i = 0; i < siteCount; i++) {
        	JSONArray tmp = hArray.getJSONArray(i);
        	for (int j = 0; j < stateCount; j++) {
        		h[i][j] = tmp.getDouble(j);
        	}
        }    	

    	JSONArray JArray = o.getJSONArray("J");
        for (int i = 0; i < siteCount; i++) {
        	JSONArray tmp = JArray.getJSONArray(i);
            for (int j = i + 1; j < siteCount; j++) {
            	JSONArray tmp2 = tmp.getJSONArray(j - i - 1);
            	for (int a = 0; a < stateCount; a++) {
                	JSONArray tmp3 = tmp2.getJSONArray(a);
                	for (int b = 0; b < stateCount; b++) {
                		J[i][j][a][b] = tmp3.getDouble(b);
                	}
            	}
        	}
        }    	
    }
    
    // --- Helper to extract contact scores (Frobenius Norm) ---
    public double getContactScore(int i, int j) {
        if (i >= j) return 0.0;
        double sumSq = 0.0;
        // In full DCA, we would perform Average Product Correction (APC) here.
        // This is the raw interaction strength.
        for (int a = 0; a < stateCount; a++) {
            for (int b = 0; b < stateCount; b++) {
                // Often we subtract the row/col means (gauge invariance), 
                // but raw J norm is a decent first approx.
                double val = J[i][j][a][b];
                sumSq += val * val;
            }
        }
        return Math.sqrt(sumSq);
    }
    
    /**
     * Converts the asymmetric parameters learned by plmDCA into 
     * a symmetric Frobenius norm score (APC corrected usually).
     */
    public double[][] getContactScores() {
        double[][] scores = new double[siteCount][siteCount];

        // 1. Symmetrise and Compute Raw Frobenius Norm
        // plmDCA learns J_ij (effect of j on i) and J_ji (effect of i on j) separately.
        // We average them: J_final = (J_ij + J_ji) / 2
        
        double[] rowSums = new double[siteCount]; // For APC
        double[] colSums = new double[siteCount];
        double totalSum = 0.0;

        for (int i = 0; i < siteCount; i++) {
            for (int j = i + 1; j < siteCount; j++) {
                double normSq = 0.0;

                for (int a = 0; a < stateCount; a++) {
                    for (int b = 0; b < stateCount; b++) {
                        // Average the asymmetric couplings
                        // J_asym[target][neighbour][target_aa][neighbour_aa]
                        double val1 = J[i][j][b][a]; // predicting i using j
                        double val2 = J[j][i][a][b]; // predicting j using i
                        
                        // Zero-sum gauge correction is usually applied here in production tools
                        // For simplicity, we just average and square
                        double avg = 0.5 * (val1 + val2);
                        normSq += avg * avg;
                    }
                }
                
                double fn = Math.sqrt(normSq);
                scores[i][j] = fn;
                scores[j][i] = fn; // Symmetric
                
                rowSums[i] += fn;
                colSums[j] += fn; // Symmetric so colSums == rowSums
                totalSum += fn;
            }
        }

        // 2. Average Product Correction (APC)
        // APC_ij = Score_ij - (Score_i. * Score_.j) / Score_..
        // This removes background phylogenetic noise.
        
        double avgTotal = totalSum / (siteCount * (siteCount - 1) / 2.0); // Mean score
        
        // Actually, the standard APC formula uses sums over the row/col (L-1)
        for (int i = 0; i < siteCount; i++) {
            rowSums[i] /= (siteCount - 1);
        }

        double[][] apcScores = new double[siteCount][siteCount];
        double globalAvg = totalSum / (siteCount * siteCount); // Approx

        for (int i = 0; i < siteCount; i++) {
            for (int j = i + 1; j < siteCount; j++) {
                double correction = (rowSums[i] * rowSums[j]) / globalAvg;
                double val = scores[i][j] - correction;
                apcScores[i][j] = val;
                apcScores[j][i] = val;
            }
        }

        return apcScores;
    }

}
