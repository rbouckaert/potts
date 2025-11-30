package potts.dca;

import java.util.Arrays;

import org.json.JSONArray;
import org.json.JSONException;
import org.json.JSONObject;

import beast.base.core.Description;

@Description("Base class for representing a direct-coupling analysis (DCA)")
public class DCA {
    protected int stateCount;    // Number of states (e.g., 21 for Amino Acids)
    protected int siteCount;     // Length of sequence
    protected int sequenceCount; // Number of sequences in MSA

    // --- Parameters to Learn ---
    // h[i][a] : Field for site i, state a
    protected double [][] h; 
    // J[i][j][a][b] : Coupling for site i (state a) and site j (state b)
    // Note: We only utilise entries where i < j to save redundant math
    protected double [][][][] J;

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
    
}
