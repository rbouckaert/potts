package potts.tools;

import org.json.JSONArray;
import org.json.JSONObject;

import beast.base.core.Description;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.FileWriter;

@Description("convert file produced by github.com/anna-pa-m/adabmDCA so it can be read as JSON file")
public class BMA2JSON {

	public static void main(String[] args) {
		// Example usage
		try {
			convertPottsToJson(args[0], args[1], 84, 21, 449);
			System.out.println("Conversion complete!");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	/**
	 * Reads a Potts model text file and exports it to the specific JSON format.
	 * 
	 * @param inputPath  Path to the text file (h and J parameters)
	 * @param outputPath Path where the JSON will be saved
	 * @param siteCount  Number of sites (L)
	 * @param stateCount Number of states (q, usually 21)
	 * @param seqCount   Number of sequences used for training
	 */
	public static void convertPottsToJson(String inputPath, String outputPath, int siteCount, int stateCount,
			int seqCount) throws Exception {

		// 1. Initialize data structures
		double[][] h = new double[siteCount][stateCount];
		// Note: siteCount x siteCount x stateCount x stateCount
		double[][][][] J = new double[siteCount][siteCount][stateCount][stateCount];

		// 2. Parse the text file
		try (BufferedReader br = new BufferedReader(new FileReader(inputPath))) {
			String line;
			while ((line = br.readLine()) != null) {
				line = line.trim();
				if (line.isEmpty())
					continue;

				String[] parts = line.split("\\s+");
				String type = parts[0];

				if (type.equals("h")) {
					// Format: h i a value
					int i = Integer.parseInt(parts[1]);
					int a = Integer.parseInt(parts[2]);
					double val = Double.parseDouble(parts[3]);
					// Adjust indices if file is 1-indexed (common in bioinformatics)
					// If your file is 0-indexed, remove the -1
					h[i][a] = val;

				} else if (type.equals("J")) {
					// Format: J i j a b value
					int i = Integer.parseInt(parts[1]);
					int j = Integer.parseInt(parts[2]);
					int a = Integer.parseInt(parts[3]);
					int b = Integer.parseInt(parts[4]);
					double val = Double.parseDouble(parts[5]);

					J[i][j][a][b] = val;
					J[j][i][a][b] = val;
				}
			}
		}

		// 3. Construct JSON following the specific nesting logic of the fromJSON method
		JSONObject root = new JSONObject();
		root.put("stateCount", stateCount);
		root.put("siteCount", siteCount);
		root.put("seqCount", seqCount);

		// Map 'h' to JSONArray
		JSONArray hArray = new JSONArray();
		for (int i = 0; i < siteCount; i++) {
			JSONArray siteArr = new JSONArray();
			for (int a = 0; a < stateCount; a++) {
				siteArr.put(h[i][a]);
			}
			hArray.put(siteArr);
		}
		root.put("h", hArray);

		// Map 'J' to JSONArray using the (j - i - 1) relative indexing logic
		JSONArray jArrayOuter = new JSONArray();
		for (int i = 0; i < siteCount; i++) {
			JSONArray iArr = new JSONArray();
			// The parser expects iArr to contain arrays for j starting from i + 1
			for (int j = i + 1; j < siteCount; j++) {
				JSONArray abArr = new JSONArray();
				for (int a = 0; a < stateCount; a++) {
					JSONArray bArr = new JSONArray();
					for (int b = 0; b < stateCount; b++) {
						bArr.put(J[i][j][a][b]);
					}
					abArr.put(bArr);
				}
				iArr.put(abArr); // This becomes index (j - i - 1) in the parser
			}
			jArrayOuter.put(iArr);
		}
		root.put("J", jArrayOuter);

		// 4. Write to file
		try (FileWriter writer = new FileWriter(outputPath)) {
			writer.write(root.toString());
		}
	}
}