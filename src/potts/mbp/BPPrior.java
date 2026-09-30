package potts.mbp;

import java.util.List;
import java.util.Map;
import static java.util.Map.entry;

import beast.base.core.Description;
import beast.base.evolution.datatype.Aminoacid;

@Description("Base pair prior: takes 2 amino acide sequences as input and for every "
		+ "site contribute p for paired and 1-p for unpaired")
public class BPPrior extends MBPPrior {
	
	
	@Override
	public void initAndValidate() {
		super.initAndValidate();

		Map<String,String> pairMap  = Map.ofEntries(
				entry("A","T"),
				entry("C","G"),
				entry("G","C"),
				entry("T","A")
				);
		
		// basepair map from amino acid to codon list
		Map<String, List<String>> mapAA2Codon = Map.ofEntries(
		        entry("A", List.of("GCT", "GCC", "GCA", "GCG")),
		        entry("R", List.of("CGT", "CGC", "CGA", "CGG", "AGA", "AGG")),
		        entry("N", List.of("AAT", "AAC")),
		        entry("D", List.of("GAT", "GAC")),
		        entry("C", List.of("TGT", "TGC")),
		        entry("Q", List.of("CAA", "CAG")),
		        entry("E", List.of("GAA", "GAG")),
		        entry("G", List.of("GGT", "GGC", "GGA", "GGG")),
		        entry("H", List.of("CAT", "CAC")),
		        entry("I", List.of("ATT", "ATC", "ATA")),
		        entry("L", List.of("TTA", "TTG", "CTT", "CTC", "CTA", "CTG")),
		        entry("K", List.of("AAA", "AAG")),
		        entry("M", List.of("ATG")),
		        entry("F", List.of("TTT", "TTC")),
		        entry("P", List.of("CCT", "CCC", "CCA", "CCG")),
		        entry("S", List.of("TCT", "TCC", "TCA", "TCG", "AGT", "AGC")),
		        entry("T", List.of("ACT", "ACC", "ACA", "ACG")),
		        entry("W", List.of("TGG")),
		        entry("Y", List.of("TAT", "TAC")),
		        entry("V", List.of("GTT", "GTC", "GTA", "GTG")),
		        entry("*", List.of("TAA", "TAG", "TGA")) // Stop codons
		    );

		
		Aminoacid aminoacid = new Aminoacid();
		match = new double[20][20];
		for (int i = 0; i < 20; i++) {
			String c1 = aminoacid.getCharacter(i);
			for (int j = 0; j < 20; j++) {
				String c2 = aminoacid.getCharacter(j);
				
				System.out.print(c1 + " " + c2 + " ");
				List<String> condons1 = mapAA2Codon.get(c1);
				List<String> condons2 = mapAA2Codon.get(c2);
				
				for (String codon1 : condons1) {
					for (String codon2 : condons2) {
						String pair2 = pairMap.get(codon2.charAt(2)+"") + pairMap.get(codon2.charAt(1)+"") + pairMap.get(codon2.charAt(0)+"");
						//System.out.print(codon1 + "," + codon2 + ">" + pair2 +" ");
						
						if (codon1.equals(pair2)) {
							match[i][j] = 1.0;
							//System.out.print("[[[match]]] ");
						}
					}
				}
				//System.out.println();
			}
		}
	}
	
	
	
}
