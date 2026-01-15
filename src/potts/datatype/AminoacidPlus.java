package potts.datatype;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType.Base;

@Description("DataType for amino acids with explicit character for gaps.")
public class AminoacidPlus extends Base {

    public AminoacidPlus() {
        stateCount = 21;
        codeLength = 1;
        codeMap = "ACDEFGHIKLMNPQRSTVWY" + GAP_CHAR + "X" + MISSING_CHAR;

        mapCodeToStateSet = new int[23][];
        for (int i = 0; i < 21; i++) {
            mapCodeToStateSet[i] = new int[1];
            mapCodeToStateSet[i][0] = i;
        }
        int[] all = new int[21];
        for (int i = 0; i < 21; i++) {
            all[i] = i;
        }
        mapCodeToStateSet[21] = all;
        mapCodeToStateSet[22] = all;
    }

    @Override
    public String getTypeDescription() {
        return "aminoacidplus";
    }

}
