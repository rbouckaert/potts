package potts.datatype;

import beast.base.core.Description;
import beast.base.evolution.datatype.DataType.Base;

@Description("DataType for amino acids with explicit character for gaps.")
public class AminoacidPlus extends Base {

    public AminoacidPlus() {
        stateCount = 21;
        codeLength = 1;
        codeMap = "ACDEFGHIKLMNPQRSTVWY" + "X" + GAP_CHAR + MISSING_CHAR;

        mapCodeToStateSet = new int[23][];
        for (int i = 0; i < 20; i++) {
            mapCodeToStateSet[i] = new int[1];
            mapCodeToStateSet[i][0] = i;
        }
        int[] all = new int[20];
        for (int i = 0; i < 20; i++) {
            all[i] = i;
        }
        mapCodeToStateSet[20] = all;
        mapCodeToStateSet[21] = new int[] {20};
        mapCodeToStateSet[22] = new int[] {20};
    }

    @Override
    public String getTypeDescription() {
        return "aminoacidplus";
    }

}
