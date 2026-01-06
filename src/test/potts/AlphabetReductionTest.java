package test.potts;

import java.io.File;
import java.io.IOException;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import beast.base.evolution.datatype.Aminoacid;
import beast.base.evolution.datatype.DataType;

public class AlphabetReductionTest {
	
	DataType dataType = new Aminoacid();
	int [] order;

	public class Coef implements Comparable<Coef> {
		public Coef(int i2, int j2, double c2) {
			i = i2;
			j = j2;
			c= c2;
		}
		int i, j;
		double c;
		
		@Override
		public int compareTo(Coef o) {
			Coef other = (Coef) o;
			if (c > other.c) {
				return 1;
			}
			if (c < other.c) {
				return -1;
			}
			return 0;
		}
		
		@Override
		public String toString() {
			return dataType.getCharacter(i) + " " + dataType.getCharacter(j) + " (" +c + ")";
		}
	}
	
	void test() {
		order = new int[20];
		List<Integer> seq = dataType.stringToEncoding(getAlphabet());
		for (int i = 0; i < 20; i++) {
			order[i] = seq.get(i);
		}
		double [][] M = getM();
		
		double [][] distanceMatrix = new double[20][20]; 
		
		// calculate correlation coefficients
		List<Coef> list = new ArrayList<>();
		for (int i = 0; i < 20; i++) {
			for (int j = 0; j < 20; j++) {
				double c = corr(M, i, j);
				list.add(new Coef(i, j, c));
				distanceMatrix[i][j] = c;
			}
		}
		
		
		try {
			PrintStream out = new PrintStream(new File("/tmp/m.csv"));
			for (int i = 0; i < 20; i++) {
				out.print(getAlphabet().charAt(i));
				if (i < 19) {
					out.print(",");
				}
			}
			out.println();
			for (int i = 0; i < 20; i++) {
				for (int j = 0; j < 20; j++) {
					out.print(1.0-distanceMatrix[order[i]][order[j]]);
					if (j < 19) {
						out.print(",");
					}
				}
				out.println();
			}
			out.close();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		// sort coefficients
		Collections.sort(list);
		
		// print mergers
		int k = list.size()-1;
		int [] cluster = new int[20];
		for (int i = 0; i < 20; i++) {
			cluster[i] = i;
		}
		int clusterCount = 20;
		printCluster(cluster);
		System.out.println();
		while (clusterCount > 1) {
			Coef c = list.get(k--);
			if (cluster[c.i] != cluster[c.j]) {
				// merge
				int newCluster = Math.min(cluster[c.i], cluster[c.j]);
				int oldCluster = Math.max(cluster[c.i], cluster[c.j]);
				
				for (int j = 0; j < 20; j++) {
					if (cluster[j] == oldCluster) {
						cluster[j] = newCluster;
					}
				}
				clusterCount--;
				printCluster(cluster);
				System.out.println(c);
			}
		}
		
		
	}
	
	
    private void printCluster(int[] cluster) {
		for (int m = 0; m < 20; m++) {
			int i = order[m];
			System.out.print("(");
			for (int k = 0; k < 20; k++) {
				int j = order[k];
				if (cluster[j] == i) {
					System.out.print(dataType.getCharacter(j));
				}
			}
			System.out.print(")");
		}
		System.out.print(" ");		
	}


	private String getAlphabet() {return  "ARNDCQEGHILKMFPSTWYV";}
    
	  private final float[][] blossom50residueScores =
	            /* A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V */
	  { /* A */ {  5                                                          },
	    /* R */ { -2, 7                                                       },
	    /* N */ { -1,-1, 7                                                    },
	    /* D */ { -2,-2, 2, 8                                                 },
	    /* C */ { -1,-4,-2,-4,13                                              },
	    /* Q */ { -1, 1, 0, 0,-3, 7                                           },
	    /* E */ { -1, 0, 0, 2,-3, 2, 6                                        },
	    /* G */ {  0,-3, 0,-1,-3,-2,-3, 8                                     },
	    /* H */ { -2, 0, 1,-1,-3, 1, 0,-2,10                                  },
	    /* I */ { -1,-4,-3,-4,-2,-3,-4,-4,-4, 5                               },
	    /* L */ { -2,-3,-4,-4,-2,-2,-3,-4,-3, 2, 5                            },
	    /* K */ { -1, 3, 0,-1,-3, 2, 1,-2, 0,-3,-3, 6                         },
	    /* M */ { -1,-2,-2,-4,-2, 0,-2,-3,-1, 2, 3,-2, 7                      },
	    /* F */ { -3,-3,-4,-5,-2,-4,-3,-4,-1, 0, 1,-4, 0, 8                   },
	    /* P */ { -1,-3,-2,-1,-4,-1,-1,-2,-2,-3,-4,-1,-3,-4,10                },
	    /* S */ {  1,-1, 1, 0,-1, 0,-1, 0,-1,-3,-3, 0,-2,-3,-1, 5             },
	    /* T */ {  0,-1, 0,-1,-1,-1,-1,-2,-2,-1,-1,-1,-1,-2,-1, 2, 5          },
	    /* W */ { -3,-3,-4,-5,-5,-1,-3,-3,-3,-3,-2,-3,-1, 1,-4,-4,-3,15       },
	    /* Y */ { -2,-1,-2,-3,-3,-1,-2,-3, 2,-1,-1,-2, 0, 4,-3,-2,-2, 2, 8    },
	    /* V */ {  0,-3,-3,-4,-1,-3,-3,-4,-4, 4, 1,-3, 1,-1,-3,-2, 0,-3,-1, 5 }
	            /* A  R  N  D  C  Q  E  G  H  I  L  K  M  F  P  S  T  W  Y  V */
	  };

	  
	  private final float[][] pam250residueScores = {

	            /*  A   R   N   D   C   Q   E   G   H   I   L   K   M   F   P   S   T   W   Y   V */
	            {   2},
	            {  -2,  6},
	            {   0,  0,  2},
	            {   0, -1,  2,  4},
	            {  -2, -4, -4, -5, 12},
	            {   0,  1,  1,  2, -5,  4},
	            {   0, -1,  1,  3, -5,  2,  4},
	            {   1, -3,  0,  1, -3, -1,  0,  5},
	            {  -1,  2,  2,  1, -3,  3,  1, -2,  6},
	            {  -1, -2, -2, -2, -2, -2, -2, -3, -2,  5},
	            {  -2, -3, -3, -4, -6, -2, -3, -4, -2,  2,  6},
	            {  -1,  3,  1,  0, -5,  1,  0, -2,  0, -2, -3,  5},
	            {  -1,  0, -2, -3, -5, -1, -2, -3, -2,  2,  4,  0,  6},
	            {  -3, -4, -3, -6, -4, -5, -5, -5, -2,  1,  2, -5,  0,  9},
	            {   1,  0,  0, -1, -3,  0, -1,  0,  0, -2, -3, -1, -2, -5,  6},
	            {   1,  0,  1,  0,  0, -1,  0,  1, -1, -1, -3,  0, -2, -3,  1,  2},
	            {   1, -1,  0,  0, -2, -1,  0,  0, -1,  0, -2,  0, -1, -3,  0,  1,  3},
	            {  -6,  2, -4, -7, -8, -5, -7, -7, -3, -5, -2, -3, -4,  0, -6, -2, -5, 17},
	            {  -3, -4, -2, -4,  0, -4, -4, -5,  0, -1, -1, -4, -2,  7, -5, -3, -3,  0, 10},
	            {   0, -2, -2, -2, -2, -2, -2, -1, -2,  4,  2, -2,  2, -1, -1, -1,  0, -6, -2,  4}};

	  private double[][] getM() {
		double [][] M = new double[20][20];
		
		double [][] score = buildScores(blossom50residueScores);
		//double [][] score = buildScores(pam250residueScores);
		
		for (int i = 0; i < 20; i++) {
			char c1 = dataType.getCharacter(i).charAt(0);
			for (int j = 0; j < 20; j++) {
				char c2 = dataType.getCharacter(j).charAt(0);
				M[i][j] = score[c1][c2];
			}
		}
		
		// shift to ensure mean = 0 in each row
		for (int i = 0; i < 20; i++) {
			double sum = 0;
			for (int j = 0; j < 20; j++) {
				sum += M[i][j];
			}
			double mean = sum/20;
			for (int j = 0; j < 20; j++) {
				M[i][j] -= mean;
			}
		}
		return M;
	}
	
    /**
     * @param scores float[][] with position [i][j] holding the score for
     *        getScore(getAlphabet().charAt(i), getAlphabet().charAt(j)).
     */
    protected double [][] buildScores(float[][] scores) {
        String states = getAlphabet().toUpperCase();
        // Allow lowercase and uppercase states (ASCII code <= 127):
        double [][] score = new double[127][127];
        for (int i=0; i<states.length(); i++) {
            char a = states.charAt(i);
            char lca = Character.toLowerCase(a);
            for (int j=0; j<=i; j++) {
                char b = states.charAt(j);
                char lcb = Character.toLowerCase(b);
                score[a][b] = score[b][a]
                    = score[a][lcb] = score[lcb][a]
                    = score[lca][b] = score[b][lca]
                    = score[lca][lcb] = score[lcb][lca]
                    = scores[i][j];
            }
        }
        return score;
    }


	private static double corr(double [][] M, int i1, int i2) {
		double sum = 0;
		double v1 = 0, v2 = 0;
		for (int i = 0; i < 20; i++) {
			sum += M[i1][i] * M[i2][i];
			v1 += M[i1][i] * M[i1][i];
			v2 += M[i2][i] * M[i2][i];
		}
		return sum / Math.sqrt(v1 * v2);
	}

	public static void main(String[] args) {
		AlphabetReductionTest t = new AlphabetReductionTest();
		t.test();
	}
		
}
