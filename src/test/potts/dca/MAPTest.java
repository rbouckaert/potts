package test.potts.dca;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;

import beast.base.util.Randomizer;
import potts.dca.DCA;
import potts.dca.DCAFromFile;
import potts.dca.DCASequenceSimulator;

public class MAPTest {
	final static double EPS = 1e-6;
	final static int N = 25;
	static double [][]  log;
	
	private static double findMAPState(int[] seq, DCA dca, int instance) {
		//initState(seq, dca);
		
		double maxLogP = calcLogP(dca, seq);
		System.err.print(maxLogP + " ");
		
		int stepCount = dca.getSiteCount() * 100;
		int stateCount = dca.getStateCount();
		
		int [] bestseq = new int[seq.length];
		
		
		double logP = maxLogP;
		double temp = 1.0;
		double delta = 5.0/stepCount;
		int deltaStep = 1 + stepCount / 10;
		double sumLogP = 0;
		int sumCount = 0;
		int burnin = stepCount / 2;
		
		if (log == null) {
			log = new double[N][stepCount];
		}
		
		for (int step = 0; step < stepCount; step++) {
            // Try to mutate every site once (Standard sweep)
            for (int i = 0; i < dca.getSiteCount(); i++) {
                int oldState = seq[i];
                int newState = Randomizer.nextInt(stateCount);
                
                if (oldState == newState) continue;

                // Calculate change in Hamiltonian (Score)
                // Delta H = H(new) - H(old)
                // We accept if Delta H > 0 (more probable) or with prob exp(Delta H)
                
                double deltaH = DCASequenceSimulator.computeDeltaHamiltonian(dca, seq, i, oldState, newState);
                
//double logP1 = calcLogP(dca, seq);                
//if (Math.abs(logP1 -logP) > EPS) {
//	int h = 3;
//	h++;
//}
                // Metropolis Criterion
                if (deltaH >= 0 || Randomizer.nextDouble() < Math.exp(deltaH) * temp) {
                    seq[i] = newState; // Accept mutation
                    logP += deltaH;
//double logP2 = calcLogP(dca, seq);
//if (Math.abs(logP2 -logP) > EPS) {
//	int h = 3;
//	h++;
//}
                    if (logP > maxLogP) {
                    	maxLogP = logP;
                    	System.arraycopy(seq, 0, bestseq, 0, seq.length);
                    }
                }
            }
//            temp = Math.exp(-step * delta);
            if (step % deltaStep == 0) {
            	System.err.print('.');
//            	System.err.println(logP);
            }
            log[instance][step] = logP;
            if (step > burnin && logP > 100) {
            	sumLogP += logP;
            	sumCount++;
            }
        }
		
		sumLogP /= sumCount;
		if (true) {
			return sumLogP;
		}
		
//		logP = calcLogP(dca, bestseq);
		return maxLogP;
	}
	
	private static double calcLogP(DCA dca, int[] seq) {
		double logP = 0;
		for (int i = 0; i < dca.getSiteCount(); i++) {
        // 1. Field Energy
			double delta = dca.getH()[i][seq[i]];
			logP += delta;
		}

        // 2. Coupling Energy with all other sites j
		for (int i = 0; i < dca.getSiteCount(); i++) {
			int newState = seq[i];
			double delta = 0;
			double [][][] J = dca.getJ()[i];
			for (int j = i + 1; j < dca.getSiteCount(); j++) {
				int neighborState = seq[j];
	            delta += J[j][newState][neighborState];
	        }
			logP += delta;
		}
        return logP;
 	}

	/** initialise with state that maximises product of DCA.h **/
	private static void initState(int[] state, DCA dca) {
		double [][] h = dca.getH();
		for (int i = 0; i < dca.getSiteCount(); i++) {
			double [] pattern = h[i];
			double max = pattern[0];
			int iMax = 0;
			for (int j = 0; j < pattern.length; j++) {
				if (pattern[j] > max) {
					max = pattern[j];
					iMax = j;
				}
			}
			state[i] = iMax;
		}
	}

	// usage: MAPTest <file.dca>
	public static void main(String[] args) {
		Randomizer.setSeed(127);
		DCAFromFile dca = new DCAFromFile();
		dca.initByName("fileName", args[0]);
		
		int [] state = new int[dca.getSiteCount()];
		for (int i = 0; i < state.length; i++) {
			state[i] = Randomizer.nextInt(dca.getStateCount());
		}

		double sum = 0;
		double sum2 = 0;
		for (int i = 0; i < N; i++) {
			
			for (int k = 0; k < state.length; k++) {
				state[k] = Randomizer.nextInt(dca.getStateCount());
			}

			double logP = findMAPState(state, dca, i);
			sum += logP;
			sum2 += logP * logP;
			System.out.println(logP);
		}
		sum /= N;
		sum2 /= N;
		System.out.println("mean = " + sum + "  stdev = " + Math.sqrt(sum2 - sum * sum));
		
		try {
			PrintStream trace= new PrintStream("/tmp/trace.log");
			trace.print("sample\t");
			for (int i = 0; i < N; i++) {
				trace.print("logP" + i + "\t");
				
				Arrays.sort(log[i]);
			}
			trace.println();
			for (int i = 0; i < log[0].length; i++) {
				trace.print(i + "\t");
				for (int j = 0; j < N; j++) {
					trace.print(log[j][i] + "\t");
				}
				trace.println();
			}
			trace.close();
			System.err.println("trace stored in /tmp/trace.log");
		} catch (IOException e) {
			e.printStackTrace();
		}
	}


}
