package test.potts.dca;

import java.io.IOException;
import java.io.PrintStream;
import java.util.Arrays;

import beast.base.util.Randomizer;
import potts.dca.DCA;
import potts.dca.DCAFromFile;
import potts.dca.DCASequenceSimulator;

public class PFTest {
	final static double EPS = 1e-6;
	final static int N = 25;

	final static int PARTICLE_COUNT = 20;
	final static int STEP_COUNT_MULTIPLIER = 9;
	final static double TEMP = 1.5;
	
	static double [][]  log;
	static int [][]  particles;
	
	private static double findMAPState(int[][] particles, DCA dca, int instance) {
		int[][] particles2 = new int[PARTICLE_COUNT][dca.getSiteCount()];
		double [] currentLogP2 = new double[PARTICLE_COUNT];
		
		double [] currentLogP = new double[PARTICLE_COUNT];
		for (int i = 0; i < PARTICLE_COUNT; i++) {
			currentLogP[i] = calcLogP(dca, particles[i]);
		}
		System.err.print(currentLogP + " ");
		
		int stepCount = dca.getSiteCount() * STEP_COUNT_MULTIPLIER;
		int stateCount = dca.getStateCount();
		
		double temp = 1.0;
		int deltaStep = 1 + stepCount / 10;
		double sumLogP = 0;
		int sumCount = 0;
		
		if (log == null) {
			log = new double[N][stepCount];
		}
		
		for (int step = 0; step < stepCount; step++) {
            // Try to mutate every site once (Standard sweep)
			for (int k = 0; k < PARTICLE_COUNT; k++) {
				int [] seq = particles[k];
				
				int [] permutation = new int[dca.getSiteCount()];
				for (int i = 0; i < permutation.length; i++) {
					permutation[i] = i;
				}
				for (int i = 0; i < permutation.length; i++) {
					int target = Randomizer.nextInt(dca.getSiteCount());
					int tmp = permutation[target];
					permutation[target] = permutation[i];
					permutation[i] = tmp;
				}
				
	            for (int j = 0; j < dca.getSiteCount(); j++) {
	            	int i = permutation[j];
	                int oldState = seq[i];
	                int newState = Randomizer.nextInt(stateCount);
	                
	                if (oldState == newState) continue;
	
	                // Calculate change in Hamiltonian (Score)
	                // Delta H = H(new) - H(old)
	                // We accept if Delta H > 0 (more probable) or with prob exp(Delta H)
	                
	                double deltaH = DCASequenceSimulator.computeDeltaHamiltonian(dca, seq, i, oldState, newState);

	                // Metropolis Criterion
	                if (deltaH >= 0 || Randomizer.nextDouble() < Math.exp(deltaH) * temp) {
	                    seq[i] = newState; // Accept mutation
	                    currentLogP[k] += deltaH;
	                }
	            }
			}
            if (step % deltaStep == 0) {
            	System.err.print('.');
            }

            // resample
    		double [] logP = currentLogP.clone();
    		// normalise & transform from log space to real space
    		double max = logP[0];
    		for (double d : logP) {
    			max = Math.max(max, d);
    		}
    		for (int i = 0; i < PARTICLE_COUNT; i++) {
    			logP[i] = Math.exp((logP[i] - max) / TEMP);
    		}
            double sum = 0;
            for (double d : logP) { 
            	sum += d;
            }
            logP[0] /= sum;
            for (int i = 1; i < PARTICLE_COUNT; i++) {
    			logP[i] = logP[i]/sum + logP[i-1];
    		}
            int [] particleIndices = new int[PARTICLE_COUNT];
            for (int i = 0; i < PARTICLE_COUNT; i++) {
            	particleIndices[i] = Randomizer.randomChoice(logP);
            }
            Arrays.sort(particleIndices);
            
            for (int i = 0; i < PARTICLE_COUNT; i++) {
            	currentLogP2[i] = currentLogP[particleIndices[i]];
            	System.arraycopy(particles[particleIndices[i]], 0, particles2[i], 0, dca.getSiteCount());
            }
            double [] tmp = currentLogP2; currentLogP2 = currentLogP; currentLogP = tmp;
            int [][] tmp2 = particles2; particles2 = particles; particles = tmp2;
            
            
            double totalLogP = 0;
            for (double d : currentLogP) {
            	totalLogP += d;
            }
            totalLogP /= PARTICLE_COUNT;
            
            log[instance][step] = totalLogP;
        	sumLogP += totalLogP;
        	sumCount++;
        }
		
		sumLogP /= sumCount;
		return sumLogP;
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


	// usage: PFTest <file.dca>
	public static void main(String[] args) {
		long start = System.currentTimeMillis();
		Randomizer.setSeed(127);
		DCAFromFile dca = new DCAFromFile();
		dca.initByName("fileName", args[0]);
		

		double sum = 0;
		double sum2 = 0;
		particles = new int[PARTICLE_COUNT][dca.getSiteCount()];
		for (int i = 0; i < N; i++) {
			

			for (int k = 0; k < particles.length; k++) {
				for (int j = 0; j < particles[0].length; j++) {
					particles[k][j] = Randomizer.nextInt(dca.getStateCount());
				}
			}

			double logP = findMAPState(particles, dca, i);
			sum += logP;
			sum2 += logP * logP;
			System.out.println(logP);
		}
		sum /= N;
		sum2 /= N;
		System.out.println("mean = " + sum + "  stdev = " + Math.sqrt(sum2 - sum * sum));
		long end = System.currentTimeMillis();
		System.out.println("Done in " + (end - start) + "ms");
		
		
		try {
			PrintStream trace= new PrintStream("/tmp/trace.log");
			trace.print("sample\t");
			for (int i = 0; i < N; i++) {
				trace.print("logP" + i + "\t");
//				Arrays.sort(log[i]);
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
