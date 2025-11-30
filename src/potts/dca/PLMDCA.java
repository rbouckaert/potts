package potts.dca;

import java.util.Arrays;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;

import beast.base.util.Randomizer;

public class PLMDCA extends DCA {

    private final double lambda; // L2 Regularization strength (e.g. 0.01)
    
    // The Data
    private final int[][] msa;
    private final double[] weights;
    private final double meff;

    // The Result: Asymmetric Coupling Matrix
    // J_asym[i][j][a][b] is the interaction calculated when predicting site i using site j.
    // Note: In plmDCA, J[i][j] is not necessarily equal to J[j][i] during training.
    // We average them at the end.
//    private double[][][][] J;
//    private double[][] h;

    public PLMDCA(int[][] msa, double[] weights, int q, double lambda) {
        this.msa = msa;
        this.weights = weights;
        this.sequenceCount = msa.length;
        this.siteCount = msa[0].length;
        this.stateCount = q;
        this.lambda = lambda;
        
        // Calculate effective sample size
        double sumW = 0.0;
        for (double w : weights) sumW += w;
        this.meff = sumW;

        // Initialise parameter storage
        this.h = new double[siteCount][q];
        this.J = new double[siteCount][siteCount][q][q];
    }

    /**
     * Main training routine.
     * Uses parallel threads to optimise each site independently.
     */
    public void train(int iterations, double learningRate) {
        System.out.println("Starting plmDCA optimization...");

        // plmDCA is "embarrassingly parallel"
        // We can solve for site r=0, r=1, ... r=L independently.
        ExecutorService executor = Executors.newFixedThreadPool(Runtime.getRuntime().availableProcessors());

        for (int r = 0; r < siteCount; r++) {
            final int siteIndex = r;
            executor.submit(() -> optimiseSite(siteIndex, iterations, learningRate));
        }

        executor.shutdown();
        try {
            executor.awaitTermination(1, TimeUnit.DAYS);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        System.out.println("Optimization complete.");
    }

    /**
     * Solves the local optimisation problem for a specific target site 'r'.
     * Maximises P(Sequence[r] | Sequence[all other sites])
     */
    private void optimiseSite(int r, int iterations, double learningRate) {
        // Local parameters for site r
        // fields[a] corresponds to h_r(a)
        double[] local_h = new double[stateCount];
        
        // couplings[j][b][a] corresponds to J_rj(b, a)
        // We map site j (state b) -> target r (state a)
        double[][][] local_J = new double[siteCount][stateCount][stateCount];

        // Momentum buffers for Gradient Descent (Nesterov/Momentum)
        double[] v_h = new double[stateCount];
        double[][][] v_J = new double[siteCount][stateCount][stateCount];
        double mu = 0.5; // Momentum coefficient

        // Temporary storage for energies and probabilities per sequence
        double[] energies = new double[stateCount];
        double[] probs = new double[stateCount];

        // Gradient Descent Loop
        for (int iter = 0; iter < iterations; iter++) {
            
            // 1. Reset Gradients
            double[] grad_h = new double[stateCount];
            double[][][] grad_J = new double[siteCount][stateCount][stateCount]; // initialised to 0.0

            // 2. Compute Gradient over all sequences
            for (int m = 0; m < sequenceCount; m++) {
                int[] seq = msa[m];
                int targetAA = seq[r]; // The true amino acid at site r
                double w = weights[m];

                // A. Calculate Local Hamiltonian (Energy) for each possible AA at site r
                // E(a) = h_r(a) + Sum_{j!=r} J_rj(seq[j], a)
                double maxE = -Double.MAX_VALUE; // For LogSumExp stability

                for (int a = 0; a < stateCount; a++) {
                    double energy = local_h[a];
                    
                    for (int j = 0; j < siteCount; j++) {
                        if (r == j) continue;
                        int neighbourAA = seq[j];
                        energy += local_J[j][neighbourAA][a];
                    }
                    
                    energies[a] = energy;
                    if (energy > maxE) maxE = energy;
                }

                // B. Softmax -> Probabilities
                double zVal = 0.0;
                for (int a = 0; a < stateCount; a++) {
                    // Subtract maxE to prevent overflow in exp()
                    probs[a] = Math.exp(energies[a] - maxE);
                    zVal += probs[a];
                }
                
                // Normalise
                double invZ = 1.0 / zVal;
                for (int a = 0; a < stateCount; a++) {
                    probs[a] *= invZ;
                }

                // C. Accumulate Gradients
                // Gradient = Weight * (Indicator(TrueAA) - PredictedProb)
                // Note: We want to MAXIMISE likelihood, so we move in direction of gradient.
                
                // Update Field Gradient
                for (int a = 0; a < stateCount; a++) {
                    double diff = ((a == targetAA) ? 1.0 : 0.0) - probs[a];
                    grad_h[a] += w * diff;
                }

                // Update Coupling Gradient
                for (int j = 0; j < siteCount; j++) {
                    if (r == j) continue;
                    int neighbourAA = seq[j];
                    
                    for (int a = 0; a < stateCount; a++) {
                        double diff = ((a == targetAA) ? 1.0 : 0.0) - probs[a];
                        grad_J[j][neighbourAA][a] += w * diff;
                    }
                }
            } // End Sequence Loop

            // 3. Apply L2 Regularization and Update Parameters
            // Param_new = Param_old + LearningRate * (Gradient - 2*lambda*Param_old)
            
            // Normalise Gradient by M_eff (standard practice in DCA)
            double norm = 1.0 / meff;

            // Update h
            for (int a = 0; a < stateCount; a++) {
                double reg = 2 * lambda * local_h[a];
                double g = (grad_h[a] * norm) - reg;
                
                v_h[a] = mu * v_h[a] + learningRate * g;
                local_h[a] += v_h[a];
            }

            // Update J
            for (int j = 0; j < siteCount; j++) {
                if (r == j) continue;
                for (int b = 0; b < stateCount; b++) {
                    for (int a = 0; a < stateCount; a++) {
                        double reg = 2 * lambda * local_J[j][b][a];
                        double g = (grad_J[j][b][a] * norm) - reg;
                        
                        v_J[j][b][a] = mu * v_J[j][b][a] + learningRate * g;
                        local_J[j][b][a] += v_J[j][b][a];
                    }
                }
            }
        } // End Iteration Loop

        // Store results in global arrays
        // Synchronised is not needed because each thread writes to a unique 'r' index
        this.h[r] = local_h;
        this.J[r] = local_J;
        
        // Progress indicator
        if (r % 10 == 0) System.out.print(".");
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

    // --- Main for testing ---
    public static void main(String[] args) {
    	
    	long start = System.currentTimeMillis();
        // Create dummy data
        int L = 50;
        int M = 500;
        int q = 21; // Simplified alphabet
        double[] weights = new double[M];
        Arrays.fill(weights, 1.0);

        // Fill random data
        int[][] mockData = new int[M][L];
        
        // Create correlated data: Sites 0 and 1 are coupled (if 0 is A, 1 is A)
        for (int m = 0; m < M; m++) {
            for(int i=0; i < L; i++) {
            	mockData[m][i] = Randomizer.nextInt(q);
            }
            
            // Artificial Correlation:
            // Force site 1 to match site 0 usually
            if(Randomizer.nextDouble() > 0.1) {
            	mockData[m][1] = mockData[m][0];
            }

            if(Randomizer.nextDouble() > 0.1) {
            	mockData[m][3] = mockData[m][1];
            }
        }

        PLMDCA plm = new PLMDCA(mockData, weights, q, 0.01);
        
        // Train
        // 100 iterations, Learning Rate 0.1
        plm.train(100, 0.1); 

        // Get Scores
        double[][] scores = plm.getContactScores();
        System.out.println("Top Score (0,1): " + scores[0][1]);
        
        System.out.println("\nInferred Coupling Scores (Frobenius Norm):");
        // We expect high score between 0 and 1
        M = Math.min(5,  M);
        for(int i=0; i<M; i++) {
            for(int j=i+1; j<M; j++) {
                double score = scores[i][j];
                System.out.printf("Site %d - %d : %.4f\n", i, j, score);
            }
        }

        
    	long end = System.currentTimeMillis();

    	System.err.println("Done in " + (end-start) + " ms");
   }
}