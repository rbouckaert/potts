package potts.dca;

import java.util.Arrays;

import beast.base.core.Log;
import beast.base.util.Randomizer;


/*


 
Given a set of sequences S_i
We assume the evolutionary data follows a Boltzmann distribution defined by the Potts Hamiltonian:
        
P(S∣h,J)=\frac{1}{Z(h,J)}exp⁡(∑_ih_i(S_i)+∑_{i<j}J_{ij}(S_i,S_j))

where Z(h,J) is the partition function defined as 

Z=∑_{all sequences S}​exp(−E(S))

is impossible to calculate due to the large number of sequences in the sum.


Frobenius Norm for pair i,j:

Score_{ij}​=\sqrt{∑_{a,b}​_Jij​(a,b)^2}


The goal is to find the parameters θ={h,J} that maximise the log-likelihood of the Multiple Sequence Alignment (MSA).
        
L(θ)=\frac{1}{M}∑_{m=1}^M\log ⁡P(S(m)∣θ)

Where S(m) are the sequences in your alignment.

 
## Implementation Details

Persistent Contrastive Divergence (PCD): It maintains a population of "fantasy particles" (simulated sequences) 
that persist between updates, rather than restarting from random noise every time. This makes convergence 
significantly faster.

Memory Efficiency: It only stores J[i][j] for i<j to save memory (symmetric matrix).

Vectorisation: Java doesn't have native vector/matrix operations like Python's NumPy, so the loops are explicit.

Sign Convention: This implementation maximises the "Hamiltonian" (Score).
* Score=∑h+∑J
* Probability ∝exp⁡(Score)
This is mathematically equivalent to minimising Energy (E=−Score).





## Key Considerations for Production Use

Sequence Weighting: Real MSAs contain clusters of very similar sequences (phylogenetic bias). 
In a production tool, you must calculate a weight wm for each sequence (usually 1/number of neighbours)
and use weighted sums when calculating fi_obs and fij_obs.

APC (Average Product Correction): The raw Frobenius norms calculated in getContactScore will 
contain background noise. You usually subtract the background using the APC formula:
 
F_{ij}^{APC}=F_{ij}−F_{i⋅}F_{.j}/F_{..}

Gauge Invariance: The Potts model is over-parameterised (shifting all energies by a constant 
doesn't change probabilities). Production codes usually apply "Zero-Sum gauge" constraints 
to the matrices to keep parameters centred around zero.

Runtime: This is O(M⋅L^2⋅Steps) For a protein of length 100, L2=10,000L2=10,000. This Java code 
is fast enough for small proteins or teaching, but for massive proteins (L > 500), you would 
likely use C++ with SIMD instructions or GPU acceleration.
 */

public class BoltzmannDCA extends DCA {

    // --- Hyperparameters ---
    private double eta;        // Learning rate (e.g., 0.05)
    private double lambda;     // L2 Regularisation penalty (e.g., 0.01)

    // --- Statistics ---
    // Observed from Input MSA (fixed target)
    private double [][] fi_obs;
    private double [][][][] fij_obs;
    
    // Calculated from Model (simulation)
    private double [][] fi_model;
    private double [][][][] fij_model;

    // --- Simulation State (PCD) ---
    // These sequences persist across training steps
    private int [][] mcmcSequences;

    /**
     * @param msa Input Multiple Sequence Alignment [num_seqs][seq_len]
     * @param stateCount Alphabet size (e.g., 21)
     * @param learningRate Step size for gradient ascent (e.g. 0.01)
     * @param regularisation L2 penalty weight (e.g. 0.001)
     */
    public BoltzmannDCA(int[][] msa, int stateCount, double learningRate, double regularisation) {
        this.sequenceCount = msa.length;
        this.siteCount = msa[0].length;
        this.stateCount = stateCount;
        this.eta = learningRate;
        this.lambda = regularisation;

        // Initialise Parameters (starts at 0.0)
        this.h = new double[siteCount][stateCount];
        this.J = new double[siteCount][siteCount][stateCount][stateCount];

        // Initialise Statistics Arrays
        this.fi_obs = new double[siteCount][stateCount];
        this.fij_obs = new double[siteCount][siteCount][stateCount][stateCount];
        this.fi_model = new double[siteCount][stateCount];
        this.fij_model = new double[siteCount][siteCount][stateCount][stateCount];

        // 1. Calculate Observed Statistics from the MSA
        calculateStatistics(msa, this.fi_obs, this.fij_obs);

        // 2. Initialise MCMC Chains
        // In PCD, we often start the chains as clones of the real data
        // to start closer to equilibrium.
        this.mcmcSequences = new int[sequenceCount][siteCount];
        for (int m = 0; m < sequenceCount; m++) {
            System.arraycopy(msa[m], 0, this.mcmcSequences[m], 0, siteCount);
        }
    }

    /**
     * Main training loop.
     * @param epochs Number of gradient updates to perform.
     * @param mcSteps Number of Metropolis sweeps per epoch.
     */
    public void train(int epochs, int mcSteps) {
        Log.info("Starting Boltzmann Learning on " + sequenceCount + " sequences, L=" + siteCount + "...");

        for (int epoch = 0; epoch < epochs; epoch++) {
            
            // Step 1: Run MCMC (Monte Carlo) to evolve the persistent chains
            // This brings the population closer to the distribution defined by current h and J
            runMetropolisHastings(mcSteps);

            // Step 2: Calculate statistics of the model's current output
            calculateStatistics(this.mcmcSequences, this.fi_model, this.fij_model);

            // Step 3: Update Parameters (Gradient Ascent)
            updateParameters();

            // Optional: Logging
            if (epoch % 100 == 0) {
                double diff = calculateAverageError();
                System.out.printf("Epoch %d: Avg Error (Obs - Model) = %.6f\n", epoch, diff);
            }
        }
    }

    /**
     * Updates h and J based on the difference between Observed and Model stats.
     */
    private void updateParameters() {
        // Update Fields (h)
        for (int i = 0; i < siteCount; i++) {
            for (int a = 0; a < stateCount; a++) {
                double grad = fi_obs[i][a] - fi_model[i][a];
                // New = Old + LearningRate * (Gradient - Regularisation)
                h[i][a] += eta * (grad - 2 * lambda * h[i][a]);
            }
        }

        // Update Couplings (J)
        for (int i = 0; i < siteCount; i++) {
            for (int j = i + 1; j < siteCount; j++) { // Only update i < j
                for (int a = 0; a < stateCount; a++) {
                    for (int b = 0; b < stateCount; b++) {
                        double grad = fij_obs[i][j][a][b] - fij_model[i][j][a][b];
                        J[i][j][a][b] += eta * (grad - 2 * lambda * J[i][j][a][b]);
                        // Enforce symmetry logic implicitly by only storing i<j
                    }
                }
            }
        }
    }

    /**
     * Runs Metropolis-Hastings sampling on the persistent chains.
     * @param steps Number of full sweeps (attempts per site)
     */
    private void runMetropolisHastings(int steps) {
        for (int step = 0; step < steps; step++) {
            for (int m = 0; m < sequenceCount; m++) {
                int[] seq = mcmcSequences[m];
                
                // Try to mutate every site once (Standard sweep)
                for (int i = 0; i < siteCount; i++) {
                    int oldState = seq[i];
                    int newState = Randomizer.nextInt(stateCount);
                    
                    if (oldState == newState) continue;

                    // Calculate change in Hamiltonian (Score)
                    // Delta H = H(new) - H(old)
                    // We accept if Delta H > 0 (more probable) or with prob exp(Delta H)
                    
                    double deltaH = computeDeltaHamiltonian(seq, i, oldState, newState);
                    
                    // Metropolis Criterion
                    if (deltaH >= 0 || Randomizer.nextDouble() < Math.exp(deltaH)) {
                        seq[i] = newState; // Accept mutation
                    }
                }
            }
        }
    }

    /**
     * Efficiently calculates the change in energy for a single mutation.
     * deltaH = (h_new - h_old) + Sum_neighbors(J_new_neighbor - J_old_neighbor)
     */
    private double computeDeltaHamiltonian(int[] seq, int i, int oldState, int newState) {
        // 1. Change in Field Energy
        double delta = h[i][newState] - h[i][oldState];

        // 2. Change in Coupling Energy with all other sites j
        for (int j = 0; j < siteCount; j++) {
            if (i == j) continue;
            
            int neighborState = seq[j];
            
            // Access J respecting i < j storage convention
            if (i < j) {
                delta += J[i][j][newState][neighborState] - J[i][j][oldState][neighborState];
            } else {
                // if i > j, we look up J[j][i][neighbour][target]
                delta += J[j][i][neighborState][newState] - J[j][i][neighborState][oldState];
            }
        }
        return delta;
    }

    /**
     * Helper to compute single and pairwise frequencies.
     * Note: In a real library, this would use sequence weighting.
     */
    private void calculateStatistics(int[][] sequences, double[][] fi, double[][][][] fij) {
        // 1. Reset counts
        for (double[] row : fi) Arrays.fill(row, 0.0);
        for (int i = 0; i < siteCount; i++) {
            for (int j = i + 1; j < siteCount; j++) {
                for (double[] row : fij[i][j]) Arrays.fill(row, 0.0);
            }
        }

        // 2. Accumulate counts
        for (int m = 0; m < sequenceCount; m++) {
            int[] seq = sequences[m];
            
            for (int i = 0; i < siteCount; i++) {
                fi[i][seq[i]]++;
                
                for (int j = i + 1; j < siteCount; j++) {
                    fij[i][j][seq[i]][seq[j]]++;
                }
            }
        }

        // 3. Normalise to probabilities
        double invM = 1.0 / sequenceCount;
        for (int i = 0; i < siteCount; i++) {
            for (int a = 0; a < stateCount; a++) {
                fi[i][a] *= invM;
            }
            for (int j = i + 1; j < siteCount; j++) {
                for (int a = 0; a < stateCount; a++) {
                    for (int b = 0; b < stateCount; b++) {
                        fij[i][j][a][b] *= invM;
                    }
                }
            }
        }
    }
    
    
    /**
     * Calculates single and pair frequencies with Sequence Weighting.
     * 
     * @param sequences The array of sequences (e.g. the MSA)
     * @param weights   The weight for each sequence. 
     *                  If null, assumes weight = 1.0 for all (Unweighted).
     * @param fi        Output: Single site frequencies
     * @param fij       Output: Pairwise frequencies
     */
    private void calculateStatistics(int[][] sequences, double[] weights, 
                                     double[][] fi, double[][][][] fij) {
        
        int numSeqs = sequences.length;
        double mEff = 0.0; // Effective number of sequences

        // 1. Reset counts to 0
        for (double[] row : fi) Arrays.fill(row, 0.0);
        for (int i = 0; i < siteCount; i++) {
            for (int j = i + 1; j < siteCount; j++) {
                for (double[] row : fij[i][j]) Arrays.fill(row, 0.0);
            }
        }

        // 2. Accumulate Weighted Counts
        for (int m = 0; m < numSeqs; m++) {
            int[] seq = sequences[m];
            
            // If weights array is provided, use it. Otherwise default to 1.0.
            double w = (weights != null) ? weights[m] : 1.0;
            
            // Track the total effective weight (denominator)
            mEff += w;

            for (int i = 0; i < siteCount; i++) {
                // Instead of ++, we add the weight
                fi[i][seq[i]] += w;
                
                for (int j = i + 1; j < siteCount; j++) {
                    fij[i][j][seq[i]][seq[j]] += w;
                }
            }
        }

        // 3. Normalise to probabilities using M_eff
        // Avoid division by zero check if empty array passed
        if (mEff == 0.0) return; 
        
        double invMeff = 1.0 / mEff;

        for (int i = 0; i < siteCount; i++) {
            for (int a = 0; a < stateCount; a++) {
                fi[i][a] *= invMeff;
            }
            for (int j = i + 1; j < siteCount; j++) {
                for (int a = 0; a < stateCount; a++) {
                    for (int b = 0; b < stateCount; b++) {
                        fij[i][j][a][b] *= invMeff;
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
    private double calculateAverageError() {
        double totalDiff = 0.0;
        int count = 0;
        for(int i=0; i<siteCount; i++) {
            for(int j=i+1; j<siteCount; j++) {
                for(int a=0; a<stateCount; a++) {
                    for(int b=0; b<stateCount; b++) {
                        totalDiff += Math.abs(fij_obs[i][j][a][b] - fij_model[i][j][a][b]);
                        count++;
                    }
                }
            }
        }
        return totalDiff / count;
    }

    // --- Main Method for Testing ---
    public static void main(String[] args) {
    	long start = System.currentTimeMillis();
        // Generate a tiny synthetic dataset (10 sequences, length 5, 4 states)
        // In reality, you would parse a FASTA file here and map AAs to 0..20
        int numSeqs = 51;
        int len = 50;
        int states = 3; 
        
        int[][] mockData = new int[numSeqs][len];
        
        // Create correlated data: Sites 0 and 1 are coupled (if 0 is A, 1 is A)
        for (int m = 0; m < numSeqs; m++) {
            for(int i=0; i<len; i++) {
            	mockData[m][i] = Randomizer.nextInt(states);
            }
            
            // Artificial Correlation:
            // Force site 1 to match site 0 usually
            if(Randomizer.nextDouble() > 0.1) {
            	mockData[m][1] = mockData[m][0];
            }
        }

        BoltzmannDCA dca = new BoltzmannDCA(mockData, states, 0.5, 0.001);
        
        // Train
        // In real usage: 2000+ epochs, 10+ mcSteps
        dca.train(2000, 2); 
        
        Log.info("\nInferred Coupling Scores (Frobenius Norm):");
        // We expect high score between 0 and 1
        for(int i=0; i<len; i++) {
            for(int j=i+1; j<len; j++) {
                double score = dca.getContactScore(i, j);
                System.out.printf("Site %d - %d : %.4f\n", i, j, score);
            }
        }
        
        Log.info("\nInferred Coupling Scores (Frobenius Norm) PCA corrected:");
        double [][] scores = dca.getContactScores();
        for (int i = 0; i < scores.length; i++) {
        	System.out.println(Arrays.toString(scores[i]));
        }
        
    	long end = System.currentTimeMillis();
        System.err.println("Done in " + (end-start) + " ms");
    }
    
    

}