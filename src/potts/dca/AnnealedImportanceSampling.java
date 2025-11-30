package potts.dca;

import beast.base.util.Randomizer;




/**
 * 
This is a Java implementation of **Annealed Importance Sampling (AIS)** tailored for the Potts Model.

### Key Technical Details
1.  **Hamiltonian Sign Convention:** Consistent with the previous `BoltzmannDCA` code, this uses 
 		**Score/Hamiltonian ($H$)**, where $P(S) \propto e^{+H(S)}$.
    *   This reverses the sign compared to standard physics textbooks (which use Energy $E$, where $P \propto e^{-E}$).
    *   Weight update: $log(w) += (\beta_{new} - \beta_{old}) \times H_{couplings}$.
2.  **LogSumExp:** Crucial for combining results from multiple runs. Since $Z$ is astronomically large 
	($10^{100}+$ for proteins), we must perform all math in log-space to avoid floating-point overflow.
3.  **Optimization:** The code calculates the Independent Partition Function ($Z_A$) analytically 
	($O(L \cdot q)$) and then uses AIS to find the correction factor.

 
 ### How to use this with your previous DCA code
If you want to calculate the probability of a specific sequence $S$ after training your model:

1.  Train your model using `BoltzmannDCA` or `PLMDCA` to get `h` and `J`.
2.  Pass those `h` and `J` matrices into this `AnnealedImportanceSampling` class.
3.  Run `estimateLogZ(20000, 20)`. This returns `lnZ`.
4.  Calculate the Hamiltonian/Score of your sequence:
    ```java
    double score = 0;
    for(int i=0; i<L; i++) score += h[i][seq[i]];
    for(int i=0; i<L; i++) 
       for(int j=i+1; j<L; j++) 
          score += J[i][j][seq[i]][seq[j]];
    ```
5.  **Final Probability:**
    `double logProb = score - lnZ;`
    
 */
public class AnnealedImportanceSampling {

    private final int siteCount;
    private final int stateCount;
    private final double[][] h;
    private final double[][][][] J; // Assumed symmetric storage or accessible via i<j

    /**
     * @param h Fields [L][q]
     * @param J Couplings [L][L][q][q]
     * @param q Number of states
     */
    public AnnealedImportanceSampling(double[][] h, double[][][][] J, int q) {
        this.siteCount = h.length;
        this.stateCount = q;
        this.h = h;
        this.J = J;
    }

    /**
     * Main method to estimate Log(Z).
     * 
     * @param nBetas Number of intermediate temperature steps (e.g., 10,000)
     * @param nRuns  Number of independent AIS paths to average (e.g., 10 to 100)
     * @return The estimated natural logarithm of the Partition Function (ln Z)
     */
    public double estimateLogZ(int nBetas, int nRuns) {
        // 1. Calculate analytic Z for the independent model (Beta = 0)
        double logZ_indep = calculateIndependentLogZ();
        
        System.out.println("Analytic LogZ (Independent): " + logZ_indep);

        // 2. Perform AIS runs to estimate the ratio Z_model / Z_indep
        double[] logWeights = new double[nRuns];

        // Parallelise runs if desired, here sequential for clarity
        for (int i = 0; i < nRuns; i++) {
            logWeights[i] = runSingleAISPath(nBetas);
            if ((i + 1) % 10 == 0) System.out.printf("Finished run %d/%d\n", i + 1, nRuns);
        }

        // 3. Average the runs using LogSumExp trick
        // Z_est = Z_indep * (1/N * sum(weights))
        // log(Z_est) = log(Z_indep) - log(N) + log(sum(exp(logWeights)))
        double logSumWeights = logSumExp(logWeights);
        double logZ_correction = logSumWeights - Math.log(nRuns);

        System.out.println("AIS Correction Factor: " + logZ_correction);
        
        return logZ_indep + logZ_correction;
    }

    /**
     * Runs one AIS particle from Beta=0 to Beta=1.
     */
    private double runSingleAISPath(int nBetas) {
        // 1. Start with a sample from the Independent Model (Beta=0)
        int[] seq = sampleFromIndependentModel();
        
        double logWeight = 0.0;
        
        // Define Beta Schedule (Linear 0 -> 1)
        // Note: For very difficult models, sigmoidal schedules are better
        double stepSize = 1.0 / nBetas;

        for (int k = 0; k < nBetas; k++) {
            double beta_old = k * stepSize;
            double beta_new = (k + 1) * stepSize;

            // 2. Calculate Coupling Energy (H_J) for current sequence
            // Note: H_total = H_h + beta * H_J
            double H_J = calculateCouplingScore(seq);

            // 3. Update Weight
            // In Hamiltonian formalism: w = exp((beta_new - beta_old) * H_J)
            logWeight += (beta_new - beta_old) * H_J;

            // 4. MCMC Transition
            // Equilibrate sequence to distribution defined by beta_new
            // Typically 1 full sweep (L attempts) is sufficient per small beta step
            runMCMCSweep(seq, beta_new);
        }

        return logWeight;
    }

    /**
     * Exact calculation of Log Z for the model with J=0.
     * Z_indep = Product_i ( Sum_a exp(h_i(a)) )
     */
    private double calculateIndependentLogZ() {
        double totalLogZ = 0.0;
        for (int i = 0; i < siteCount; i++) {
            double[] localH = h[i];
            // Compute log(sum(exp(h))) safely
            totalLogZ += logSumExp(localH);
        }
        return totalLogZ;
    }

    /**
     * Samples a sequence from the independent model (J=0).
     * Probability of AA 'a' at site 'i' is softmax(h[i]).
     */
    private int[] sampleFromIndependentModel() {
        int[] seq = new int[siteCount];
        for (int i = 0; i < siteCount; i++) {
            double[] localH = h[i];
            
            // Standard trick to sample from logits
            // 1. Compute max for stability
            double maxH = -Double.MAX_VALUE;
            for(double v : localH) maxH = Math.max(maxH, v);
            
            // 2. Compute probabilities
            double sumP = 0.0;
            double[] probs = new double[stateCount];
            for (int a = 0; a < stateCount; a++) {
                probs[a] = Math.exp(localH[a] - maxH);
                sumP += probs[a];
            }
            
            // 3. Sample
            double r = Randomizer.nextDouble() * sumP;
            double cumsum = 0.0;
            for (int a = 0; a < stateCount; a++) {
                cumsum += probs[a];
                if (r <= cumsum) {
                    seq[i] = a;
                    break;
                }
            }
        }
        return seq;
    }

    /**
     * Calculates ONLY the coupling component of the Hamiltonian.
     * Score = Sum_{i<j} J_ij(a, b)
     */
    private double calculateCouplingScore(int[] seq) {
        double score = 0.0;
        for (int i = 0; i < siteCount; i++) {
            for (int j = i + 1; j < siteCount; j++) {
                // Access J assuming J[i][j] is valid for i < j
                score += J[i][j][seq[i]][seq[j]];
            }
        }
        return score;
    }

    /**
     * Runs one pass of Metropolis-Hastings to update sequence.
     * Target Hamiltonian: H(S) = H_field(S) + beta * H_coupling(S)
     */
    private void runMCMCSweep(int[] seq, double beta) {
        for (int i = 0; i < siteCount; i++) {
            // Pick random site and new state
            int site = Randomizer.nextInt(siteCount);
            int oldState = seq[site];
            int newState = Randomizer.nextInt(stateCount);
            
            if (oldState == newState) continue;

            // Calculate Delta Hamiltonian
            // 1. Field contribution
            double deltaH = h[site][newState] - h[site][oldState];

            // 2. Coupling contribution (scaled by Beta)
            double deltaJ = 0.0;
            for (int j = 0; j < siteCount; j++) {
                if (site == j) continue;
                int neighborState = seq[j];
                
                if (site < j) {
                    deltaJ += J[site][j][newState][neighborState] - J[site][j][oldState][neighborState];
                } else {
                    deltaJ += J[j][site][neighborState][newState] - J[j][site][neighborState][oldState];
                }
            }
            
            deltaH += beta * deltaJ;

            // Metropolis Criterion (Maximizing Score)
            // Accept if dH > 0 OR random < exp(dH)
            if (deltaH >= 0 || Randomizer.nextDouble() < Math.exp(deltaH)) {
                seq[site] = newState;
            }
        }
    }

    /**
     * Helper: LogSumExp for an array.
     * computes log(sum(exp(x_i))) safely.
     */
    public static double logSumExp(double[] values) {
        if (values.length == 0) return Double.NEGATIVE_INFINITY;
        double maxVal = values[0];
        for (double v : values) maxVal = Math.max(maxVal, v);
        
        double sumExp = 0.0;
        for (double v : values) {
            sumExp += Math.exp(v - maxVal);
        }
        return maxVal + Math.log(sumExp);
    }
    
    // -- Main for Testing --
    public static void main(String[] args) {
        // Setup tiny model
        int L = 10;
        int q = 2; // Binary (Spin glass like)
        
        double[][] h = new double[L][q];
        double[][][][] J = new double[L][L][q][q];
        
        // Make field prefer state 0
        for(int i=0; i<L; i++) h[i][0] = 1.0; 
        
        // Make neighbors want to be different (antiferromagnetic)
        for(int i=0; i<L-1; i++) {
             J[i][i+1][0][1] = 2.0;
             J[i][i+1][1][0] = 2.0;
        }

        AnnealedImportanceSampling ais = new AnnealedImportanceSampling(h, J, q);
        
        // Run Estimation
        // 10,000 steps, average over 10 runs
        double logZ = ais.estimateLogZ(10000, 10);
        
        System.out.println("Final Estimated Log Z: " + logZ);
        // Z should be roughly L*1.0 (fields) + contributions from J
    }
}
