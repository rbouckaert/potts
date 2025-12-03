package potts.likelihood;

import beagle.Beagle;
import beast.base.core.Description;
import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.core.Log;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.sitemodel.SiteModelInterface;
import beast.base.evolution.sitemodel.SiteModelInterface.Base;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.inference.parameter.IntegerParameter;
import beast.base.util.Randomizer;

@Description("Tree-likelihood that allows site specific root states that are allowed to change")
public class TreeLikelihoodWithRootStates extends beast.base.evolution.likelihood.TreeLikelihood {

	final public Input<IntegerParameter> rootFrequenciesSequenceInput = new Input<>("rootstates", 
			  "specifies site specific root states instead of root frequencies. ", Validate.REQUIRED);

	// private double [][] rootFrequenciesSequence;
	private int [] rootFrequenciesSequence;
	private int [] storedRootFrequenciesSequence;
	
	private int siteCount, stateCount, categoryCount, patternCount;


	@Override
	public void initAndValidate() {
		super.initAndValidate();
			
		if (dataInput.get().siteWeightsInput.get() != null) {
			throw new IllegalArgumentException("siteWeightsInput not implemented yet in " + this.getClass().getName());
		}
		
		stateCount = dataInput.get().getMaxStateCount();
		siteCount = dataInput.get().getSiteCount();
        // account for alignment weights here
		patternCount = dataInput.get().getPatternCount();
		categoryCount = ((SiteModelInterface.Base)siteModelInput.get()).getCategoryCount();
		if (categoryCount <= 0) {
			categoryCount = 1;
		}
		
		if (rootFrequenciesSequenceInput.get() != null) {
			
			rootFrequenciesSequence = new int[siteCount];
			storedRootFrequenciesSequence = new int[siteCount];
			
			// sanity check
			if (rootFrequenciesInput.get() != null) {
				throw new IllegalArgumentException("Either rootFrequencies or rootfreqseq can be specified, but not both");
			}
			
			
			initRootFrequencies();
			
			// sanity check
			if (siteCount != rootFrequenciesSequence.length) {
				throw new IllegalArgumentException("root sequence length (" + rootFrequenciesSequence.length + ") differs from alignment length("+ siteCount + ")");
			}
			
			patternLogLikelihoods = new double[siteCount];
			m_siteModel = (Base) siteModelInput.get();
			substitutionModel = m_siteModel.getSubstitutionModel();

		} else {
			rootFrequenciesSequence = null;
		}
	}
	
	

	protected void initRootFrequencies() {
		IntegerParameter seq = rootFrequenciesSequenceInput.get();
		
		// initialising root sequence based on frequencies occurring at each site
		Log.warning("initialising root sequence based on frequencies occurring at each site");
		
		Alignment data = dataInput.get();
		seq.setDimension(siteCount);
		for (int i = 0; i < siteCount; i++) {
			int [] pattern = data.getPattern(data.getPatternIndex(i));
			double [] probs = new double[stateCount];
			for (int k : pattern) {
				probs[k]++;
			}
			int newState = Randomizer.randomChoicePDF(probs);
			seq.setValue(i, newState);
		}
		
		Integer [] values = seq.getValues();
		for (int i = 0; i < siteCount; i++) {
			rootFrequenciesSequence[i] = values[i];
		}
	}

	
	@Override
	public double calculateLogP() {
        if (beagle != null) {
            logP = beagle.calculateLogP();
            if (rootFrequenciesSequence != null) {
            	logP = recalculateBeagleLogPWithRootFrequences();
            }
            return logP;
        }
        logP = super.calculateLogP();;
		return logP;
	}	


	private double [] rootpartials2 = null;
	private double [] rootpartials = null;
	
	private double recalculateBeagleLogPWithRootFrequences() {
		if (rootpartials == null) {
			rootpartials = new double[patternCount * stateCount];
			if (categoryCount > 1) {
				rootpartials2 = new double[patternCount * stateCount * categoryCount];
			}
		} else if (rootpartials.length != patternCount * stateCount) {
            System.out.println("root partials not equal to num patterns " + rootpartials.length + " != " + (patternCount * stateCount));
			rootpartials = new double[patternCount * stateCount];
			if (categoryCount > 1) {
				rootpartials2 = new double[patternCount * stateCount * categoryCount];
			}
		}
				
		int number = treeInput.get().getRoot().getNr();
		int node = beagle.getPartialBufferHelper().getOffsetIndex(number);
		if (categoryCount <= 1) {
			beagle.getBeagle().getPartials(node, Beagle.NONE, rootpartials);
		} else {
			beagle.getBeagle().getPartials(node, Beagle.NONE, rootpartials2);
			// integrate categories first
            final double[] proportions = ((SiteModelInterface.Base)siteModelInput.get()).getCategoryProportions(treeInput.get().getRoot());
            calculateIntegratePartials(rootpartials2, proportions, rootpartials);
		}
		
		double [] beaglePatternLogLikelihoods = beagle.getPatternLogLikelihoods();
		double [] scaleFactor = new double[patternCount];
		double [] freqs = substitutionModel.getFrequencies();
        int u = 0;
		for (int i = 0; i < patternCount; i++) {
            double sum = 0.0;
            for (int j = 0; j < stateCount; j++) {
                sum += freqs[j] * rootpartials[u];
                u++;
            }
            scaleFactor[i] = beaglePatternLogLikelihoods[i] - Math.log(sum);
		}
		
        for (int k = 0; k < siteCount; k++) {
        	int j = dataInput.get().getPatternIndex(k);
        	int v = j * stateCount;
            double sum = 0.0;
            sum += rootpartials[v + rootFrequenciesSequence[k]];
            patternLogLikelihoods[k] = Math.log(sum) + scaleFactor[j];
        }
		calcLogP();
		return logP;
	}

	
	protected void calculateIntegratePartials(double[] inPartials, double[] proportions, double[] outPartials) {

        int u = 0;
        int v = 0;
        for (int k = 0; k < patternCount; k++) {

            for (int i = 0; i < stateCount; i++) {

                outPartials[u] = inPartials[v] * proportions[0];
                u++;
                v++;
            }
        }


        for (int l = 1; l < categoryCount; l++) {
            u = 0;

            for (int k = 0; k < patternCount; k++) {

                for (int i = 0; i < stateCount; i++) {

                    outPartials[u] += inPartials[v] * proportions[l];
                    u++;
                    v++;
                }
            }
        }
    }
	
    protected void calcLogP() {
    	if (rootFrequenciesSequence != null) {
            logP = 0.0;

            if (dataInput.get().siteWeightsInput.get() != null) {
                // handle alignment with weights
                for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                    logP += patternLogLikelihoods[i] * dataInput.get().getPatternWeight(i);
                }
            } else {
                // full alignment without weights
                for (int i = 0; i < siteCount; i++) {
                    logP += patternLogLikelihoods[i];
                }
            }
            if (useAscertainedSitePatterns) {
                // at this point, logP contains contributions of ascertained sites, which it should not have.
                final double ascertainmentCorrection = dataInput.get().getAscertainmentCorrection(patternLogLikelihoods);
                
                int totalPatternWeight = 0;
                for (int i = 0; i < dataInput.get().getPatternCount(); i++) {
                	totalPatternWeight += dataInput.get().getPatternWeight(i);
                }
            	// TODO: test this
                // here, we subtract ascertainmentCorrection once for each site that should contribute to the likelihood (i.e. non-ascertaining sites)
                // and subtract the ascertainmentCorrection for the ascertaining sites already accumulated in logP above 
                logP +=  -ascertainmentCorrection * (totalPatternWeight + 1);
            }
    		
    	} else {
    		super.calcLogP();
    	}
    }

	

	@Override
	protected int traverse(Node node) {
        int update = (node.isDirty() | hasDirt);

        final int nodeIndex = node.getNr();

        final double branchRate = branchRateModel.getRateForBranch(node);
        final double branchTime = node.getLength() * branchRate;

        // First update the transition probability matrix(ices) for this branch
        //if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_StoredBranchLengths[nodeIndex])) {
        if (!node.isRoot() && (update != Tree.IS_CLEAN || branchTime != m_branchLengths[nodeIndex])) {
            m_branchLengths[nodeIndex] = branchTime;
            final Node parent = node.getParent();
            likelihoodCore.setNodeMatrixForUpdate(nodeIndex);
            for (int i = 0; i < m_siteModel.getCategoryCount(); i++) {
                final double jointBranchRate = m_siteModel.getRateForCategory(i, node) * branchRate;
                substitutionModel.getTransitionProbabilities(node, parent.getHeight(), node.getHeight(), jointBranchRate, probabilities);
                //System.out.println(node.getNr() + " " + Arrays.toString(m_fProbabilities));
                likelihoodCore.setNodeMatrix(nodeIndex, i, probabilities);
            }
            update |= Tree.IS_DIRTY;
        }

        // If the node is internal, update the partial likelihoods.
        if (!node.isLeaf()) {

            // Traverse down the two child nodes
            final Node child1 = node.getLeft(); //Two children
            final int update1 = traverse(child1);

            final Node child2 = node.getRight();
            final int update2 = traverse(child2);

            // If either child node was updated then update this node too
            if (update1 != Tree.IS_CLEAN || update2 != Tree.IS_CLEAN) {

                final int childNum1 = child1.getNr();
                final int childNum2 = child2.getNr();

                likelihoodCore.setNodePartialsForUpdate(nodeIndex);
                update |= (update1 | update2);
                if (update >= Tree.IS_FILTHY) {
                    likelihoodCore.setNodeStatesForUpdate(nodeIndex);
                }

                if (m_siteModel.integrateAcrossCategories()) {
                    likelihoodCore.calculatePartials(childNum1, childNum2, nodeIndex);
                } else {
                    throw new RuntimeException("Error TreeLikelihood 201: Site categories not supported");
                    //m_pLikelihoodCore->calculatePartials(childNum1, childNum2, nodeNum, siteCategories);
                }

                if (node.isRoot()) {
                    // No parent this is the root of the beast.tree -
                    // calculate the pattern likelihoods

                    final double[] proportions = m_siteModel.getCategoryProportions(node);
                    likelihoodCore.integratePartials(node.getNr(), proportions, m_fRootPartials);

                    if (getConstantPattern() != null) { // && !SiteModel.g_bUseOriginal) {
                        proportionInvariant = m_siteModel.getProportionInvariant();
                        // some portion of sites is invariant, so adjust root partials for this
                        for (final int i : getConstantPattern()) {
                            m_fRootPartials[i] += proportionInvariant;
                        }
                    }

                    if (this.rootFrequenciesSequence != null) {
                    	// use site specific root frequencies
                    	calculateLogLikelihoods(m_fRootPartials, this.rootFrequenciesSequence, patternLogLikelihoods);
                    } else {
                    	double[] rootFrequencies = substitutionModel.getFrequencies();
                    	if (rootFrequenciesInput.get() != null) {
                    		rootFrequencies = rootFrequenciesInput.get().getFreqs();
                    	}
                    	likelihoodCore.calculateLogLikelihoods(m_fRootPartials, rootFrequencies, patternLogLikelihoods);
                    }
                }

            }
        }
        return update;
    } // traverse
	

	private void calculateLogLikelihoods(
			double[] partials, 
			int[] rootFrequencies,
			double[] outLogLikelihoods) {

		Alignment data = dataInput.get();
//        if (data.siteWeightsInput.get() != null) {
//            // handle site weights
//            for (int k = 0; k < siteCount; k++) {
//                double sum = 0.0;
//                int i = data.getPatternIndex(k);
//                sum += partials[stateCount * i + rootFrequencies[i]];
//                outLogLikelihoods[k] = Math.log(sum) + getLikelihoodCore().getLogScalingFactor(k);
//            }
//        } else {
            // alignment and root are unweighted
            for (int k = 0; k < siteCount; k++) {
                double sum = 0.0;
                int i = data.getPatternIndex(k);
                sum += partials[stateCount * i + rootFrequencies[k]];
                outLogLikelihoods[k] = Math.log(sum) + getLikelihoodCore().getLogScalingFactor(i);
            }
//        }
	}

	
	@Override
	public void store() {
		System.arraycopy(rootFrequenciesSequence, 0, storedRootFrequenciesSequence, 0, siteCount);
		
		super.store();
	}
	
	@Override
	public void restore() {
		int [] tmp = storedRootFrequenciesSequence;
		storedRootFrequenciesSequence = rootFrequenciesSequence;
		rootFrequenciesSequence = tmp;
		
		super.restore();
	}


	@Override
	protected boolean requiresRecalculation() {
		boolean isDirty = super.requiresRecalculation();
		if (rootFrequenciesSequenceInput.get().isDirtyCalculation()) {
			initRootFrequencies();
			isDirty = true;
		}
		return isDirty;
	}
}
