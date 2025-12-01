package potts.likelihood;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.json.JSONException;
import org.json.JSONObject;

import beast.base.core.Input;
import beast.base.core.Input.Validate;
import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.branchratemodel.BranchRateModel;
import beast.base.evolution.branchratemodel.StrictClockModel;
import beast.base.evolution.likelihood.GenericTreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.tree.Node;
import beast.base.evolution.tree.Tree;
import beast.base.evolution.tree.TreeInterface;
import beast.base.util.Randomizer;
import beastfx.app.inputeditor.BeautiDoc;
import potts.dca.DCA;

public class PottsTreeLikelihood extends GenericTreeLikelihood {
	final public Input<File> dcaInput = new Input<>("dca", "DCA json file previously trained on an alignment", Validate.REQUIRED);
	final public Input<Double> temperatureInput = new Input<>("temperature", "temperature balancing the effect of Potts model Pt and substitution model Ps. "
			+ "Mutations are chosen proportional to Pt^1/t * Ps (thus 1/t", 1.0);

	int [][][] sequences;
	int siteCount;
	int stateCount;
	
    /**
     * BEASTObject associated with inputs. Since none of the inputs are StateNodes, it
     * is safe to link to them only once, during initAndValidate.
     */
    protected SubstitutionModel substitutionModel;
    protected SiteModel.Base m_siteModel;
    protected BranchRateModel.Base branchRateModel;

    public SubstitutionModel getSubstitutionModel() {return substitutionModel;}
    /**
     * flag to indicate the
     * // when CLEAN=0, nothing needs to be recalculated for the node
     * // when DIRTY=1 indicates a node partial needs to be recalculated
     * // when FILTHY=2 indicates the indices for the node need to be recalculated
     * // (often not necessary while node partial recalculation is required)
     */
    protected int hasDirt;

    /**
     * Lengths of the branches in the tree associated with each of the nodes
     * in the tree through their node  numbers. By comparing whether the
     * current branch length differs from stored branch lengths, it is tested
     * whether a node is dirty and needs to be recomputed (there may be other
     * reasons as well...).
     * These lengths take branch rate models in account.
     */
    protected double[] m_branchLengths;
    protected double[] storedBranchLengths;

    /**
     * memory allocation for probability tables obtained from the SiteModel *
     */
    protected double[] logProbabilities;

    protected int matrixSize;

    /**
     * flag to indicate ascertainment correction should be applied *
     */
    protected boolean useAscertainedSitePatterns = false;

    /**
     * alias for the data 
     */
    protected Alignment alignment;


    /** direct coupling analysis **/
    DCA dca;
    
    /**
     * dealing with proportion of site being invariant *
     */
    protected double proportionInvariant = 0;
    public double getProportionInvariant() {
    	return proportionInvariant;
    }
	public void setProportionInvariant(double proportionInvariant) {
		this.proportionInvariant = proportionInvariant;
	}

    @Override
    public void initAndValidate() {
        // sanity check: make sure data is an Alignment 
    	if (!(dataInput.get() instanceof Alignment)) {
    		throw new RuntimeException("Expected Alignment as data, not " + dataInput.get().getClass().getName());
    	}
        alignment = (Alignment) dataInput.get();

        int nodeCount = treeInput.get().getNodeCount();
        if (!(siteModelInput.get() instanceof SiteModel.Base)) {
        	throw new IllegalArgumentException("siteModel input should be of type SiteModel.Base");
        }
        m_siteModel = (SiteModel.Base) siteModelInput.get();
        m_siteModel.setDataType(alignment.getDataType());
        substitutionModel = m_siteModel.substModelInput.get();

        if (branchRateModelInput.get() != null) {
            branchRateModel = branchRateModelInput.get();
        } else {
            branchRateModel = new StrictClockModel();
        }
        m_branchLengths = new double[nodeCount];
        storedBranchLengths = new double[nodeCount];

        stateCount = alignment.getMaxStateCount();

        proportionInvariant = m_siteModel.getProportionInvariant();
        m_siteModel.setPropInvariantIsCategory(false);

        siteCount = alignment.getSiteCount();
        sequences = new int[2][nodeCount][siteCount];

        setStates(treeInput.get().getRoot(), siteCount);
        hasDirt = Tree.IS_FILTHY;

        matrixSize = (stateCount + 1) * (stateCount + 1);
        logProbabilities = new double[(stateCount + 1) * (stateCount + 1)];
        Arrays.fill(logProbabilities, 1.0);

        if (alignment.isAscertained) {
            useAscertainedSitePatterns = true;
        }
        
        try {
        	dca = new DCA();
			String json = BeautiDoc.load(dcaInput.get());
			dca.fromJSON(new JSONObject(json));
			
			if (dca.getSiteCount() != siteCount) {
				throw new IllegalArgumentException ("Site count of DCA (" + dca.getSiteCount() + ") and alignment (" + siteCount + ") should match");
			}
			if (dca.getStateCount() != stateCount+1) {
				throw new IllegalArgumentException ("State count of DCA (" + dca.getStateCount() + ") and alignment (" + stateCount + ") should match");
			}
        } catch (IOException | JSONException  e) {
			e.printStackTrace();
		}

    }

    
    /**
     * set leaf and internal states *
     */
    protected void setStates(Node node, int siteCount) {
        if (node.isLeaf()) {
            int i;
            int[] states = new int[siteCount];
            int taxonIndex = getTaxonIndex(node.getID(), alignment);
            for (i = 0; i < siteCount; i++) {
                int code = alignment.getPattern(taxonIndex, alignment.getPatternIndex(i));
                int[] statesForCode = alignment.getDataType().getStatesForCode(code);
                if (statesForCode.length==1)
                    states[i] = statesForCode[0];
                else
                    states[i] = code; // Causes ambiguous states to be ignored.
            }
            System.arraycopy(states, 0, sequences[0][taxonIndex], 0, siteCount);
        } else {
            setStates(node.getLeft(), siteCount);
            setStates(node.getRight(), siteCount);

            // internal nodes get random mixture of left and right sequences
        	int [] currentSeq = sequences[0][node.getNr()];
        	int [] leftSeq = sequences[0][node.getLeft().getNr()];
        	int [] righttSeq = sequences[0][node.getRight().getNr()];
        	for (int i = 0; i < siteCount; i++) {
        		currentSeq[i] = Randomizer.nextBoolean() ? leftSeq[i] : righttSeq[i];
        	}

        }
    }

    /**
    *
    * @param taxon the taxon name as a string
    * @param data the alignment
    * @return the taxon index of the given taxon name for accessing its sequence data in the given alignment,
    *         or -1 if the taxon is not in the alignment.
    */
   private int getTaxonIndex(String taxon, Alignment data) {
       int taxonIndex = data.getTaxonIndex(taxon);
       if (taxonIndex == -1) {
       	if (taxon.startsWith("'") || taxon.startsWith("\"")) {
               taxonIndex = data.getTaxonIndex(taxon.substring(1, taxon.length() - 1));
           }
           if (taxonIndex == -1) {
           	throw new RuntimeException("Could not find sequence " + taxon + " in the alignment");
           }
       }
       return taxonIndex;
	}

   
   
	double temperatureFactor;

   @Override
   public double calculateLogP() {
       logP = 0;
		temperatureFactor = 1.0 / temperatureInput.get();
       TreeInterface tree = treeInput.get();
       for (int i = tree.getLeafNodeCount(); i < tree.getNodeCount(); i++) {
    	   Node node = tree.getNode(i);
    	   if (!node.isRoot()) {
    		   logP += substmodelContribution(node); 
    	   }
		   logP += dcaModelContribution(node); 
       }
       return logP;
   }
   
   private double dcaModelContribution(Node node) {
	   double logP = 0;
       int [] seq = sequences[0][node.getNr()];
       double [][] h = dca.getH();
       double [][][][] J = dca.getJ();
       for (int i = 0; i < siteCount; i++) {
    	   int a = seq[i];
    	   logP += h[i][a];
           for (int j = i + 1; j < siteCount; j++) {
        	   int b = seq[j];
        	   logP += J[i][j][a][b];
           }
       }
	   return logP;
   }
   
   private double substmodelContribution(Node node) {
       final double branchRate = branchRateModel.getRateForBranch(node);
       final double jointBranchRate = m_siteModel.getRateForCategory(0, node) * branchRate;
       substitutionModel.getTransitionProbabilities(node, node.getParent().getHeight(), node.getHeight(), jointBranchRate, logProbabilities);
       
       for (int i = 0; i < logProbabilities.length; i++) {
    	   logProbabilities[i] = Math.log(logProbabilities[i]);
       }
       
       double logP = 0;
       int [] seq = sequences[0][node.getNr()];
       int [] parentSeq = sequences[0][node.getParent().getNr()];
       for (int i = 0; i < siteCount; i++) {
    	   logP += logProbabilities[parentSeq[i] * stateCount + seq[i]];
       }
	   return logP;
   }


}
