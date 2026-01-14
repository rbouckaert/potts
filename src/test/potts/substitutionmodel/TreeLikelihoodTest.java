package test.potts.substitutionmodel;


import org.junit.jupiter.api.Test;

import beast.base.evolution.alignment.Alignment;
import beast.base.evolution.alignment.Sequence;
import beast.base.evolution.likelihood.TreeLikelihood;
import beast.base.evolution.sitemodel.SiteModel;
import beast.base.evolution.substitutionmodel.Blosum62;
import beast.base.evolution.substitutionmodel.CPREV;
import beast.base.evolution.substitutionmodel.Dayhoff;
import beast.base.evolution.substitutionmodel.EmpiricalSubstitutionModel;
import beast.base.evolution.substitutionmodel.JTT;
import beast.base.evolution.substitutionmodel.MTREV;
import beast.base.evolution.substitutionmodel.SubstitutionModel;
import beast.base.evolution.substitutionmodel.WAG;
import beast.base.evolution.tree.Tree;
import potts.datatype.AminoacidPlus;
import potts.substitutionmodel.EmpiricalModelPlus;

import static org.junit.jupiter.api.Assertions.assertEquals;
import test.beast.BEASTTestCase;

/**
 * This test mimics the testLikelihood.xml file from Beast 1, which compares Beast 1 results to PAUP results.
 * So, it these tests succeed, then Beast II calculates the same for these simple models as Beast 1 and PAUP.
 * *
 */
public class TreeLikelihoodTest  {

    public TreeLikelihoodTest() {
        super();
    }

    protected TreeLikelihood newTreeLikelihood() {
    	System.setProperty("java.only","true");
        return new TreeLikelihood();
    }
    
    
    static public Alignment getAminoAcidPlusAlignment() throws Exception {
        Sequence Struthio_camelus = new Sequence("Struthio_camelus", "V-YPNTNEEGKEVVLPKILSPIGSDGVYSNELANIEYTNVSKAAAAAFATVDDYKPVPLDYMLDSKTSNKNNVVESSGTLRHFGK");
        Sequence Rhea_americana = new Sequence("Rhea_americana", "V-YPNTNEEGKEVLLPEILNPVGTDGVYSNELANIEYTNVAKDAAAAFATVDDHKPVSLEYMLDSKTSNKDNVVESNGTLSHFGK");
        Sequence Pterocnemia_pennata = new Sequence("Pterocnemia_pennata", "VKYPNTNEEGKEVLLPEILNPVGADGVYSNELANIEYTNVSKDHDEVFATVDDHKPVSLEYMLDSKTSNKDNVVESNGTLSHFGK");
        Alignment data = new Alignment();
        data.initByName("sequence", Struthio_camelus, "sequence", Rhea_americana, "sequence", Pterocnemia_pennata,
                "userDataType", new AminoacidPlus()
        );
        return data;
    }


    void aminoacidModelTest(EmpiricalSubstitutionModel eSubstModel, double expectedValue) throws Exception {
    	
    	SubstitutionModel.Base substModel = new EmpiricalModelPlus();
    	substModel.initByName("substModel", eSubstModel);
    	
        Alignment data = getAminoAcidPlusAlignment();
        Tree tree = BEASTTestCase.getAminoAcidTree(data);
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", substModel);

        TreeLikelihood likelihood = newTreeLikelihood();
        likelihood.initByName("data", data, "tree", tree, "siteModel", siteModel);
        double logP = 0;
        logP = likelihood.calculateLogP();
        assertEquals(expectedValue, logP, BEASTTestCase.PRECISION);
    }

    @Test
    public void testAminoAcidLikelihoodWAG() throws Exception {
        // Set up WAG model
        WAG wag = new WAG();
        wag.initAndValidate();
        aminoacidModelTest(wag, -378.8259637262415);
    }



    @Test
    public void testAminoAcidLikelihoodJTT() throws Exception {
        // JTT
        JTT jtt = new JTT();
        jtt.initAndValidate();
        aminoacidModelTest(jtt, -379.24972167046724);

    }

    @Test
    public void testAminoAcidLikelihoodBlosum62() throws Exception {
        // Blosum62
        Blosum62 blosum62 = new Blosum62();
        blosum62.initAndValidate();
        aminoacidModelTest(blosum62, -383.9021300914031);

    }

    @Test
    public void testAminoAcidLikelihoodDayhoff() throws Exception {
        // Dayhoff
        Dayhoff dayhoff = new Dayhoff();
        dayhoff.initAndValidate();
        aminoacidModelTest(dayhoff, -383.29937320519474);
    }

    @Test
    public void testAminoAcidLikelihoodcpRev() throws Exception {
        // cpRev
        CPREV cpRev = new CPREV();
        cpRev.initAndValidate();
        aminoacidModelTest(cpRev, -390.06748619568515);
    }

    @Test
    public void testAminoAcidLikelihoodMTRev() throws Exception {
        // MTRev
        MTREV mtRev = new MTREV();
        mtRev.initAndValidate();
        aminoacidModelTest(mtRev, -401.1500209849424);

    }


    /**
     * Test only effective when BEAGLE installed - otherwise it will always
     * pass since BEAGLE code will never be run.
     *
     * @throws Exception
     */
    @Test
    public void testBeagleRNALikelihood() throws Exception {

        Alignment data = getAminoAcidPlusAlignment();

        Tree tree = BEASTTestCase.getAminoAcidTree(data);

    	SubstitutionModel.Base substModel = new EmpiricalModelPlus();
    	WAG wag = new WAG();
    	wag.initAndValidate();
    	substModel.initByName("substModel", wag);
    	
        SiteModel siteModel = new SiteModel();
        siteModel.initByName("mutationRate", "1.0", "gammaCategoryCount", 1, "substModel", substModel);

        TreeLikelihood likelihoodNoBeagle = newTreeLikelihood();
        likelihoodNoBeagle.initByName("data", data, "tree", tree, "siteModel", siteModel);
        double logLnoBeagle = likelihoodNoBeagle.calculateLogP();

        System.setProperty("java.only", "false");
        TreeLikelihood likelihoodBeagle = new TreeLikelihood();
        likelihoodBeagle.initByName("data", data, "tree", tree, "siteModel", siteModel);
        double logLBeagle = likelihoodBeagle.calculateLogP();

        assertEquals(logLBeagle, logLnoBeagle, BEASTTestCase.PRECISION);
        assertEquals(logLBeagle, -378.8259637262415, BEASTTestCase.PRECISION);
    }
    
} // class TreeLikelihoodTest
