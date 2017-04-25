package uk.ac.ebi.pride.tools.protein_details_fetcher;

import java.util.ArrayList;
import java.util.Map;

import uk.ac.ebi.pride.tools.protein_details_fetcher.model.Protein;
import uk.ac.ebi.pride.tools.protein_details_fetcher.model.Protein.STATUS;


import junit.framework.TestCase;

public class ProteinDetailFetcherTest extends TestCase {

  public void testGetProteinDetails() {
    ProteinDetailFetcher fetcher = new ProteinDetailFetcher();
    ArrayList<String> accessions = new ArrayList<>(3);
    accessions.add("P12345");
    accessions.add("P12346");
    accessions.add("P12347");
    accessions.add("TRFE_RAT"); // id for P12346
    accessions.add("ENSP00000263100");
    accessions.add("NP_004788");
    try {
      Map<String, Protein> proteins = fetcher.getProteinDetails(accessions);
      assertEquals(6, proteins.size());
      Protein p = proteins.get("P12345");
      assertNotNull(p);
      assertEquals("P12345", p.getAccession());
      assertEquals("Aspartate aminotransferase, mitochondrial (mAspAT) (EC 2.6.1.1) (EC 2.6.1.7) (Fatty acid-binding protein) (FABP-1) (Glutamate oxaloacetate transaminase 2) (Kynurenine aminotransferase 4) (Kynurenine aminotransferase IV) (Kynurenine--oxoglutarate transaminase 4) (Kynurenine--oxoglutarate transaminase IV) (Plasma membrane-associated fatty acid-binding protein) (FABPpm) (Transaminase A)", p.getName());
      assertEquals("MALLHSARVLSGVASAFHPGLAAAASARASSWWAHVEMGPPDPILGVTEAYKRDTNSKKMNLGVGAYRDDNGKPYVLPSVRKAEAQIAAKGLDKEYLPIGGLAEFCRASAELALGENSEVVKSGRFVTVQTISGTGALRIGASFLQRFFKFSRDVFLPKPSWGNHTPIFRDAGMQLQSYRYYDPKTCGFDFTGALEDISKIPEQSVLLLHACAHNPTGVDPRPEQWKEIATVVKKRNLFAFFDMAYQGFASGDGDKDAWAVRHFIEQGINVCLCQSYAKNMGLYGERVGAFTVICKDADEAKRVESQLKILIRPMYSNPPIHGARIASTILTSPDLRKQWLQEVKGMADRIIGMRTQLVSNLKKEGSTHSWQHITDQIGMFCFTGLKPEQVERLTKEFSIYMTKDGRISVAGVTSGNVGYLAHAIHQVTK", p.getSequenceString());
      p = proteins.get("TRFE_RAT");
      assertNotNull(p);
      assertEquals("TRFE_RAT", p.getAccession());
      assertEquals("Serotransferrin (Transferrin) (Beta-1 metal-binding globulin) (Liver regeneration-related protein LRRG03) (Siderophilin)", p.getName());
      assertEquals("MRFAVGALLACAALGLCLAVPDKTVKWCAVSEHENTKCISFRDHMKTVLPADGPRLACVKKTSYQDCIKAISGGEADAITLDGGWVYDAGLTPNNLKPVAAEFYGSLEHPQTHYLAVAVVKKGTDFQLNQLQGKKSCHTGLGRSAGWIIPIGLLFCNLPEPRKPLEKAVASFFSGSCVPCADPVAFPQLCQLCPGCGCSPTQPFFGYVGAFKCLRDGGGDVAFVKHTTIFEVLPQKADRDQYELLCLDNTRKPVDQYEDCYLARIPSHAVVARNGDGKEDLIWEILKVAQEHFGKGKSKDFQLFGSPLGKDLLFKDSAFGLLRVPPRMDYRLYLGHSYVTAIRNQREGVCPEGSIDSAPVKWCALSHQERAKCDEWSVSSNGQIECESAESTEDCIDKIVNGEADAMSLDGGHAYIAGQCGLVPVMAENYDISSCTNPQSDVFPKGYYAVAVVKASDSSINWNNLKGKKSCHTGVDRTAGWNIPMGLLFSRINHCKFDEFFSQGCAPGYKKNSTLCDLCIGPAKCAPNNREGYNGYTGAFQCLVEKGDVAFVKHQTVLENTNGKNTAAWAKDLKQEDFQLLCPDGTKKPVTEFATCHLAQAPNHVVVSRKEKAARVSTVLTAQKDLFWKGDKDCTGNFCLFRSSTKDLLFRDDTKCLTKLPEGTTYEEYLGAEYLQAVGNIRKCSTSRLLEACTFHKS", p.getSequenceString());
      p = proteins.get("P12346");
      assertNotNull(p);
      assertEquals("P12346", p.getAccession());
      assertEquals("Serotransferrin (Transferrin) (Beta-1 metal-binding globulin) (Liver regeneration-related protein LRRG03) (Siderophilin)", p.getName());
      assertEquals("MRFAVGALLACAALGLCLAVPDKTVKWCAVSEHENTKCISFRDHMKTVLPADGPRLACVKKTSYQDCIKAISGGEADAITLDGGWVYDAGLTPNNLKPVAAEFYGSLEHPQTHYLAVAVVKKGTDFQLNQLQGKKSCHTGLGRSAGWIIPIGLLFCNLPEPRKPLEKAVASFFSGSCVPCADPVAFPQLCQLCPGCGCSPTQPFFGYVGAFKCLRDGGGDVAFVKHTTIFEVLPQKADRDQYELLCLDNTRKPVDQYEDCYLARIPSHAVVARNGDGKEDLIWEILKVAQEHFGKGKSKDFQLFGSPLGKDLLFKDSAFGLLRVPPRMDYRLYLGHSYVTAIRNQREGVCPEGSIDSAPVKWCALSHQERAKCDEWSVSSNGQIECESAESTEDCIDKIVNGEADAMSLDGGHAYIAGQCGLVPVMAENYDISSCTNPQSDVFPKGYYAVAVVKASDSSINWNNLKGKKSCHTGVDRTAGWNIPMGLLFSRINHCKFDEFFSQGCAPGYKKNSTLCDLCIGPAKCAPNNREGYNGYTGAFQCLVEKGDVAFVKHQTVLENTNGKNTAAWAKDLKQEDFQLLCPDGTKKPVTEFATCHLAQAPNHVVVSRKEKAARVSTVLTAQKDLFWKGDKDCTGNFCLFRSSTKDLLFRDDTKCLTKLPEGTTYEEYLGAEYLQAVGNIRKCSTSRLLEACTFHKS", p.getSequenceString());
    }
    catch(Exception e) {
      e.printStackTrace();
      fail(e.getMessage());
    }
  }

  public void testDemerged() {
    ProteinDetailFetcher fetcher = new ProteinDetailFetcher();
    ArrayList<String> accessions = new ArrayList<>(3);
    accessions.add("P02551");
    accessions.add("P17073");
    try {
      Map<String, Protein> proteins = fetcher.getProteinDetails(accessions);
      Protein p = proteins.get("P02551");
      assertNotNull(p);
      assertEquals("P02551", p.getAccession());
      assertEquals(3, p.getReplacingProteins().size());
      assertEquals(STATUS.DEMERGED, p.getStatus());
      p = proteins.get("P17073");
      assertNotNull(p);
      assertEquals("P17073", p.getAccession());
      assertEquals(2, p.getReplacingProteins().size());
      assertEquals(STATUS.DEMERGED, p.getStatus());
    }
    catch(Exception e) {
      e.printStackTrace();
      fail(e.getMessage());
    }
  }

  public void testMerged() {
    ProteinDetailFetcher fetcher = new ProteinDetailFetcher();
    ArrayList<String> accessions = new ArrayList<>(3);
    accessions.add("O88541");
    accessions.add("Q63332");
    try {
      Map<String, Protein> proteins = fetcher.getProteinDetails(accessions);
      Protein p = proteins.get("O88541");
      assertNotNull(p);
      assertEquals("O88541", p.getAccession());
      assertEquals(1, p.getReplacingProteins().size());
      assertEquals("P24368", p.getReplacingProteins().get(0).getAccession());
      assertEquals(STATUS.MERGED, p.getStatus());
      p = proteins.get("Q63332");
      assertNotNull(p);
      assertEquals("Q63332", p.getAccession());
      assertEquals(1, p.getReplacingProteins().size());
      assertEquals("Q63041", p.getReplacingProteins().get(0).getAccession());
      assertEquals(STATUS.MERGED, p.getStatus());
    }
    catch(Exception e) {
      e.printStackTrace();
      fail(e.getMessage());
    }
  }

  public void testNcbi() {
    ProteinDetailFetcher fetcher = new ProteinDetailFetcher();
    ArrayList<String> accessions = new ArrayList<>(3);
    accessions.add("XP_216558");
    accessions.add("XP_217055");
    accessions.add("XP_224609");
    accessions.add("XP_218966.2"); // deleted
    try {
      Map<String, Protein> proteins = fetcher.getProteinDetails(accessions);
      Protein p;
      p = proteins.get("XP_216558");
      assertNotNull(p);
      assertEquals("XP_216558", p.getAccession());
      p = proteins.get("XP_217055");
      assertNotNull(p);
      assertEquals("XP_217055", p.getAccession());
      assertEquals(STATUS.UNKNOWN, p.getStatus());
      assertEquals(0, p.getReplacingProteins().size());
      p = proteins.get("XP_218966");
      assertNotNull(p);
      assertEquals("XP_218966", p.getAccession());
      assertEquals(STATUS.UNKNOWN, p.getStatus());
    }
    catch(Exception e) {
      e.printStackTrace();
      fail(e.getMessage());
    }
  }

  public void testUniprotIds() {
    ProteinDetailFetcher fetcher = new ProteinDetailFetcher();
    ArrayList<String> accessions = new ArrayList<>(3);
    accessions.add("CH601_ECOK1");
    accessions.add("PYGB_BOVIN");
    try {
      Map<String, Protein> proteins = fetcher.getProteinDetails(accessions);
      Protein p;
      p = proteins.get("CH601_ECOK1");
      assertNotNull(p);
      assertEquals("CH601_ECOK1", p.getAccession());
      assertEquals(STATUS.ACTIVE, p.getStatus());
      assertEquals("60 kDa chaperonin 1 (GroEL protein 1) (Protein Cpn60 1)", p.getName());
      assertEquals("MAAKDVKFGNDARVKMLRGVNVLADAVKVTLGPKGRNVVLDKSFGAPTITKDGVSVAREIELEDKFENMGAQMVKEVASKANDAAGDGTTTATVLAQAIITEGLKAVAAGMNPMDLKRGIDKAVTAAVEELKALSVPCSDSKAIAQVGTISANSDETVGKLIAEAMDKVGKEGVITVEDGTGLQDELDVVEGMQFDRGYLSPYFINKPETGAVELESPFILLADKKISNIREMLPVLEAVAKAGKPLLIIAEDVEGEALATLVVNTMRGIVKVAAVKAPGFGDRRKAMLQDIATLTGGTVISEEIGMELEKATLEDLGQAKRVVINKDTTTIIDGVGEEAAIQGRVAQIRQQIEEATSDYDREKLQERVAKLAGGVAVIKVGAATEVEMKEKKARVEDALHATRAAVEEGVVAGGGVALIRVASKLADLRGQNEDQNVGIKVALRAMEAPLRQIVLNCGEEPSVVANTVKGGDGNYGYNAATEEYGNMIDMGILDPTKVTRSALQYAASVAGLMITTECMVTDLPKNDAADLGAAGGMGGMGGMGGMM", p.getSequenceString());
      p = proteins.get("PYGB_BOVIN");
      assertNotNull(p);
      assertEquals(STATUS.ACTIVE, p.getStatus());
      assertEquals("Glycogen phosphorylase, brain form (EC 2.4.1.1)", p.getName());
      assertEquals("MAKPLTDGERRKQISVRGLAGLGDVAEVRKSFNRHLHFTLVKDRNVATRRDYYLALAHTVRDHLVGRWIRTQQRYYERDPKRIYYLSLEFYMGRTLQNTMVNLGLQNACDEAIYQLGLDLEELEEIEEDAGLGNGGLGRLAACFLDSMATLGLAAYGYGIRYEFGIFNQKIVNGWQVEEADDWLRYGNPWEKARPEYMLPVHFYGRVEHSPEGVRWLDTQVVLAMPYDTPVPGYKNDTVNTMRLWSAKAPNDFKLHDFNVGGYIEAVLDRNLAENISRVLYPNDNFFEGKELRLKQEYFVVAATLQDIIRRFKSSKFGCRDPVRTSFETFPDKVAIQLNDTHPALAIPELMRILVDVEKVDWDKAWEITKKTCAYTNHTVLPEALERWPVSMFEKLLPRHLDIIYAINQRHLDHVAALFPGDVDRLRRMSVIEEGDCKRINMAHLCVIGSHAVNGVARIHSEIVRQSVFKDFYELEPEKFQNKTNGITPRRWLLLCNPGLAETIVERIGEGFLTDLSQLKKLLPLVGDEALIRDVAQVKQENKVKFSAFLEKQYGVKVNPSSMFDVHVKRIHEYKRQLLNCLHVVTLYNRIKKDPTQAFVPRTVMIGGKAAPGYHMAKKIIKLVTSIGNIVNHDPIVGDRLKVIFLENYRVSLAEKVIPAADLSQQISTAGTEASGTGNMKFMLNGALTIGTMDGANVEMAEEAGAENLFIFGLRVEDVEALDRKGYNAHEYYDRLPELRQAVDQINGGFFSPREPDCFKDVVNMLLNHDRFKVFADYEAYVACQARVDQLYRNPKEWTKKVIRNIACSGKFSSDRTITEYAHDIWGAEPPALQTPPPSLPRD", p.getSequenceString());
    }
    catch(Exception e) {
      e.printStackTrace();
      fail(e.getMessage());
    }
  }

  public void testDeletedUniprot() {
    ProteinDetailFetcher fetcher = new ProteinDetailFetcher();
    ArrayList<String> accessions = new ArrayList<>(3);
    accessions.add("Q91WS9");
    try {
      Map<String, Protein> proteins = fetcher.getProteinDetails(accessions);
      Protein p;
      p = proteins.get("Q91WS9");
      assertNotNull(p);
      assertEquals("Q91WS9", p.getAccession());
      assertEquals(STATUS.DELETED, p.getStatus());
    }
    catch(Exception e) {
      e.printStackTrace();
      fail(e.getMessage());
    }
  }
}
