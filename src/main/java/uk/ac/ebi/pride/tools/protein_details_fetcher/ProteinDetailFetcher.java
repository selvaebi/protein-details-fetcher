package uk.ac.ebi.pride.tools.protein_details_fetcher;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.*;
import java.util.concurrent.CompletableFuture;
import java.util.concurrent.ExecutionException;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.apache.commons.lang3.StringUtils;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.pride.tools.protein_details_fetcher.model.Protein;
import uk.ac.ebi.pride.tools.protein_details_fetcher.model.Protein.PROPERTY;
import uk.ac.ebi.pride.tools.protein_details_fetcher.model.Protein.STATUS;
import uk.ac.ebi.pride.tools.protein_details_fetcher.util.ProteinAccessionPattern;
import uk.ac.ebi.uniprot.dataservice.client.Client;
import uk.ac.ebi.uniprot.dataservice.client.ServiceFactory;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.*;
import uk.ac.ebi.uniprot.dataservice.client.alignment.blast.input.DatabaseOption;

public class ProteinDetailFetcher {

  private static final Logger logger = LoggerFactory.getLogger(ProteinDetailFetcher.class);

  private static final String TAB = "\t";
  private static final String EOL = "\n";
  private static final String UNIPROT_PROTEIN = "P31946L";

  /**
   * Queries used to fetch the actual protein
   * details.
   * @author jg
   *
   */
  public enum DETAILS_QUERY{
    NCBI_FASTA("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=%s&rettype=fasta&tool=protein_details_fetcher"),
    UNIPROT("http://www.uniprot.org/uniprot/?query=%s&format=tab&columns=id,protein%%20names,sequence,reviewed,entry%%20name,organism-id"),
    UNIPROT_FASTA("http://www.uniprot.org/uniprot/%s");

    private String queryString;

    DETAILS_QUERY(String formatString) {
      this.queryString = formatString;
    }

    public String getQueryString() {
      return queryString;
    }
  }

  /**
   * Queries used to map accession systems.
   * @author jg
   *
   */
  public enum MAPPING_QUERY {
    ENSEMBL_TO_UNIPROT("http://www.uniprot.org/mapping/?from=ENSEMBL_PRO_ID&to=ACC&query=%s&format=tab");


    private String queryString;

    MAPPING_QUERY(String queryString) {
      this.queryString = queryString;
    }

    public String getQueryString() {
      return queryString;
    }
  }

  /**
   * Supported types of protein accessions.
   * @author jg
   *
   */
  public enum AccessionType{UNIPROT_ACC, UNIPROT_ID, UNIPARC, REFSEQ, ENSEMBL, UNKNOWN}

  /**
   * Returns the (guessed) accession type for the passed
   * accession. In case the accession is not recognized
   * UNKNOWN is returned.
   * @param accession The accession to guess the type for.
   * @return The accession's type.
   */
  public static AccessionType getAccessionType(String accession) {
    AccessionType result = AccessionType.UNKNOWN;
    if (ProteinAccessionPattern.isSwissprotAccession(accession)) {
      result =  AccessionType.UNIPROT_ACC; // swissprot accession
    } else if (ProteinAccessionPattern.isSwissprotEntryName(accession)) {
      result = AccessionType.UNIPROT_ID; // swissprot entry
    } else if (ProteinAccessionPattern.isUniparcAccession(accession)) {
      result = AccessionType.UNIPARC; // uniparc
    } else if (ProteinAccessionPattern.isEnsemblAccession(accession)) {
      result = AccessionType.ENSEMBL; // ENSEMBL
    } else  if (ProteinAccessionPattern.isRefseqAccession(accession)) {
      result = AccessionType.REFSEQ; // NCBI
    }
    return result;
  }

  /**
   * Returns various details for the given protein (f.e. name,
   * sequence).
   * @param accessions The protein's accession.
   * @return A Protein object containing the additional information.
   * @throws Exception error when retrieving protein accession
   */
  public HashMap<String, Protein> getProteinDetails(Collection<String> accessions) throws Exception {
    HashMap<AccessionType, ArrayList<String>> sortedAccessions = new HashMap<>();
    for (String accession : accessions) {
      AccessionType accType = getAccessionType(accession);
      if (!sortedAccessions.containsKey(accType)) {
        sortedAccessions.put(accType, new ArrayList<>());
      }
      accession = accession.replaceAll("\\.\\d+$", "");
      sortedAccessions.get(accType).add(accession);
    }
    HashMap<String, Protein> proteins = new HashMap<>();
    for (AccessionType accType : sortedAccessions.keySet()) {
      switch(accType) {
        case UNIPROT_ACC:
          proteins.putAll(getUniProtDetails(sortedAccessions.get(accType)));
          break;
        case UNIPROT_ID:
          proteins.putAll(getUniProtDetails(sortedAccessions.get(accType)));
          break;
        case ENSEMBL:
          proteins.putAll(getEnsemblDetails(sortedAccessions.get(accType)));
          break;
        case REFSEQ:
          proteins.putAll(getNcbiDetails(sortedAccessions.get(accType)));
          break;
      }
    }
    for (String accession : accessions) {
      accession = accession.replaceAll("\\.\\d+$", ""); // remove any version information from the accession
      if (!proteins.containsKey(accession)) {
        Protein p = new Protein(accession); // add empty protein objects for all proteins that could not be retrieved
        p.setStatus(STATUS.UNKNOWN);
        proteins.put(accession, p);
      }
    }
    return proteins;
  }

  public String getUniprotProteinByBlast(String sequence, float thershold){
    // Create UniProt blast service
    ServiceFactory serviceFactoryInstance = Client.getServiceFactoryInstance();
    UniProtBlastService uniProtBlastService = serviceFactoryInstance.getUniProtBlastService();
    uniProtBlastService.start();
    // Create a blast input with a Database and sequence
    BlastInput input = new BlastInput.Builder(DatabaseOption.TREMBL, sequence).build();
    // Submitting the input to the service will return a completable future
    CompletableFuture<BlastResult<UniProtHit>> resultFuture = uniProtBlastService.runBlast(input);

    BlastResult blastResult = null;
    String resultId = null;
    try {
      blastResult = resultFuture.get();
      logger.info("Number of blast hits: " + blastResult.getNumberOfHits());
      for (Object obejct : blastResult.hits()) {
        UniProtHit hit = (UniProtHit) obejct;
        System.out.println(((UniProtHit)hit).getSummary().getEntryAc() + "\t" +
                ((UniProtHit)hit).getEntry().getPrimaryUniProtAccession().getValue());
        Optional<Alignment> hitResult = hit.getSummary().getAlignments().stream().filter(x -> (x.getIdentity() >= thershold)).findFirst();
        if(hitResult.isPresent()){
          resultId = hit.getEntry().getUniProtId().toString();
        }
      }
    } catch (InterruptedException e) {

    } catch (ExecutionException e) {
      e.printStackTrace();
    } finally {
      uniProtBlastService.stop();
    }
    return resultId;
  }

  /**
   * Returns the details for the given identifiers in
   * a HashMap with the identifier's key or GI number as key and a Protein
   * object holding the details as value.
   *
   * @param accessions The accession to get the details for.
   * @return A Protein object containing the protein's details.
   * @throws Exception any problems getting NCBI details
   */
  private HashMap<String, Protein> getNcbiDetails(Collection<String> accessions) throws Exception {
    StringBuilder query = new StringBuilder();
    for (String accession : accessions) {
      query.append((query.length() > 0) ? "," : "").append(accession);
    }
    String fastas = getPage(String.format(DETAILS_QUERY.NCBI_FASTA.getQueryString(), query.toString()));
    String[] lines = fastas.split(EOL);
    StringBuilder fasta = new StringBuilder(); // parse the fasta entries
    HashMap<String, Protein> proteins = new HashMap<>();
    for (String line : lines) {
      if (fasta.length() > 0 && line.startsWith(">")) {
        Protein protein = convertNcbiFastaToProtein(fasta.toString());
        if (protein != null) {
          String accessionNoVersion = protein.getAccession().replaceAll("`.\\d+$", "");
          proteins.put(accessionNoVersion, protein);
        }
        fasta = new StringBuilder();
      }
      fasta.append(line).append(EOL);
    }
    if (!StringUtils.isEmpty(fasta.toString())) {
      Protein protein = convertNcbiFastaToProtein(fasta.toString());
      if (protein != null)
        proteins.put(protein.getAccession(), protein);
    }
    return proteins;
  }

  /**
   * Converts an NCBI gi fasta entry to a Protein object.
   * @param fasta The fasta to convert.
   * @return Protein Returns the converter Protein object or null in case nothing was found.
   * @throws Exception any problems mapping to a protein
   */
  private Protein convertNcbiFastaToProtein(String fasta) throws Exception {
    if (fasta.trim().length()==0 || (!fasta.startsWith(">") && fasta.contains("Nothing has been found"))) {
      return null;
    }
    String header = fasta.substring(0, fasta.indexOf(EOL)); // only use the first line
    String sequence = fasta.substring(fasta.indexOf(EOL) + 1);
    Pattern pat = Pattern.compile(">(\\S+)\\s([^\\n]+)"); // extract the protein name
    Matcher matcher = pat.matcher(header);
    if (!matcher.find()) {
      throw new Exception("Unexpected fasta format encountered:\n" + fasta);
    }
    String accession = matcher.group(1);
    String name = matcher.group(2).trim();
    String accessionVersion = null;
    if (accession != null) {
      int index = accession.lastIndexOf('.');
      if (index != -1) {
        accessionVersion = accession.substring(index + 1);
      }
      accession = accession.replaceAll("\\.\\d+", ""); // remove the version info from the accession
    }
    Protein protein = new Protein(accession);
    protein.setName(name);
    protein.setSequenceString(sequence);
    protein.setProperty(PROPERTY.SOURCE, "NCBI accession");
    if (accessionVersion != null) {
      protein.setProperty(PROPERTY.ACCESSION_VERSION, accessionVersion);
    }
    logger.info("Protein: " + protein.getName() + " " + protein.getAccession() + " " + protein.getSequenceString());
    return protein;
  }

  /**
   * Retrieves the protein details from UniProt from the
   * given UniProt accessions. Returns null in case nothing
   * was retrieved. <br />
   * <b>Warning:</b> The function first tries to retrieve the protein details
   * expecting them to be accessions. Only if no results are retrieved that way
   * a second request is send interpreting the passed strings as ids. In case
   * accessions and ids are mixed, only the proteins identified through accessions
   * will be returned.
   * @param accessions The UniProt accessions of the proteins.
   * @return A Collection of Protein objects containing the proteins' details or null if the accession doesn't exist
   * @throws Exception In case something went wrong
   */
  private Map<String, Protein> getUniProtDetails(Collection<String> accessions) throws Exception {
    StringBuilder query = new StringBuilder();
    Boolean usingAccession = true;
    List<String> isoforms = new ArrayList<>();
    for (String accession : accessions) {
      if(accession.contains("-") || accession.contains("-")){
        query.append((query.length() > 1) ? "%20or%20" : "").append("sequence:").append(accession);
        isoforms.add(accession);
      }else{
        query.append((query.length() > 1) ? "%20or%20" : "").append("accession:").append(accession);
      }
    }
    String url = String.format(DETAILS_QUERY.UNIPROT.getQueryString(), query.toString());
    String page = getPage(url);
    if ("".equals(page.trim())) {
      usingAccession = false;
      page = getPage(String.format(DETAILS_QUERY.UNIPROT.getQueryString(), query.toString().replace("accession:", "mnemonic:")));
    }
    if(page.equals("")){
      return Collections.emptyMap();
    }
    String[] lines = page.split(EOL);
    if (lines.length < 2) {
      return Collections.emptyMap();
    }
    Map<String,String> isoformMapSequence = new HashMap<>();
    if(isoforms.size() > 0)
      isoformMapSequence = getFastaSequencePage(isoforms);
    HashMap<String, Integer> fieldIndex = new HashMap<>();
    HashMap<String, Protein> proteins = new HashMap<>();
    for (int lineIndex = 0; lineIndex < lines.length; lineIndex++) {
      String[] fields = lines[lineIndex].split(TAB);
      if (lineIndex == 0) {  // if it's the first line build the field index
        fieldIndex.put("Entry", 0); //hard coding due to Uniprot header bug [help #132525]
        fieldIndex.put("Protein names", 1);
        fieldIndex.put("Sequence", 2);
        fieldIndex.put("Status", 3);
        fieldIndex.put("Entry name", 4);
        fieldIndex.put("Organism ID", 5);
        /*for (int i = 0; i < fields.length; i++) {
          fieldIndex.put(fields[i], i); // re-instate when Uniprot fixed [help #132525]
        }*/
        if (!fieldIndex.containsKey("Entry") || !fieldIndex.containsKey("Protein names") ||
            !fieldIndex.containsKey("Sequence") || !fieldIndex.containsKey("Status") ||
            !fieldIndex.containsKey("Entry name")) {
          throw new Exception("Unexpected UniProt response retrieved."); // ensure the required fields were found
        }
        continue;
      }
      if (fields.length < 1) {
        continue;
      }
      if (fields[1].startsWith("Merged into")) { // check if the protein was demerged
        String newAccession = fields[1].substring(12, 18);
        Protein replacingP = proteins.get(newAccession);
        Protein p = new Protein(fields[0]);
        p.setStatus(STATUS.MERGED);
        p.setProperty(PROPERTY.STATUS_INFO, fields[1]);
        p.getReplacingProteins().add(replacingP);
        proteins.put(fields[0], p);
      } else if (fields[1].startsWith("Demerged into")) {
        String accession = fields[ fieldIndex.get("Entry") ];
        if (!accessions.contains(accession) && fields.length >= 5) {
          accession = fields[4];
        }
        Protein p = new Protein(accession);
        p.setStatus(STATUS.DEMERGED);
        p.setProperty(PROPERTY.STATUS_INFO, fields[1]);
        ArrayList<String> replacingAccessions = extractUniprotReplacingAccessions(fields[1]);
        for (String acc : replacingAccessions) {
          if (proteins.containsKey(acc)) {
            p.getReplacingProteins().add(proteins.get(acc));
          }
        }
        if (!proteins.containsKey(accession)) {
          proteins.put(accession, p);
        }
      }
      else if (fields[1].startsWith("Deleted")) {
        Protein p = new Protein(fields[fieldIndex.get("Entry")]);
        p.setStatus(STATUS.DELETED);
        proteins.put(fields[0], p);
      }
      else {
        if (fields.length != fieldIndex.size()) {   // make sure the line is in the expected format
          throw new Exception("Unexpected UniProt answer retrieved. Line has a different number of fields than defined in the header: <" + lines[lineIndex] + ">");
        }
        String accession = (usingAccession) ? fields[fieldIndex.get("Entry")] : fields[fieldIndex.get("Entry name")];
        String accessionIsoform = checkIsoformSequence(accessions, accession);
        Protein p = new Protein(accession);
        p.setName(fields[fieldIndex.get("Protein names")]);
        p.setProperty(PROPERTY.SOURCE, (fields[fieldIndex.get("Status")].equals("reviewed")) ? "UniProt/Swiss-Prot" : "UniProt/TrEMBL");
        p.setStatus(STATUS.ACTIVE);
        String sequence = fields[fieldIndex.get("Sequence")];
        sequence = sequence.replaceAll("\\s", "");
        p.setSequenceString(sequence);
        String organismID = fields[fieldIndex.get("Organism ID")];
        p.setOrganismId(organismID);
        proteins.put(accession, p);
        if(accessionIsoform != null){
          p.setSequenceString((isoformMapSequence.getOrDefault(accessionIsoform, sequence)));
          p.setAccession(accessionIsoform);
          proteins.put(accessionIsoform,p);
        }
      }
    }
    return proteins;
  }

  /**
   * Extracts the replacing accession from a UniProt status
   * line Demerged into XXXXX, XXXXX and XXXXX.
   * @param string
   * @return
   */
  private ArrayList<String> extractUniprotReplacingAccessions(String string) {
    ArrayList<String> extractedAccessions = new ArrayList<>(2);
    String accessions = string.substring(14); // get the accessions
    char nextChar;
    do {
      String accession = accessions.substring(0, 6);
      nextChar = accessions.charAt(6);
      extractedAccessions.add(accession);
      if (nextChar == ' ') {
        accessions = accessions.substring(11);
      }
      else if (nextChar != '.') {
        accessions = accessions.substring(8);
      }
    } while (nextChar != '.' && accessions.length() >= 7);
    return extractedAccessions;
  }

  /**
   * Returns the details for the given Ensembl identifiers in
   * a HashMap with the Ensembl identifier as key and a Protein
   * object holding the details as value. Uses the UniProt
   * mapping service to map ensembl entries to UniProt accessions.
   * Then the call is passed on to getUniProtDetails.
   *
   * @param accessions The Ensembl accession to get the name for.
   * @return A Protein object containing the protein's details.
   * @throws Exception
   */
  private HashMap<String, Protein> getEnsemblDetails(Collection<String> accessions) throws Exception {
    // use the uniprot mapping service to convert the ensembl accessions
    String query = "";
    for (String acc : accessions) {
      query += (query.length() > 0 ? "," : "") + acc;
    }
    String page = getPage(String.format(MAPPING_QUERY.ENSEMBL_TO_UNIPROT.getQueryString(), query));
    String[] lines = page.split(EOL);
    HashMap<String, String> uniprotToEnsemblMapping = new HashMap<>();
    HashMap<String, Integer> fieldMapping = new HashMap<>();
    for (int i = 0; i < lines.length; i++) {
      String[] fields = lines[i].split(TAB);
      if (i == 0) {
        for (int j = 0; j < fields.length; j++)
          fieldMapping.put(fields[j], j);
        if (!fieldMapping.containsKey("To") || !fieldMapping.containsKey("From")) {
          throw new Exception("Unexpected response retrieved from UniProt mapping service.");
        }
        continue;
      }
      uniprotToEnsemblMapping.put(fields[fieldMapping.get("To")], fields[fieldMapping.get("From")]);
    }
    Map<String, Protein> proteins = new HashMap<>();
    if(uniprotToEnsemblMapping.size() > 0 ) {
      proteins = getUniProtDetails(uniprotToEnsemblMapping.keySet());
    }
    HashMap<String, Protein> ensemblProteins = new HashMap<>();
    for (Protein p : proteins.values()) {
      String ensemblAcc = uniprotToEnsemblMapping.get(p.getAccession());
      if (ensemblAcc == null) {
        continue;
      }
      p.setAccession(ensemblAcc);
      ensemblProteins.put(ensemblAcc, p);
    }
    return ensemblProteins;
  }

  /**
   * Gets the page from the given address. Returns the
   * retrieved page as a string.
   *
   * @param urlString The address of the resource to retrieve.
   * @return The page as a String
   * @throws Exception Thrown on any problem.
   */
  private String getPage(String urlString) throws Exception {
    URL url = new URL(urlString);
    HttpURLConnection connection = (HttpURLConnection) url.openConnection();
    connection.connect();
    BufferedReader in = null;
    StringBuilder page = new StringBuilder();
    try {
      in = new BufferedReader(new InputStreamReader(connection.getInputStream()));
      String line;
      while ((line = in.readLine()) != null) {
        page.append(line);
        page.append("\n");
      }
    } catch (IOException ioe) {
      logger.warn("Failed to read web page");
    } finally {
      if (in != null) {
        in.close();
      }
    }
    return page.toString();
  }

  /**
   * The only way to retrieve the isoforms sequences is one-by-one. The current method
   * get a Map of isoform accessions and retrieve
   * @param isoforms a list of isoforms
   * @return a map of isoform accessions and its sequence
   */
  private Map<String, String> getFastaSequencePage(List<String> isoforms) throws Exception {
    Map<String, String> isoformSequences = new HashMap<>();
    for(String accession : isoforms){
      String query = accession + ".fasta";
      String url = String.format(DETAILS_QUERY.UNIPROT_FASTA.getQueryString(), query);
      String page = getPage(url);
      String[] lines = page.split(EOL);
      if (!page.equals("") && lines.length > 2) { // if there's only one line or the page was empty no protein names were retrieved
        StringBuilder sequence = new StringBuilder();
        for(int i = 1; i < lines.length; i++){
          sequence.append(lines[i]);
        }
        isoformSequences.put(accession, sequence.toString());
      }
    }
    return isoformSequences;
  }

  private String checkIsoformSequence(Collection<String> accessions, String accession){
    String isoform = null;
    for(String accessionString : accessions){
      if(accessionString.contains(accession)){
        isoform = accessionString;
        break;
      }
    }
    return isoform;
  }

  /**
   * From a protein Ids check all the services (Uniprot, EntrezGene) to know if the web services
   * are available.
   *
   * @return true if Uniprot is available, false otherwise
   */
  public static boolean checkUniprotService(){
    boolean servicesAvailable = false;
    Collection<String> uniprotIds = new ArrayList<>();
    uniprotIds.add(UNIPROT_PROTEIN);
    ProteinDetailFetcher proteinDetailFetcher = new ProteinDetailFetcher();
    try{
      HashMap<String, Protein> uniprotDetails = proteinDetailFetcher.getProteinDetails(uniprotIds);
      if( uniprotDetails != null && !uniprotDetails.isEmpty()) {
        servicesAvailable = true;
      }
    } catch (Exception e) {
      logger.error("There is no internet connection", e);
    }
    return servicesAvailable;
  }
}