package uk.ac.ebi.pride.tools.protein_details_fetcher;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.net.HttpURLConnection;
import java.net.URL;
import java.util.*;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import org.jdom.Document;
import org.jdom.Element;
import org.jdom.input.SAXBuilder;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import uk.ac.ebi.pride.tools.protein_details_fetcher.model.Protein;
import uk.ac.ebi.pride.tools.protein_details_fetcher.model.Protein.PROPERTY;
import uk.ac.ebi.pride.tools.protein_details_fetcher.model.Protein.STATUS;
import uk.ac.ebi.pride.tools.protein_details_fetcher.util.ProteinAccessionPattern;

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
		NCBI_FASTA("http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=%s&rettype=fasta&tool=protein_details_fetcher"),
		UNIPROT("http://www.uniprot.org/uniprot/?query=%s&format=tab&columns=id,protein%%20names,sequence,reviewed,entry%%20name,organism-id"),
        UNIPROT_FASTA("http://www.uniprot.org/uniprot/%s");
		
		private String queryString;
		
		private DETAILS_QUERY(String formatString) {
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
		ENSEMBL_TO_UNIPROT("http://www.uniprot.org/mapping/?from=ENSEMBL_PRO_ID&to=ACC&query=%s&format=tab"),
		IPI_TO_UNIPROT("http://www.uniprot.org/mapping/?from=P_IPI&to=&query=%s&format=tab");
		
		private String queryString;
		
		private MAPPING_QUERY(String queryString) {
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
	public enum AccessionType{UNIPROT_ACC, UNIPROT_ID, UNIPARC, IPI, REFSEQ, ENSEMBL, GI, ENSEMBL_TRANSCRIPT, accessionType, UNKNOWN}
	
    // query string for the NCBI esummary tool
    public final String ESUMMARY_QUERY_STRING = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?db=protein&id=%s";
    // query string for the NCBI esearch tool
    public final String ESEARCH_QUERY_STRING = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=protein&term=%s";
    // query string to fetch the protein sequence from NCBI
    public final String EFETCH_FORMAT_STRING = "http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=%s&rettype=fasta&tool=pride_inspector";

    private HashMap<String, String> pageBuffer = new HashMap<String, String>();
    // TODO: in case this class should be used in a multi-thread environment, this member variable should be changed to static + concurrentHashMap

    /**
     * Returns the (guessed) accession type for the passed
     * accession. In case the accession is not recognized
     * UNKNOWN is returned.
     * @param accession The accession to guess the type for.
     * @return The accession's type.
     */
    public AccessionType getAccessionType(String accession) {
    	// swissprot accession
    	if (ProteinAccessionPattern.isSwissprotAccession(accession))
            return AccessionType.UNIPROT_ACC;
    	
    	if (ProteinAccessionPattern.isSwissprotEntryName(accession))
    		return AccessionType.UNIPROT_ID;

        // uniparc
        if (ProteinAccessionPattern.isUniparcAccession(accession))
            return AccessionType.UNIPARC;

        // IPI
        if (ProteinAccessionPattern.isIPIAccession(accession))
        	return AccessionType.IPI;

        // ENSEMBL
        if (ProteinAccessionPattern.isEnsemblAccession(accession))
        	return AccessionType.ENSEMBL;

		if(ProteinAccessionPattern.isEnsemblTranscriptAccession(accession))
			return AccessionType.ENSEMBL_TRANSCRIPT;

        // NCBI
        if (ProteinAccessionPattern.isRefseqAccession(accession))
        	return AccessionType.REFSEQ;

        // GI
        if (ProteinAccessionPattern.isGIAccession(accession))
        	return AccessionType.GI;
        
        return AccessionType.UNKNOWN;
    }
    
    /**
     * Returns various details for the given protein (f.e. name,
     * sequence).
     * @param accessions The protein's accession.
     * @return A Protein object containing the additional information.
     * @throws Exception error when retrieving protein accession
     */
    public HashMap<String, Protein> getProteinDetails(Collection<String> accessions) throws Exception {
    	// sort the passed accessions into Lists based on the (guessed) identifier system
    	HashMap<AccessionType, ArrayList<String>> sortedAccessions = new HashMap<AccessionType, ArrayList<String>>();
    	
    	for (String accession : accessions) {
    		// get the accessions type
    		AccessionType accType = getAccessionType(accession);
    		
    		// put the accession in the respective ArrayList
    		if (!sortedAccessions.containsKey(accType))
    			sortedAccessions.put(accType, new ArrayList<String>());
    		
    		// remove any version information from the accession
    		accession = accession.replaceAll("\\.\\d+$", "");
    		sortedAccessions.get(accType).add(accession);
    	}
    	
    	// map the accessions
    	HashMap<String, Protein> proteins = new HashMap<String, Protein>();
    	for (AccessionType accType : sortedAccessions.keySet()) {
    		switch(accType) {
    			case UNIPROT_ACC:
    				proteins.putAll(getUniProtDetails(sortedAccessions.get(accType)));
    				break;
    			case UNIPROT_ID:
    				proteins.putAll(getUniProtDetails(sortedAccessions.get(accType)));
    				break;
    			case IPI:
    				proteins.putAll(getIpiDetails(sortedAccessions.get(accType)));
    				break;
    			case ENSEMBL:
    				proteins.putAll(getEnsemblDetails(sortedAccessions.get(accType)));
    				break;
    			case GI:
    				proteins.putAll(getNcbiDetails(sortedAccessions.get(accType), true));
    				break;
    			case REFSEQ:
    				proteins.putAll(getNcbiDetails(sortedAccessions.get(accType), false));
    				break;
    		}
    	}
    	
    	// add empty protein objects for all proteins that could not be retrieved
    	// and set the status to DELETED
    	for (String accession : accessions) {
    		// remove any version information from the accession
    		accession = accession.replaceAll("\\.\\d+$", "");
    		
    		if (!proteins.containsKey(accession)) {
    			Protein p = new Protein(accession);
    			p.setStatus(STATUS.UNKNOWN);
    			proteins.put(accession, p);
    		}
    	}
    	
        return proteins;
    }
    
    /**
     * Returns the details for the given identifiers in
     * a HashMap with the identifier's key or GI number as key and a Protein
     * object holding the details as value.
     *
     * @param accessions The accession to get the details for.
     * @return A Protein object containing the protein's details.
     * @throws Exception
     */
    private HashMap<String, Protein> getNcbiDetails(Collection<String> accessions, boolean useGiNumber) throws Exception {
    	// create the query string
    	String query = "";
    	
    	for (String accession : accessions)
    		query += ((query.length() > 0) ? "," : "") + accession;
    	
    	// get the IPI fasta entry
        String fastas = getPage(String.format(DETAILS_QUERY.NCBI_FASTA.getQueryString(), query));
        String[] lines = fastas.split(EOL);
        
        // parse the fasta entries
        String fasta = "";
        
        HashMap<String, Protein> proteins = new HashMap<String, Protein>();
        
        for (int lineIndex = 0; lineIndex < lines.length; lineIndex++) {
        	String line = lines[lineIndex];
        	
        	// process the current fasta
        	if (fasta.length() > 0 && line.startsWith(">")) {
        		Protein protein = convertNcbiFastaToProtein(fasta, useGiNumber);
        		if (protein != null) {
        			String accessionNoVersion = protein.getAccession().replaceAll("\\.\\d+$", "");
        			proteins.put(accessionNoVersion, protein);
        		}
                
                fasta = "";
        	}
        	
        	fasta += line + EOL;
        }
        
        if (!"".equals(fasta)) {
        	Protein protein = convertNcbiFastaToProtein(fasta, useGiNumber);
        	
        	if (protein != null)
        		proteins.put(protein.getAccession(), protein);
        }
        
        // TODO: add the protein's status
       proteins = enrichNCBIProteins(proteins);
        
       return proteins;
    }
    
        private HashMap<String, Protein> enrichNCBIProteins(HashMap<String, Protein> proteins) throws Exception {
		// build the query
    	String query = "";
    	
    	for (String accession : proteins.keySet()) {
    		Protein p = proteins.get(accession);    		
    		query += (query.length() > 0 ? "," : "") + p.getProperty(PROPERTY.GI_NUMBER);
    	}
    	
    	// get the result
    	String page = getPage(String.format(ESUMMARY_QUERY_STRING, query));
    	
    	// get the properties
    	HashMap<Integer, HashMap<String, String>> proteinProperties = getNcbiProperties(page);
    	
    	// enrich the existing proteins
    	for (Protein p : proteins.values()) {
    		// get the properties
    		HashMap<String, String> properties = proteinProperties.get(Integer.parseInt(p.getProperty(PROPERTY.GI_NUMBER)));
    		
    		if (properties != null) {
    			// check the status
    			String status = properties.get("Status");
    			
    			if ("live".equals(status))
    				p.setStatus(STATUS.ACTIVE);
    			else if (!"live".equals(status) && (properties.get("ReplacedBy") == null || "".equals(properties.get("ReplacedBy"))))
    				p.setStatus(STATUS.DELETED);
    			else {
    				p.setStatus(STATUS.CHANGED); 
    				
    				String replacedBy = properties.get("ReplacedBy");
    				p.getReplacingProteins().add(new Protein(replacedBy));
    			}
    		}
    	}
    	
    	return proteins;
	}
    


    /**
     * Returns the properties of the given identifier fetched
     * using the NCBI esummary tool and returns them as a HashMap.
     *
     * @param xml The accession to retrieve the properties for.
     * @return A HashMap with the property's name as key and its value as value.
     * @throws Exception Thrown in case something went wrong.
     */
    private HashMap<Integer, HashMap<String, String>> getNcbiProperties(String xml) throws Exception {
        // create the xml object
        SAXBuilder builder = new SAXBuilder();
        Document doc = builder.build(new StringReader(xml));

        // get the document summary element
        Element root = doc.getRootElement();

        if (root == null)
            throw new Exception("Failed to parse NCBI XML snipplet");

        @SuppressWarnings("unchecked")
		List<Element> docSums = root.getChildren("DocSum");
        HashMap<Integer, HashMap<String, String>> elementProperties = new HashMap<Integer, HashMap<String,String>>();
        
        for (Element docSum : docSums) {
	        if (docSum == null)
	            throw new Exception("Failed to parse NCBI XML snipplet");
	
	        // get all the items
	        @SuppressWarnings("unchecked")
			List<Element> items = docSum.getChildren("Item");
	
	        // initialize the return variable
	        HashMap<String, String> properties = new HashMap<String, String>();
	
	        // parse the items
	        for (Element item : items) {
	            properties.put(item.getAttributeValue("Name"), item.getValue());
	        }
	        
	        // get the id
	        String gi = docSum.getChildText("Id");
	        
	        elementProperties.put(Integer.parseInt(gi), properties);
        }

        return elementProperties;
    }

	/**
     * Converts an NCBI gi fasta entry to a Protein object.
     * @param fasta The fasta to convert.
     * @param useGi Indicates whether the GI number or the source accession should be set as the Protein's accession.
     * @return Protein Returns the converter Protein object or null in case nothing was found.
     * @throws Exception
     */
    private Protein convertNcbiFastaToProtein(String fasta, boolean useGi) throws Exception {
    	// make sure something was found
    	if (fasta.trim().length()==0 || (!fasta.startsWith(">") && fasta.contains("Nothing has been found"))) {
    		return null;
    	}    		
    	
    	// only use the first line
        String header = fasta.substring(0, fasta.indexOf(EOL));
        
        // get the sequence
        String sequence = fasta.substring(fasta.indexOf(EOL) + 1);
        // remove all whitespaces
        sequence = sequence.replaceAll("\\s", "");

        // extract the protein name
        Pattern pat = Pattern.compile(">[^|]+\\|([^|]+)\\|([^|]+)\\|([^|]*)\\|(.*)");

        Matcher matcher = pat.matcher(header);

        // make sure it matches
        if (!matcher.find())
            throw new Exception("Unexpected fasta format encountered:\n" + fasta);

        String gi = matcher.group(1);
        String source = matcher.group(2);
        String accession = matcher.group(3);
        String name = matcher.group(4).trim();
        
        String accessionVersion = null;
        
        if ("ref".equals(source))
        	source = "RefSeq";
        if (accession != null) {
        	int index = accession.lastIndexOf('.');
        	if (index != -1)
        		accessionVersion = accession.substring(index + 1);
        	// remove the version info from the accession
        	accession = accession.replaceAll("\\.\\d+", "");
        }
        
        
        // create the protein object
        Protein protein = new Protein((useGi || accession == null) ? gi : accession);
        protein.setName(name);
        protein.setSequenceString(sequence);
        protein.setProperty(PROPERTY.SOURCE, useGi ? "NCBI gi" : source);
        protein.setProperty(PROPERTY.GI_NUMBER, gi);
        if (accessionVersion != null)
        	protein.setProperty(PROPERTY.ACCESSION_VERSION, accessionVersion);
        
        return protein;
    }
    
    /**
     * Returns the details for the given IPI identifiers in
     * a HashMap with the IPI identifier as key and a Protein
     * object holding the details as value.  Uses the UniProt
     * mapping service to map IPI entries to UniProt accessions.
     * Then the call is passed on to getUniProtDetails.
     *
     * @param accessions The IPI accession to get the name for.
     * @return A Protein object containing the protein's details.
     * @throws Exception
     */
    private HashMap<String, Protein> getIpiDetails(Collection<String> accessions) throws Exception {
        // use the uniprot mapping service to convert the ipi accessions
        String query = "";

        for (String acc : accessions)
            query += (query.length() > 0 ? "," : "") + acc;

        // map the accessions
        String page = getPage(String.format(MAPPING_QUERY.IPI_TO_UNIPROT.getQueryString(), query));
        String[] lines = page.split(EOL);

        HashMap<String, String> uniprotToIpiMapping = new HashMap<String, String>();
        HashMap<String, Integer> fieldMapping = new HashMap<String, Integer>();

        if (lines.length<2) {
            // no mappings found
            return new HashMap<String, Protein>();
        }

        for (int i = 0; i < lines.length; i++) {
            String[] fields = lines[i].split(TAB);

            if (i == 0) {
                for (int j = 0; j < fields.length; j++)
                    fieldMapping.put(fields[j], j);

                // make sure the required fields are there
                if (!fieldMapping.containsKey("To") || !fieldMapping.containsKey("From"))
                    throw new Exception("Unexpected response retrieved from UniProt mapping service.");

                continue;
            }

            // save the mapping
            uniprotToIpiMapping.put(fields[fieldMapping.get("To")], fields[fieldMapping.get("From")]);
        }

        // get the UniProt mappings
        Map<String, Protein> proteins = new HashMap<String, Protein>();
        if(uniprotToIpiMapping.size() > 0)
            proteins = getUniProtDetails(uniprotToIpiMapping.keySet());

        // create a new HashMap changing the UniProt accessions to IPI accessions
        HashMap<String, Protein> ipiProteins = new HashMap<String, Protein>();

        for (Protein p : proteins.values()) {
            String ipiAcc = uniprotToIpiMapping.get(p.getAccession());
            if (ipiAcc != null) {
                p.setAccession(ipiAcc);
                ipiProteins.put(ipiAcc, p);
            }
        }
        return ipiProteins;
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
    	// build the query string for the accessions

    	String query = "";
    	Boolean usingAccession = true;

        List<String> isoforms = new ArrayList<String>();

    	for (String accession : accessions) {
    		if(accession.contains("-") || accession.contains("-")){
                query += ((query.length() > 1) ? "%20or%20" : "") + "sequence:" + accession;
                isoforms.add(accession);
            }else{
                query += ((query.length() > 1) ? "%20or%20" : "") + "accession:" + accession;
            }
        }

    	String url = String.format(DETAILS_QUERY.UNIPROT.getQueryString(), query);
    	//System.out.println(url);

    	// get the page
        String page = getPage(url);
        if ("".equals(page.trim())) {
        	// try with uniprot ids
        	usingAccession = false;
            page = getPage(String.format(DETAILS_QUERY.UNIPROT.getQueryString(), query.replace("accession:", "mnemonic:")));
        }
        if(page.equals("")){
            return Collections.emptyMap();
        }

        String[] lines = page.split(EOL);

        if (lines.length < 2) {
            return Collections.emptyMap();
        }

        //retrieve the isoform sequences for each isoforms ID
        Map<String,String> isoformMapSequence = new HashMap<String, String>();
        if(isoforms.size() > 0)
            isoformMapSequence = getFastaSequencePage(isoforms);

        HashMap<String, Integer> fieldIndex = new HashMap<String, Integer>();
        HashMap<String, Protein> proteins = new HashMap<String, Protein>();
        
        // create the retrieved proteins
        for (int lineIndex = 0; lineIndex < lines.length; lineIndex++) {
        	 String[] fields = lines[lineIndex].split(TAB);
        	// if it's the first line build the field index
        	if (lineIndex == 0) {
        		for (int i = 0; i < fields.length; i++)
        			fieldIndex.put(fields[i], i);
        		
        		// make sure the required fields were found
        		if (!fieldIndex.containsKey("Entry") || !fieldIndex.containsKey("Protein names") ||
        			!fieldIndex.containsKey("Sequence") || !fieldIndex.containsKey("Status") ||
        			!fieldIndex.containsKey("Entry name"))
        			throw new Exception("Unexpected UniProt response retrieved.");
        		
        		continue;
        	}

        	if (fields.length < 1)
        		continue;
        	
        	// check if the protein was demerged
        	if (fields[1].startsWith("Merged into")) {
        		// extract the new accession
        		String newAccession = fields[1].substring(12, 18);
        		// get the protein
        		Protein replacingP = proteins.get(newAccession);
        		// create the current protein
        		Protein p = new Protein(fields[0]);
        		
        		// set the status
        		p.setStatus(STATUS.MERGED);
        		p.setProperty(PROPERTY.STATUS_INFO, fields[1]);
        		p.getReplacingProteins().add(replacingP);
        		
        		proteins.put(fields[0], p);
        	}
        	else if (fields[1].startsWith("Demerged into")) {
        		// check if the accession or id should be used (not clear for ids)
        		String accession = fields[ fieldIndex.get("Entry") ];
        		
        		if (!accessions.contains(accession) && fields.length >= 5) {
        			accession = fields[4];
        		}
        		
        		// create a new protein
        		Protein p = new Protein(accession);
        		
        		// set the status
        		p.setStatus(STATUS.DEMERGED);
        		p.setProperty(PROPERTY.STATUS_INFO, fields[1]);
        		
        		// add the replacing proteins
        		ArrayList<String> replacingAccessions = extractUniprotReplacingAccessions(fields[1]);
        		for (String acc : replacingAccessions) {
        			if (proteins.containsKey(acc))
        				p.getReplacingProteins().add(proteins.get(acc));
        		}	
        		
        		if (!proteins.containsKey(accession))
        			proteins.put(accession, p);
        	}
        	else if (fields[1].startsWith("Deleted")) {
        		// create a new protein
        		Protein p = new Protein(fields[fieldIndex.get("Entry")]);
        		
        		p.setStatus(STATUS.DELETED);
        		
        		proteins.put(fields[0], p);
        	}
        	else {
        		// make sure the line is in the expected format
            	if (fields.length != fieldIndex.size())
            		throw new Exception("Unexpected UniProt answer retrieved. Line has a different number of fields than defined in the header: <" + lines[lineIndex] + ">");
            	
            	// create the protein object
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
                    p.setSequenceString((isoformMapSequence.containsKey(accessionIsoform)?isoformMapSequence.get(accessionIsoform):sequence));
                    p.setAccession(accessionIsoform);
                    proteins.put(accessionIsoform,p);
                }
        	}
        }
        
        // return the proteins
        return proteins;
    }
    
    /**
     * Extracts the replacing accession from a UniProt status
     * line Demerged into XXXXX, XXXXX and XXXXX.
     * @param string
     * @return
     */
    private ArrayList<String> extractUniprotReplacingAccessions(String string) {
    	ArrayList<String> extractedAccessions = new ArrayList<String>(2);
    	
    	// get the accessions
    	String accessions = string.substring(14);
    	char nextChar;
    	
    	do {
    		String accession = accessions.substring(0, 6);
    		nextChar = accessions.charAt(6);
    		
    		extractedAccessions.add(accession);
    		
    		if (nextChar == ' ')
    			accessions = accessions.substring(11);
    		else if (nextChar != '.')
    			accessions = accessions.substring(8);
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
    	
    	for (String acc : accessions)
    		query += (query.length() > 0 ? "," : "") + acc;
    	
    	// map the accessions
    	String page = getPage(String.format(MAPPING_QUERY.ENSEMBL_TO_UNIPROT.getQueryString(), query));
    	String[] lines = page.split(EOL);
    	
    	HashMap<String, String> uniprotToEnsemblMapping = new HashMap<String, String>();
    	HashMap<String, Integer> fieldMapping = new HashMap<String, Integer>();
    	
    	for (int i = 0; i < lines.length; i++) {
    		String[] fields = lines[i].split(TAB);
    		
    		if (i == 0) {
    			for (int j = 0; j < fields.length; j++)
    				fieldMapping.put(fields[j], j);
    			
    			// make sure the required fields are there
    			if (!fieldMapping.containsKey("To") || !fieldMapping.containsKey("From"))
    				throw new Exception("Unexpected response retrieved from UniProt mapping service.");
    			
    			continue;
    		}
    		
    		// save the mapping
    		uniprotToEnsemblMapping.put(fields[fieldMapping.get("To")], fields[fieldMapping.get("From")]);
    	}
    	
    	// get the UniProt mappings
        Map<String, Protein> proteins = new HashMap<String, Protein>();
        if(uniprotToEnsemblMapping.size() > 0 )
    	   proteins = getUniProtDetails(uniprotToEnsemblMapping.keySet());
    	
    	// create a new HashMap changing the UniProt accessions to ENSEMBL accessions
    	HashMap<String, Protein> ensemblProteins = new HashMap<String, Protein>();
    	
    	for (Protein p : proteins.values()) {
    		String ensemblAcc = uniprotToEnsemblMapping.get(p.getAccession());
    		
    		if (ensemblAcc == null)
    			continue;
    		
    		p.setAccession(ensemblAcc);
    		
    		ensemblProteins.put(ensemblAcc, p);
    	}
    	
    	return ensemblProteins;
    }
    
    /**
     * Retrieves the page at the given location and parses the
     * retrieved tab-delimited table.
     * @param location
     * @return
     * @throws Exception
     */
    private ArrayList<HashMap<String, String>> getTable(String location) throws Exception {
    	// get the page
    	String page = getPage(location);
    	
    	HashMap<String, Integer> header = new HashMap<String, Integer>();
    	ArrayList<HashMap<String, String>> table = new ArrayList<HashMap<String,String>>();
    	
    	String[] lines = page.split(EOL);
    	
    	for (int i = 0; i < lines.length; i++) {
    		String[] fields = lines[i].split(TAB);
    		
    		// build the header if it's the first line
    		if (i == 0) {
    			for (int j = 0; j < fields.length; j++)
    				header.put(fields[j], j);
    			
    			continue;
    		}
    		
    		// parse the line
    		HashMap<String, String> line = new HashMap<String, String>(header.size());
    		
    		for (String headerField : header.keySet()) {
    			int nFieldPot = header.get(headerField);
    			
    			if (nFieldPot >= 0 && nFieldPot < fields.length)
    				line.put(headerField, fields[nFieldPot]);
    		}
    		
    		table.add(line);
    	}
    	
    	return table;
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
        // check if the page is cached
        if (pageBuffer.containsKey(urlString))
            return pageBuffer.get(urlString);

        // create the url
        URL url = new URL(urlString);

        // send the request
        HttpURLConnection connection = (HttpURLConnection) url.openConnection();

        connection.connect();

        // get the page
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
     * @param isoforms
     * @return
     */
    private Map<String, String> getFastaSequencePage(List<String> isoforms) throws Exception {

        Map<String, String> isoformSequences = new HashMap<String, String>();

        for(String accession: isoforms){
           String query = accession + ".fasta";
            String url = String.format(DETAILS_QUERY.UNIPROT_FASTA.getQueryString(), query);
            // System.out.println(url);
            // get the page
            String page = getPage(url);
            String[] lines = page.split(EOL);

            // if there's only one line or the page was empty no protein names were retrieved
            if (!page.equals("") && lines.length > 2) {
                String sequence = "";
                for(int i = 1; i < lines.length; i++){
                    sequence = sequence + lines[i];
                }
               isoformSequences.put(accession,sequence);
            }
        }
        return isoformSequences;
    }

    private String checkIsoformSequence(Collection<String> accessions, String accession){
        String isoform = null;
        for(String accessionString: accessions){
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
     * @return boolean
     */
    public static boolean checkUniprotService(){
        boolean servicesAvailable = false;
        Collection<String> uniprotIds = new ArrayList<String>();
        uniprotIds.add(UNIPROT_PROTEIN);
        ProteinDetailFetcher proteinDetailFetcher = new ProteinDetailFetcher();
        try{
        HashMap<String, Protein> uniprotDetails = proteinDetailFetcher.getProteinDetails(uniprotIds);
        if( uniprotDetails != null && !uniprotDetails.isEmpty())
            return true;
        }catch (Exception e) {
            logger.error("There is no internet connection", e);
            return false;
        }
        return servicesAvailable;
    }
}