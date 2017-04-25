package uk.ac.ebi.pride.tools.protein_details_fetcher.model;

import org.apache.commons.lang3.StringUtils;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

/**
 * <p>Object to store protein related details</p>
 * User: rwang
 * Date: 08/06/11
 * Time: 16:42
 */

public class Protein implements Serializable {
    public enum STATUS {UNKNOWN, ACTIVE, DELETED, CHANGED, DEMERGED, MERGED, ERROR};
    /**
     * Describes the various properties a protein
     * object might have.
     * @author jg
     *
     */
    public enum PROPERTY{STATUS_INFO, SOURCE, ERROR_STRING, ACCESSION_VERSION};

    /**
     * The protein's properties.
     */
    private HashMap<PROPERTY, String> properties = new HashMap<>();
    /**
     * Human readable name for the protein
     */
    private String name = null;
    /**
     * Protein accession
     */
    private String accession = null;
    /**
     * Protein sequence
     */
    private String sequenceString = null;
    /**
     * The protein's status.
     */
    private STATUS status = STATUS.UNKNOWN;
    /**
     * If a protein was DEMERGED or MERGED this
     * can hold the replacing proteins.
     */

    private String organismId = null;

    private List<Protein> replacingProteins = new ArrayList<>();


    public Protein(String accession) {
        if (accession == null || "".equals(accession.trim())) {
            throw new IllegalArgumentException("Protein accession cannot be NULL");
        }
        this.accession = accession;
    }

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    public String getAccession() {
        return accession;
    }

    public void setAccession(String accession) {
        this.accession = accession;
    }

    public String getSequenceString() {
        return sequenceString;
    }

    public void setSequenceString(String sequence) {
        if (!StringUtils.isEmpty(sequence)) {
            this.sequenceString = sequence.toUpperCase();
        } else {
            this.sequenceString = null;
        }
    }

    public STATUS getStatus() {
        return status;
    }

    public void setStatus(STATUS status) {
        this.status = status;
    }

    public List<Protein> getReplacingProteins() {
        return replacingProteins;
    }

    public String getSubSequenceString(int start, int stop) {
        String result = null;
        if (sequenceString != null && sequenceString.length() >= stop && start >= 1 && start <= stop) {
            result =  this.sequenceString.substring(start - 1, stop);
        }
        return result;
    }

    /**
     * To check whether a given sequence is a sub sequence of the protein
     * @param subSeq    given sequence
     * @param start start position
     * @param stop  stop position
     * @return  boolean true means exist
     */
    public boolean hasSubSequenceString(String subSeq, int start, int stop) {
        String targetSeq = getSubSequenceString(start, stop);
        return targetSeq != null && subSeq != null && targetSeq.equals(subSeq.toUpperCase());
    }

    /**
     * To check whether a given sequence is a sub sequence of the protein
     * @param subSeq    given sequence
     * @return  boolean true means exist
     */
    public boolean hasSubSequenceString(String subSeq) {
        return sequenceString != null && subSeq != null && sequenceString.contains(subSeq);
    }


    /**
     * Search a given sub sequence within the protein sequence
     * return a set of starting positions which matches the given sub sequence
     *
     * @param subSeq    given sub sequence
     * @return Return a Set of positions
     */
    public Set<Integer> searchStartingPosition(String subSeq) {
        Set<Integer> position = new HashSet<>();
        if (sequenceString != null && subSeq != null) {
            int previousIndex = -1;
            int index;
            while((index = (previousIndex == -1 ?
                sequenceString.indexOf(subSeq) :
                sequenceString.indexOf(subSeq, previousIndex + 1)))
                > -1) {
                position.add(index);
                previousIndex = index;
            }
        }
        return position;
    }

    /**
     * Sets the given property.
     * @param property PROPERTY with details of the protein
     * @param value The value of the PROPERTY
     */
    public void setProperty(PROPERTY property, String value) {
        properties.put(property, value);
    }

    /**
     * Returns the value of the specified property
     * or null in case it wasn't set.
     * @param property The property to fetch.
     * @return The property's value as String or null if it wasn't set.
     */
    public String getProperty(PROPERTY property) {
        return properties.get(property);
    }

    public boolean hasProperty(PROPERTY property) {
        return properties.containsKey(property);
    }

    /**
     * Return the Id of the organism
     * @return The id of the organism
     */
    public String getOrganismId() {
        return organismId;
    }

    /**
     * Set the organism
     *
     * @param organismId Organism ID
     */
    public void setOrganismId(String organismId) {
        this.organismId = organismId;
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Protein)) return false;

        Protein protein = (Protein) o;

        if (!accession.equals(protein.accession)) return false;
        if (name != null ? !name.equals(protein.name) : protein.name != null) return false;
        if (sequenceString != null ? !sequenceString.equals(protein.sequenceString) : protein.sequenceString != null)
            return false;
        if (properties.get(PROPERTY.SOURCE) != null ? !properties.get(PROPERTY.SOURCE).equals(protein.getProperty(PROPERTY.SOURCE)) : protein.getProperty(PROPERTY.SOURCE) != null) return false;

        return true;
    }

    @Override
    public int hashCode() {
        int result = name != null ? name.hashCode() : 0;
        result = 31 * result + accession.hashCode();
        result = 31 * result + (properties.get(PROPERTY.SOURCE) != null ? properties.get(PROPERTY.SOURCE).hashCode() : 0);
        result = 31 * result + (sequenceString != null ? sequenceString.hashCode() : 0);
        return result;
    }
}
