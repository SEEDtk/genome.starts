/**
 *
 */
package org.theseed.genome.starts;

import java.util.Arrays;
import java.util.HashMap;
import java.util.Iterator;

import org.apache.commons.lang3.StringUtils;

/**
 * This object contains numerous utilities for analyzing data relating to ORFs.  It is used
 * to iterate through a contig, producing start codon information.  For each start codon, we
 * need to know
 *
 *  1. GC content of the ORF
 *  2. Length of the ORF
 *  3. Distance from the start to the end of the ORF
 *  4. The length of the ribosomal binding site
 *  5. The gap to the ribosomal binding site
 *  6. The identity of the start codon (ATG, GTG, TTG)
 *  7. The confidence of stop codon
 *  8. The identity of the stop codon (TAA, TAG, TGA)
 *  9. The amino acid profile of the region between the start and the stop
 *
 *  The two codon identifiers are coded as a one-hot of three values.  There are 21 amino acids.  So, the output will
 *  be an array of 33 double-precision values.
 *
 * @author Bruce Parrello
 *
 */
public class StartCodonFinder implements Iterator<StartCodonFinder.StartCodon>, Iterable<StartCodonFinder.StartCodon> {


    /** the four Shine-Delgarno consensus sequences */
    private static final String[] RBS =
        { "aggaggtgat", "agggggtgat", "gggaggtgat", "ggggggtgat" };

    /** the three start codons */
    private static final String[] STARTS = { "atg", "gtg", "ttg" };

    /** the three stop codons */
    private static final String[] STOPS = { "taa", "tag", "tga" };

    /** genetic code 11 translation table */
    @SuppressWarnings("serial")
    private static final HashMap<String, String> GENETIC_CODE_11 = new HashMap<String, String>() {{
        put("aaa","K"); put("aac","N"); put("aag","K"); put("aat","N"); put("aca","T");
        put("acc","T"); put("acg","T"); put("act","T"); put("aga","R"); put("agc","S");
        put("agg","R"); put("agt","S"); put("ata","I"); put("atc","I"); put("atg","M");
        put("att","I"); put("caa","Q"); put("cac","H"); put("cag","Q"); put("cat","H");
        put("cca","P"); put("ccc","P"); put("ccg","P"); put("cct","P"); put("cga","R");
        put("cgc","R"); put("cgg","R"); put("cgt","R"); put("cta","L"); put("ctc","L");
        put("ctg","L"); put("ctt","L"); put("gaa","E"); put("gac","D"); put("gag","E");
        put("gat","D"); put("gca","A"); put("gcc","A"); put("gcg","A"); put("gct","A");
        put("gga","G"); put("ggc","G"); put("ggg","G"); put("ggt","G"); put("gta","V");
        put("gtc","V"); put("gtg","V"); put("gtt","V"); put("taa","*"); put("tac","Y");
        put("tag","*"); put("tat","Y"); put("tca","S"); put("tcc","S"); put("tcg","S");
        put("tct","S"); put("tga","*"); put("tgc","C"); put("tgg","W"); put("tgt","C");
        put("tta","L"); put("ttc","F"); put("ttg","L"); put("ttt","F");
    }};
    /** conversion table from amino acid codes to numbers */
    private static final String AA_LIST = "ACDEFGHIKLMNPQRSTVWXY";
    /** epsilon for keeping counts away from zero */
    private static final double EPSILON = 1e-10;
    /** number of amino acids */
    private static final int AA_COUNT = AA_LIST.length();

    // MAP OF OUTPUT ARRAY
    public static final int GC_IDX = 0;
    public static final int ORF_LEN_IDX = 1;
    public static final int REGION_LEN_IDX = 2;
    public static final int RBS_LEN_IDX = 3;
    public static final int RBS_GAP_IDX = 4;
    public static final int START_CODON_IDX = 5;		// 3 entries
    public static final int STOP_CONFIDENCE_IDX = 8;
    public static final int STOP_CODON_IDX = 9; 		// 3 entries
    public static final int AA_PROFILE_IDX = 12; 		// 21 entries
    public static final int OUTPUT_LEN = 12 + AA_COUNT;

    /**
     * This class encapsulates information about a start codon we find.
     *
     * @author Bruce Parrello
     *
     */
    public static class StartCodon {

        /** position (1-based) of the codon */
        private int loc;
        /** array of values describing the codon's situation */
        private double[] data;

        /**
         * Initialize a new start descriptor.
         *
         * @param loc	location of the start
         * @param data	array of double-precision values describing the environment of the start (will be cloned)
         */
        protected StartCodon(int loc, double[] data) {
            this.loc = loc;
            this.data = data.clone();
        }

        /**
         * @return the location (1-based) of the start
         */
        public int getLoc() {
            return loc;
        }

        /**
         * @return the data array describing the start and the ORF that contains it
         */
        public double[] getData() {
            return data;
        }

    }


    /**
     * @return 	the index of the start codon at the current position, or -1 if this
     * 			position is not a start codon
     *
     * @param sequence	the contig DNA sequence
     * @param loc		the position (1-based) in the contig
     */
    public static int checkStart(String sequence, int loc) {
        String codon = StringUtils.substring(sequence, loc - 1, loc + 2);
        int retVal = STARTS.length - 1;
        while (retVal >= 0 && ! codon.contentEquals(STARTS[retVal])) retVal--;
        return retVal;
    }

    /**
     * @return 	the length of the ribosomal binding sequence at the current position, or
     * 			0 if there is no ribosomal binding sequence here
     *
     * @param sequence	the contig DNA sequence
     * @param loc		the position (1-based) in the contig
     */
    public static int checkRBS(String sequence, int loc) {
        String tetramer = StringUtils.substring(sequence, loc - 1, loc + 3);
        int retVal = 0;
        if (tetramer.length() >= 4) {
            // Loop through the possible RBS sequences trying to find a match.
            for (int i = 0; i < RBS.length; i++) {
                int pos = RBS[i].indexOf(tetramer);
                while (pos >= 0) {
                    // Here the tetramer matches part of the RBS.  We need to figure
                    // out how long we can extend the match.
                    String subSequence = StringUtils.substring(sequence, loc - 1, loc + 9 - pos);
                    String rbsSection = RBS[i].substring(pos);
                    int length = StringUtils.indexOfDifference(subSequence, rbsSection);
                    if (length < 0) length = subSequence.length();
                    if (retVal < length) retVal = length;
                    pos = RBS[i].indexOf(tetramer, pos + 1);
                }
            }
        }
        return retVal;
    }

    /**
     * @return the fractional GC content of a region in the specified sequence
     *
     * @param sequence		the sequence in question
     * @param endLocation	the location (1-based) past the end of the target region
     * @param length		the number of base pairs in the target region
     */
    public static double gcContent(String sequence, int endLocation, int length) {
        int gcCount = 0;
        int endPosition = endLocation - 1;
        for (int i = endPosition - length; i < endPosition; i++) {
            switch (sequence.charAt(i)) {
            case 'g':
            case 'c':
                gcCount++;
                break;
            }
        }
        return ((double) gcCount / length);
    }

    /**
     * @return the amino acid profile of a region in the specified sequence, in the form of log fraction
     *
     * @param sequence		the sequence in question
     * @param startLocation	the location (1-based) of the start of the target region
     * @param endLocation	the location (1-based) past the end of the target region
     */
    public static double[] aaContent(String sequence, int startLocation, int endLocation) {
        // Start with an array of zero counts.
        int[] counts = new int[AA_COUNT];
        Arrays.fill(counts, 0);
        // Run through the appropriate section of the sequence, 3 letters at a time.
        int endPosition = endLocation - 1;
        int total = 0;
        for (int i = startLocation - 1; i < endPosition; i += 3) {
            String codon = GENETIC_CODE_11.get(sequence.substring(i, i+3));
            int idx = AA_LIST.indexOf(codon);
            if (idx >= 0) counts[idx]++;
            total++;
        }
        // Now we take the log of the count + epsilon / total.
        double[] retVal = new double[AA_COUNT];
        for (int i = 0; i < AA_COUNT; i++) {
            retVal[i] = Math.log((counts[i] + EPSILON) / total);
        }
        return retVal;
    }

    // FIELDS
    /** the contig being traversed */
    private String contig;
    /** position of the last possible start */
    private int lastPos;
    /** position (1-based) of the next start codon, or 0 if none exists */
    private int nextPos;
    /** length of the last known ribosomal binding sequence (0 if none) */
    private int rbsLen;
    /** position (1-based) of the last known ribosomal binding sequence (0 if none) */
    private int rbsPos;
    /** output array for the next start codon */
    private double[] outputs;
    /** tracker for finding the ORFs */
    ContigOrfTracker orfTracker;

    /**
     * Initialize this object for iterating through the potential starts in a contig.
     *
     * @param sequence	the DNA sequence of the contig
     * @param tracker	the ORF tracker for the contig
     */
    public StartCodonFinder(String sequence, ContigOrfTracker tracker) {
        this.contig = sequence.toLowerCase();
        this.orfTracker = tracker;
        // The last position is 6 before the end, but it is 1-based.
        this.lastPos = sequence.length() - 5;
        // Create the output array.
        this.outputs = new double[33];
        // Clear the RBS data.
        this.rbsLen = 0;
        this.rbsPos = 0;
        // Find the first start.
        this.findNextStart(0);
    }

    /**
     * Find the next potential start codon and compute its data.
     *
     * @param startLoc	position from which to start looking
     */
    private void findNextStart(int startLoc) {
        // We begin by sliding through the sequence, looking for RBS sections and start codons.  We stop at the
        // first start codon that is inside an ORF.
        boolean done = false;
        // This variable slides through the contig.
        int pos = startLoc;
        while (pos < this.lastPos && ! done) {
            pos++;
            int rbs = checkRBS(this.contig, pos);
            // Note we reject the RBS if it's inside another one.
            if (rbs > 0 && pos > this.rbsLen + this.rbsPos) {
                // Here we found an RBS, so we need to save its information.
                this.rbsLen = rbs;
                this.rbsPos = pos;
            }
            // If this says we have a start, we need to check it.
            int startType = checkStart(this.contig, pos);
            if (startType >= 0) {
                OrfTracker.StopCodon orf = this.orfTracker.findOrf(pos);
                if (orf != null) {
                    // We are inside an ORF.  This is potentially a real start.
                    double[] data = new double[OUTPUT_LEN];
                    data[GC_IDX] = gcContent(this.contig, orf.getLocation(), this.orfTracker.getOrfLength());
                    data[ORF_LEN_IDX] = this.orfTracker.getOrfLength();
                    data[REGION_LEN_IDX] = orf.getLocation() - pos;
                    data[RBS_LEN_IDX] = this.rbsLen;
                    data[RBS_GAP_IDX] = pos - (this.rbsPos + this.rbsLen);
                    // Fill in the one-hot for the start.
                    for (int i = 0; i < 3; i++)
                        data[START_CODON_IDX + i] = (i == startType ? 1 : 0);
                    data[STOP_CONFIDENCE_IDX] = orf.getConfidence();
                    // Fill in the one-hot for the stop.
                    String stop = orf.getCodon();
                    for (int i = 0; i < 3; i++)
                        data[STOP_CODON_IDX + i] = (stop.contentEquals(STOPS[i]) ? 1 : 0);
                    // Finally, fill in the amino acid profile.
                    double[] aaProfile = aaContent(this.contig, pos, orf.getLocation());
                    System.arraycopy(aaProfile, 0, data, AA_PROFILE_IDX, aaProfile.length);
                    // Store the output.
                    this.outputs = data;
                    this.nextPos = pos;
                    done = true;
                }
            }
        }
        // Did we find something? If not, flag end-of-file.
        if (! done) this.nextPos = 0;
    }

    @Override
    public boolean hasNext() {
        return (this.nextPos != 0);
    }

    @Override
    public StartCodon next() {
        StartCodon retVal = null;
        if (this.nextPos > 0) {
            retVal = new StartCodon(this.nextPos, this.outputs);
            this.findNextStart(this.nextPos);
        }
        return retVal;
    }

    @Override
    public Iterator<StartCodon> iterator() {
        return this;
    }


}
