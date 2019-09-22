/**
 *
 */
package org.theseed.genome.starts;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import org.apache.commons.lang3.StringUtils;
import org.theseed.io.TabbedLineReader;

/**
 * This object tracks all the ORFs for a contig.  For each frame in the contig, it keeps an ordered list of
 * the potential stop codons.  For each such codon, it records the identity of the codon (TAA, TAG, TGA), its
 * location (which is the basis of the ordering), and the confidence that the preceding ORF contains a coding
 * region.
 *
 * A static method is provided to create an object instance for each contig in a genome by
 * processing the stop prediction file produced by the stop-prediction neural network.
 * Once the file is read, the map returned can be used to extract the ORFs in each frame of each contig.
 *
 * The key columns in the input file are "location", which contains the contig ID and the codon position separated
 * by a semicolon, "codon", which contains the three-letter codon itself, "predicted", and "confidence".  We want
 * our confidence to be the confidence that the stop is a true stop, so if the "predicted" is "other", we invert
 * the confidence by subtracting it from 1.  The input file is a tabbed dataset.  It is also possible the input
 * will contain START codons.  If this is the case, we track the confidence of those, too.
 *
 * @author Bruce Parrello
 *
 */
public class ContigOrfTracker {


    // FIELDS
    /** array of orf trackers, one per frame */
    private OrfTracker[] frames;
    /** length of the last ORF requested */
    private int orfLength;
    /** map of start codon locations to confidences */
    private HashMap<Integer, Double> startMap;

    /**
     * Construct a blank, empty, stop-prediction file.
     */
    public ContigOrfTracker() {
        clear();
    }


    /**
     * Read a start/stop prediction file and return a mapping that provides access to the
     * ContigOrfTracker for each contig by contig ID.
     *
     * @param inFile	input file containing the stop predictions for an entire genome
     *
     * @return a map keyed by contig ID that returns each contig's ORF tracker
     * @throws IOException
     */
    public static Map<String, ContigOrfTracker> readPredictionFile(File inFile) throws IOException {
        // Create the master hash map.
        HashMap<String, ContigOrfTracker> retVal = new HashMap<String, ContigOrfTracker>();
        // Prepare to read the file.
        TabbedLineReader tabStream = new TabbedLineReader(inFile);
        try {
            int locCol = tabStream.findField("location");
            int codonCol = tabStream.findField("codon");
            int typeCol = tabStream.findField("predicted");
            int confCol = tabStream.findField("confidence");
            // Loop through the input.
            for (TabbedLineReader.Line line : tabStream) {
                // Get the location string and parse it into the contig ID and position.
                String fullLocation = line.get(locCol);
                String[] parts = StringUtils.split(fullLocation, ';');
                int location = Integer.parseInt(parts[1]);
                // Get the contig tracker for this contig.
                ContigOrfTracker contigTracker;
                if (! retVal.containsKey(parts[0])) {
                    contigTracker = new ContigOrfTracker();
                    retVal.put(parts[0], contigTracker);
                } else {
                    contigTracker = retVal.get(parts[0]);
                }
                // Determine the codon type.
                String codon = line.get(codonCol).toLowerCase();
                if (StartCodonFinder.testStart(codon) >= 0) {
                    // Here we have a start codon.
                    contigTracker.addStart(location, line.get(typeCol), line.getDouble(confCol));
                } else {
                    // Here we have a stop codon.
                    contigTracker.addStop(location, codon, line.get(typeCol),
                            line.getDouble(confCol));
                }
            }
        } finally {
            tabStream.close();
        }
        return retVal;
    }

    /**
     * Store a start codon in the start map.
     *
     * @param location		location of the codon in the contig (1-based)
     * @param type			type of prediction ("other" or "start")
     * @param confidence	confidence in prediction
     */
    public void addStart(int location, String type, double confidence) {
        Double realConf = (type.contentEquals("start") ? confidence : 1.0 - confidence);
        this.startMap.put(location, realConf);
    }

    /**
     * Add a stop codon prediction to the appropriate frame.
     *
     * @param location		location of the stop codon (1-based)
     * @param codon			three-letter DNA string for the codon
     * @param type			type of prediction ("other" or "stop")
     * @param confidence	confidence of the prediction
     */
    public void addStop(int location, String codon, String type, double confidence) {
        int frame = location % 3;
        this.frames[frame].addStop(location, codon, type, confidence);
    }

    /**
     * Initialize the frame array.
     */
    private void clear() {
        this.frames = new OrfTracker[3];
        for (int i = 0; i < 3; i++)
            this.frames[i] = new OrfTracker();
        this.orfLength = 0;
        this.startMap = new HashMap<Integer, Double>();
    }

    /**
     * Find the ORF containing the specified location.
     *
     * @param location		location whose ORF is desired
     *
     * @return 	the terminating stop codon for the ORF, or NULL if the location is
     * 			in the tail of the contig (after the last stop)
     */
    public OrfTracker.StopCodon findOrf(int location) {
        int frame = location % 3;
        OrfTracker.StopCodon retVal = this.frames[frame].findOrf(location);
        this.orfLength = this.frames[frame].getOrfLength();
        return retVal;
    }

    /**
     * @return the length of the last ORF returned to the user
     */
    public int getOrfLength() {
        return this.orfLength;
    }

    /**
     * @return the confidence in a start codon at the specified location
     *
     * @param location	location of the start in question
     */
    public double getStartConfidence(int location) {
        Integer searcher = Integer.valueOf(location);
        double retVal = 0.0;
        Double found = this.startMap.get(searcher);
        if (found != null) {
            retVal = found;
        }
        return retVal;
    }

}
