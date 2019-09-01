/**
 *
 */
package org.theseed.genome.starts;

import java.util.Iterator;
import java.util.NavigableSet;
import java.util.TreeSet;

/**
 * This class keeps a sorted list of the stop codons in a single frame of a contig.  The codons
 * are kept sorted by location so that we can easily compute the open reading frame surrounding
 * a particular position.
 *
 * @author Bruce Parrello
 *
 */
public class OrfTracker implements Iterable<OrfTracker.StopCodon> {

    /**
     * This class represents a record in the file and indicates a possible stop codon.
     *
     * @author Bruce Parrello
     *
     */
    public static class StopCodon implements Comparable<StopCodon> {

        // FIELDS
        /** the position (1-based) of the stop codon in the contig */
        private int location;
        /** the three-letter codon text (lower case) */
        private String codon;
        /** the confidence that this codon is the stop for a coding region */
        private double confidence;

        /**
         * Create a new stop codon representative.
         *
         * @param location		the location (1-based) of the codon in the contig
         * @param codon			the three letters representing the nucleons of the codon
         * @param type			the predicted status of the codon ("other" or "stop")
         * @param confidence	the confidence of the prediction
         */
        protected StopCodon(int location, String codon, String type, double confidence) {
            this.location = location;
            this.codon = codon.toLowerCase();
            if (type.contentEquals("stop"))
                this.confidence = confidence;
            else
                this.confidence = 1.0 - confidence;
        }

        /**
         * We compare two stop codons by location, so they are stored in location order within
         * the frame.
         */
        @Override
        public int compareTo(StopCodon o) {
            return this.location - o.location;
        }

        /**
         * @return the location of this stop codon
         */
        public int getLocation() {
            return location;
        }

        /**
         * @return the DNA value of this stop codon
         */
        public String getCodon() {
            return codon;
        }

        /**
         * @return the probability that this is a true stop
         */
        public double getConfidence() {
            return confidence;
        }

        @Override
        public int hashCode() {
            final int prime = 31;
            int result = 1;
            result = prime * result + location;
            return result;
        }

        @Override
        public boolean equals(Object obj) {
            if (this == obj)
                return true;
            if (obj == null)
                return false;
            if (getClass() != obj.getClass())
                return false;
            StopCodon other = (StopCodon) obj;
            if (location != other.location)
                return false;
            return true;
        }

    }

    /**
     * This is an iterator through this object, returning all the stop codons in order.
     * Its main job is to keep the current-stop and orf-length members of the enclosing
     * class in sync.
     *
     * @author Bruce Parrello
     *
     */
    public class OrfIterator implements Iterator<StopCodon> {

        /** iterator through the underlying tree set */
        private Iterator<StopCodon> treeIter;

        /**
         * Initialize the iterator.
         */
        public OrfIterator() {
            treeIter = OrfTracker.this.stopList.iterator();
            OrfTracker.this.orfStop = null;
        }

        @Override
        public boolean hasNext() {
            return treeIter.hasNext();
        }

        @Override
        public StopCodon next() {
            StopCodon retVal = this.treeIter.next();
            if (retVal != null) {
                OrfTracker.this.setOrfLength(retVal, OrfTracker.this.orfStop);
                OrfTracker.this.orfStop = retVal;
            }
            return retVal;
        }

    }

    // FIELDS
    /** tree containing the stop codons (including both the true and false stops) */
    private NavigableSet<StopCodon> stopList;
    /** length of the current ORF, in base pairs (this is the interior length, between the stops) */
    private int orfLength;
    /** ending stop codon of the current ORF */
    private StopCodon orfStop;
    /** dummy codon for searching */
    private StopCodon searchCodon;

    /**
     * Construct a new, blank ORF tracker.
     */
    public OrfTracker() {
        this.stopList = new TreeSet<StopCodon>();
        this.orfLength = 0;
        this.orfStop = null;
        this.searchCodon = new StopCodon(0, "atg", "stop", 0);
    }

    /**
     * Store the new ORF length and closing stop codon.
     *
     * @param newStop	closing stop codon
     * @param oldStop	opening stop codon (or NULL, if there is none)
     */
    protected void setOrfLength(StopCodon newStop, StopCodon oldStop) {
        if (oldStop == null) {
            this.orfLength = (newStop.location - 1) / 3 * 3;
        } else {
            this.orfLength = newStop.location - oldStop.location - 3;
        }
        this.orfStop = newStop;
    }

    /**
     * Add a location to the tracker.
     *
     * @param location		the location (1-based) of the codon in the contig
     * @param codon			the three letters representing the nucleons of the codon
     * @param type			the predicted status of the codon ("other" or "stop")
     * @param confidence	the confidence of the prediction
     */
    public void addStop(int location, String codon, String type, double confidence) {
        StopCodon stop = new StopCodon(location, codon, type, confidence);
        this.stopList.add(stop);
    }

    /**
     * Position on the ORF containing a specified location.
     *
     * @param location	location to contain the ORF
     *
     * @return the stop codon at the end of the ORF
     */
    public StopCodon findOrf(int location) {
        this.searchCodon.location = location;
        StopCodon retVal = this.stopList.ceiling(searchCodon);
        if (retVal != null) {
            // Get the opening stop, so we can compute the length.  The length is the interior length, between the stops.
            StopCodon opener = this.stopList.lower(retVal);
            // Store the length and save the closing stop.
            this.setOrfLength(retVal, opener);
        }
        // Return the closing stop.
        return retVal;
    }

    /**
     * @return the length of the ORF at the current position
     */
    public int getOrfLength() {
        return orfLength;
    }

    /**
     * @return the stop codon for the ORF at the current position
     */
    public StopCodon getOrfStop() {
        return orfStop;
    }

    @Override
    public Iterator<StopCodon> iterator() {
        return this.new OrfIterator();
    }



}
