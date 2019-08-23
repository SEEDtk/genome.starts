/**
 *
 */
package org.theseed.genome.starts;

/**
 * This object reads a file of stop predictions.  For each frame in each contig, it keeps an ordered list of
 * the potential stop codons.  For each such codon, it records the identity of the codon (TAA, TAG, TGA), its
 * location (which is the basis of the ordering), and the confidence that the preceding ORF contains a coding
 * region.  Once the file is read, this object can be used to extract the ORFs in each frame of each contig.
 *
 * The key columns in the input file are "location", which contains the contig ID and the codon position separated
 * by a semicolon, "codon", which contains the three-letter codon itself, "predicted", and "confidence".  We want
 * our confidence to be the confidence that the stop is a true stop, so if the "predicted" is "other", we invert
 * the confidence by subtracting it from 1.  The input file is a tabbed dataset.
 *
 * @author Bruce Parrello
 *
 */
public class StopPredictionFile {

    //TODO body of StopPredictionFile

}
