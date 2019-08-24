/**
 *
 */
package org.theseed.genome.starts;

import org.apache.commons.lang3.StringUtils;

/**
 * This object contains utilities for analyzing start-related data for a location in a
 * contig.  This includes trying to find a ribosomal binding sequence as well as determining
 * if we have an actual start codon.
 *
 * @author Bruce Parrello
 *
 */
public class ContigNucleonAnalyzer {

	/** the four Shine-Delgarno consensus sequences */
	private static final String[] RBS =
		{ "aggaggtgat", "agggggtgat", "gggaggtgat", "ggggggtgat" };

	/** the three start codons */
	private static final String[] STARTS = { "atg", "gtg", "ttg" };

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

}
