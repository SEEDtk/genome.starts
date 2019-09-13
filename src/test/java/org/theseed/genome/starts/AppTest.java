package org.theseed.genome.starts;

import junit.framework.Test;

import junit.framework.TestCase;
import junit.framework.TestSuite;

import static org.hamcrest.MatcherAssert.assertThat;
import static org.hamcrest.Matchers.*;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;

import org.theseed.genome.Genome;
import org.theseed.io.TabbedLineReader;
import org.theseed.proteins.RoleMap;

/**
 * Unit test for simple App.
 */
public class AppTest extends TestCase
{
    /**
     * Create the test case
     *
     * @param testName name of the test case
     */
    public AppTest( String testName )
    {
        super( testName );
    }

    /**
     * @return the suite of tests being tested
     */
    public static Test suite()
    {
        return new TestSuite( AppTest.class );
    }

    /**
     * test the ORF tracker
     */
    public void testOrfTracker() {
        OrfTracker tracker = new OrfTracker();
        tracker.addStop(30, "ATG", "other", 0.75);
        tracker.addStop(69, "tag", "stop", 0.36);
        tracker.addStop(102, "TAA", "other", 0.2);
        OrfTracker.StopCodon stop1 = tracker.findOrf(72);
        assertNotNull(stop1);
        assertThat(stop1.getLocation(), equalTo(102));
        assertThat(stop1.getCodon(), equalTo("taa"));
        assertThat(stop1.getConfidence(), closeTo(0.8, 0.0001));
        assertThat(tracker.getOrfLength(), equalTo(30));
        assertThat(tracker.getOrfStop(), sameInstance(stop1));
        OrfTracker.StopCodon stop2 = tracker.findOrf(21);
        assertNotNull(stop2);
        assertThat(stop2.getLocation(), equalTo(30));
        assertThat(tracker.getOrfLength(), equalTo(27));
        assertThat(stop1.compareTo(stop2), greaterThan(0));
        stop1 = tracker.findOrf(110);
        assertNull(stop1);
        OrfTracker tracker2 = new OrfTracker();
        tracker2.addStop(30, "taa", "stop", 0.7);
        OrfTracker.StopCodon stop3 = tracker2.findOrf(30);
        assertThat(stop3.getLocation(), equalTo(30));
        assertThat(stop3, equalTo(stop2));
        assertThat(stop3.getCodon(), equalTo("taa"));
        assertThat(stop3.getConfidence(), closeTo(0.7, 0.0001));
        Iterator<OrfTracker.StopCodon> iter = tracker.iterator();
        assertTrue(iter.hasNext());
        stop1 = iter.next();
        assertThat(stop1.getLocation(), equalTo(30));
        assertThat(tracker.getOrfLength(), equalTo(27));
        assertTrue(iter.hasNext());
        stop1 = iter.next();
        assertThat(stop1.getLocation(), equalTo(69));
        assertThat(stop1, sameInstance(tracker.getOrfStop()));
        assertThat(tracker.getOrfLength(), equalTo(36));
        assertTrue(iter.hasNext());
        stop1 = iter.next();
        assertThat(stop1.getLocation(), equalTo(102));
        assertThat(tracker.getOrfLength(), equalTo(30));
        assertFalse(iter.hasNext());
    }

    /**
     * Test the basic functionality of the StopPredictionFile.
     */
    public void testContigOrfTracker() {
        ContigOrfTracker tracker = new ContigOrfTracker();
        tracker.addStop(112, "tag", "stop", 0.1);
        tracker.addStop(102, "tga", "stop", 0.2);
        tracker.addStop(36, "tag", "other", 0.3);
        tracker.addStop(333, "tga", "other", 0.4);
        tracker.addStop(200, "tga", "stop", 0.5);
        tracker.addStop(95, "taa", "other", 0.6);
        tracker.addStop(38, "taa", "stop", 0.7);
        tracker.addStop(40, "TGA", "other", 0.8);
        tracker.addStop(61, "taa", "other", 0.9);
        OrfTracker.StopCodon stop1 = tracker.findOrf(100);
        assertThat(stop1.getLocation(), equalTo(112));
        assertThat(tracker.getOrfLength(), equalTo(48));
        assertThat(stop1.getCodon(), equalTo("tag"));
        assertThat(stop1.getConfidence(), closeTo(0.1, 0.0001));
        stop1 = tracker.findOrf(300);
        assertThat(stop1.getLocation(), equalTo(333));
        assertThat(stop1.getConfidence(), equalTo(0.6));
        stop1 = tracker.findOrf(301);
        assertNull(stop1);
        stop1 = tracker.findOrf(1);
        assertThat(stop1.getLocation(), equalTo(40));
        assertThat(stop1.getCodon(), equalTo("tga"));
        assertThat(tracker.getOrfLength(), equalTo(39));
        stop1 = tracker.findOrf(60);
        assertThat(stop1.getLocation(), equalTo(102));
        assertThat(tracker.getOrfLength(), equalTo(63));
        tracker.addStart(10, "other", 0.6);
        tracker.addStart(22, "start", 0.9);
        assertThat(tracker.getStartConfidence(10), closeTo(0.4, 1e-5));
        assertThat(tracker.getStartConfidence(22), closeTo(0.9, 1e-5));
        assertThat(tracker.getStartConfidence(102), equalTo(0.0));
    }

    /**
     * test a full genome prediction file
     *
     * @throws IOException
     */
    public void testPredictionFile() throws IOException {
        File inFile = new File("src/test", "test.tbl");
        Map<String, ContigOrfTracker> contigMap = ContigOrfTracker.readPredictionFile(inFile);
        assertThat(contigMap.keySet().size(), equalTo(2));
        ContigOrfTracker tracker = contigMap.get("NC002112");
        assertNotNull(tracker);
        OrfTracker.StopCodon stop1 = tracker.findOrf(100);
        assertThat(stop1.getLocation(), equalTo(112));
        assertThat(tracker.getOrfLength(), equalTo(48));
        assertThat(stop1.getCodon(), equalTo("tag"));
        assertThat(stop1.getConfidence(), closeTo(0.1, 0.0001));
        stop1 = tracker.findOrf(300);
        assertThat(stop1.getLocation(), equalTo(333));
        assertThat(stop1.getConfidence(), equalTo(0.6));
        stop1 = tracker.findOrf(301);
        assertNull(stop1);
        stop1 = tracker.findOrf(1);
        assertThat(stop1.getLocation(), equalTo(40));
        assertThat(stop1.getCodon(), equalTo("tga"));
        assertThat(tracker.getOrfLength(), equalTo(39));
        stop1 = tracker.findOrf(60);
        assertThat(stop1.getLocation(), equalTo(102));
        assertThat(tracker.getOrfLength(), equalTo(63));
        // test the starts
        assertThat(tracker.getStartConfidence(203), closeTo(0.85, 1e-5));
        assertThat(tracker.getStartConfidence(303), closeTo(0.25, 1e-5));
        assertThat(tracker.getStartConfidence(403), closeTo(0.35, 1e-5));
        assertThat(tracker.getStartConfidence(503), equalTo(0.0));
        tracker = contigMap.get("NC_004347");
        assertNotNull(tracker);
        stop1 = tracker.findOrf(231);
        assertThat(stop1.getLocation(), equalTo(237));
        assertThat(stop1.getCodon(), equalTo("taa"));
        assertThat(stop1.getConfidence(), closeTo(0.0154, 0.00001));
        assertThat(tracker.getStartConfidence(203), equalTo(0.0));
    }

    /**
     * test the nucleon analyzer for RBS
     */
    public void testRBSs() {
        // Start off with the tail-end tetramers.
        //             400004000040000654000054000054000000
        String seq1 = "ggtgcgtgactgatcggtgatcggtgacgtgatccc";
        int[] expected = new int[] { 4, 0, 0, 0, 0, 4, 0, 0, 0, 0, 4, 0, 0, 0,
                                     0, 6, 5, 4, 0, 0, 0, 0, 5, 4, 0, 0, 0, 0,
                                     5, 4, 0, 0, 0, 0, 0, 0 };
        for (int i = 0; i < seq1.length(); i++) {
            assertThat("Error at position " + i, StartCodonFinder.checkRBS(seq1, i+1),
                    equalTo(expected[i]));
        }
        // Try some of the longer ones.
        //		7654000076540000A987654000065400007654000
        seq1 = "ggaggtgcaggaggtcggggggtgatcggggtgcaggtgat";
        expected = new int[] { 7, 6, 5, 4, 0, 0, 0, 0, 7, 6, 5, 4, 0, 0, 0, 0,
                              10, 9, 8, 7, 6, 5, 4, 0, 0, 0, 0, 6, 5, 4, 0, 0,
                               0, 0, 7, 6, 5, 4, 0, 0, 0 };
        for (int i = 0; i < seq1.length(); i++) {
            assertThat("Error at position " + i, StartCodonFinder.checkRBS(seq1, i+1),
                    equalTo(expected[i]));
        }
        // A final list to get the ones we missed.
        seq1 = "aggaggtgatccc";
        expected = new int[] { 10, 9, 8, 7, 6, 5, 4, 0, 0, 0, 0, 0, 0 };
        for (int i = 0; i < seq1.length(); i++) {
            assertThat("Error at position " + i, StartCodonFinder.checkRBS(seq1, i+1),
                    equalTo(expected[i]));
        }
        // Now we try all possible combinations.
        String[] seqs = new String[] {"aggaggtgat", "gggaggtgat", "agggggtgat", "ggggggtgat"};
        for (String seq : seqs)
            for (int i = 0; i < 7; i++)
                for (int j = 10 - i; j >= 4; j--) {
                    String chunk = seq.substring(i, i + j);
                    assertThat("Error checking " + chunk, StartCodonFinder.checkRBS(chunk, 1),
                            equalTo(j));
                }

    }

    /**
     * test Pribnow detection
     */
    public void testPribnow() {
        String seq = "tatataataattctaatcaacct";
        assertThat(StartCodonFinder.pribnowBox(seq, 1), closeTo(2.82, 1e-4));
        assertThat(StartCodonFinder.pribnowBox(seq, 2), closeTo(0.49, 1e-4));
        assertThat(StartCodonFinder.pribnowBox(seq, 3), closeTo(4.20, 1e-4));
        assertThat(StartCodonFinder.pribnowBox(seq, 4), closeTo(0.59, 1e-4));
        assertThat(StartCodonFinder.pribnowBox(seq, 5), closeTo(2.20, 1e-4));
    }

    /**
     * test start detection, AA counting, and GC content
     */
    public void testCoding() {
        //            ----2--0----1--
        String seq = "atttttgatgccgtg";
        int[] expected = new int[] { -1, -1, -1, -1, 2, -1, -1, 0, -1, -1, -1, -1, 1, -1, -1 };
        for (int i = 0; i < seq.length(); i++) {
            assertThat("Error at position " + i, StartCodonFinder.checkStart(seq, i+1),
                    equalTo(expected[i]));
        }
        seq = "tgatctgggctgatcactttacaaccggttcttcagcctgaagaacacccatgaggctacattgggacagatcaaggttgctccgtgaagaaacgtctcgatatgatagagcaactccttggtg";
        double gc = StartCodonFinder.gcContent(seq, 86, 72);
        assertThat(gc, equalTo(0.50));
        gc = StartCodonFinder.gcContent(seq, 11, 10);
        assertThat(gc, closeTo(0.60, 1e-5));
        seq = "acgtactgtgatctgggctgtatcactttacaaccggttcttcagcctgaagaacacccatgaggctacattgggacagatcaaggttgctctgaagaaacgtctcga" +
              "tatgatagagcaactccttggtgcagtaaaatcctatgcggcgtttgagaaggataccctggagcgggtaactgccatgcgtgcatcggtgggaactgca" +
              "ggtgccggggatcttgcaggtattgaggctgaatcacgttcgatgcttggaagactttttgcagtgatggaaaactatcctgacctcaaaacctcccaga" +
              "cggtaacaaccctgatggactcggtaaaacaggttgaagatgagattgcaaggcagcgctacacagcaaacaatgttgccgagcagttcaacacaatgat" +
              "ggatacaataccctcaaatattattggcaaatttgcaggcctcataaaactccagtaccttgagtttgctggagaagatattgagaagagacctgaaatc" +
              "tcattttaattatttttatctagagtctttaaggaccagatagtagttacatatctagtagaggaacctgaaggatatccccattccagttttccactaa" +
              "ttagaatatcccgataactgaatagaatccgtaatgtatctatacagttctctcatggaatgtatctgaacaatccatctgaaatacttacaatcgtgta" +
              "gcaattgggtggaattgttatcagtaagtcagaccagatcacgtcgtctctacgatgacttttccatctgcaaagatcctgaccaccccaccactttctg" +
              "atacaacaatccctaaagcctgggtctcctgagtaattgcagcgatcgaaacatggcgagttccaaacccttttggcaccgtaactctgctggtatcaac" +
              "acgtacgtacgt";
        double[] counts = new double[] { 2.000000, 2.000000, 3.333333,  2.666667, 2.666667, 2.333333,
                                         3.333333, 5.000000, 3.000000, 13.333333, 2.333333, 5.666667,
                                         5.666667, 4.666667, 9.000000, 11.333333, 5.000000, 3.333333,
                                         2.666667, 0.000000, 3.333333 };
        double[] aaVector = StartCodonFinder.aaContent(seq, 9, 909);
        for (int i = 0; i < 21; i++)
            assertThat("Error in count " + i, aaVector[i], closeTo(counts[i], 1e-5));
    }

    /**
     * test the start codon finder
     *
     * @throws IOException
     */
    public void testStartFinder() throws IOException {
        File stopFile = new File("src/test", "testStops.tbl");
        Map<String, ContigOrfTracker> contigMap = ContigOrfTracker.readPredictionFile(stopFile);
        ContigOrfTracker contigTracker = contigMap.get("NC1301");
        String contigSeq =
                "tgatctgggctgtatcactttacaaccggttcttcagcctgaagaacacccatgaggctacattgggacagatcaaggttgctctgaagaaacgtctcga" +
                "tatgatagagcaactccttggtgcagtaaaatcctatgcggcgtttgagaaggataccctggagcgggtaactgccatgcgtgcatcggtgggaactgca" +
                "ggtgccggggatcttgcaggtattgaggctgaatcacgttcgatgcttggaagactttttgcagtgatggaaaactatcctgacctcaaaacctcccaga" +
                "cggtaacaaccctgatggactcggtaaaacaggttgaagatgagattgcaaggcagcgctacacagcaaacaatgttgccgagcagttcaacacaatgat" +
                "ggatacaataccctcaaatattattggcaaatttgcaggcctcataaaactccagtaccttgagtttgctggagaagatattgagaagagacctgaaatc" +
                "tcattttaattatttttatctagagtctttaaggaccagatagtagttacatatctagtagaggaacctgaaggatatccccattccagttttccactaa" +
                "ttagaatatcccgataactgaatagaatccgtaatgtatctatacagttctctcatggaatgtatctgaacaatccatctgaaatacttacaatcgtgta" +
                "gcaattgggtggaattgttatcagtaagtcagaccagatcacgtcgtctctacgatgacttttccatctgcaaagatcctgaccaccccaccactttctg" +
                "atacaacaatccctaaagcctgggtctcctgagtaattgcagcgatcgaaacatggcgagttccaaacccttttggcaccgtaactctgctggtatcaac";
        File startFile = new File("src/test", "testStarts.tbl");
        TabbedLineReader startReader = new TabbedLineReader(startFile);
        int locIdx = startReader.findField("location");
        int startCodonIdx = startReader.findField("codonI");
        int orfLengthIdx = startReader.findField("orfLength");
        int regionLengthIdx = startReader.findField("regionLength");
        int stopCodonIdx = startReader.findField("stopI");
        int confIdx = startReader.findField("confidence");
        int rbsIdx = startReader.findField("rbsLen");
        int rbsGapIdx = startReader.findField("rbsGap");
        for (StartCodonFinder.StartCodon startFound : new StartCodonFinder("NC1301", contigSeq, contigTracker)) {
            TabbedLineReader.Line line = startReader.next();
            int location = line.getInt(locIdx);
            assertThat("Location error", startFound.getLoc(), equalTo(line.getInt(locIdx)));
            double[] data = startFound.getData();
            int codon = line.getInt(startCodonIdx);
            int codon0 = line.getInt(stopCodonIdx);
            for (int i = 0; i < 3; i++) {
                double datum = data[StartCodonFinder.START_CODON_IDX + i];
                if (i != codon) {
                    assertThat("Error expecting " + location, datum, equalTo(0.0));
                } else {
                    assertThat("Error expecting " + location, datum, equalTo(1.0));
                }
                datum = data[StartCodonFinder.STOP_CODON_IDX + i];
                if (i != codon0) {
                    assertThat("Error expecting " + location, datum, equalTo(0.0));
                } else {
                    assertThat("Error expecting " + location, datum, equalTo(1.0));
                }
            }
            double pribScore = 0.0;
            for (int i = 0; i < 5; i++) {
                int pPos = location - 12 + i;
                pribScore += StartCodonFinder.pribnowBox(contigSeq, pPos);
            }
            assertThat("Error expecting " + location, data[StartCodonFinder.PRIBNOW_SCORE_IDX],
                    equalTo(pribScore));
            double regionLen = data[StartCodonFinder.REGION_LEN_IDX];
            double orfLen = data[StartCodonFinder.ORF_LEN_IDX];
            int orfEndPoint = startFound.getLoc() + (int) regionLen;
            assertThat("Error expecting " + location, data[StartCodonFinder.ORF_LEN_IDX], equalTo(line.getDouble(orfLengthIdx)));
            assertThat("Error expecting " + location, data[StartCodonFinder.REGION_LEN_IDX], equalTo(line.getDouble(regionLengthIdx)));
            assertThat("Error expecting " + location, data[StartCodonFinder.STOP_CONFIDENCE_IDX], closeTo(line.getDouble(confIdx), 1e-5));
            assertThat("Error expecting " + location, data[StartCodonFinder.RBS_GAP_IDX], equalTo(line.getDouble(rbsGapIdx)));
            assertThat("Error expecting " + location, data[StartCodonFinder.RBS_LEN_IDX], equalTo(line.getDouble(rbsIdx)));
            double gc = StartCodonFinder.gcContent(contigSeq, orfEndPoint, (int) orfLen);
            assertThat("Error expecting " + location, data[StartCodonFinder.GC_IDX], equalTo(gc));
            double[] aaVec = StartCodonFinder.aaContent(contigSeq, startFound.getLoc(), orfEndPoint);
            for (int i = 0; i < aaVec.length; i++)
                assertThat("Error expecting " + location + " aa[" + i + "]", data[StartCodonFinder.AA_PROFILE_IDX + i], equalTo(aaVec[i]));
        }
        TabbedLineReader.Line line = startReader.next();
        assertThat(line.getInt(locIdx), equalTo(-1));
        startReader.close();
    }

    /**
     * Test the contig starts object.
     *
     * @throws IOException
     * @throws NumberFormatException
     */
    public void testContigStarts() throws NumberFormatException, IOException {
    	File gFile = new File("src/test", "372461.17.gto");
    	Genome genome = new Genome(gFile);
    	Map<String, ContigStarts> contigStartMap = ContigStarts.contigStartsMap(genome, '+', null);
    	File dFile = new File("src/test", "starts.tbl");
    	TabbedLineReader dataFile = new TabbedLineReader(dFile);
    	for (TabbedLineReader.Line line : dataFile) {
    		ContigStarts startMap = contigStartMap.get(line.get(0));
    		String roles = startMap.getRoles(line.getInt(1));
    		String expect = line.get(2);
    		if (expect.contentEquals("x")) expect = null;
    		assertThat(roles, equalTo(expect));
    	}
    	dataFile.close();
    	// Now with mapped roles.
    	RoleMap roleMap = RoleMap.load(new File("src/test", "roles.tbl"));
    	contigStartMap = ContigStarts.contigStartsMap(genome, '+', roleMap);
    	dataFile = new TabbedLineReader(new File("src/test", "starts2.tbl"));
    	for (TabbedLineReader.Line line : dataFile) {
    		ContigStarts startMap = contigStartMap.get(line.get(0));
    		String roles = startMap.getRoles(line.getInt(1));
    		String expect = line.get(2);
    		if (expect.contentEquals("x")) expect = null;
    		assertThat(roles, equalTo(expect));
    	}
    	dataFile.close();
    }

}
