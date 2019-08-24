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
        assertThat(tracker.getOrfLength(), equalTo(30));
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
        assertThat(tracker.getOrfLength(), equalTo(30));
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
    	assertThat(tracker.getOrfLength(), equalTo(40));
    	stop1 = tracker.findOrf(60);
    	assertThat(stop1.getLocation(), equalTo(102));
    	assertThat(tracker.getOrfLength(), equalTo(63));
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
    	assertThat(tracker.getOrfLength(), equalTo(40));
    	stop1 = tracker.findOrf(60);
    	assertThat(stop1.getLocation(), equalTo(102));
    	assertThat(tracker.getOrfLength(), equalTo(63));
    	tracker = contigMap.get("NC_004347");
    	assertNotNull(tracker);
    	stop1 = tracker.findOrf(231);
    	assertThat(stop1.getLocation(), equalTo(237));
    	assertThat(stop1.getCodon(), equalTo("taa"));
    	assertThat(stop1.getConfidence(), closeTo(0.0154, 0.00001));
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
    		assertThat("Error at position " + i, ContigNucleonAnalyzer.checkRBS(seq1, i+1),
    				equalTo(expected[i]));
    	}
    	// Try some of the longer ones.
    	//		7654000076540000A987654000065400007654000
    	seq1 = "ggaggtgcaggaggtcggggggtgatcggggtgcaggtgat";
    	expected = new int[] { 7, 6, 5, 4, 0, 0, 0, 0, 7, 6, 5, 4, 0, 0, 0, 0,
    						  10, 9, 8, 7, 6, 5, 4, 0, 0, 0, 0, 6, 5, 4, 0, 0,
    						   0, 0, 7, 6, 5, 4, 0, 0, 0 };
    	for (int i = 0; i < seq1.length(); i++) {
    		assertThat("Error at position " + i, ContigNucleonAnalyzer.checkRBS(seq1, i+1),
    				equalTo(expected[i]));
    	}
    	// A final list to get the ones we missed.
    	seq1 = "aggaggtgatccc";
    	expected = new int[] { 10, 9, 8, 7, 6, 5, 4, 0, 0, 0, 0, 0, 0 };
    	for (int i = 0; i < seq1.length(); i++) {
    		assertThat("Error at position " + i, ContigNucleonAnalyzer.checkRBS(seq1, i+1),
    				equalTo(expected[i]));
    	}
    	// Now we try all possible combinations.
    	String[] seqs = new String[] {"aggaggtgat", "gggaggtgat", "agggggtgat", "ggggggtgat"};
    	for (String seq : seqs)
    		for (int i = 0; i < 7; i++)
    			for (int j = 10 - i; j >= 4; j--) {
    				String chunk = seq.substring(i, i + j);
    				assertThat("Error checking " + chunk, ContigNucleonAnalyzer.checkRBS(chunk, 1),
    						equalTo(j));
    			}

    }

    /**
     * test start detection
     */
    public void testStarts() {
    	//            ----2--0----1--
    	String seq = "atttttgatgccgtg";
    	int[] expected = new int[] { -1, -1, -1, -1, 2, -1, -1, 0, -1, -1, -1, -1, 1, -1, -1 };
    	for (int i = 0; i < seq.length(); i++) {
    		assertThat("Error at position " + i, ContigNucleonAnalyzer.checkStart(seq, i+1),
    				equalTo(expected[i]));
    	}
    }

}
