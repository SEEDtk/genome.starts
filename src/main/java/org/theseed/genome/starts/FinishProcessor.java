/**
 *
 */
package org.theseed.genome.starts;

import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.io.BalancedOutputStream;
import org.theseed.io.TabbedLineReader;
import org.theseed.locations.Location;
import org.theseed.utils.ICommand;

// TODO test FinishProcessor

/**
 * This process reads the output from the Combo neural net (rank starts and stops based on neighborhood) and the FancyStart neural net
 * (rank starts based on statistical characteristics) and picks the best-ranked start for each ORF.
 *
 * The two positional parameters are the Combo output file and the FancyStart output file, respectively.  The FancyStart output file
 * is echoed output with an prefixed column ("final") indicating if the start is the best one or not.
 *
 * @author Bruce Parrello
 *
 */
public class FinishProcessor implements ICommand {

    /** This class represents the best start found so far for an ORF */
    private class GoodStart {
        /** input line for the start */
        private TabbedLineReader.Line startLine;
        /** confidence in the start */
        private double confidence;
        /** location of the start in the contig */
        private int location;

        /**
         * Create a new start record.
         *
         * @param line	original input line for the start
         * @param loc	location of the start in the contig
         */
        public GoodStart(TabbedLineReader.Line line, int loc) {
            this.startLine = line;
            this.confidence = line.getDouble(confColumn);
            this.location = loc;
        }

        /**
         * Compare a new start to this one.  If it is better, store it.
         *
         * @param line	input line for the new start
         * @param loc	location of the new start in the contig
         *
         * @return the input line for the worse start.
         */
        public String merge(TabbedLineReader.Line line, int loc) {
            // We default to returning the new location.
            TabbedLineReader.Line retVal = line;
            double conf = line.getDouble(confColumn);
            if (conf > this.confidence || conf == this.confidence && loc < this.location) {
                // Here the new location is better. Return the old one and store the new one.
                retVal = this.startLine;
                this.startLine = line;
                this.confidence = conf;
                this.location = loc;
            }
            return retVal.getAll();
        }

        /**
         * @return the original input line for the start
         */
        public String getStartLine() {
            return startLine.getAll();
        }



    }
    // FIELDS
    /** map of contig IDs to ORFs */
    private Map<String, ContigOrfTracker> contigMap;
    /** input file for called starts */
    private TabbedLineReader startStream;
    /** input column containing the start confidence */
    private int confColumn;
    /** input column containing the start location */
    private int locColumn;
    /** input column containing the prediction */
    private int predColumn;
    /** output stream */
    private BalancedOutputStream outStream;
    /** map of ORFs to best starts */
    private Map<Location, GoodStart> startMap;


    // COMMAND LINE

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** TRUE if we want progress messages */
    @Option(name="-v", aliases={"--verbose", "--debug"}, usage="display progress on STDERR")
    private boolean debug;

    @Argument(index=0, metaVar="comboOutFile", usage="output file from start/stop caller")
    private File startStopFile;

    @Argument(index=1, metaVar="fancyOutFile", usage="output file from advanced start caller")
    private File startFile;


    /**
     * Parse command-line options to specify the parameters of this object.
     *
     * @param args	an array of the command-line parameters and options
     *
     * @return TRUE if successful, FALSE if the parameters are invalid
     */
    @Override
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.help = false;
        this.debug = false;
        // Parse the command line.
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
                // Read the start/stop file to compute the ORFs.
                if (this.debug) System.err.println("Reading ORFs from " + this.startStopFile);
                this.contigMap = ContigOrfTracker.readPredictionFile(this.startStopFile);
                // Open the start file to begin input and initialize the output file.
                this.startStream = new TabbedLineReader(this.startFile);
                this.outStream = new BalancedOutputStream(0.0, System.out);
                this.outStream.writeImmediate("final", this.startStream.header());
                // Find the key input columns.
                this.confColumn = this.startStream.findField("confidence");
                this.locColumn = this.startStream.findField("location");
                this.predColumn = this.startStream.findField("predicted");
                // Create the start map.
                this.startMap = new HashMap<Location, GoodStart>();
                // Denote we're ready.
                retVal = true;
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            // For parameter errors, we display the command usage.
            parser.printUsage(System.err);
        } catch (IOException e) {
            System.err.println(e.getMessage());
        }
        return retVal;
    }

    @Override
    public void run() {
        // First we read all the starts, remembering the good ones.
        for (TabbedLineReader.Line startLine : this.startStream) {
            // Is this a good start?
            String prediction = startLine.get(predColumn);
            if (prediction.contentEquals("other")) {
                // No. Write it out.
                this.outStream.write("other", startLine.getAll());
            } else {
                // Parse the location for this start.  The first part (0) is the contig, the second
                // (1) is the position.
                String[] locParts = StringUtils.split(startLine.get(locColumn), ';');
                // Find the ORF containing this start and form a region from it.
                ContigOrfTracker contigOrfs = this.contigMap.get(locParts[0]);
                int location = Integer.parseInt(locParts[1]);
                OrfTracker.StopCodon orfStop = contigOrfs.findOrf(location);
                if (orfStop != null) {
                    // Form a location from the ORF.
                    int origin = orfStop.getLocation() - contigOrfs.getOrfLength() + 1;
                    Location orfLoc = Location.create(locParts[0], "+", origin, orfStop.getLocation());
                    // Get the best start for this ORF.
                    GoodStart bestStart = this.startMap.get(orfLoc);
                    if (bestStart == null) {
                        // This is our first appearance of the ORF.  The new start is automatically the best.
                        this.startMap.put(orfLoc, new GoodStart(startLine, location));
                    } else {
                        // Pick the better start and output the other.
                        String worseLine = bestStart.merge(startLine, location);
                        this.outStream.write("other", worseLine);
                    }
                }
            }
        }
        // All the starts have been processed.  The ones remaining are the good ones.
        for (GoodStart bestStart : this.startMap.values()) {
            this.outStream.write("start", bestStart.getStartLine());
        }
    }
}
