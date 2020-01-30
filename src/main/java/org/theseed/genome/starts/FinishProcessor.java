/**
 *
 */
package org.theseed.genome.starts;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.HashMap;
import java.util.Map;

import org.apache.commons.lang3.StringUtils;
import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.io.TabbedLineReader;
import org.theseed.locations.Location;
import org.theseed.utils.ICommand;

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
        /** ORF location */
        private Location orfLoc;

        /**
         * Create a new start record.
         *
         * @param line	original input line for the start
         * @param loc	location of the start in the contig
         */
        public GoodStart(TabbedLineReader.Line line, int loc, Location orfLoc) {
            this.startLine = line;
            this.confidence = line.getDouble(confColumn);
            this.location = loc;
            this.orfLoc = orfLoc;
        }

        /**
         * Compare a new start to this one.  If it is better, store it.
         *
         * @param line	input line for the new start
         * @param loc	location of the new start in the contig
         *
         * @return the input line for the worse start.
         */
        public GoodStart merge(TabbedLineReader.Line line, int loc, Location orfLoc) {
            // We default to returning the new location.
            GoodStart retVal;
            GoodStart buffer = new GoodStart(line, loc, orfLoc);
            double conf = line.getDouble(confColumn);
            if (conf > this.confidence || conf == this.confidence && loc < this.location) {
                // Here the new location is better. Return the old one and store the new one.
                buffer.startLine = this.startLine;
                buffer.confidence = this.confidence;
                buffer.location = this.location;
                buffer.orfLoc = this.orfLoc;
                retVal = buffer;
                this.startLine = line;
                this.confidence = conf;
                this.location = loc;
                this.orfLoc = orfLoc;
            } else {
                // Here the old location is better.  Return the new one and keep this one in place.
                retVal = buffer;
            }
            return retVal;
        }

        /**
         * @return the original input line for the start
         */
        public String getStartLine() {
            return startLine.getAll();
        }

        /**
         * @return the confidence
         */
        public double getConfidence() {
            return confidence;
        }

        /**
         * @return the start location
         */
        public int getLocation() {
            return location;
        }

        /**
         * @return the contig ID
         */
        public String getContig() {
            return orfLoc.getContigId();
        }

        /**
         * @return the stop location
         */
        public int getStopLoc() {
            return orfLoc.getEnd();
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
    private PrintWriter outStream;
    /** map of ORFs to best starts */
    private Map<Location, GoodStart> startMap;
    /** writer for output */
    private Write writer;

    // OUTPUT WRITERS

    private abstract class Write {

        public abstract void header();

        public abstract void dataLine(String mode, GoodStart start);

    }

    private class NormalWrite extends Write {

        @Override
        public void header() {
            outStream.println("final\t" + startStream.header());
        }

        @Override
        public void dataLine(String mode, GoodStart start) {
            outStream.println(mode + "\t" + start.getStartLine());
        }

    }

    private class AltWriter extends Write {

        @Override
        public void header() {
            outStream.println("contig\tstart\tstop\tconfidence\tstrand\ttype");
        }

        @Override
        public void dataLine(String mode, GoodStart start) {
            if (mode.contentEquals("start")) {
                outStream.format("%s\t%d\t%d\t%8.6f\t+\tCDS%n", start.getContig(),
                        start.getLocation(), start.getStopLoc(), start.getConfidence());
            }
        }

    }

    // COMMAND LINE

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** TRUE if we want progress messages */
    @Option(name="-v", aliases={"--verbose", "--debug"}, usage="display progress on STDERR")
    private boolean debug;

    /** TRUE if we want the alternate output format */
    @Option(name="-a", aliases={"--alt"}, usage="alternate output format with ORF data")
    private boolean altFormat;

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
                // Create the writer.
                this.writer = (this.altFormat ? this.new AltWriter() : this.new NormalWrite());
                // Open the start file to begin input and initialize the output file.
                this.startStream = new TabbedLineReader(this.startFile);
                this.outStream = new PrintWriter(System.out);
                this.writer.header();
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
            // Parse the location for this start.  The first part (0) is the contig, the second
            // (1) is the position.
            String[] locParts = StringUtils.split(startLine.get(locColumn), ';');
            // Find the ORF containing this start and form a region from it.
            ContigOrfTracker contigOrfs = this.contigMap.get(locParts[0]);
            int location = Integer.parseInt(locParts[1]);
            OrfTracker.StopCodon orfStop = contigOrfs.findOrf(location);
            Location orfLoc = null;
            if (orfStop != null) {
                // Form a location from the ORF.
                int origin = orfStop.getLocation() - contigOrfs.getOrfLength() + 1;
                orfLoc = Location.create(locParts[0], "+", origin, orfStop.getLocation());
            }
            // Is this a good start?
            String prediction = startLine.get(predColumn);
            if (prediction.contentEquals("other") || orfLoc == null) {
                // No. Write it out.
                this.writer.dataLine("other", new GoodStart(startLine, location, orfLoc));
            } else {
                // Get the best start for this ORF.
                GoodStart bestStart = this.startMap.get(orfLoc);
                if (bestStart == null) {
                    // This is our first appearance of the ORF.  The new start is automatically the best.
                    this.startMap.put(orfLoc, new GoodStart(startLine, location, orfLoc));
                } else {
                    // Pick the better start and output the other.
                    GoodStart worseStart = bestStart.merge(startLine, location, orfLoc);
                    this.writer.dataLine("other", worseStart);
                }
            }
        }
        // All the starts have been processed.  The ones remaining are the good ones.
        for (GoodStart bestStart : this.startMap.values()) {
            this.writer.dataLine("start", bestStart);
        }
        this.outStream.close();
    }

}
