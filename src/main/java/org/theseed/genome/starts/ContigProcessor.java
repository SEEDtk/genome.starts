/**
 *
 */
package org.theseed.genome.starts;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
import java.util.Random;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.genome.Contig;
import org.theseed.genome.Genome;
import org.theseed.genome.GenomeDirectory;
import org.theseed.io.BalancedOutputStream;
import org.theseed.io.Shuffler;
import org.theseed.locations.LocationList;
import org.theseed.utils.ICommand;

/**
 * This class generates training data for the start-caller.  The positional parameter is the name of a directory containing
 * a GTO file and a predicted-stop file for each genome.  The GTO files will have the suffix ".gto", and the predicted-stop
 * files the suffix ".stops.tbl", with both files using the genome ID as the base part of the name.
 *
 * The objective is to produce a training set for a deep learning module to call starts.
 *
 * The command-line options are as follows.
 *
 * -b	indicates the output should be balanced; the records are held in memory and then
 * 		a distributed, balanced subset is output; the value should be a number from 1.0
 * 		to 2.0, indicating the maximum number of output records per class as a fraction of
 * 		the smallest class's size
 * -m	maximum number of starts to output per genome; the default is 400
 *
 * @author Bruce Parrello
 *
 */
public class ContigProcessor implements ICommand {

    // FIELDS
    /** genome directory object */
    private GenomeDirectory inputDir;


    // COMMAND LINE

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** TRUE if we want progress messages */
    @Option(name="-v", aliases={"--verbose", "--debug"}, usage="display progress on STDERR")
    private boolean debug;

    /** balanced output fuzz factor */
    @Option(name="-b", aliases={"--balance", "--fuzz"}, metaVar="1.2", usage="specify class-balanced output")
    private double fuzzFactor;

    /** output size control */
    @Option(name="-m", aliases={"--max", "--maximum"}, metaVar="1000",
            usage="maximum number of potential starts to write per genome")
    private int maxOut;

    /** name of the directory containing the input files */
    @Argument(metaVar="gtoDir", usage="directory containing .gto and .stops.tbl files")
    private File genomeDir;

    /**
     * Parse command-line options to specify the parameters of this object.
     *
     * @param args	an array of the command-line parameters and options
     *
     * @return TRUE if successful, FALSE if the parameters are invalid
     */
    public boolean parseCommand(String[] args) {
        boolean retVal = false;
        // Set the defaults.
        this.help = false;
        this.debug = false;
        this.maxOut = 400;
        // Parse the command line.
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
                // Validate the genome directory.
                if (! this.genomeDir.isDirectory()) {
                    throw new FileNotFoundException(genomeDir.getPath() + " is not a valid directory.");
                }
                // Validate the fuzz factor.
                if (this.fuzzFactor != 0 && (this.fuzzFactor < 1.0 || this.fuzzFactor > 2.0)) {
                    throw new IllegalArgumentException("Balance factor must be 0 (off) or between 1.0 and 2.0 inclusive.");
                }
                // Load the genome directory.
                this.inputDir = new GenomeDirectory(this.genomeDir);
                new Random();
                // We made it this far, we can run the application.
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

    public void run() {
        try {
            // Start the output with a header line.
            BalancedOutputStream outStream = new BalancedOutputStream(this.fuzzFactor, System.out);
            outStream.writeImmediate("frame", StartCodonFinder.StartCodon.header());
            // Loop through the genomes.
            for (Genome genome : this.inputDir) {
                if (this.debug) System.err.format("Processing genome %s: %s.%n", genome.getId(), genome.getName());
                // Get the stop prediction file.
                File stopFile = new File(this.genomeDir, genome.getId() + ".stops.tbl");
                Map<String, ContigOrfTracker> contigMap = ContigOrfTracker.readPredictionFile(stopFile);
                Map<String, LocationList> locMap = LocationList.createGenomeCodingMap(genome);
                // Now loop through the contigs, accumulating starts.
                Shuffler<StartCodonFinder.StartCodon> buffer = new Shuffler<StartCodonFinder.StartCodon>(200000);
                for (Contig contig : genome.getContigs()) {
                    ContigOrfTracker orfTracker = contigMap.get(contig.getId());
                    // Need to check for NULL in case the contig had no potential stops in it.
                    if (orfTracker == null) {
                        if (this.debug) System.err.format("Skipping contig %s in genome %s-- no stops.%n", contig.getId(), genome.getId());
                    } else {
                        if (this.debug) System.err.format("Processing contig %s in genome %s.%n", contig.getId(), genome.getId());
                        // Loop through the start codons, writing them out.
                        StartCodonFinder startFinder = new StartCodonFinder(contig.getId(), contig.getSequence(), orfTracker);
                        buffer.addSequence(startFinder);
                    }
                }
                // Next, we pick random starts from the genome.
                if (this.debug) System.err.format("Shuffling %d starts from genome.%n", buffer.size());
                buffer.shuffle(this.maxOut);
                // Finally, write the chosen starts.
                if (this.debug)System.err.println("Writing the chosen starts.");
                Iterator<StartCodonFinder.StartCodon> outStarts = buffer.limitedIter(this.maxOut);
                while (outStarts.hasNext()) {
                    StartCodonFinder.StartCodon start = outStarts.next();
                    String contigId = start.getContigId();
                    LocationList contigLocs = locMap.get(contigId);
                    LocationList.Edge edge = contigLocs.isEdge(start.getLoc(), false);
                    outStream.write(edge.toString(), start.toString());
                }
            }
            // Close and flush the output.
            outStream.close();
        } catch (Exception e) {
            e.printStackTrace(System.err);
        }
    }

}
