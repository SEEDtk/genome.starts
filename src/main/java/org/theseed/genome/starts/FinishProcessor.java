/**
 *
 */
package org.theseed.genome.starts;

import java.io.File;
import java.io.IOException;
import java.util.Map;

import org.kohsuke.args4j.Argument;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.io.TabbedLineReader;
import org.theseed.utils.ICommand;

/**
 * This process reads the output from the Combo neural net (rank starts and stops based on neighborhood) and the FancyStart neural net
 * (rank starts based on statistical characteristics) and picks the best-ranked start for each ORF.
 *
 * The two positional parameters are the Combo output file and the FancyStart output file, respectively.  The FancyStart output file
 * is echoed output with an additional column ("final") indicating if the start is the best one or not.
 *
 * @author Bruce Parrello
 *
 */
public class FinishProcessor implements ICommand {

    // FIELDS
    /** map of contig IDs to ORFs */
    Map<String, ContigOrfTracker> contigMap;
    /** input file for called starts */
    TabbedLineReader startStream;

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
                // Open the stop file to begin input.
                this.startStream = new TabbedLineReader(this.startFile);
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
        // TODO Auto-generated method stub

    }

}
