/**
 *
 */
package org.theseed.genome.starts;

import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;
import org.theseed.utils.ICommand;

/**
 * This class produces a prediction file for input to the start-calling neural net.  The prediction file contains the
 * input parameters expected by the neural net, but it also includes the expected prediction and the functional assignment
 * for each start's protein.  This helps us to determine how useful the start-caller is.
 *
 * The positional parameters are the name of a GTO file containing the genome, and the name of the file containing the
 * start and stop predictions from the start/stop predictor.
 *
 * @author Bruce Parrello
 *
 */
public class GenomeProcessor implements ICommand {

    // COMMAND LINE

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** TRUE if we want progress messages */
    @Option(name="-v", aliases={"--verbose", "--debug"}, usage="display progress on STDERR")
    private boolean debug;




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
        // Parse the command line.
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
                // TODO parsing body for TEST case
                // We made it this far, we can run the application.
                retVal = true;
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            // For parameter errors, we display the command usage.
            parser.printUsage(System.err);
//        } catch (IOException e) {
//            System.err.println(e.getMessage());
        }
        return retVal;
    }

    public void run() {
        // TODO execute body for TEST case
    }
}
