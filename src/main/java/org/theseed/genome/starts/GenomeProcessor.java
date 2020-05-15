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
import org.theseed.genome.Contig;
import org.theseed.genome.Genome;
import org.theseed.proteins.RoleMap;
import org.theseed.utils.ICommand;

/**
 * This class produces a prediction file for input to the start-calling neural net.  The prediction file contains the
 * input parameters expected by the neural net, but it also includes the expected prediction and the functional assignment
 * for each start's protein.  This helps us to determine how useful the start-caller is.
 *
 * The positional parameters are the name of a GTO file containing the genome, and the name of the file containing the
 * start and stop predictions from the start/stop predictor.
 *
 * The command-line options are as follows:
 *
 *  -v	display progress messages on STDERR
 *
 *  -r	name of a file containing a mapping from role IDs to names; if this option is specified,
 *  	the output will contain role IDs instead of functions for each true start
 *
 * @author Bruce Parrello
 *
 */
public class GenomeProcessor implements ICommand {

	// FIELDS
	/** map of role IDs to role names, or NULL if the output is just to contain functions */
	private RoleMap roleMap;
	/** genome loaded from the GTO file */
	private Genome genome;
	/** map from contig IDs to ORF trackers */
	private Map<String, ContigOrfTracker> orfTrackerMap;
	/** map from starts to roles */
	private Map<String, ContigStarts> startRoleMap;

    // COMMAND LINE

    /** help option */
    @Option(name="-h", aliases={"--help"}, help=true)
    private boolean help;

    /** TRUE if we want progress messages */
    @Option(name="-v", aliases={"--verbose", "--debug"}, usage="display progress on STDERR")
    private boolean debug;

    /** role ID mapping file */
    @Option(name="-r", aliases={"--roles"}, usage="name of file containing the role ID map")
    private File roleFile;

    /** name of the genome file */
    @Argument(index=0, metaVar="genomeFile", usage="name of genome GTO file")
    private File genomeFile;

    /** name of the start/stop call file */
    @Argument(index=1, metaVar="startStopFile", usage="file of predictions from start/stop caller")
    private File startStopFile;

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
        this.roleFile = null;
        // Parse the command line.
        CmdLineParser parser = new CmdLineParser(this);
        try {
            parser.parseArgument(args);
            if (this.help) {
                parser.printUsage(System.err);
            } else {
            	// Create the role map.
            	if (this.roleFile == null) {
            		this.roleMap = null;
            		if (debug) System.err.println("Functional assignments will be output.");
            	} else {
            		this.roleMap = RoleMap.load(this.roleFile);
            		if (debug) System.err.format("%d roles loaded from %s%n", this.roleMap.fullSize(),
            				this.roleFile);
            	}
                // Load the genome.
            	this.genome = new Genome(this.genomeFile);
            	if (debug) System.err.format("%s loaded from file %s%n", this.genome,
            			this.genomeFile);
            	// Load the stop prediction file.
            	this.orfTrackerMap = ContigOrfTracker.readPredictionFile(this.startStopFile);
            	if (debug) System.err.println("ORF predictions loaded from " +
            			this.startStopFile + ".");
            	// Load the role maps.
            	this.startRoleMap = ContigStarts.contigStartsMap(this.genome, '+', this.roleMap);
            	if (debug) System.err.println("Start roles computed.");
                // We made it this far, we can run the application.
                retVal = true;
            }
        } catch (CmdLineException e) {
            System.err.println(e.getMessage());
            // For parameter errors, we display the command usage.
            parser.printUsage(System.err);
        } catch (IOException e) {
            e.printStackTrace(System.err);
        }
        return retVal;
    }

    public void run() {
    	// Start the output file.
    	System.out.println(StartCodonFinder.StartCodon.header() + "\texpect\troles");
    	// Loop through the contigs.
    	for (Contig contig : this.genome.getContigs()) {
    		String contigId = contig.getId();
    		if (debug) System.err.println("Processing contig " + contigId);
    		ContigStarts roleMapper = this.startRoleMap.get(contigId);
    		ContigOrfTracker orfTracker = this.orfTrackerMap.get(contigId);
    		String contigSeq = contig.getSequence();
    		StartCodonFinder startFinder = new StartCodonFinder(contigId, contigSeq, orfTracker);
    		// Loop through the starts in the contig.
    		for (StartCodonFinder.StartCodon start : startFinder) {
    			String expect;
    			String roles = roleMapper.getRoles(start.getLoc());
    			if (roles == null) {
    				expect = "other";
    				roles = "";
    			} else {
    				expect = "start";
    			}
    			System.out.println(start.toString() + "\t" + expect + "\t" + roles);
    		}
    	}
    }
}
