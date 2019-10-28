package org.theseed.genome.starts;

import java.util.Arrays;

import org.theseed.utils.ICommand;

/**
 * This is an application to generate input and training sets for a gene start detector.
 *
 * The supported commands are
 *
 * 	train	generate a training set
 * 	predict	generate an input set for finding starts
 * 	test	generate an input set for verification
 *  finish  compute final starts for each ORF
 *
 */
public class App
{
    public static void main( String[] args )
    {
        // Get the control parameter.
        String command = args[0];
        String[] newArgs = Arrays.copyOfRange(args, 1, args.length);
        ICommand processor;
        // Parse the parameters.
        switch (command) {
        case "train" :
            processor = new ContigProcessor();
            break;
        case "predict" :
            processor = new FastaProcessor();
            break;
        case "test" :
            processor = new GenomeProcessor();
            break;
        case "finish" :
            processor = new FinishProcessor();
            break;
        default :
            throw new RuntimeException("Invalid command " + command + ": must be \"train\", \"test\", \"predict\", or \"finish\".");
        }
        boolean ok = processor.parseCommand(newArgs);
        if (ok) {
            processor.run();
        }
    }
}
