package deletionsAnalysis;

import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

/**
 * Created by german on 09.10.14.
 */
public class OptionsParse {

    public Options OptionParse() {
        /*
         * Create command line options
         * @param void
         * @return options
         */
        Options option = new Options();

        option.addOption(OptionBuilder.withLongOpt("qc")
                .withDescription("chi squared quantile")
                .hasArg()
                .withArgName("qc")
                .create("qc"));

        option.addOption(OptionBuilder.withLongOpt("b")
                .withDescription("path to BED file")
                .hasArg()
                .withArgName("b")
                .create("bed"));

        option.addOption(OptionBuilder.withLongOpt("d")
                .withDescription("path to data (tab delimited) file")
                .hasArg()
                .withArgName("d")
                .create("data"));

        option.addOption(OptionBuilder.withLongOpt("mc")
                .withDescription("minimum correlation threshold ")
                .hasArg()
                .withArgName("mc")
                .create("mc"));

        option.addOption(OptionBuilder.withLongOpt("mnm")
                .withDescription("maximum number of linear models ")
                .hasArg()
                .withArgName("mnm")
                .create("mnm"));

        option.addOption(OptionBuilder.withLongOpt("nne")
                .withDescription("number of models for non-efficiency ")
                .hasArg()
                .withArgName("nne")
                .create("nne"));

        option.addOption(OptionBuilder.withLongOpt("nod")
                .withDescription("number of models for outlier detection ")
                .hasArg()
                .withArgName("nod")
                .create("nod"));

        option.addOption(OptionBuilder.withLongOpt("lb")
                .withDescription("lower bound")
                .hasArg()
                .withArgName("lb")
                .create("lb"));

        option.addOption(OptionBuilder.withLongOpt("ub")
                .withDescription("upper bound")
                .hasArg()
                .withArgName("ub")
                .create("ub"));

        option.addOption(OptionBuilder.withLongOpt("dist")
                .withDescription("minimum distance between amplicons ")
                .hasArg()
                .withArgName("dist")
                .create("dist"));

        option.addOption(OptionBuilder.withLongOpt("lcb")
                .withDescription("minimum bound for low covered samples.")
                .hasArg()
                .withArgName("lcb")
                .create("lcb"));

        option.addOption(OptionBuilder.withLongOpt("lca")
                .withDescription("minimum bound for low covered amplicon")
                .hasArg()
                .withArgName("lca")
                .create("lca"));

        option.addOption(OptionBuilder.withLongOpt("f")
                .withDescription("output file")
                .hasArg()
                .withArgName("f")
                .create("f"));

        option.addOption(OptionBuilder.withLongOpt("lss")
                .withDescription("learning sample size")
                .hasArg()
                .withArgName("lss")
                .create("lss"));

        option.addOption(OptionBuilder.withLongOpt("h")
                .withDescription("help message")
                .hasArg()
                .withArgName("help")
                .create());

        return option;
    }

}
