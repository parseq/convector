package deletionsAnalysis;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by german on 16.10.14.
 */
public class ReaderData {
    /*
    class for parsing files with coverages (tab delimited). Example:
        Gene    Target	Sample1	Sample2	Sample3
        N/A     AMP0001	    174	    146	    122
        N/A	    AMPL002	    405 	344	    245
     */
    private Map<String, ArrayList<Integer>> coverages = new HashMap<String, ArrayList<Integer>>();
    private ArrayList<String> samplesNames = new ArrayList<String>();
    private static final Logger log = Logger.getLogger( Solver.class.getName() );


    ReaderData(String dataFileName) {
        try {
            BufferedReader br = new BufferedReader(new FileReader(dataFileName));
            for (String line; (line = br.readLine()) != null; ) {
                if (line.startsWith("Gene\tTarget")) {
                    String[] tmpListForMap = line.split("\\t");
                    for (int i = 2; i < tmpListForMap.length; i++) {
                        samplesNames.add(tmpListForMap[i]);
                    }
                } else {
                    String[] tmpListForMap = line.split("\\t");
                    ArrayList values = new ArrayList<String>();

                    String tmpAmplID = tmpListForMap[1];

                    for (int i = 2; i != tmpListForMap.length; i++)
                        values.add(Integer.parseInt(tmpListForMap[i]));
                    coverages.put(tmpAmplID, values);
                }
            }
            log.log(Level.FINE, "We have found " + samplesNames.size() + " samples.");

        } catch (FileNotFoundException e) {
            log.log(Level.SEVERE, "\nSorry, but your file with information about coverages were not found.");
            log.log(Level.SEVERE, e.getMessage());
            log.log(Level.SEVERE, "Please, be sure that your file exists.");
            System.exit(1);
        } catch (IOException e) {
            log.log(Level.SEVERE, "\nSorry, but we had some problems with reading from your bed file.");
            log.log(Level.SEVERE, e.getMessage());
            log.log(Level.SEVERE, "Please, be sure that this tool has access to your input files.");
            System.exit(1);
        }
    }

    public Map<String, ArrayList<Integer>> getCoverages() {
        return coverages;
    }

    public ArrayList<String> getSamplesNames() {
        return samplesNames;
    }
}