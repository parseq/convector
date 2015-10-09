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
public class ReaderBED {
    /*
    class for parsing BED-formatted files (tab delimited).
    Example:
        track name=The_Best_Panel
        chr0	111111111	111111222	AMPL001DUMMY	TMIG exon 1	N/A
        chr1	2222222	333333333	AMPL002DUMMY	AMIG exon 1	N/A
     */
    private Map<String, ArrayList<String>> infoAboutAmplicon = new HashMap<String, ArrayList<String>>();
    public String bedFileName = "Unknown BED file";
    private static final Logger log = Logger.getLogger( Solver.class.getName() );


    ReaderBED(String bedFileName) {
        this.bedFileName = bedFileName;
        try {
            BufferedReader br = new BufferedReader(new FileReader(bedFileName));
            for (String line; (line = br.readLine()) != null; ) {
                if (line.startsWith("chr") || line.startsWith("clean") || line.startsWith("del")) {
                    String[] tmpListForMap = line.split("\\t");
                    ArrayList values = new ArrayList<String>();
                    for (int i = 0; i != tmpListForMap.length; i++)
                        if (i != 3) {
                            values.add(tmpListForMap[i]);
                        }
                    infoAboutAmplicon.put(tmpListForMap[3], values);
                }
            }
        } catch (FileNotFoundException e) {
            log.log(Level.SEVERE, "\nSorry, but your file with information about amplicons were not found.");
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

    public Map<String, ArrayList<String>> getInfoAboutAmplicon() {
        return infoAboutAmplicon;
    }


    public void removeAmplFromBed(String amplName) {
        if (infoAboutAmplicon.containsKey(amplName)) {
            infoAboutAmplicon.remove(amplName);
        }
    }
}
