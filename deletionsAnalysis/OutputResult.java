package deletionsAnalysis;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.TreeMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by german on 16.10.14.
 */
public class OutputResult {
    /*
    just a class for useful and handy output.
    @params all the data (except data about homozygous amplicons and low covered samples).
    @return file with results for further analysis.
     */
    private ArrayList<String> dirtySamples = new ArrayList<String>();
    private static final Logger log = Logger.getLogger( Solver.class.getName() );


    protected OutputResult(String filename, ArrayList<String> listOfExcludedAmpls, ArrayList<String> nonEffecitveAmpls,
                           HashMap<String, HashMap<String, Integer>> sampleToAmplicons, String versionAndParams,
                           ArrayList<String> samplesNames, ReaderBED infoAboutAmplicons,
                           Integer minNumOfModelsForOutlierDetection,
                           String totalCoverages, String datasetFile) throws IOException {
        String outFileName = filename + ".xls";
        Path file = Paths.get(outFileName);
        try {
            // Create the empty file with default permissions, etc.
            Files.createFile(file);
        } catch (FileAlreadyExistsException x) {
            log.log(Level.WARNING, "File named %s" +
                            " already exists. For accurate analysis, please, name your files differently.\n" +
                            "Now, please, repeat the analysis and choose the name for output file wisely. Now - overwriting.", file
            );
        } catch (IOException x) {
            // Some other sort of failure, such as permissions.
            log.log(Level.SEVERE, "createFile error: %s%n", x);
        }

        FileWriter fw = new FileWriter(outFileName);
        BufferedWriter bw = new BufferedWriter(fw);

        bw.write("Description\t" + versionAndParams + "\n");
        bw.write("Dataset was taken from file\t" + datasetFile + "\n");
        bw.write("BED file used from file\t" + infoAboutAmplicons.bedFileName + "\n");

        bw.newLine();
        String topString = "Chr\tStartpos\tStoppos\tAmplicon name\tDescription\t";
        for (String str : samplesNames) {
            topString += str;
            topString += "\t";
        }
        bw.write(topString + "\n");
        bw.write("\t\t\t\tTotal reads:\t" + totalCoverages + "\n");

        AmpliconSorter bvc = new AmpliconSorter(infoAboutAmplicons.getInfoAboutAmplicon());

        TreeMap<String, ArrayList<String>> sorted_map = new TreeMap<String, ArrayList<String>>(bvc);
        sorted_map.putAll(infoAboutAmplicons.getInfoAboutAmplicon());
        int delCounter = 0;
        int dupCounter = 0;
        for (Map.Entry<String, ArrayList<String>> e : sorted_map.entrySet()) {
            String tmpString = "";
            tmpString += e.getValue().get(0) + "\t";
            tmpString += e.getValue().get(1) + "\t";
            tmpString += e.getValue().get(2) + "\t";
            tmpString += e.getKey() + "\t";
            tmpString += e.getValue().get(3) + "\t";


            if (listOfExcludedAmpls.contains(e.getKey())) {
                for (String name : samplesNames) {
                    tmpString += "EX" + "\t";
                }
            } else if (nonEffecitveAmpls.contains(e.getKey())) {
                for (String name : samplesNames) {
                    tmpString += "NE" + "\t";
                }
            } else {
                for (String name : samplesNames) {
                    if (sampleToAmplicons.containsKey(name)) {
                        if (sampleToAmplicons.get(name).containsKey(e.getKey())) {
                            if (sampleToAmplicons.get(name).get(e.getKey()) >= minNumOfModelsForOutlierDetection) {
                                tmpString += "DEL" + "\t";
                                if (!dirtySamples.contains(name)) {
                                    dirtySamples.add(name);
                                }
                                delCounter++;
                            } else if (sampleToAmplicons.get(name).get(e.getKey()) < -minNumOfModelsForOutlierDetection) {
                                tmpString += "^\t";
                                if (!dirtySamples.contains(name)) {
                                    dirtySamples.add(name);
                                }
                                dupCounter++;
                            } else {
                                tmpString += " " + "\t";
                            }
                        }
                    }
                }
            }
            bw.write(tmpString + "\n");
            bw.flush();
        }
        bw.write("\n");
        bw.write("Total deletions: " + delCounter + "\n");
        bw.write("Total duplications: " + dupCounter + "\n");

        bw.close();
    }

    public ArrayList<String> getDirtySamples() {
        return dirtySamples;
    }
}

