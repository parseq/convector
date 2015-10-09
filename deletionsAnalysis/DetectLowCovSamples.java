package deletionsAnalysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by german on 15.08.14.
 */
public class DetectLowCovSamples {
    /*
    Class for detection and report of low covered samples
     */
    HashMap<Integer, Double> totals = new HashMap<Integer, Double>();
    ArrayList<Integer> listToExclude = new ArrayList<Integer>();
    int lowCovBound;
    private static final Logger log = Logger.getLogger( Solver.class.getName() );


    public DetectLowCovSamples(ArrayList<String> samplesNames, HashMap<String, Amplicon> dataAboutAmplicons, Integer lowCovBound) {
        this.lowCovBound = lowCovBound;
        findLowCoveredSample(samplesNames, dataAboutAmplicons);
    }

    public String returnTotalCoverage(ArrayList<String> samplesNames) {
        /*
        @params list of samples in dataset
        @return string for the top of the file (coverages)
         */
        String totalCoverages = "";
        for (String name : samplesNames) {
            totalCoverages += (totals.get(samplesNames.indexOf(name))).intValue() + "\t";
        }
        return totalCoverages;
    }

    private void findLowCoveredSample(ArrayList<String> samplesNames, HashMap<String, Amplicon> dataAboutAmplicons) {
        for (int i = 0; i < samplesNames.size(); i++) {
            totals.put(i, 0.0);
        }

        for (Map.Entry<String, Amplicon> entry : dataAboutAmplicons.entrySet()) {
            for (int i = 0; i < entry.getValue().getInitialCoverages().size(); i++) {
                Double x = totals.get(i);
                x += ((entry.getValue().getInitialCoverages().get(i)));
                totals.put(i, x);
            }
        }

        for (Map.Entry<Integer, Double> entry : totals.entrySet()) {
            if (entry.getValue() < lowCovBound) {
                listToExclude.add(entry.getKey());
                log.log(Level.WARNING,  "Sample " + samplesNames.get(entry.getKey()) + " too low covered " + entry.getValue());
            }
        }
    }


}
