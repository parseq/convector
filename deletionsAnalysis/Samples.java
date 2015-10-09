package deletionsAnalysis;

import java.io.Serializable;
import java.util.HashMap;
import java.util.Map;
import java.util.ArrayList;

/**
 * Created by german on 15.11.14.
 */
public class Samples implements Serializable {
    /*
    Class for containing coverages (transformed and intiial)
     */
    private HashMap<String, HashMap<String, Double>> coverageOfAmpliconsInSample = new HashMap<String, HashMap<String, Double> >();
    private HashMap<String, HashMap<String, Double>> initialCoverageOfAmpliconsInSample = new HashMap<String, HashMap<String, Double> >();
    private ArrayList<String> finallyAmplNames;

    public Samples(ArrayList<String> samplesNames, Map<String, Amplicon> dataAboutAmplicons, ArrayList<String> nonEffecitveAmpls) {
        /*
        Constructor that creates object with coverages
        @param list of names of samples, map of {Amplicon name : Amplicon object}, list of amplicons that were non effective (NE)
         */
        for (int i = 0; i < samplesNames.size(); i++) {
            coverageOfAmpliconsInSample.put(samplesNames.get(i), new HashMap<String, Double>());
            initialCoverageOfAmpliconsInSample.put(samplesNames.get(i), new HashMap<String, Double>());
        }

        finallyAmplNames = new ArrayList<String>();
        for (Map.Entry en : dataAboutAmplicons.entrySet()) {
            if (!nonEffecitveAmpls.contains((String)en.getKey())) {
                Amplicon tmpAmpl = (Amplicon) en.getValue();
                finallyAmplNames.add(tmpAmpl.returnID());
            }
        }
        for (int i = 0; i < samplesNames.size(); i++) {
            Map<String, Double> covOfAmplInSample = coverageOfAmpliconsInSample.get(samplesNames.get(i));
            Map<String, Double> initCovOfAmplInSample = initialCoverageOfAmpliconsInSample.get(samplesNames.get(i));

            for (Map.Entry en : dataAboutAmplicons.entrySet()) {
                if (!nonEffecitveAmpls.contains((String)en.getKey())) {
                    Amplicon tmpAmpl = (Amplicon) en.getValue();
                    covOfAmplInSample.put((String) en.getKey(), tmpAmpl.getCoverages().get(i));
                    initCovOfAmplInSample.put((String) en.getKey(), tmpAmpl.getCoverages().get(i));
                }
            }
        }
    }

    public Double returnCov(String sample, String ampl) {
        /*
        @param name of sample and name of amplicon
        @return coverage (transformed)
         */
        return coverageOfAmpliconsInSample.get(sample).get(ampl);
    }

    public Double returnInitCov(String sample, String ampl) {
        /*
        @param name of sample and name of amplicon
        @return coverage (initial raw coverage)
         */
        return initialCoverageOfAmpliconsInSample.get(sample).get(ampl);
    }


    public void normalizeCov(String sample, String ampl, double factor) {
        /*
        @param name of sample and name of amplicon and normalization value
        @return transform coverages and put it in coverageOfAmpliconsInSample
         */
        double value = (coverageOfAmpliconsInSample.get(sample).get(ampl));

        coverageOfAmpliconsInSample.get(sample).put(ampl, (value - factor));
    }

    public ArrayList<String> getFinallyAmplNames() {
        /*
        @return name of amplicons that are suitable for analysis
         */
        return finallyAmplNames;
    }

}
