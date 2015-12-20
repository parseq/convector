package deletionsAnalysis;

import java.io.Serializable;
import java.util.ArrayList;
import java.util.Map;

/**
 * Created by german on 07.08.14.
 */
public class Amplicon implements Serializable{
    private ArrayList<Double> coverages;
    private ArrayList<Double> initialCoverages;
    private ArrayList<Double> predictedCoverages = new ArrayList<Double>();
    private String ampliconID;
    private String chr;
    private int start;
    private int end;
    private int thresholdInefficiency;
    private boolean excluded;
    private String exon;
    private static final Integer threshold_homo = 10;


    public Amplicon(String ampliconID, Map<String, ArrayList<String>> infoAboutAmplicon, Map<String, ArrayList<Integer>> coveragesNumbers,
                    Integer thresholdInefficiency) {
        excluded = false;
        this.ampliconID = ampliconID;
        coverages = new ArrayList<Double>();
        this.thresholdInefficiency = thresholdInefficiency;
        initialCoverages = new ArrayList<Double>();

        initializeWithBedMap(infoAboutAmplicon);

        initializeWithData(coveragesNumbers);
    }

    private void initializeWithBedMap(Map<String, ArrayList<String>> infoAboutAmplicon) {
        /*
        @params hashMap with information about amplicon from bed file
         */
        chr = infoAboutAmplicon.get(ampliconID).get(0);
        start = Integer.parseInt(infoAboutAmplicon.get(ampliconID).get(1));
        end = Integer.parseInt(infoAboutAmplicon.get(ampliconID).get(2));
        exon = infoAboutAmplicon.get(ampliconID).get(3);
    }

    private void initializeWithData(Map<String, ArrayList<Integer>> coveragesNumbers) {
        /*
        @params coverages from file with coverages
         */
        ArrayList<Integer> cov = coveragesNumbers.get(ampliconID);
        if (Statistics.sm(cov) < thresholdInefficiency) {
            excluded = true;
        } else {
            for (Integer num : cov) {
                initialCoverages.add(1.0 * num);
                if (num >= threshold_homo) {
                    double rawValue = (Math.log(1.0 * num) / Math.log(2));
                    coverages.add(rawValue);
                } else {
                    coverages.add(0.0);
                }
            }
        }
    }

    public Integer distance(Amplicon otherAmplicon) {
        /*
        @params other amplicon to find distance between them
        @return If amplicon came from other chromosome: return INFTY, else: return distance in BP between their locations
         */
        if (chr.equals(otherAmplicon.chr))
            return Math.min(Math.abs(start - otherAmplicon.end), Math.abs(end - otherAmplicon.start));
        return Integer.MAX_VALUE;
    }

    public boolean isExcluded() {
        /*
        @return if this amplicon was excluded
         */
        return excluded;
    }

    public ArrayList<Double> getCoverages() {
        /*
        @return coverages of amplicon
         */
        return coverages;
    }

    public ArrayList<Double> getShiftedCoverages(double shift) {
        /*
        @return coverages of amplicon shifted by specified value
         */
        ArrayList<Double> shiftedCoverages = new ArrayList<Double>();
        for (int i = 0; i < coverages.size(); i++) shiftedCoverages.add(coverages.get(i) + shift);
        return shiftedCoverages;
    }

    public String returnID() {
        /*
        @return name of the chromosome
         */
        return ampliconID;
    }

    public void modifyCoverages(ArrayList<Double> newCoverages) {
        coverages = newCoverages;
    }

    public void normalize(int i, Double factor) {
        /*
        function for some sorts of normalizations (can be useful for 2nd stage of algorithm)
         */
        Double x = coverages.get(i) / factor;
        coverages.set(i, x);
    }

    public ArrayList<Double> getInitialCoverages() {
        return initialCoverages;
    }

    public String getExon() {
        return exon;
    }


    public void addElemToListOfPredictedCoverages(double elem) {
        predictedCoverages.add(elem);
    }
}