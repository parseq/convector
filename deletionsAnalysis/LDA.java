package deletionsAnalysis;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Serializable;
import java.nio.file.FileAlreadyExistsException;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.apache.commons.math3.distribution.ChiSquaredDistribution;

/**
 * Created by german on 13.11.14.
 */
public class LDA implements Serializable{
    private Samples samplesToLDA; // object of Sample class, contains coverages of samples (log transformed and initial)
    private ArrayList<String> dirtySamples; // names of samples that were detected as containing CNVs
    private Map<String, Double> diagMatrixOfCovariancesMap = new HashMap<String, Double>(); // robust regularized matrix of covariances, scale estimations
    private Map<String, Double> medians = new HashMap<String, Double>(); // estimations of location for coverages
    private ArrayList<String> namesOfSamples; // string with names of all samples
    private ArrayList<String> finalAmpls; // list of amplicons that are suitable for the analysis
    private ArrayList<String> newDirtySamples = new ArrayList<String>(); // new names of samples that were detected as containing CNVs
    public boolean numberOfDirtySamplesChanged = false; // check if the total number of dirty samples decreased on each step
    private double qChisq = 0.0; // quantile of chi squared distribution with 1 degree of freedom
    private HashMap<String, ArrayList<String>> ampliconNamesFromExons; // map, {exon : [amplicons' names]}
    private double logDup = Math.log(1.5) / Math.log(2); // precalculated value for distance between normal and duplicated amplicons
    private double logDel = 1.0; // Math.log(2) / Math.log(2)
    private HashMap<String, Amplicon> dataAboutAmplicons; // {amplicon name : object of class Amplicon}
    private ArrayList<String> nonEffectiveAmpls; // list of non effective amplicons
    private static final Logger log = Logger.getLogger( Solver.class.getName() );
    private int numberOfTests = 0;
    private double alpha = 0.001;
    private ArrayList<Double> pValues = new ArrayList<Double>();
    private HashMap<String, HashMap<String, Double>> distanceToDel = new HashMap<String, HashMap<String, Double>>();
    private HashMap<String, HashMap<String, Double>> distanceToDup = new HashMap<String, HashMap<String, Double>>();
    private HashMap<String, HashMap<String, Double>> distanceToNorm = new HashMap<String, HashMap<String, Double>>();




    public LDA(Samples samples, HashMap<String, Amplicon> dataAboutAmplicons, ArrayList<String> namesOfSamples,
               ArrayList<String> dirtySamples, ArrayList<String> nonEffectiveAmpls,
               HashMap<String, PriorityQueue<Pair>> amplCorrelationPriorityForLDA, Double minCorThreshold, Integer minDistance,
               int minNumOfModels, int numberOfTestsPerSample, HashMap<String, HashMap<String, Double>> predictedMedians, Integer numberOfTests) {
        this.numberOfTests = numberOfTests;
        org.apache.commons.math3.distribution.ChiSquaredDistribution chiSq = new org.apache.commons.math3.distribution.ChiSquaredDistribution(1);
        this.qChisq = chiSq.inverseCumulativeProbability(1 - alpha / numberOfTests);
        this.dirtySamples = dirtySamples;
        this.finalAmpls = samples.getFinallyAmplNames();
        this.namesOfSamples = namesOfSamples;
        this.ampliconNamesFromExons = new HashMap<String, ArrayList<String>>();
        this.dataAboutAmplicons = dataAboutAmplicons;
        this.nonEffectiveAmpls = nonEffectiveAmpls;


        for (Map.Entry key : dataAboutAmplicons.entrySet()) {
            Amplicon amplToReadInfoAboutExons = (Amplicon) key.getValue();
            if (!ampliconNamesFromExons.containsKey(amplToReadInfoAboutExons.getExon()) && !nonEffectiveAmpls.contains(amplToReadInfoAboutExons.returnID())) {
                ampliconNamesFromExons.put(amplToReadInfoAboutExons.getExon(), new ArrayList<String>());
            }
        }

        for (Map.Entry key : dataAboutAmplicons.entrySet()) {
            Amplicon amplToReadInfoAboutExons = (Amplicon) key.getValue();
            if (!nonEffectiveAmpls.contains(amplToReadInfoAboutExons.returnID())) {
                ampliconNamesFromExons.get(amplToReadInfoAboutExons.getExon()).add(amplToReadInfoAboutExons.returnID());
            }
        }

        for (Map.Entry key : amplCorrelationPriorityForLDA.entrySet()) {
            String amplName = (String) key.getKey();

            if (!nonEffectiveAmpls.contains(amplName)) {

                PriorityQueue<Pair> priorityByCorrelation = (PriorityQueue<Pair>) key.getValue();

                double tmpCor = 1.0;

                ArrayList<Double> factorsForSampleAndAmpl = new ArrayList<Double>();
                for (int i = 0; i < namesOfSamples.size(); i++) {
                    factorsForSampleAndAmpl.add(0.0);
                }

                Pair e = priorityByCorrelation.poll();

                ArrayList<ArrayList<Double>> toCalculateMean = new ArrayList<ArrayList<Double>>();
                for (int i = 0; i < namesOfSamples.size(); i++) {
                    toCalculateMean.add(new ArrayList<Double>());
                }

                int counter_of_values_to_compare_with = 0;

                while ((e.getCorrelation() > minCorThreshold && counter_of_values_to_compare_with < minNumOfModels) && !priorityByCorrelation.isEmpty()) {
                    e = priorityByCorrelation.poll();
                    tmpCor = e.getCorrelation();
                    if (dataAboutAmplicons.get(e.getAmplName()).distance(dataAboutAmplicons.get(amplName)) > minDistance) {// && counter < minNumOfModels) {
                        Amplicon correlatedAmpl = dataAboutAmplicons.get(e.getAmplName());
                        for (int i = 0; i < namesOfSamples.size(); i++) {
                            //factorsForSampleAndAmpl.set(i, factorsForSampleAndAmpl.get(i) + (correlatedAmpl.getCoverages().get(i)));
                            //toCalculateMean.get(i).add(correlatedAmpl.getCoverages().get(i));
                            toCalculateMean.get(i).add(predictedMedians.get(amplName).get(namesOfSamples.get(i)));

                        }
                        counter_of_values_to_compare_with++;
                    }
                }

                for (int i = 0; i < namesOfSamples.size(); i++) {
                    //double sum = 0.0;
                    //for (int j = 0; j < toCalculateMean.get(i).size(); ++j) sum += toCalculateMean.get(i).get(j);
                    factorsForSampleAndAmpl.set(i, Statistics.sm(toCalculateMean.get(i)));
                }
                for (int i = 0; i < namesOfSamples.size(); i++) {
                    samples.normalizeCov(namesOfSamples.get(i), amplName, factorsForSampleAndAmpl.get(i));
                }
            }
        }
        samplesToLDA = samples;

        fillDiagMatrixOfCovs();
    }


    public void fillDiagMatrixOfCovs() {
        ArrayList<String> amplNames = samplesToLDA.getFinallyAmplNames();

        for (String amplName : amplNames) {

            ArrayList<Double> toCalcVar = new ArrayList<Double>();
            for (String sampleName : namesOfSamples) {
                if (!dirtySamples.contains(sampleName)) {
                    toCalcVar.add(samplesToLDA.returnCov(sampleName, amplName));
                }
            }
            double robustVarianceOfCurrentAmpl = 0.0;
            double robustLocationEstimator = 0.0;
            robustVarianceOfCurrentAmpl = Math.pow(Statistics.snEstimator(toCalcVar), 2);

            robustLocationEstimator = Statistics.medW(toCalcVar);

            medians.put(amplName, robustLocationEstimator);
            diagMatrixOfCovariancesMap.put(amplName, 1.0 / robustVarianceOfCurrentAmpl);
        }
    }

    public double getDist(double val, String amplName) {
        return Math.min(Math.pow(val - medians.get(amplName), 2) * diagMatrixOfCovariancesMap.get(amplName), qChisq);
    }

    public double getNormPlot(double val, String amplName) {
        return (val - medians.get(amplName)) * Math.sqrt(diagMatrixOfCovariancesMap.get(amplName));
    }

    public double getDelPlot(double val, String amplName) {
        return (val - medians.get(amplName) + logDel) * Math.sqrt(0.5 * diagMatrixOfCovariancesMap.get(amplName));
    }

    public double getDupPlot(double val, String amplName) {
        return (val - medians.get(amplName) - logDup) * Math.sqrt(1.5 * diagMatrixOfCovariancesMap.get(amplName));    }

    public double getDistDel(double val, String amplName) {
        return Math.min(Math.pow(val - medians.get(amplName) + logDel, 2) * diagMatrixOfCovariancesMap.get(amplName) * 0.5, qChisq);
    }

    public double getDistDup(double val, String amplName) {
        return Math.min(Math.pow(val - medians.get(amplName) - logDup, 2) * diagMatrixOfCovariancesMap.get(amplName) * 1.5, qChisq);
    }


    public void outRes(
            String filename, String versionAndParams,
            ArrayList<String> samplesNames, ReaderBED infoAboutAmplicons,
            String datasetFile) throws IOException {
        String outFileName = filename + "_after_LDA.xls";
        Path file = Paths.get(outFileName);



        ArrayList<String> samplesToCheck = null;

        samplesToCheck = dirtySamples;

        try {
            // Create the empty file with default permissions, etc.
            Files.createFile(file);
        } catch (FileAlreadyExistsException x) {
            log.log(Level.FINE, "File named %s" +
                            " already exists. For accurate analysis, please, name your files differently.\n" +
                            "Now, please, repeat the analysis and choose the name for output file wisely. Overwriting of the file. ", file
            );
        } catch (IOException x) {
            // Some other sort of failure, such as permissions.
            log.log(Level.SEVERE, "createFile error: %s%n", x);
        }
        log.log(Level.FINE, outFileName);
        FileWriter fw = new FileWriter(outFileName);
        BufferedWriter bw = new BufferedWriter(fw);

        bw.write("Description\t" + versionAndParams + "\n");
        bw.write("The results after LDA analysis" + "\n");
        bw.write("Dataset was taken from file\t" + datasetFile + "\n");
        bw.write("BED file used from file\t" + infoAboutAmplicons.bedFileName + "\n");

        bw.newLine();
        String topString = "Chr\tStartpos\tStoppos\tAmplicon name\tDescription\t";
        for (String str : samplesToCheck) {
            topString += str;
            topString += "\t";
        }
        bw.write(topString + "\n");

        AmpliconSorter bvc = new AmpliconSorter(infoAboutAmplicons.getInfoAboutAmplicon());

        TreeMap<String, ArrayList<String>> sorted_map = new TreeMap<String, ArrayList<String>>(bvc);
        sorted_map.putAll(infoAboutAmplicons.getInfoAboutAmplicon());
        int delCounter = 0;
        int dupCounter = 0;

        HashMap<String, HashMap<String, Integer>> candidateExonsAndSamplesDel = new HashMap<String, HashMap<String, Integer>>();
        HashMap<String, HashMap<String, Integer>> candidateExonsAndSamplesDup = new HashMap<String, HashMap<String, Integer>>();

        HashMap<String, HashMap<String, Integer>> candidateExonsAndSamplesHomoDel = new HashMap<String, HashMap<String, Integer>>();
        HashMap<String, HashMap<String, Double>> candidateExonsAndSamplesDistToWt = new HashMap<String, HashMap<String, Double>>();
        HashMap<String, HashMap<String, Double>> candidateExonsAndSamplesDistToDel = new HashMap<String, HashMap<String, Double>>();
        HashMap<String, HashMap<String, Double>> candidateExonsAndSamplesDistToDup = new HashMap<String, HashMap<String, Double>>();
        for (String name : samplesToCheck) {
            distanceToDel.put(name, new HashMap<String, Double>());
            distanceToDup.put(name, new HashMap<String, Double>());
            distanceToNorm.put(name, new HashMap<String, Double>());


            for (Map.Entry<String, ArrayList<String>> e : sorted_map.entrySet()) {
                if (!nonEffectiveAmpls.contains(e.getKey()) && dataAboutAmplicons.containsKey(e.getKey())) {
                    double distToWt = getNormPlot(samplesToLDA.returnCov(name, e.getKey()), e.getKey());
                    double distToDel = getDelPlot(samplesToLDA.returnCov(name, e.getKey()), e.getKey());
                    double distToDup = getDupPlot(samplesToLDA.returnCov(name, e.getKey()), e.getKey());
                    distanceToNorm.get(name).put(e.getKey(), distToWt);
                    distanceToDel.get(name).put(e.getKey(), distToDel);
                    distanceToDup.get(name).put(e.getKey(), distToDup);
                }
            }


        }

        for (String name : samplesToCheck) {
            candidateExonsAndSamplesDel.put(name, new HashMap<String, Integer>());
            candidateExonsAndSamplesDup.put(name, new HashMap<String, Integer>());

            candidateExonsAndSamplesHomoDel.put(name, new HashMap<String, Integer>());
            candidateExonsAndSamplesDistToWt.put(name, new HashMap<String, Double>());
            candidateExonsAndSamplesDistToDup.put(name, new HashMap<String, Double>());
            candidateExonsAndSamplesDistToDel.put(name, new HashMap<String, Double>());

            for (Map.Entry<String, ArrayList<String>> exonAndListOfAmplicons : ampliconNamesFromExons.entrySet()) {
                candidateExonsAndSamplesDel.get(name).put(exonAndListOfAmplicons.getKey(), 0);
                candidateExonsAndSamplesDup.get(name).put(exonAndListOfAmplicons.getKey(), 0);
                candidateExonsAndSamplesHomoDel.get(name).put(exonAndListOfAmplicons.getKey(), 0);
                candidateExonsAndSamplesDistToWt.get(name).put(exonAndListOfAmplicons.getKey(), 0.0);
                candidateExonsAndSamplesDistToDel.get(name).put(exonAndListOfAmplicons.getKey(), 0.0);
                candidateExonsAndSamplesDistToDup.get(name).put(exonAndListOfAmplicons.getKey(), 0.0);
            }
        }

        for (Map.Entry<String, ArrayList<String>> e : sorted_map.entrySet()) {
            if (!nonEffectiveAmpls.contains(e.getKey()) && dataAboutAmplicons.containsKey(e.getKey())) {
                for (String name : samplesToCheck) {
                    String currentExon = dataAboutAmplicons.get(e.getKey()).getExon();

                    double distToWt = getDist(samplesToLDA.returnCov(name, e.getKey()), e.getKey());
                    double distToDel = getDistDel(samplesToLDA.returnCov(name, e.getKey()), e.getKey());
                    double distToDup = getDistDup(samplesToLDA.returnCov(name, e.getKey()), e.getKey());

                    candidateExonsAndSamplesDistToWt.get(name).put(currentExon, candidateExonsAndSamplesDistToWt.get(name).get(currentExon) + distToWt);
                    candidateExonsAndSamplesDistToDel.get(name).put(currentExon, candidateExonsAndSamplesDistToDel.get(name).get(currentExon) + distToDel);
                    candidateExonsAndSamplesDistToDup.get(name).put(currentExon, candidateExonsAndSamplesDistToDup.get(name).get(currentExon) + distToDup);

                    if ( // distance to Deletion model is not very big
                            (distToDel < qChisq) || distToDel < distToWt) {
                        candidateExonsAndSamplesDel.get(name).put(currentExon, candidateExonsAndSamplesDel.get(name).get(currentExon) + 1);
                    }
                    if ( // distance to Duplication model is not very big
                            (distToDup < qChisq) || distToDup < distToWt) {
                        candidateExonsAndSamplesDup.get(name).put(currentExon, candidateExonsAndSamplesDup.get(name).get(currentExon) - 1);
                    }
                    if (samplesToLDA.returnInitCov(name, e.getKey()) < 0.01) { // create a mark for determining Homozygous Deletion
                            candidateExonsAndSamplesHomoDel.get(name).put(currentExon, candidateExonsAndSamplesHomoDel.get(name).get(currentExon) + 1);
                        }
                    }
            }
        }

        HashMap<String, ArrayList<String>> deleteriousExons = new HashMap<String, ArrayList<String>>();
        HashMap<String, ArrayList<String>> duplicatedExons = new HashMap<String, ArrayList<String>>();
        HashMap<String, ArrayList<String>> homoDelExons = new HashMap<String, ArrayList<String>>();

        //int totalNumberOfTests = 0;
        for (String name : samplesToCheck) {
            deleteriousExons.put(name, new ArrayList<String>());
            duplicatedExons.put(name, new ArrayList<String>());
            homoDelExons.put(name, new ArrayList<String>());
            for (Map.Entry<String, ArrayList<String>> exonAndListOfAmplicons : ampliconNamesFromExons.entrySet()) {
                if (candidateExonsAndSamplesDel.get(name).get(exonAndListOfAmplicons.getKey()) == exonAndListOfAmplicons.getValue().size()) {
                    deleteriousExons.get(name).add(exonAndListOfAmplicons.getKey());
                    //totalNumberOfTests++;
                } else if ((candidateExonsAndSamplesDup.get(name).get(exonAndListOfAmplicons.getKey()) == -exonAndListOfAmplicons.getValue().size())) {
                    duplicatedExons.get(name).add(exonAndListOfAmplicons.getKey());
                    //totalNumberOfTests++;
                }
                if ((candidateExonsAndSamplesHomoDel.get(name).get(exonAndListOfAmplicons.getKey()) == exonAndListOfAmplicons.getValue().size())) {
                    homoDelExons.get(name).add(exonAndListOfAmplicons.getKey());
                }
            }
        }
        //numberOfTests = totalNumberOfTests;
        log.log(Level.INFO, "Overall number of suspicious exons: " + numberOfTests);

        for (Map.Entry<String, ArrayList<String>> e : sorted_map.entrySet()) {
            if (finalAmpls.contains(e.getKey())) {
                String currentExon = dataAboutAmplicons.get(e.getKey()).getExon();

                for (String name : samplesToCheck) {
                    double distExonToWt = candidateExonsAndSamplesDistToWt.get(name).get(currentExon);
                    org.apache.commons.math3.distribution.ChiSquaredDistribution chiSq = new org.apache.commons.math3.distribution.ChiSquaredDistribution(ampliconNamesFromExons.get(currentExon).size());
                    pValues.add(1 - chiSq.cumulativeProbability(distExonToWt));
                }
            }
        }
        Collections.sort(pValues);
        int k = 1;
        int m = numberOfTests;

        for (int i = 0; i < pValues.size(); i++) {
           if (pValues.get(i) < (k * alpha / m)) {
               k++;
           } else {
               break;
           }
        }
        alpha *= (1.0 * k)/m;

        for (Map.Entry<String, ArrayList<String>> e : sorted_map.entrySet()) {
            if (finalAmpls.contains(e.getKey())) {
                String currentExon = dataAboutAmplicons.get(e.getKey()).getExon();
                String tmpString = "";
                tmpString += e.getValue().get(0) + "\t";
                tmpString += e.getValue().get(1) + "\t";
                tmpString += e.getValue().get(2) + "\t";
                tmpString += e.getKey() + "\t";
                tmpString += e.getValue().get(3) + "\t";

                for (String name : samplesToCheck) {
                    double distExonToDel = candidateExonsAndSamplesDistToDel.get(name).get(currentExon);
                    double distExonToWt = candidateExonsAndSamplesDistToWt.get(name).get(currentExon);
                    double distExonToDup = candidateExonsAndSamplesDistToDup.get(name).get(currentExon);
                    org.apache.commons.math3.distribution.ChiSquaredDistribution chiSq = new org.apache.commons.math3.distribution.ChiSquaredDistribution(ampliconNamesFromExons.get(currentExon).size());


                    if (homoDelExons.get(name).contains(currentExon)) {
                        tmpString += "HOMO DEL\t";
                        if (!newDirtySamples.contains(name)) {
                            newDirtySamples.add(name);
                        }
                        delCounter++;
                    } else if (deleteriousExons.get(name).contains(currentExon) &&
                            strongEvidenceDeletion(currentExon, distExonToDel, distExonToWt)
                            && notSimilarToNormalButSimilarToDeletion(currentExon, distExonToDel, distExonToWt)
                            )
                    {
                        tmpString += "DEL                " + "bf: " + chiSq.density(distExonToDel) / chiSq.density(distExonToWt) + "\t";
                        if (!newDirtySamples.contains(name)) {
                            newDirtySamples.add(name);
                        }
                        delCounter++;
                    } else if (duplicatedExons.get(name).contains(currentExon) &&
                            strongEvidenceDuplication(currentExon, distExonToDup, distExonToWt)
                            && notSimilarToNormalButSimilarToDuplication(currentExon, distExonToDup, distExonToWt)
                            )
                    {
                        tmpString += "^                     " + "bf: " + chiSq.density(distExonToDup) / chiSq.density(distExonToWt) + "\t";
                        if (!newDirtySamples.contains(name)) {
                            newDirtySamples.add(name);
                        }
                        dupCounter++;
                    }
                    else {
                        tmpString += " " + "\t";
                    }
                }

                bw.write(tmpString + "\n");
                bw.flush();
            }
        }
        bw.write("\n");
        bw.write("Total deletions: " + delCounter + "\n");
        bw.write("Total duplications: " + dupCounter + "\n");

        bw.close();
        if (dirtySamples.size() > newDirtySamples.size()) {
            numberOfDirtySamplesChanged = true;
        } else if (dirtySamples.size() < newDirtySamples.size()) {
            log.log(Level.SEVERE, "The number of dirty samples if less than number of new dirty samples. It is normal only if it is the final step.");
        }
    }

    public ArrayList<String> getNewDirtySamples() {
        return newDirtySamples;
    }

    private boolean notSimilarToNormalButSimilarToDeletion(String currentExon, double distExonToDel, double distExonToWt) {
        /*
        determines if exon is more closer to Deletion CNV state, but clearly not normal
         */
        org.apache.commons.math3.distribution.ChiSquaredDistribution chiSq = new org.apache.commons.math3.distribution.ChiSquaredDistribution(ampliconNamesFromExons.get(currentExon).size());

        return (
                        distExonToWt > chiSq.inverseCumulativeProbability(1 - alpha)
                )
        ;
    }

    private boolean notSimilarToNormalButSimilarToDuplication(String currentExon, double distExonToDup, double distExonToWt) {
        /*
        determines if exon is more closer to Duplication CNV state, but clearly not normal
         */
        org.apache.commons.math3.distribution.ChiSquaredDistribution chiSq = new org.apache.commons.math3.distribution.ChiSquaredDistribution(ampliconNamesFromExons.get(currentExon).size());
        return  (
                       distExonToWt > chiSq.inverseCumulativeProbability(1 - alpha)
        );
    }

    private boolean strongEvidenceDeletion(String currentExon, double distExonToDel, double distExonToWt) {
        /*
        determines if exon is more closer to Deletion CNV state, but clearly not normal
         */
        org.apache.commons.math3.distribution.ChiSquaredDistribution chiSq = new org.apache.commons.math3.distribution.ChiSquaredDistribution(ampliconNamesFromExons.get(currentExon).size());

        return (
                chiSq.density(distExonToDel) / chiSq.density(distExonToWt) > 100

        )
                ;
    }


    private boolean strongEvidenceDuplication(String currentExon, double distExonToDup, double distExonToWt) {
        /*
        determines if exon is more closer to Deletion CNV state, but clearly not normal
         */
        org.apache.commons.math3.distribution.ChiSquaredDistribution chiSq = new org.apache.commons.math3.distribution.ChiSquaredDistribution(ampliconNamesFromExons.get(currentExon).size());

        return (
                chiSq.density(distExonToDup) / chiSq.density(distExonToWt) > 100
        )
                ;
    }

    public void printDistances(String distanceFileName) {
        Path file = Paths.get(distanceFileName);
        try {
            FileWriter fw = new FileWriter(distanceFileName);
            BufferedWriter bw = new BufferedWriter(fw);
            String samples = "\t";
            for (String sample : newDirtySamples) {
                samples += sample + "\t";
                samples += sample + "\t";
                samples += sample + "\t";
            }
            samples += "\n";
            bw.write(samples);
            for (String ampl : finalAmpls) {
                String stringDistances = ampl + "\t";
                for (String sample : newDirtySamples) {
                    stringDistances += distanceToNorm.get(sample).get(ampl) + "\t";
                    stringDistances += distanceToDel.get(sample).get(ampl) + "\t";
                    stringDistances += distanceToDup.get(sample).get(ampl) + "\t";
                }
                bw.write(stringDistances + "\n");
            }
            bw.flush();
            bw.close();
        } catch (Exception e){
            log.log(Level.SEVERE, "Problem with creation of file with distances.");
            e.printStackTrace();
        }
    }
}