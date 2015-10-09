package deletionsAnalysis;

import org.apache.commons.cli.*;

import java.io.IOException;
import java.io.Serializable;
import java.util.*;
import java.io.File;
import java.util.logging.*;
/**
 * Created by german on 11.07.14.
 */

class PairComparator implements Comparator<Pair>, Serializable {
    @Override
    public int compare(Pair x, Pair y)
            /*
             * Simple Double comparator. Analogue of Double.compare().
             */ {
        if (x.getCorrelation() > y.getCorrelation()) {
            return -1;
        }
        if (x.getCorrelation() < y.getCorrelation()) {
            return 1;
        }
        return 0;
    }
}

class AmpliconSorter implements Comparator<String> {
    Map<String, ArrayList<String>> base;

    public AmpliconSorter(Map<String, ArrayList<String>> base) {
        this.base = base;
    }
    /*
    sort amplicons by position on chromosome.
    we compare chromosome names at first and then we compare positions on chromosome.
     */
    @Override
    public int compare(String x, String y) {
        /*
        simple comparator method: compare names of chromosomes and then positions on it
         */
        if (base.get(x).get(0).equals(base.get(y).get(0))) {
            if (Integer.parseInt(base.get(x).get(1)) > Integer.parseInt(base.get(y).get(1))) {
                return 1;
            }
            if (Integer.parseInt(base.get(x).get(1)) < Integer.parseInt(base.get(y).get(1))) {
                return -1;
            }
        } else {
            return (base.get(x).get(0).compareToIgnoreCase(base.get(y).get(0)));
        }
        return 0;
    }
}



class Solver {
    private static final Logger log = Logger.getLogger( Solver.class.getName() );


    private void createDirectoryIfNeeded(String directoryName) {
        File theDir = new File(directoryName);

        // if the directory does not exist, create it
        if (!theDir.exists())
        {
            log.log(Level.INFO, "creating directory: " + directoryName);
            theDir.mkdir();
        }
    }

    HashMap<String, PriorityQueue<Pair>> getAmplCorrelationPriority(HashMap<String, Amplicon> dataAboutAmplicons,
                                                                    ArrayList<String> samplesNames, Comparator<Pair> comparator,
                                                                    ArrayList<String> dirtySamples, Boolean robustCor, Homozygothe homo
    ) {
        HashMap<String, PriorityQueue<Pair>> amplCorrelationPriority = new HashMap<String, PriorityQueue<Pair>>();
        HashMap<String, Amplicon> dataAboutAmpliconsWithoutDirty = new HashMap<String, Amplicon>();
        log.log(Level.FINE, "In training set for LDA: " + (samplesNames.size() - dirtySamples.size()) + " samples.");

        for ( Map.Entry<String, Amplicon> idOfAmpl : dataAboutAmplicons.entrySet()) {
            try {
                Amplicon newAmplicon = (Amplicon) ObjectCloner.deepCopy(idOfAmpl.getValue());
                ArrayList<Double> newCoverages = new ArrayList<Double>();
                for (int i = 0; i < samplesNames.size(); i++) {
                    String str = samplesNames.get(i);
                    if (robustCor) {
                        newCoverages.add(newAmplicon.getCoverages().get(i));
                    }
                    else if (!dirtySamples.contains(str)) {
                        newCoverages.add(newAmplicon.getCoverages().get(i));
                    }
                }
                newAmplicon.modifyCoverages(newCoverages);
                dataAboutAmpliconsWithoutDirty.put(idOfAmpl.getKey(), newAmplicon);
            } catch (Exception e) {
                log.log(Level.SEVERE, e.getMessage());
            }
        }

        for (Map.Entry<String, Amplicon> entry1 : dataAboutAmpliconsWithoutDirty.entrySet()) {
            PriorityQueue<Pair> priorityByCorrelation =
                    new PriorityQueue<Pair>(dataAboutAmpliconsWithoutDirty.size(), comparator);
            for (Map.Entry<String, Amplicon> entry2 : dataAboutAmpliconsWithoutDirty.entrySet()) {
                Double cor = 0.0;
                if (robustCor) {
                    cor = Statistics.robustCor(entry1.getValue().getCoverages(), entry2.getValue().getCoverages());
                } else cor = Statistics.cor(entry1.getValue().getCoverages(), entry2.getValue().getCoverages());
                Pair pairOfAmpls =
                        new Pair(cor, entry2.getKey());
                priorityByCorrelation.add(pairOfAmpls);
            }
            amplCorrelationPriority.put(entry1.getKey(), priorityByCorrelation);
        }
        return amplCorrelationPriority;
    }

    public Solver(String[] args) {
        /*
         * Initialize and solve all the script.
         *
         * @param arguments from command line
         * @return starts functions and creates a xls file with the results of the analysis
         */

        OptionsParse setOfOptions = new OptionsParse();

        Options option = setOfOptions.OptionParse();


        CommandLineParser parser = new GnuParser();

        HelpFormatter formatter = new HelpFormatter();
        CommandLine cmd = null;

        try {
            cmd = parser.parse(option, args);
        } catch (ParseException e) {
            log.log(Level.SEVERE, "Parsing of command line args failed.  Reason: " + e.getMessage());
            log.log(Level.SEVERE, e.getMessage());
        }

        if (cmd.hasOption("h")) {
            formatter.printHelp("usage", option);
            System.exit(0);
        }
        // quantiles for one degree of freedom
        double[] quantilesChisq = {23.92813, 20.83729, 19.51142, 16.44811, 15.13671, 12.11567, 10.82757, 8.807468, 7.879439, 6.634897};
        // quantiles for two degrees of freedom
        double[] quantilesChisq2df = {29.01732, 27.63102, 24.41215, 23.02585, 19.80698, 18.42068, 15.2018, 13.81551, 10.59663, 9.21034};

        double qChisq = quantilesChisq[3];
        double qChisqForDouble = quantilesChisq2df[3];

        int quantChisq = 0;
        try {
            quantChisq = Integer.parseInt(cmd.getOptionValue("qc"));
            qChisq = quantilesChisq[quantChisq];
            qChisqForDouble = quantilesChisq2df[quantChisq];
        } catch (Exception e) {

            log.log(Level.WARNING, "WARNING: minimum chisq for LDA quantile should be integer. Default value used.");
            qChisq = quantilesChisq[3];
            qChisqForDouble = quantilesChisq2df[3];
        }

        String filename = cmd.getOptionValue("f");
        double minCorThreshold = 0.75;
        try {
            minCorThreshold = Double.parseDouble(cmd.getOptionValue("mc"));
        } catch (Exception e) {
            log.log(Level.WARNING, "WARNING: minimum correlation threshold should be double. Default value used.");
        }

        int maxNumOfModels = 5;
        try {
            maxNumOfModels = Integer.parseInt(cmd.getOptionValue("mnm"));
        } catch (Exception e) {
            log.log(Level.WARNING, "WARNING: max number of models should be integer. Default value used.");
        }

        int numOfModelForNonEfficiency = 3;
        try {
            numOfModelForNonEfficiency = Integer.parseInt(cmd.getOptionValue("nne"));
        } catch (Exception e) {
            log.log(Level.WARNING, "WARNING: number of models for non efficiency should be integer. Default value used.");
        }

        int minNumOfModelsForOutlierDetection = 3;
        try {
            minNumOfModelsForOutlierDetection = Integer.parseInt(cmd.getOptionValue("nod"));
        } catch (Exception e) {
            log.log(Level.WARNING, "WARNING: number of models for outlier detection should be integer. Default value used.");
        }

        double lowerBoundOutlier = -2.59;
        try {
            lowerBoundOutlier = Double.parseDouble(cmd.getOptionValue("lb"));
        } catch (Exception e) {
            log.log(Level.WARNING, "WARNING: lower bound should be floating point number. Default value used.");
        }


        double upperBoundOutlier = 2.59;
        try {
            upperBoundOutlier = Double.parseDouble(cmd.getOptionValue("ub"));
        } catch (Exception e) {
            log.log(Level.WARNING, "WARNING: upper bound should be floating point number. Default value used.");
        }

        int minDistanceBetween = 1000000;
        try {
            minDistanceBetween = Integer.parseInt(cmd.getOptionValue("dist"));
        } catch (Exception e) {
            log.log(Level.WARNING, "WARNING: distance between paired amplicons should be integer. Default value used.");
        }

        int lowCovBound = 25000;
        try {
            lowCovBound = Integer.parseInt(cmd.getOptionValue("lcb"));
        } catch (Exception e) {
            log.log(Level.WARNING, "WARNING: bound for low covered samples should be integer. Default value used.");
        }

        double learningSampleSize = 0.8;
        try {
            learningSampleSize = Double.parseDouble(cmd.getOptionValue("lss"));
        } catch (Exception e) {
            log.log(Level.WARNING, "WARNING: bound for learning sample size was specified in a wrong format. Default value used.");
        }

        int thresholdForInefficiency = 100;
        try {
            thresholdForInefficiency = Integer.parseInt(cmd.getOptionValue("lca"));
        } catch (Exception e) {
            log.log(Level.WARNING, "WARNING: threshold for non efficiency should be integer. Default value used.");
        }

        VersionAndParams topStrings = new VersionAndParams(minCorThreshold, maxNumOfModels, numOfModelForNonEfficiency,
                minNumOfModelsForOutlierDetection, lowerBoundOutlier, upperBoundOutlier,
                minDistanceBetween, thresholdForInefficiency);
        String versionAndParams = topStrings.getVersAndParams();

        ReaderBED readerBed = new ReaderBED(cmd.getOptionValue("b"));
        ReaderData readerData = new ReaderData(cmd.getOptionValue("d"));
        ArrayList<String> namesOfMissedAmpl = new ArrayList<String>();
        for (Map.Entry<String, ArrayList<String>> entry : readerBed.getInfoAboutAmplicon().entrySet()) {
            if (!readerData.getCoverages().containsKey(entry.getKey())) {
                namesOfMissedAmpl.add(entry.getKey());
            }
        }
        for (String name : namesOfMissedAmpl) {
            readerBed.removeAmplFromBed(name);
        }
        HashMap<String, Amplicon> dataAboutAmplicons = new HashMap<String, Amplicon>();
        final ArrayList<String> samplesNames = readerData.getSamplesNames();


        ArrayList<String> listOfExcluded = new ArrayList<String>();
        int total_size = 0;

        for (Map.Entry<String, ArrayList<String>> entry : readerBed.getInfoAboutAmplicon().entrySet()) {
            boolean b = !(entry.getKey().isEmpty());
            if (b && readerData.getCoverages().containsKey(entry.getKey())) {
                Amplicon tmpAmplicon = new Amplicon(entry.getKey(), readerBed.getInfoAboutAmplicon(), readerData.getCoverages(), thresholdForInefficiency);

                if (tmpAmplicon.isExcluded()) {
                    listOfExcluded.add(entry.getKey());
                } else {
                    total_size++;
                    dataAboutAmplicons.put(entry.getKey(), tmpAmplicon);
                }
            }
        }


        log.log(Level.FINE, "PREPROCESSING...");
        log.log(Level.FINE,"\n\tNOTE: It is important to analyse the results of preprocessing and to decide what you want to do with low\n" +
                "\tcovered samples and homozygous amplicons.\n" +
                "\tRemember: sometimes homozygous deletions is just an artifacts of PCR or mapping quality\n" +
                "\tof reads in this amplicon can be equal to 0. Check this manually and remove \n" +
                "\tsuspicious samples from the further analysis\n");
        log.log(Level.FINE, "Step 1. Searching for low covered amplicons.");
        // Report about low covered amplicons, this will explain 'EX' labels
        log.log(Level.FINE, "------------------------------------------------------------");
        log.log(Level.FINE, "These amplicons won't work (their mathematical expectation < " + thresholdForInefficiency + "): ");
        for (String elem : listOfExcluded) {
            log.log(Level.FINE, "\t" + elem);
        }
        total_size -= listOfExcluded.size();
        log.log(Level.FINE, "--------------- THESE AMPLICONS ARE EXCLUDED ---------------");
        log.log(Level.FINE, "");

        log.log(Level.FINE, "Step 2. Searching for low covered samples.Value used: " + lowCovBound);
        log.log(Level.FINE, "------------------------------------------------------------");
        String totalCoverages = "";
        DetectLowCovSamples delLowCov = new DetectLowCovSamples(samplesNames, dataAboutAmplicons, lowCovBound);
        log.log(Level.FINE, "-------------- THESE SAMPLES ARE TOO LOW COVERED ---------------");

        log.log(Level.FINE, "Step 3. Regression analysis.");

        HashMap<Integer, Double> totals = delLowCov.totals;
        totalCoverages = delLowCov.returnTotalCoverage(samplesNames);


        Comparator<Pair> comparator = new PairComparator();

        HashMap<String, PriorityQueue<Pair>> amplCorrelationPriority = new HashMap<String, PriorityQueue<Pair>>();
        HashMap<String, PriorityQueue<Pair>> amplCorrelationPriorityForLDA = new HashMap<String, PriorityQueue<Pair>>();
        HashMap<String, HashMap<String, ArrayList<Double>>> predictedValues = new
                HashMap<String, HashMap<String, ArrayList<Double>>>();



        HashMap<String, HashMap<String, Integer>> sampleToAmplicons = new HashMap<String, HashMap<String, Integer>>();
        ArrayList<String> nonEffecitveAmpls = new ArrayList<String>();

        for (String elem : samplesNames) {
            HashMap<String, Integer> tmpHash = new HashMap<String, Integer>();
            HashMap<String, ArrayList<Double>> tmpPredicted = new HashMap<String, ArrayList<Double>>();
            for (String ampl : dataAboutAmplicons.keySet()) {
                tmpHash.put(ampl, 0);
                tmpPredicted.put(ampl, new ArrayList<Double>());
            }
            sampleToAmplicons.put(elem, tmpHash);
            predictedValues.put(elem, tmpPredicted);
        }
        HashMap<String, Double> correlationsBetweenAmplicons = new HashMap<String, Double>();
        for (Map.Entry<String, Amplicon> entry1 : dataAboutAmplicons.entrySet()) {
            for (Map.Entry<String, Amplicon> entry2 : dataAboutAmplicons.entrySet()) {
                if (!correlationsBetweenAmplicons.containsKey(entry1.getKey() + entry2.getKey())) {
                    Double tmpRobustCorrelation = Statistics.robustCor(entry1.getValue().getCoverages(), entry2.getValue().getCoverages());
                    correlationsBetweenAmplicons.put(entry1.getKey() + entry2.getKey(), tmpRobustCorrelation);
                    correlationsBetweenAmplicons.put(entry2.getKey() + entry1.getKey(), tmpRobustCorrelation);
                }
            }
        }

        for (Map.Entry<String, Amplicon> entry1 : dataAboutAmplicons.entrySet()) {
            PriorityQueue<Pair> priorityByCorrelation =
                    new PriorityQueue<Pair>(dataAboutAmplicons.size(), comparator);
            for (Map.Entry<String, Amplicon> entry2 : dataAboutAmplicons.entrySet()) {
                Double robCor = correlationsBetweenAmplicons.get(entry1.getKey() + entry2.getKey());
                Pair pairOfAmpls =
                        new Pair(robCor, entry2.getKey());
                priorityByCorrelation.add(pairOfAmpls);

            }
            amplCorrelationPriority.put(entry1.getKey(), priorityByCorrelation);

            double tmpCor = 1.0;
            int counter = 0;
            Pair e = priorityByCorrelation.poll();
            ArrayList<String> potentiallyDuplicated = new ArrayList<String>();

            HashSet<String> setOfDeleterious = new HashSet<String>(samplesNames);
            HashSet<String> setOfDuplicated = new HashSet<String>(samplesNames);

            while (tmpCor > minCorThreshold && !priorityByCorrelation.isEmpty() && counter < maxNumOfModels) {

                ArrayList<String> potentiallyDeleterious = new ArrayList<String>();
                e = priorityByCorrelation.poll();

                if (entry1.getValue().distance(dataAboutAmplicons.get(e.getAmplName())) > minDistanceBetween) {
                    {
                        tmpCor = e.getCorrelation();
                        LinearModel lm = new LinearModel(dataAboutAmplicons.get(e.getAmplName()).getCoverages(), entry1.getValue().getCoverages());
                        LinearModel lmDel = new LinearModel(dataAboutAmplicons.get(e.getAmplName()).getCoverages(), entry1.getValue().getShiftedCoverages(-1.0));
                        lmDel.setStudentizedResidualsShift(entry1.getValue().getCoverages());
                        LinearModel lmDup = new LinearModel(dataAboutAmplicons.get(e.getAmplName()).getCoverages(), entry1.getValue().getShiftedCoverages(0.5849625));
                        lmDup.setStudentizedResidualsShift(entry1.getValue().getCoverages());

                        ArrayList<Double> jackknifedResiduals = lm.getInternallyStudentizedResiduals();
                        ArrayList<Double> jackknifedResidualsDel = lmDel.getInternallyStudentizedResidualsShift();
                        ArrayList<Double> jackknifedResidualsDup = lmDup.getInternallyStudentizedResidualsShift();

                        for (int i = 0; i < jackknifedResiduals.size(); i++) {
                            if (jackknifedResiduals.get(i) < lowerBoundOutlier &&
                                    (jackknifedResidualsDel.get(i) < -1.0 * lowerBoundOutlier || Math.abs(jackknifedResidualsDel.get(i)) < Math.abs(jackknifedResiduals.get(i)))) {

                                potentiallyDeleterious.add(samplesNames.get(i));
                                Integer x = sampleToAmplicons.get(samplesNames.get(i)).get(entry1.getKey()) + 1;
                                sampleToAmplicons.get(samplesNames.get(i)).put(entry1.getKey(), x);

                            } else if (jackknifedResiduals.get(i) > upperBoundOutlier &&
                                    (jackknifedResidualsDup.get(i) > -1.0 * upperBoundOutlier || Math.abs(jackknifedResidualsDup.get(i)) < Math.abs(jackknifedResiduals.get(i)))) {
                                potentiallyDuplicated.add(samplesNames.get(i));
                                Integer x = sampleToAmplicons.get(samplesNames.get(i)).get(entry1.getKey()) - 1;
                                sampleToAmplicons.get(samplesNames.get(i)).put(entry1.getKey(), x);

                            }
                        }
                        HashSet<String> tmpDeleterious = new HashSet<String>(potentiallyDeleterious);
                        HashSet<String> tmpDuplicated = new HashSet<String>(potentiallyDuplicated);

                        setOfDeleterious.retainAll(tmpDeleterious);
                        setOfDuplicated.retainAll(tmpDuplicated);

                        counter++;
                    }
                }
            }

            if (counter < numOfModelForNonEfficiency) {
                log.log(Level.FINE, "Amplicon " + entry1.getValue().returnID() + " does not have enough located far enough and correlated amplicons for analysis");
                total_size--;
                nonEffecitveAmpls.add(entry1.getValue().returnID());
            } else {
                log.log(Level.FINE, "Analysis of " + entry1.getValue().returnID() + " was completed succesfully.");
            }
        }
        for (String sample : samplesNames) {
            for (Map.Entry<String, Amplicon> entry :dataAboutAmplicons.entrySet()) {
                entry.getValue().addElemToListOfPredictedCoverages(Statistics.med(predictedValues.get(sample).get(entry.getKey())));
            }
        }

        log.log(Level.FINE, "-------ANALYSIS WITH THE HELP OF SCORING FUNCTION-------");
        for (Map.Entry<String, HashMap<String, Integer>> entry1 : sampleToAmplicons.entrySet()) {
            ArrayList<String> potentiallyDeleteriousWithScoring = new ArrayList<String>();
            for (String key : entry1.getValue().keySet()) {
                if (entry1.getValue().get(key) >= minNumOfModelsForOutlierDetection) {
                    potentiallyDeleteriousWithScoring.add(key);
                }
            }
        }


        // Report about samples that contain homozygous amplicons
        log.log(Level.FINE, "Step 4. Seach for homozygous amplicons.");
        log.log(Level.FINE, "------------------------------------------------------------");
        Homozygothe homo = new Homozygothe(dataAboutAmplicons, nonEffecitveAmpls);
        homo.reportHomozygotheSample(samplesNames);

        log.log(Level.FINE, "-------------- THESE SAMPLES ARE HOMOZYGOTES ---------------");
        log.log(Level.FINE, "");
        Samples sam = new Samples(samplesNames, dataAboutAmplicons, nonEffecitveAmpls);
        ArrayList<String> dirtySamples;
        try {
            OutputResult outRes = null;
            outRes = new OutputResult(filename, listOfExcluded, nonEffecitveAmpls,
                    sampleToAmplicons, versionAndParams,
                    samplesNames, readerBed, minNumOfModelsForOutlierDetection,
                    totalCoverages, cmd.getOptionValue("d"));


            final Samples samClone = (Samples) (ObjectCloner.deepCopy(sam));

            dirtySamples = outRes.getDirtySamples();
            final ArrayList<String> dirtySamplesAfterFirstStep = (ArrayList<String>) (ObjectCloner.deepCopy(dirtySamples));

            int maximum_cycles_of_algorithm = 10;
            LDA lda = null;

            createDirectoryIfNeeded("./tmpResultsCNV");
            File file = new File("./tmpResultsCNV");
            String[] myFiles;
            myFiles = file.list();
            for (int k = 0; k < myFiles.length; k++) {
                File myFile = new File(file, myFiles[k]);
                myFile.delete();
            }

            if (samplesNames.size() <= 20) {
                log.log(Level.WARNING, "The second stage of the alogrithm can not be succesfull due the small number of samples.");
                qChisq = quantilesChisq[quantChisq];
                for (int i = 0; i < maximum_cycles_of_algorithm; i++) {
                    log.log(Level.FINE, "LDA Step " + i);
                    amplCorrelationPriorityForLDA = getAmplCorrelationPriority(dataAboutAmplicons,
                            samplesNames, comparator, dirtySamples, true, homo);
                    sam = (Samples) (ObjectCloner.deepCopy(samClone));
                    lda = new LDA(sam, dataAboutAmplicons, samplesNames,
                            dirtySamples, nonEffecitveAmpls, amplCorrelationPriorityForLDA, minCorThreshold,
                            minDistanceBetween, maxNumOfModels, qChisq, qChisqForDouble);
                    lda.outRes("./tmpResultsCNV/" + filename + "_" + quantChisq + "_step_" + i, versionAndParams, samplesNames, readerBed, cmd.getOptionValue("d"));
                    dirtySamples = lda.getNewDirtySamples();
                    if (lda.numberOfDirtySamplesChanged == false || i == maximum_cycles_of_algorithm) {
                        lda.outRes(filename, versionAndParams, samplesNames, readerBed, cmd.getOptionValue("d"));
                        break;
                    }
                }
            } else {
                LDA theBestLDA = null;

                int q = quantChisq;
                qChisq = quantilesChisq[q];
                dirtySamples = (ArrayList<String>) (ObjectCloner.deepCopy(dirtySamplesAfterFirstStep));
                for (int i = 0; i < maximum_cycles_of_algorithm; i++) {
                    log.log(Level.FINE, "LDA Step " + i + " with quantile " + qChisq);
                    if (samplesNames.size() - dirtySamples.size() > Math.max(samplesNames.size() / 2 + 1, 20)) {
                        amplCorrelationPriorityForLDA = getAmplCorrelationPriority(dataAboutAmplicons,
                                samplesNames, comparator, dirtySamples, false, homo);
                    } else {
                        amplCorrelationPriorityForLDA = getAmplCorrelationPriority(dataAboutAmplicons,
                                samplesNames, comparator, dirtySamples, true, homo);
                    }

                    sam = (Samples) (ObjectCloner.deepCopy(samClone));
                    lda = new LDA(sam, dataAboutAmplicons, samplesNames,
                            dirtySamples, nonEffecitveAmpls, amplCorrelationPriorityForLDA, minCorThreshold,
                            minDistanceBetween, maxNumOfModels, qChisq, qChisqForDouble);

                    dirtySamples = lda.getNewDirtySamples();

                    lda.outRes("./tmpResultsCNV/" + filename + "_" + dirtySamples.size() + "_" + q + "_step_" + i, versionAndParams, samplesNames, readerBed, cmd.getOptionValue("d"));

                    if (lda.numberOfDirtySamplesChanged == false) {
                        theBestLDA = (LDA)ObjectCloner.deepCopy(lda);
                        log.log(Level.INFO, "The sufficient results were obtained using qChisq: " + quantChisq + " with sensitivity option " + q);
                        break;
                    }
                }
                theBestLDA.outRes(filename , versionAndParams, samplesNames, readerBed, cmd.getOptionValue("d"));
            }
        } catch (IOException e) {
            log.log(Level.SEVERE, "Input / output exception! See Stack Trace. \n" + e.getMessage());
        } catch (Exception e) {
            log.log(Level.SEVERE, "Unknown error. See Stack Trace. \n" + e.getMessage());
        }

    }

}

public class Main {

    public static void main(String[] args) {
        /*
         * The main method of the main class
         *
         * @param arguments of command line
         * @return nothing
         */

        try {
            try {
                LogManager.getLogManager().readConfiguration(
                        Solver.class.getResourceAsStream("./logging.properties"));
            } catch (IOException e) {
                System.err.println("Could not setup logger configuration: " + e.toString());
            }
            Solver solve = new Solver(args);
        } catch (Exception e) {
            System.err.println("\nSmth wrong with the tool! All we can do is to output a stack trace.");
            e.printStackTrace();
            System.err.println("Please, show this message to your system administrator or, if you do not have one, to your mom.");

        }

    }
}

