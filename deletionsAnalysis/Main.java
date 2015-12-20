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

    Random rand = new Random(1990);
    private void createDirectoryIfNeeded(String directoryName) {
        File theDir = new File(directoryName);

        // if the directory does not exist, create it
        if (!theDir.exists())
        {
            log.log(Level.INFO, "creating directory: " + directoryName);
            theDir.mkdir();
        }
    }

    HashMap<String, PriorityQueue<Pair>> getAmplCorrelationPriorityRob(HashMap<String, Amplicon> dataAboutAmplicons,
                                                                       ArrayList<String> samplesNames, Comparator<Pair> comparator
    ) {
        HashMap<String, PriorityQueue<Pair>> amplCorrelationPriority = new HashMap<String, PriorityQueue<Pair>>();
        HashMap<String, Amplicon> dataAboutAmpliconsWithoutDirty = new HashMap<String, Amplicon>();

        for ( Map.Entry<String, Amplicon> idOfAmpl : dataAboutAmplicons.entrySet()) {
            try {
                Amplicon newAmplicon = (Amplicon) ObjectCloner.deepCopy(idOfAmpl.getValue());
                ArrayList<Double> newCoverages = new ArrayList<Double>();
                for (int i = 0; i < samplesNames.size(); i++) {
                    String str = samplesNames.get(i);
                    newCoverages.add(newAmplicon.getCoverages().get(i));
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
                cor = Statistics.robustCor(entry1.getValue().getCoverages(), entry2.getValue().getCoverages());
                Pair pairOfAmpls =
                        new Pair(cor, entry2.getKey());
                priorityByCorrelation.add(pairOfAmpls);
            }
            amplCorrelationPriority.put(entry1.getKey(), priorityByCorrelation);
        }
        return amplCorrelationPriority;
    }



    HashMap<String, PriorityQueue<Pair>> getAmplCorrelationPrioritySimple(HashMap<String, Amplicon> dataAboutAmplicons,
                                                                          ArrayList<String> samplesNames, Comparator<Pair> comparator,
                                                                          ArrayList<String> dirtySamples
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
                    if (!dirtySamples.contains(str)) {
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
                cor = Statistics.cor(entry1.getValue().getCoverages(), entry2.getValue().getCoverages());
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

        int numOfModelForNonEfficiency = 4;
        try {
            numOfModelForNonEfficiency = Integer.parseInt(cmd.getOptionValue("nne"));
        } catch (Exception e) {
            log.log(Level.WARNING, "WARNING: number of models for non efficiency should be integer. Default value used.");
        }

        int minNumOfModelsForOutlierDetection = 4;
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

        int thresholdForInefficiency = 100;
        try {
            thresholdForInefficiency = Integer.parseInt(cmd.getOptionValue("lca"));
        } catch (Exception e) {
            log.log(Level.WARNING, "WARNING: threshold for non efficiency should be integer. Default value used.");
        }

        boolean control = false;
        try {
            control = cmd.hasOption("c");
        } catch (Exception e) {
            log.log(Level.WARNING, "Control dataset was not specified. Whole pipeline will be used.");
        }
        if (control) {
            log.log(Level.WARNING, "Control dataset was specified. Only second algorithm will be used.");
        }

        VersionAndParams topStrings = new VersionAndParams(minCorThreshold, maxNumOfModels, numOfModelForNonEfficiency,
                minNumOfModelsForOutlierDetection, lowerBoundOutlier, upperBoundOutlier,
                minDistanceBetween, thresholdForInefficiency);
        String versionAndParams = topStrings.getVersAndParams();

        ReaderBED readerBed = new ReaderBED(cmd.getOptionValue("b"));
        ReaderData readerData = new ReaderData(cmd.getOptionValue("d"));
        int numberOfTestsPerSample = readerBed.getNumberOfExones();
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
        int numberOfSamples = samplesNames.size();

        Comparator<Pair> comparator = new PairComparator();

        HashMap<String, PriorityQueue<Pair>> amplCorrelationPriority = new HashMap<String, PriorityQueue<Pair>>();
        HashMap<String, PriorityQueue<Pair>> amplCorrelationPriorityForLDA = new HashMap<String, PriorityQueue<Pair>>();



        HashMap<String, HashMap<String, Integer>> sampleToAmplicons = new HashMap<String, HashMap<String, Integer>>();
        ArrayList<String> nonEffecitveAmpls = new ArrayList<String>();

        for (String elem : samplesNames) {
            HashMap<String, Integer> tmpHash = new HashMap<String, Integer>();
            for (String ampl : dataAboutAmplicons.keySet()) {
                tmpHash.put(ampl, 0);
            }
            sampleToAmplicons.put(elem, tmpHash);
        }


        amplCorrelationPriority = getAmplCorrelationPriorityRob(dataAboutAmplicons,
                samplesNames, comparator);
        HashMap<String, HashMap<String, Double> > predictedMedians = new HashMap<String, HashMap<String, Double>>();
        for (Map.Entry<String, Amplicon> entry1 : dataAboutAmplicons.entrySet()) {
            PriorityQueue<Pair> priorityByCorrelation = null;
            try {
                priorityByCorrelation = (PriorityQueue<Pair>)ObjectCloner.deepCopy(amplCorrelationPriority.get(entry1.getKey()));
            } catch (Exception e) {
                log.log(Level.SEVERE, "Prbolems with deep copy. \n" + e.getMessage());
            }
            predictedMedians.put(entry1.getKey(), new HashMap<String, Double>());

            HashMap<String, ArrayList<Double>> predictedValues = new HashMap<String, ArrayList<Double>>();
            for (int i = 0; i < entry1.getValue().getCoverages().size(); i++) {
                predictedValues.put(samplesNames.get(i), new ArrayList<Double>());
            }
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
                        lmDel.setStudentizedResidualsShift(entry1.getValue().getCoverages(), 1);
                        LinearModel lmDup = new LinearModel(dataAboutAmplicons.get(e.getAmplName()).getCoverages(), entry1.getValue().getShiftedCoverages(0.5849625));
                        lmDup.setStudentizedResidualsShift(entry1.getValue().getCoverages(), 1);

                        ArrayList<Double> predictions = lm.getPredictions();

                        for (int i = 0; i < entry1.getValue().getCoverages().size(); i++) {
                            predictedValues.get(samplesNames.get(i)).add(predictions.get(i));
                        }

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
            for (int i = 0; i < entry1.getValue().getCoverages().size(); i++) {
                predictedMedians.get(entry1.getKey()).put(samplesNames.get(i), Statistics.med(predictedValues.get(samplesNames.get(i))));
            }


            if (counter < numOfModelForNonEfficiency) {
                log.log(Level.FINE, "Amplicon " + entry1.getValue().returnID() + " does not have enough located far enough and correlated amplicons for analysis");
                total_size--;
                nonEffecitveAmpls.add(entry1.getValue().returnID());
            } else {
                log.log(Level.FINE, "Analysis of " + entry1.getValue().returnID() + " was completed succesfully.");
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
            if (control) {
                dirtySamples.clear();
                for (String name : samplesNames) {
                    if (name.startsWith("Case_")) {
                        dirtySamples.add(name);
                    }
                }
            }
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
                for (int i = 0; i < maximum_cycles_of_algorithm; i++) {
                    log.log(Level.FINE, "LDA Step " + i);
                    try {
                        amplCorrelationPriorityForLDA = (HashMap<String, PriorityQueue<Pair>>)ObjectCloner.deepCopy(amplCorrelationPriority);
                    } catch (Exception e) {
                        log.log(Level.SEVERE, "Prbolems with deep copy. \n" + e.getMessage());
                    }
                    sam = (Samples) (ObjectCloner.deepCopy(samClone));
                    lda = new LDA(sam, dataAboutAmplicons, samplesNames,
                            dirtySamples, nonEffecitveAmpls, amplCorrelationPriorityForLDA, minCorThreshold,
                            minDistanceBetween, maxNumOfModels, numberOfTestsPerSample, predictedMedians);
                    lda.outRes("./tmpResultsCNV/" + filename + "_step_" + i, versionAndParams, samplesNames, readerBed, cmd.getOptionValue("d"));
                    dirtySamples = lda.getNewDirtySamples();
                    if (lda.numberOfDirtySamplesChanged == false || i == maximum_cycles_of_algorithm) {
                        lda.outRes(filename, versionAndParams, samplesNames, readerBed, cmd.getOptionValue("d"));
                        break;
                    }
                }
            } else {
                LDA theBestLDA = null;

                dirtySamples = (ArrayList<String>) (ObjectCloner.deepCopy(dirtySamplesAfterFirstStep));
                for (int i = 0; i < maximum_cycles_of_algorithm; i++) {
                    log.log(Level.FINE, "LDA Step " + i);
                    if (samplesNames.size() - dirtySamples.size() > Math.max(samplesNames.size() / 2 + 1, 20)) {
                        try {
                            amplCorrelationPriorityForLDA = (HashMap<String, PriorityQueue<Pair>>)ObjectCloner.deepCopy(amplCorrelationPriority);
                        } catch (Exception e) {
                            log.log(Level.SEVERE, "Prbolems with deep copy. \n" + e.getMessage());
                        }
                    } else {
                        try {
                            amplCorrelationPriorityForLDA = (HashMap<String, PriorityQueue<Pair>>)ObjectCloner.deepCopy(amplCorrelationPriority);
                        } catch (Exception e) {
                            log.log(Level.SEVERE, "Prbolems with deep copy. \n" + e.getMessage());
                        }
                    }

                    sam = (Samples) (ObjectCloner.deepCopy(samClone));
                    lda = new LDA(sam, dataAboutAmplicons, samplesNames,
                            dirtySamples, nonEffecitveAmpls, amplCorrelationPriorityForLDA, minCorThreshold,
                            minDistanceBetween, maxNumOfModels, numberOfTestsPerSample, predictedMedians);

                    dirtySamples = lda.getNewDirtySamples();

                    lda.outRes("./tmpResultsCNV/" + filename + "_" + dirtySamples.size() + "_step_" + i, versionAndParams, samplesNames, readerBed, cmd.getOptionValue("d"));

                    if (lda.numberOfDirtySamplesChanged == false) {
                        theBestLDA = (LDA)ObjectCloner.deepCopy(lda);
                        theBestLDA.printDistances("distance.xls");
                        log.log(Level.INFO, "The analysis is completed");
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

