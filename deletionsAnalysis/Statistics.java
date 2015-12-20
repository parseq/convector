package deletionsAnalysis;

import java.util.ArrayList;
import java.util.Collections;
import java.util.Random;
import java.util.concurrent.*;

/**
 * Created by german on 11.07.14.
 */
public class Statistics {
    /*
    Do not touch anything here.
    */
    public static <T extends Number> Double sm(ArrayList<T> elements) {
        Double sum = 0.0;
        if (!elements.isEmpty()) {
            for (T elem : elements) {
                sum += elem.doubleValue();
            }
            return sum.doubleValue() / elements.size();
        }
        return sum;
    }

    public static double sv(ArrayList<Double> elements) {
        double sampleVariance = 0.0;
        if (!elements.isEmpty()) {
            double sampleMean = Statistics.sm(elements);
            double sampleMeanSquared = Math.pow(sampleMean, 2);
            for (double number : elements) {
                sampleVariance += (Math.pow(number, 2) - sampleMeanSquared);
            }
            sampleVariance /= (elements.size() - 1);
        }
        return sampleVariance;
    }

    public static double ssd(ArrayList<Double> elements) {
        return Math.sqrt(Statistics.sv(elements));
    }

    public static double med(ArrayList<Double> elements) {
        double median = 0.0;
        ArrayList<Double> copyList = new ArrayList<Double>();
        if (!elements.isEmpty()) {
            for (Double elem : elements) {
                copyList.add(elem);
            }
            Integer lenOfArray = copyList.size();
            Collections.sort(copyList);
            if (lenOfArray % 2 == 0) {
                Integer firstIndex = lenOfArray / 2;
                median = (copyList.get(firstIndex) + copyList.get(firstIndex - 1)) / 2;
            } else
                median = copyList.get((lenOfArray - 1) / 2);
        }
        return median;
    }


    public static double medW(ArrayList<Double> x) {
        double medW = 0.0;
        if (!x.isEmpty()) {
            ArrayList<Double> tmpX = new ArrayList<Double>();
            for (int i = 0; i < x.size(); i++) {
                for (int j = 0; j <= i; j++) {
                    tmpX.add(0.5 * (x.get(i) + x.get(j)));
                }
            }
            medW = med(tmpX);
        }
        return medW;
    }

    public static double MAD(ArrayList<Double> x) {
        double mad = 0.0;
        double medianOfX = 0.0;
        if (!x.isEmpty()) {
            medianOfX = Statistics.med(x);
            ArrayList<Double> tmpX = new ArrayList<Double>();
            for (Double elem : x) {
                tmpX.add(Math.abs(elem - medianOfX));
            }
            mad = Statistics.med(tmpX);
        }
        return mad;
    }

    public static ArrayList<Double> wavedArray(ArrayList<Double> x) {
        ArrayList<Double> result = new ArrayList<Double>();
        if (!x.isEmpty()) {

            Double madX = Statistics.MAD(x);
            Double medX = Statistics.med(x);
            for (int i = 0; i < x.size(); i++) {
                result.add((x.get(i) - medX) / (Math.sqrt(2) * madX));
            }
        }
        return result;
    }

    public static ArrayList<Double> diffOfArrays(ArrayList<Double> x, ArrayList<Double> y) {
        ArrayList<Double> result = new ArrayList<Double>();
        if (!x.isEmpty() && x.size() == y.size()) {
            for (int i = 0; i < x.size(); i++) {
                result.add(x.get(i) - y.get(i));
            }
        }
        return result;
    }

    public static ArrayList<Double> sumOfArrays(ArrayList<Double> x, ArrayList<Double> y) {
        ArrayList<Double> result = new ArrayList<Double>();
        if (!x.isEmpty() && x.size() == y.size()) {
            for (int i = 0; i < x.size(); i++) {
                result.add(x.get(i) + y.get(i));
            }
        }
        return result;
    }

    public static ArrayList<Double> trimmedSumOfResids(ArrayList<Double> lst, Integer r) {
        Double scaleEstimation = medW(lst);
        Double sum = 0.0;
        ArrayList<Double> copy = new ArrayList<Double>();
        for (Double elem : lst) {
            copy.add(elem - scaleEstimation);
        }
        Collections.sort(copy);
        ArrayList<Double> result = new ArrayList<Double>();
        for (int i = r; i < lst.size() - r; i++) {
            result.add(copy.get(i));
        }
        return result;
    }

    public static Double robustCor(ArrayList<Double> first, ArrayList<Double> second) {
        Double robCor = 0.0;
        if (!first.isEmpty() && first.size() == second.size()) {
            ArrayList<Double> sumOfWaved = Statistics.sumOfArrays(Statistics.wavedArray(first), Statistics.wavedArray(second));
            ArrayList<Double> diffOfWaved = Statistics.diffOfArrays(Statistics.wavedArray(first), Statistics.wavedArray(second));

            Double numerator = Math.pow(Statistics.snEstimator(sumOfWaved), 2) - Math.pow(Statistics.snEstimator(diffOfWaved), 2);
            Double denominator = Math.pow(Statistics.snEstimator(sumOfWaved), 2) + Math.pow(Statistics.snEstimator(diffOfWaved), 2);

            robCor = numerator / denominator;
        }
        return robCor;
    }

    public static Double robustMADCor(ArrayList<Double> first, ArrayList<Double> second) {
        Double robCor = 0.0;
        if (!first.isEmpty() && first.size() == second.size()) {
            ArrayList<Double> sumOfWaved = Statistics.sumOfArrays(Statistics.wavedArray(first), Statistics.wavedArray(second));
            ArrayList<Double> diffOfWaved = Statistics.diffOfArrays(Statistics.wavedArray(first), Statistics.wavedArray(second));

            Double numerator = Math.pow(Statistics.MAD(sumOfWaved), 2) - Math.pow(Statistics.MAD(diffOfWaved), 2);
            Double denominator = Math.pow(Statistics.MAD(sumOfWaved), 2) + Math.pow(Statistics.MAD(diffOfWaved), 2);

            robCor = numerator / denominator;
        }
        return robCor;
    }

    public static Double cor(ArrayList<Double> first, ArrayList<Double> second) {
        Double cor = 0.0;

        Double numerator = 0.0;
        Double denominatorFirstSum = 0.0;
        Double denominatorSecondSum = 0.0;

        if (!first.isEmpty() && first.size() == second.size()) {
            Double firstMean = sm(first);
            Double secondMean = sm(second);
            for (int i = 0; i < first.size(); i++) {
                numerator += (first.get(i) - firstMean) * (second.get(i) - secondMean);
                denominatorFirstSum += Math.pow(first.get(i) - firstMean, 2);
                denominatorSecondSum += Math.pow(second.get(i) - secondMean, 2);
            }
        }
        cor = numerator / Math.sqrt(denominatorFirstSum * denominatorSecondSum);
        return cor;
    }


    public static Double iqrVarEstimation(ArrayList<Double> toEstimate) {
        Double varEst = 0.0;
        ArrayList<Double> copyList = new ArrayList<Double>();
        if (!toEstimate.isEmpty()) {
            for (Double elem : toEstimate) {
                copyList.add(elem);
            }
            Collections.sort(copyList);
            if (toEstimate.size() % 2 == 0) {
                ArrayList<Double> toEstimateRight = new ArrayList<Double>(copyList.subList(0, copyList.size() / 2));
                ArrayList<Double> toEstimateLeft = new ArrayList<Double>(copyList.subList(copyList.size() / 2, copyList.size()));
                varEst = (Statistics.med(toEstimateLeft) - Statistics.med(toEstimateRight)) / 1.349;
            } else {
                ArrayList<Double> toEstimateRight = new ArrayList<Double>(copyList.subList(0, (copyList.size() - 1) / 2));
                ArrayList<Double> toEstimateLeft = new ArrayList<Double>(copyList.subList((copyList.size() + 1) / 2, copyList.size()));
                varEst = (Statistics.med(toEstimateLeft) - Statistics.med(toEstimateRight)) / 1.349;
            }
        }
        return varEst;
    }

    public static Double marVarResidEstimation(ArrayList<Double> toEstimate) {
        Double varEst = 0.0;
        ArrayList<Double> copyList = new ArrayList<Double>();
        if (!toEstimate.isEmpty()) {
            for (Double elem : toEstimate) {
                copyList.add(Math.abs(elem));
            }
            Collections.sort(copyList);
            varEst = Statistics.med(copyList) / 0.6745;
        }
        return varEst;
    }


    public static double himed(ArrayList<Double> elements) {
        Double median = 0.0;
        ArrayList<Double> copyList = new ArrayList<Double>();
        if (!elements.isEmpty()) {
            for (Double elem : elements) {
                copyList.add(elem);
            }
            Integer lenOfArray = copyList.size();
            Collections.sort(copyList);
            Integer firstIndex = ((lenOfArray) / 2) + 1;
            median = (copyList.get(firstIndex - 1));
        }
        return median;
    }

    public static double lomed(ArrayList<Double> elements) {
        Double median = 0.0;
        ArrayList<Double> copyList = new ArrayList<Double>();
        if (!elements.isEmpty()) {
            for (Double elem : elements) {
                copyList.add(elem);
            }
            Integer lenOfArray = copyList.size();
            Collections.sort(copyList);
            Integer firstIndex = (((lenOfArray + 1) / 2));
            median = (copyList.get(firstIndex - 1));
        }
        return median;
    }

    public static double getT(double alpha) {
        return Math.log(1 / alpha);
    }

    public static double getQuantile(int df, double alpha) {
        double ans = 0.0;
        if (df >= 3 && df <= 250) {
            double t = getT(alpha);
            ans = ans + df + 2 * t + 1.62 * Math.sqrt(df * t) + 0.63012 * Math.sqrt(df) * Math.log(t) - 1.12032 * Math.sqrt(df) -
                    2.48 * Math.sqrt(t) - 0.65381 * Math.log(t) - 0.22872;
            return ans;
        } else {
            if (df == 2) {

            } else if (df == 1) {

            } else if (df > 250) {
                System.err.println("Calculation of quantile with more DFs!");
            }
        }
        return ans;
    }


    public static Double snEstimator(ArrayList<Double> toEstimate) {
        Double varEst = 0.0;
        ArrayList<Double> resultList = new ArrayList<Double>();

        if (!toEstimate.isEmpty()) {
            for (int i = 0; i < toEstimate.size(); i++) {
                ArrayList<Double> tmpList = new ArrayList<Double>();
                for (int j = 0; j < toEstimate.size(); j++) {
                    tmpList.add(Math.abs(toEstimate.get(i) - toEstimate.get(j)));
                }
                resultList.add(Statistics.himed(tmpList));
            }
            Double multiplier = 1.0;
            if (toEstimate.size() < 10) {
                Double[] multipliers = {0.743, 1.851, 0.954, 1.351, 0.993, 1.198, 1.005, 1.131};
                multiplier = multipliers[toEstimate.size() - 2];
            } else  if (toEstimate.size() % 2 != 0) {
                multiplier = 1.0 * toEstimate.size() / (1.0 * toEstimate.size() - 0.9);
            }
            varEst = 1.1926 * multiplier * Statistics.lomed(resultList);
        }

        return varEst;
    }



    public static Double snEstimatorRob(ArrayList<Double> toEstimate, int iterations, Random rand) {
        Double varEst = 0.0;
        ArrayList<Double> resultList = new ArrayList<Double>();
        if (!toEstimate.isEmpty()) {
            ArrayList<Double> sns = new ArrayList<Double>();


            for (int j = 0; j < iterations; j++) {

                ArrayList<Double> first = new ArrayList<Double>();

                for (int k = 0; k < toEstimate.size(); k++) {
                    int randomNum = rand.nextInt((toEstimate.size()));
                    first.add(toEstimate.get(randomNum));
                }

                sns.add(snEstimator(first));
            }
            varEst = medW(sns);
        }

        return varEst;
    }

    public static Double snEstimatorSim(ArrayList<Double> toEstimate) {
        Random rand = new Random();
        int iterations = 10;
        Double varEst = 0.0;
        ArrayList<Double> resultList = new ArrayList<Double>();
        if (!toEstimate.isEmpty()) {
            ArrayList<Double> sns = new ArrayList<Double>();


            for (int j = 0; j < iterations; j++) {

                ArrayList<Double> first = new ArrayList<Double>();

                for (int k = 0; k < toEstimate.size(); k++) {
                    int randomNum = rand.nextInt((toEstimate.size()));
                    first.add(toEstimate.get(randomNum));
                }

                sns.add(snEstimator(first));
            }
            varEst = medW(sns);
        }

        return varEst;
    }


    /*public static Double qnEstimator(ArrayList<Double> toEstimate) {
        Double varEst = 0.0;
        ArrayList<Double> resultList = new ArrayList<Double>();
        if (!toEstimate.isEmpty()) {
            for (int i = 0; i < toEstimate.size() - 1; i++) {
                for (int j = i + 1; j < toEstimate.size(); j++) {
                    resultList.add(Math.abs(toEstimate.get(i) - toEstimate.get(j)));
                }
            }

            Collections.sort(resultList);
            int index = (int) binomial((int) Math.round(0.5 * toEstimate.size()) + 1, 2);

            varEst = 2.21914446 * resultList.get(index);
        }

        return varEst;
    }

    private static long binomial(long n, long k)
    {
        if (k>n-k)
            k=n-k;

        long b=1;
        for (long i=1, m=n; i<=k; i++, m--)
            b=b*m/i;
        return b;
    }*/
}
