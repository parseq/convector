package deletionsAnalysis;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by german on 08.08.14.
 */

class FindDuplicateInList {
    /*
    Simple procedure for finding duplicates in xVals.
    @params List to find duplicates in
    @return set of duplicates in list
     */
    public static <T extends Number> Set<T> findDuplicates(ArrayList<T> listContainingDuplicates) {

        final Set<T> setToReturn = new HashSet<T>();
        final Set<T> set1 = new HashSet<T>();

        for (T number : listContainingDuplicates) {
            if (!set1.add(number)) {
                setToReturn.add(number);
            }
        }
        return setToReturn;
    }
}

public class LinearModel {
    private double slope;
    private double intercept;
    private double sigmaEstimated;
    private ArrayList<Double> xVals;
    private ArrayList<Double> yVals;
    private ArrayList<Double> residuals;

    private ArrayList<Double> studentizedResiduals;
    private ArrayList<Double> studentizedResidualsShift;
    private static final Logger log = Logger.getLogger( Solver.class.getName() );


    public LinearModel(ArrayList<Double> x, ArrayList<Double> y) {
        /*
        constructor of Theil-Sen estimator
        @params two ArrayLists, one is predictors and one is response
        @return construct linear model
         */
        if (x.size() == y.size()) {
            xVals = new ArrayList<Double>(x);
            yVals = new ArrayList<Double>(y);
            slope = 0.0;
            intercept = 0.0;
            determineSlope();
            determineIntercept();
            residuals = new ArrayList<Double>();

            studentizedResiduals = new ArrayList<Double>();
            studentizedResidualsShift = new ArrayList<Double>();

            findResiduals();
            setStudentizedResiduals();
        } else {
            log.log(Level.WARNING, "An array of x-coords and y-coords differ. It may happen in case you have missing values of coverages.");
        }
    }

    private void determineSlope() {
        /*
        determines slope of the linear regression model
         */
        Set<Double> duplicates = FindDuplicateInList.findDuplicates(xVals);
        ArrayList<Double> slopes = new ArrayList<Double>();
        for (int i = 0; i < xVals.size(); i++) {
            for (int j = 0; j < xVals.size(); j++) {
                if ((i < j) && !((duplicates.contains(xVals.get(i))) || (duplicates.contains(xVals.get(j))))) {
                    slopes.add((yVals.get(i) - yVals.get(j)) / (xVals.get(i) - xVals.get(j)));
                }
            }
        }
        slope = Statistics.med(slopes);
    }

    private void determineIntercept() {
        /*
        determines intercept of the linear regression model
         */
        ArrayList<Double> intercepts = new ArrayList<Double>();

        for (int i = 0; i < xVals.size(); i++) {
            intercepts.add(yVals.get(i) - slope * xVals.get(i));
        }
        intercept = Statistics.med(intercepts);
    }

    private void findResiduals() {
        /*
        find residuals of linear regression model
         */
        for (int i = 0; i < xVals.size(); i++) {
            residuals.add(yVals.get(i) - xVals.get(i) * slope - intercept);
        }
    }

    private void setStudentizedResiduals() {
        /*
        find robustly Studentized residuals for wild type amplicons
         */
        findResiduals();
        this.sigmaEstimated = Statistics.snEstimator(residuals);

        double sumXSquared = 0.0;
        double squaredSumX = 0.0;
        for (int i = 0; i < xVals.size(); i++) {
            sumXSquared += Math.pow(xVals.get(i), 2);
            squaredSumX += xVals.get(i);
        }
        double mean = Statistics.sm(xVals);
        double denomWikiCalc = 0.0;
        for (int i = 0; i < xVals.size(); i++) {
            denomWikiCalc += Math.pow(xVals.get(i) - mean, 2);
        }
        for (int i = 0; i < xVals.size(); i++) {
            double denominatorWithoutSigma = 1 - 1.0 / xVals.size() - (xVals.get(i) - mean) / denomWikiCalc;
            studentizedResiduals.add(residuals.get(i) / (sigmaEstimated * denominatorWithoutSigma));
        }
    }

    public void setStudentizedResidualsShift(ArrayList<Double> shiftedY, double multiplier) {
        /*
        find robustly Studentized residuals for deleterious and duplicated amplicons
         */
        yVals = shiftedY;
        residuals = new ArrayList<Double>();
        findResiduals();
        this.sigmaEstimated = Statistics.snEstimator(residuals);

        double sumXSquared = 0.0;
        double squaredSumX = 0.0;
        for (int i = 0; i < xVals.size(); i++) {
            sumXSquared += Math.pow(xVals.get(i), 2);
            squaredSumX += xVals.get(i);
        }
        double mean = Statistics.sm(xVals);
        double denomWikiCalc = 0.0;
        for (int i = 0; i < xVals.size(); i++) {
            denomWikiCalc += Math.pow(xVals.get(i) - mean, 2);
        }
        sigmaEstimated *= multiplier;
        for (int i = 0; i < xVals.size(); i++) {
            double denominatorWithoutSigma = 1 - 1.0 / xVals.size() - (xVals.get(i) - mean) / denomWikiCalc;
            studentizedResidualsShift.add(residuals.get(i) / (sigmaEstimated * denominatorWithoutSigma));
        }
    }

    public ArrayList<Double> getInternallyStudentizedResiduals() {
        /*
        @return studentized residuals with snEstimator of standard deviation
         */
        return studentizedResiduals;
    }

    public ArrayList<Double> getInternallyStudentizedResidualsShift() {
        /*
        @return studentized residuals with snEstimator of standard deviation
         */
        return studentizedResidualsShift;
    }

    public ArrayList<Double> getPredictions() {
        ArrayList<Double> predicted = new ArrayList<Double>();
        for (Double val: xVals) {
            predicted.add(val * slope + intercept);
        }
        return predicted;
    }
}
