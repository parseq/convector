package deletionsAnalysis;

/**
 * Created by german on 09.10.14.
 */
public class VersionAndParams {
    /*
    Just a generator of the top of file.
    @params parameters from command line (or their default values, recommended by us)
     */
    private String versAndParams;

    public VersionAndParams(Double minCorThreshold, Integer maxNumOfModels, Integer numOfModelForNonEfficiency,
                            Integer minNumOfModelsForOutlierDetection, Double lowerBoundOutlier, Double upperBoundOutlier,
                            Integer minDistanceBetween, Integer thresholdForInefficiency) {
        String versionAndParams = "version 2.0 ";
        versionAndParams = versionAndParams + ", minimum correlation for models: " + minCorThreshold;
        versionAndParams = versionAndParams + ", \nmaximum number of models to test against: " + maxNumOfModels;
        versionAndParams = versionAndParams + ", \nminimum number of models for efficient estimation: " + numOfModelForNonEfficiency;
        versionAndParams = versionAndParams + ", \nnumber of models where outliers where detected: " + minNumOfModelsForOutlierDetection;
        versionAndParams = versionAndParams + ", \nlower bound for outlier detection: " + lowerBoundOutlier;
        versionAndParams = versionAndParams + ", \nupper bound for outlier detection: " + upperBoundOutlier;
        versionAndParams = versionAndParams + ", \nmin distance between correlated amplicons: " + minDistanceBetween;
        versionAndParams = versionAndParams + ", \nthreshold (SM) for inefficiency: " + thresholdForInefficiency;
        versionAndParams = versionAndParams + ", method: data came from different models";

        versAndParams = versionAndParams;

    }

    public String getVersAndParams() {
        return versAndParams;
    }
}
