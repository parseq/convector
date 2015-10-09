package deletionsAnalysis;

import java.io.Serializable;

/**
 * Created by german on 11.08.14.
 */
public class Pair implements Serializable {
    /*
    pair amplicon name - correlation for model choosing.
     */
    private double correlation;
    private String amplName;

    public Pair(double cor, String str) {
        correlation = cor;
        amplName = str;
    }

    public double getCorrelation() {
        return correlation;
    }

    public String getAmplName() {
        return amplName;
    }

}
