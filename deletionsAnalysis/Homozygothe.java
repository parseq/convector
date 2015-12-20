package deletionsAnalysis;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 * Created by german on 08.08.14.
 */

public class Homozygothe {
    private static DefaultDict<Integer, ArrayList<String>> deleteriousAmpls = new DefaultDict<Integer, ArrayList<String>>(ArrayList.class);
    private static final Logger log = Logger.getLogger( Solver.class.getName() );

    public Homozygothe(HashMap<String, Amplicon> dataAboutAmplicons, ArrayList<String> nonEffecitveAmpls) {
        /*
        @params info about amplicons and list of non effective amplicons
        @return starts a function that prints amplicons/samples with extremely low coverages
         */
        findHomozygotheSample(dataAboutAmplicons, nonEffecitveAmpls);
    }

    public static void reportHomozygotheSample(ArrayList<String> samplesNames) {
        /*
        @params names of samples
         */
        ArrayList<Integer> listToExclude = new ArrayList<Integer>();
        for (Map.Entry<Integer, ArrayList<String>> entry : deleteriousAmpls.entrySet()) {
            listToExclude.add(entry.getKey());
            log.log(Level.FINE, samplesNames.get(entry.getKey()));
            for (int i = 0; i < entry.getValue().size(); i++) {
                log.log(Level.FINE, entry.getValue().get(i));
            }
        }
    }

    private static void findHomozygotheSample(HashMap<String, Amplicon> dataAboutAmplicons, ArrayList<String> nonEffectiveAmpls) {
        /*
        @param Map, {Amplicon name : Amplicon object}, list with names of non effective amplicons
         */
        for (Map.Entry<String, Amplicon> entry : dataAboutAmplicons.entrySet()) {
            for (int i = 0; i < entry.getValue().getCoverages().size(); i++) {
                if (entry.getValue().getCoverages().get(i) == 0 && !nonEffectiveAmpls.contains(entry.getKey())) {
                    deleteriousAmpls.get(i).add(entry.getValue().returnID());
                }
            }
        }
    }
}
