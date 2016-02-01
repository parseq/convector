package parseq.chimeric_solver;

import java.util.ArrayList;
import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.util.concurrent.*;

final class SmithWaterman implements Callable<Boolean> {
    private Integer maxPercOfErrors = 20;
    private Double gapCost = -1.0;
    private Double matchCost = 1.0;
    private Double misMatchCost = -1.0;
    private String chimera;
    private String reference;

    public SmithWaterman(String chim, String ref) {
        chimera = chim;
        reference = ref;
    }

    public Boolean call() {
        if (chimera.equals(reference)) {
            return false;
        }

        int maxNumOfErrors = (int) Math.ceil(0.01 * chimera.length() * maxPercOfErrors) + 1;
        int minScoreToAlign = (chimera.length() - maxNumOfErrors);

            Double maxScore = 0.0;
            Double[][] matrixOfDp = new Double[reference.length() + 1][chimera.length() + 1];
            for (int i = 0; i < reference.length() + 1; i++) {
                for (int j = 0; j < chimera.length() + 1; j++) {
                    matrixOfDp[i][j] = 0.0;
                }
            }
            for (int i = 1; i < 1 + reference.length(); i++) {
                    for (int j = 1; j < chimera.length() + 1; j++) {
                        Double leftGap = gapCost;
                        Double rightGap = gapCost;
                        Double match = 0.0;
                        if (reference.charAt(i - 1) == chimera.charAt(j - 1)) {
                            match = matchCost;
                        } else if (reference.charAt(i - 1) == '*') {
                            match = 0.0;
                            leftGap = 0.0;
                            rightGap = 0.0;
                        } else match = misMatchCost;
                        matrixOfDp[i][j] = Math.max(Math.max(matrixOfDp[i][j - 1] + leftGap, matrixOfDp[i - 1][j]) + rightGap, Math.max(matrixOfDp[i - 1][j - 1] + match, 0));
                        if (matrixOfDp[i][j] > maxScore) {
                            maxScore = matrixOfDp[i][j];
                        }
                        if (maxScore > minScoreToAlign) {
                            return true;
                        }
                    }
                }


        return false;
    }
}

class Solver{

    private ArrayList<Boolean> canBeAligned = new ArrayList<Boolean>();
    private ArrayList<Integer> coverageIncrease = new ArrayList<Integer>();

    private ArrayList<String> references = new ArrayList<String>();
    private ArrayList<String> chimeras = new ArrayList<String>();

    private String fileWithChimeras = "tmp_chimeras.txt";
    private String fileWithReferences = "tmp_references.txt";
    private String tmpOutput = "tmp_output.txt";

    private void parseFiles() {
        try {
            BufferedReader br = new BufferedReader(new FileReader(fileWithChimeras));
            for (String line; (line = br.readLine()) != null; ) {
                String[] tmpListForMap = line.split("\\t");
                chimeras.add(tmpListForMap[0]);
                canBeAligned.add(false);
            }
        } catch (FileNotFoundException e) {
            System.err.println("\nSorry, but file with chimeras was not found.");
            System.err.println(e.getMessage());
            System.exit(1);
        } catch (IOException e) {
            System.err.println("\nSorry, but we had some problems with reading from your file with chimeras.");
            System.err.println(e.getMessage());
            System.err.println("Please, be sure that this tool has access to your input files.");
            System.exit(1);
        }

        try {
            BufferedReader br = new BufferedReader(new FileReader(fileWithReferences));
            for (String line; (line = br.readLine()) != null; ) {
                String[] tmpListForMap = line.split("\\t");
                references.add(tmpListForMap[0]);
                coverageIncrease.add(0);
            }
        } catch (FileNotFoundException e) {
            System.err.println("\nSorry, but file with references was not found.");
            System.err.println(e.getMessage());
            System.exit(1);
        } catch (IOException e) {
            System.err.println("\nSorry, but we had some problems with reading from your file with chimeras.");
            System.err.println(e.getMessage());
            System.err.println("Please, be sure that this tool has access to your input files.");
            System.exit(1);
        }
    }

    private void writeOutput() {
        Path file = Paths.get(tmpOutput);
        try {
            // Create the empty file with default permissions, etc.
            Files.delete(file);
            Files.createFile(file);
        } catch (IOException e) {
            // Some other sort of failure, such as permissions.
            System.err.format("createFile error: %s%n", e);
        }

        try {
            FileWriter fw = new FileWriter(tmpOutput);
            BufferedWriter bw = new BufferedWriter(fw);
            for (Integer elem : coverageIncrease) {
                bw.write(elem + "\n");
            }
            bw.close();
        } catch (IOException e) {
            System.err.format("smth wrong with output file");
        }


    }
    public Solver() {
        parseFiles();
        try {
            for (int i = 0; i < chimeras.size(); i++) {
                ArrayList<Future<Boolean>> futures = new ArrayList<Future<Boolean>>(references.size());

                ExecutorService service = Executors.newFixedThreadPool(16);
                CompletionService<Boolean> pool = new ExecutorCompletionService<Boolean>(service);

                for (int j = 0; j < references.size(); j++) {
                    futures.add(pool.submit(new SmithWaterman(chimeras.get(i), references.get(j))));
                }

                for (int j = 0; j < references.size(); j++) {
                    Future<Boolean> future = futures.get(j);
                    Boolean result = future.get();
                        if (result) {
                            canBeAligned.set(i, true);
                            coverageIncrease.set(j, coverageIncrease.get(j) + 1);
                            service.shutdown();
                            break;
                        }
                    }
                service.shutdown();
            }
        }
        catch (Exception e) {
            System.err.println("Interrupted exception\n" + e.getStackTrace());
        }
        writeOutput();
    }
}

public class Main {


    public static void main(String[] args) {
        Solver solve = new Solver();
    }
}
