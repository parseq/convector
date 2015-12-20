#!/usr/bin/env python2
__author__ = 'german'

import sys
import logging
from logging.config import fileConfig
import os
import argparse
import time

loginipath = ('logging_config.ini')
fileConfig(loginipath, defaults={'logfilename': 'pipeline.log'})
logger = logging.getLogger('CONVector_logger')

def get_coverages(bam_folder, bed_file, id_of_run):
    """Launch chimeric_solver and counts coverages. Save the output .xls file to the directory ./result by default"""
    string_to_cmd = ("").join(["python2 chimeric_solver.py -d ", bam_folder, " -b ", bed_file, " -r ", id_of_run, ".xls", " --conv"])
    os.system(string_to_cmd)

def get_quantiles():
    """Reads files with hard coded quantiles and gets quantiles, returns dictionaries, {degree of freedom : quantile} """
    quant_0995 = {}
    quant_099 = {}
    quant_098 = {}
    quant_095 = {}
    with open(os.path.abspath("./quantiles/quantiles0995.txt")) as f:
        table = f.readline()
        array_of_quantiles = table.split()
        for i in xrange(2, 500):
            quant_0995[i] = float(array_of_quantiles[i - 2])
    with open(os.path.abspath("./quantiles/quantiles099.txt")) as f:
        table = f.readline()
        array_of_quantiles = table.split()
        for i in xrange(2, 500):
            quant_099[i] = float(array_of_quantiles[i - 2])
    with open(os.path.abspath("./quantiles/quantiles098.txt")) as f:
        table = f.readline()
        array_of_quantiles = table.split()
        for i in xrange(2, 500):
            quant_098[i] = float(array_of_quantiles[i - 2])
    with open(os.path.abspath("./quantiles/quantiles095.txt")) as f:
        table = f.readline()
        array_of_quantiles = table.split()
        for i in xrange(2, 500):
            quant_095[i] = float(array_of_quantiles[i - 2])

    return quant_0995, quant_099, quant_098, quant_095

def get_tmp_del_files(bed_file, quant_dict_del, quant_dict_dup, id_of_run, min_cor, mode):
    """Launch CONVector using predefined parameters, output - files with CNVs
    (outputId_of_task_before_LDA and outputId_of_task_after_LDA)"""
    num_of_samples = 0
    with open("./result/" + id_of_run + "_qc.xls") as f:
        num_of_samples = len(f.readline().split()) - 2
    x = str(quant_dict_del[num_of_samples])
    y = str(quant_dict_dup[num_of_samples])
    string_to_cmd = "javac ./deletionsAnalysis/Main.java"
    os.system(string_to_cmd)
    string_to_cmd = ("").join(["java  deletionsAnalysis/Main ", " -d ./result/", id_of_run,
                               "_qc.xls", " -b ", bed_file, " -f output", id_of_run, " -mc ", min_cor, " -mnm 5 -nne 4 -nod 4 -lb -",
                               x, " -ub ", y , " -dist 1000000 -lcb 25000 -lca 50"])
    if mode:
        string_to_cmd += " -c"
    logger.warn(string_to_cmd)
    os.system(string_to_cmd)

def get_results(id_of_run, output_folder):
    """Launch finalizer and creates a report about CNVs, .xls file with name
    result_before_LDA_id_of_task.xls and result_after_LDA_id_of_task.xls"""
    string_to_get_final_results_unsupervised = ("").join(["python2 finalizer.py -i output", id_of_run, ".xls -o result_before_LDA_", id_of_run, ".xls"])
    string_to_get_final_results_supervised = ("").join(["python2 finalizer.py -i output", id_of_run, "_after_LDA.xls -o result_after_LDA_", id_of_run, ".xls"])
    os.system(string_to_get_final_results_unsupervised)
    os.system(string_to_get_final_results_supervised)

def main():
    parser = argparse.ArgumentParser(description="""This is a pipeline for the detection of CNV in the data
    obtained with parallel target sequencing and AmpliSeq. For changing the parameters you can look through this Python script""")
    parser.add_argument('--dir','-d', action="store", dest = "directory",
                        required=False, default="",
                        help = 'Directory with BAM files.')
    parser.add_argument('--bed','-b', action="store", dest = "bed",
                        required=True,
                        help = 'BED file')
    parser.add_argument('--mode','-m', action="store", dest = "mode",
                        default = "normal", required=False,
                        help = 'Mode (hard, normal, soft)')
    parser.add_argument('--id','-id', action="store", dest = "id_of_run",
                        default = "result", required=False,
                        help = 'Id of run (to output results)')
    parser.add_argument('--control','-cs', action="store", dest = "control_dataset",
                        default = "", required=False,
                        help = 'Id of control dataset (to perform quality control and to merge them)')
    parser.add_argument('--qc','-qc', action="store", dest = "QCpercent",
                        default = "0.8", required=False,
                        help = 'QC percentage')
    parser.add_argument('--mincor','-mc', action="store", dest = "minimum_correlation",
                        default = "0.7", required=False,
                        help = 'Minimum correlation threshold')
    parser.add_argument('--out','-o', action="store", dest = "output_folder",
                        default = "", required=False,
                        help = 'Folder to output files')
    parser.add_argument('--cleanControl','-cc', action="store_true", dest = "control_dataset_free_of_CNVs",
                        default = False, required=False,
                        help = 'Control dataset is free of CNVs')
    args = parser.parse_args()
    bam_folder = args.directory
    bed_file = args.bed
    bam_folder = os.path.abspath(bam_folder)
    bed_file = os.path.abspath(bed_file)
    percent = args.QCpercent
    min_cor = args.minimum_correlation
    output_folder = args.output_folder
    clean_control = args.control_dataset_free_of_CNVs


    mode = args.mode
    id_of_run = args.id_of_run
    if id_of_run == "result":
        id_of_run += bam_folder
    control_dataset = id_of_run
    if args.control_dataset:
        control_dataset = args.control_dataset
    list_of_chromosomes = ("chr7", "chr12", "chr9")

    logger.info((" ").join(["Task with id", id_of_run, "started"]))
    quant_dict = {}
    quant_0995, quant_099, quant_098, quant_095 = get_quantiles()
    if mode == "hard":
        quant_dict_del = quant_099
    elif mode == "normal":
        quant_dict_del = quant_098
    elif mode == "soft":
        quant_dict_del = quant_095
    quant_dict_dup = quant_095

    if not os.path.exists("tmp_chimeras.txt"):
        f = open('tmp_chimeras.txt', 'w')
        f.close()
    if not os.path.exists("tmp_output.txt"):
        f = open('tmp_output.txt', 'w')
        f.close()
    if not os.path.exists("tmp_references.txt"):
        f = open('tmp_references.txt', 'w')
        f.close()

    if args.directory:
        get_coverages(bam_folder, bed_file, id_of_run)
    output_file_name = id_of_run
    if not id_of_run == control_dataset:
        output_file_name = id_of_run + "_" + control_dataset
        if float(percent) < 0.9:
            percent = "0.9"
            logger.info("Too low quality control for merging test and control datasets. Value at least 0.9 should be used.")

    string_to_cmd = ("").join(["python2 qc.py ", bed_file, " ./result/", id_of_run, ".xls ./result/", control_dataset, ".xls ",  output_file_name, " ", percent])
    os.system(string_to_cmd)

    get_tmp_del_files(bed_file, quant_dict_del, quant_dict_dup, output_file_name, min_cor, clean_control)
    get_results(output_file_name, output_folder)
    string_to_cmd = (" ").join(["./visualisation.R", bed_file, "./visualisation", id_of_run])
    os.system(string_to_cmd)
    logger.info(" ".join(["Task with id", id_of_run, "finished! ;-)"]))


if __name__ == "__main__":
    main()



