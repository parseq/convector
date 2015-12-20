#!/usr/bin/env python2
__author__ = 'german'

import chimeric_solver
import sys
from collections import defaultdict
from math import log
import statistics
import os
import copy
import logging
from logging.config import fileConfig

loginipath = ('logging_config.ini')
fileConfig(loginipath, defaults={'logfilename': 'pipeline.log'})
logger = logging.getLogger('CONVector_logger')

def return_quantiles_chisq():
    """
    Reads files with quantiles of chi squared distribution, degrees of freedom: from 1 to 500 (depending on number
    of samples).
    :return: tuple of dictionaries, {degree of freedom: corresponding quantile}
    """
    quantiles99 = {}
    quantiles95 = {}
    with open("./quantiles/quantile_chisq_99.txt") as f:
        counter = 1
        for line in f:
            quantiles99[counter] = float(line.split()[0])
            counter += 1
    with open("./quantiles/quantile_chisq_95.txt") as f:
        counter = 1
        for line in f:
            quantiles95[counter] = float(line.split()[0])
            counter += 1
    return quantiles99, quantiles95

def parse_file_with_coverages(filename):
    """
    :param filename: name of file with coverages, .xls, format described in manual
    :return: dictionaries 'samples' {sample : {amplicon : coverages}};
             list with low covered amplicons;
             coverages of clean samples (without low covered amplicons, for statistical analysis);
             coverages of all samples (for output)
    """
    samples = defaultdict(dict)
    true_coverages_clean = defaultdict(dict)
    true_coverages = defaultdict(dict)
    total_number_of_samples = 0
    low_covered_ampls = []
    with open(filename) as f:
        array_of_lines = f.readlines()
        samples_names = []
    for i in xrange(len(array_of_lines)):
        line = array_of_lines[i]
        splitted_line = line.split()
        if not i:
            samples_names = splitted_line[2:]
            total_number_of_samples = len(samples_names)
        else:
            summa_of_coverages = 0
            ampl_name = splitted_line[1]
            for i in xrange(len(samples_names)):
                value = 0
                try:
                    value = (int(splitted_line[i + 2]))
                except ValueError:
                    break
                true_coverages_clean[samples_names[i]][ampl_name] = value
                true_coverages[samples_names[i]][ampl_name] = value
                # 10 - threshold for homozygous deletion
                if value > 10:
                    samples[samples_names[i]][ampl_name] = log(float(value))
                else:
                    samples[samples_names[i]][ampl_name] = 0.01
                    # 0.01 - just very low value, indicates a homo del
                summa_of_coverages += value
            if summa_of_coverages / total_number_of_samples < 50:
                for name in samples_names:
                    samples[name].pop(ampl_name, None)
                    true_coverages_clean[name].pop(ampl_name, None)
                # exclude low covered (in all samples) amplicons using threshold 50
                low_covered_ampls.append(ampl_name)
    return samples, low_covered_ampls, true_coverages_clean, true_coverages

def normalize_data_inside_chromosomes(clean_chromosomes_amplicons, samples):
    """
    :param clean_chromosomes_amplicons: coverages of amplicons inside one chromosome
    :param samples: list with names
    :return: normalized coverages
    """
    total_amount_sample_chromosome = defaultdict(dict)
    for sample, data in samples.iteritems():
        for chromosome, amplicons in clean_chromosomes_amplicons.iteritems():
            try:
                summa_of_coverages = [data[amplicon] for amplicon in amplicons]
            except KeyError:
                logger.warn("Smth went wrong...\nThe most probable thing - you are using the wrong file with coverages!")
                sys.exit(1)
            for amplicon in amplicons:
                data[amplicon] -= statistics.mean(summa_of_coverages)
            total_amount_sample_chromosome[sample][chromosome] = sum(summa_of_coverages)
    return total_amount_sample_chromosome

def form_ellipsoid(samples_to_train, amplicons_from_chromosome):
    """
    :param samples_to_train: dictionary {amplicon name : coverages in samples from training samples}
    :param amplicons_from_chromosome: coverages of samples from control dataset, list
    :return:
    """
    ellipsoid = {}
    for amplicon in amplicons_from_chromosome:
        amplicon_values = []
        for name, info in samples_to_train.iteritems():
            amplicon_values.append(info[amplicon])
        ellipsoid[amplicon] = (statistics.medianW(amplicon_values), statistics.sn_estimator(amplicon_values) ** 2)
    return ellipsoid

def diagnose_chromosome_ellipsoid(samples_to_test, ellipsoid, list_of_amplicons_to_test, qChisq, num_of_accepted):
    """

    :param samples_to_test: coverages in samples and amplicons, dict {sample name : {ampl name : coverage} }
    :param ellipsoid: pairs for each amplicon, (estimation of mean, estimation of standard deviation)
    :param list_of_amplicons_to_test: list of amplicons from one chromosome without low covered amplicons
    :param qChisq: corresponding chi square quantile
    :param num_of_accepted: number of amplicons that will be taken into account
    :return: dict of lists with normal or not coverages inside one chromosome ( 1 = normal, 0 = irregular);
             average robust residuals for calculation of ARV
    """
    normal_or_not = defaultdict(int)
    statistic_for_sample_and_chromosome = defaultdict(list)
    avtc_residuals_for_amplicons = defaultdict(list)

    for sample, coverages_of_amplicons in samples_to_test.iteritems():
        for ampl in list_of_amplicons_to_test:
            distance_to_mean = (coverages_of_amplicons[ampl] - ellipsoid[ampl][0])
            dist = (distance_to_mean ** 2) / (ellipsoid[ampl][1])
            statistic_for_sample_and_chromosome[sample].append(dist)
            avtc_residuals_for_amplicons[ampl].append(distance_to_mean ** 2)
    for sample, statistic_values in statistic_for_sample_and_chromosome.iteritems():
        if sum(sorted(statistic_values)[:num_of_accepted]) < qChisq:
            normal_or_not[sample] = 1
        else:
            normal_or_not[sample] = 0
    return normal_or_not, avtc_residuals_for_amplicons

def form_list_of_qc_negative(list_of_normal_or_not_dicts, samples_to_test_qc):
    """
    :param list_of_normal_or_not_dicts: dict {sample : list of 0 and 1, determining the irregularity of coverage inside
           chromomsome}
    :param samples_to_test_qc: list of sample in test dataset, names
    :return: list of samples that did not passed our QC algorithm
    """
    dict_of_sums = {}
    num_of_required_good_chromosomes = len(list_of_normal_or_not_dicts) - 1
    list_of_negatives = []
    for sample in samples_to_test_qc:
        dict_of_sums[sample] = 0
        for dictionary in list_of_normal_or_not_dicts:
            dict_of_sums[sample] += dictionary[sample]
    counter_of_negative = 0
    counter_of_positive = 0
    for sample in sorted(dict_of_sums.iterkeys()):
        if dict_of_sums[sample] < num_of_required_good_chromosomes:
            logger.info("QC Negative " + sample)
            counter_of_negative += 1
            list_of_negatives.append(sample)
        else:
            counter_of_positive += 1
            logger.info("QC Positive " + sample)
    logger.warn(" ".join(["Overall:", str(counter_of_negative), "samples was filtered out and", str(counter_of_positive), "were accepted"]))
    return list_of_negatives

def output_result_file(true_coverages_of_samples_to_test, true_coverages_of_samples_to_train, qc_negative_list, clean_chromosomes_amplicons, run_id, mode):
    """
    :param true_coverages_of_samples_to_test: coverages before normalization
    :param true_coverages_of_samples_to_train: coverages before normalization
    :param qc_negative_list: list of samples that did not passed QC
    :param clean_chromosomes_amplicons: amplicons from clean samples
    :param run_id: name of task
    :param mode: merging control and test sample or not
    :return: output file with coverages and QC report
    """
    ordered_list_of_samples_test = list(sorted(true_coverages_of_samples_to_test.iterkeys()))
    ordered_list_of_samples_control = list(sorted(true_coverages_of_samples_to_train.iterkeys()))
    with open("./result/" + run_id + "_qc.xls", "wb") as f:
        top_string = "Gene\tTarget\t"
        for sample in ordered_list_of_samples_test:
            if sample not in qc_negative_list:
                top_string += "Case_" + sample + "\t"
        if mode == 1:
            for sample in ordered_list_of_samples_control:
                top_string += "Control_" + sample + "\t"
        top_string += '\n'
        f.write(top_string)
        for chr, ampls in clean_chromosomes_amplicons.iteritems():
            for ampl in ampls:
                new_result_string = "N/A\t" + ampl + "\t"
                for sample in ordered_list_of_samples_test:
                    if sample not in qc_negative_list:
                        new_result_string += str(true_coverages_of_samples_to_test[sample][ampl]) + "\t"
                if mode == 1:
                    for sample in ordered_list_of_samples_control:
                        new_result_string += str(true_coverages_of_samples_to_train[sample][ampl]) + "\t"

                new_result_string += "\n"
                f.write(new_result_string)


def main():
    bed_file = sys.argv[1]
    file_with_coverages = sys.argv[2]
    file_with_training_samples = sys.argv[3]
    run_id = sys.argv[4]
    percent = float(sys.argv[5])
    directory = "./"
    panel_of_amplicons, counter_of_beds = chimeric_solver.parse_bed_file(directory, bed_file)
    samples_to_test, low_covered_ampls_test, clean_coverages_of_samples_to_test, true_coverages_of_samples_to_test = parse_file_with_coverages(file_with_coverages)
    samples_to_train, low_covered_ampls_train, clean_coverages_of_samples_to_train, true_coverages_of_samples_to_train = parse_file_with_coverages(file_with_training_samples)
    set_of_low_covered_ampls = set(low_covered_ampls_test + low_covered_ampls_train)
    mode = 0 # mode is equal to 1 if the Control (Train) and Test samples differ. Otherwise, mode = 0.
    if file_with_coverages != file_with_training_samples:
        mode = 1


    qc_by_chromosomes = defaultdict(int)
    clean_chromosomes_amplicons = defaultdict(list)
    all_amplicon_names = defaultdict(list)

    for key in panel_of_amplicons:
        for amplicon in panel_of_amplicons[key]:
            if amplicon.ID not in set_of_low_covered_ampls:
                clean_chromosomes_amplicons[key].append(amplicon.ID)
            all_amplicon_names[key].append(amplicon.ID)

    quantiles99, quantiles95 = return_quantiles_chisq()
    total_amount_sample_chromosome_test = normalize_data_inside_chromosomes(clean_chromosomes_amplicons, samples_to_test)
    total_amount_sample_chromosome_train = normalize_data_inside_chromosomes(clean_chromosomes_amplicons, samples_to_train)
    ellipsoids_for_chromosomes = {}
    ellipsoids_for_chromosomes_test = {}
    for chr, amplicons_from_chromosome in clean_chromosomes_amplicons.iteritems():
        if chr.startswith("chr"):
            ellipsoids_for_chromosomes[chr] = form_ellipsoid(samples_to_train, amplicons_from_chromosome)
            ellipsoids_for_chromosomes_test[chr] = form_ellipsoid(samples_to_test, amplicons_from_chromosome)

    list_of_normal_or_not_dicts = []
    list_of_robust_variances_control_against_control = []
    list_of_robust_variances_test_against_test = []

    for chr, ellipsoid_test in ellipsoids_for_chromosomes_test.iteritems():
        if chr in ("chr7", "chr12", "chr9"):
            list_of_amplicons_to_test = clean_chromosomes_amplicons[chr]
            for amplicon, element in ellipsoid_test.iteritems():
                if amplicon in list_of_amplicons_to_test:
                    list_of_robust_variances_test_against_test.append(element[1])

    for chr, ellipsoid in ellipsoids_for_chromosomes.iteritems():
        if chr in ("chr7", "chr12", "chr9"):
            list_of_amplicons_to_test = clean_chromosomes_amplicons[chr]
            for amplicon, element in ellipsoid.iteritems():
                if amplicon in list_of_amplicons_to_test:
                    list_of_robust_variances_control_against_control.append(element[1])
            num_of_accepted = int(len(list_of_amplicons_to_test) * percent)
            qChisq = quantiles99[(num_of_accepted)]
            normal_or_not, avtc_residuals_for_amplicons = diagnose_chromosome_ellipsoid(samples_to_test, ellipsoid, list_of_amplicons_to_test, qChisq, num_of_accepted)
            list_of_normal_or_not_dicts.append(normal_or_not)


    samples_to_test_qc = sorted(list(samples_to_test.iterkeys()))
    qc_negative_list = form_list_of_qc_negative(list_of_normal_or_not_dicts, samples_to_test_qc)


    avrcc = statistics.mean(list_of_robust_variances_control_against_control)
    avrtt = statistics.mean(list_of_robust_variances_test_against_test)

    with open("qc_control_log.txt","wb") as qc_report:
        arvc = (" ").join(["Average Robust Variance Of Control Dataset with filename", file_with_training_samples, ":", str(avrcc), "\n"])
        qc_report.write(arvc + "\n")
        arvt = (" ").join(["Average Robust Variance Of Test Dataset with filename", file_with_coverages, ":", str(avrtt), "\n"])
        qc_report.write(arvt + "\n")
        logger.info(" ".join(["ARVc =", str(arvc)]))
        logger.info(" ".join(["ARVt =", str(arvt)]))

        qc_report.write((" ").join(["Total amount of samples filtered out using Quality Control algorithm:", str(len(qc_negative_list)), "\n\n"]) )

        for negative_sample in qc_negative_list:
            qc_report.write(("").join([negative_sample, " \n - did not passed QC Control", "\n\n"]))

    output_result_file(true_coverages_of_samples_to_test, true_coverages_of_samples_to_train, qc_negative_list, all_amplicon_names, run_id, mode)


main()
