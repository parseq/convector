#!/usr/bin/env python2

import argparse
import os
import sys
from collections import defaultdict
from math import log10, floor
import time
import logging
from logging.config import fileConfig

loginipath = ('logging_config.ini')
fileConfig(loginipath, defaults={'logfilename': 'pipeline.log'})
logger = logging.getLogger('CONVector_logger')

class Amplicon:
    """
    Class for amplicon. Contains method for determining intersecion between read and amplicon.
    """
    def __init__(self, string):
        """
        :param string: string from .bed file, chr \t start coord \t end coord \t name of amplicon \t exon \t comments
        :return: object of Amplicon class
        """
        info = string.split("\t")
        assert info[1].isdigit() and info[2].isdigit(), ('''Check format of .bed file.
                                                         Coordinates of amplicons in panel
                                                         are not int.''')
        assert len(info) >= 4, ('''Check format of .bed. Some amplicons written in wrong format.''')
        self.chromosome = info[0]
        self.start_pos = int(info[1]) + 1
        self.end_pos = int(info[2]) + 1
        self.ID = info[3]
        self.length = self.end_pos - self.start_pos

    def __repr__(self):
        return " ".join(["<<<", self.ID, self.chromosome, "start:",
                         str(self.start_pos),"end: ", str(self.end_pos), ">>>"])

    def __str__(self):
        return " ".join(["<<<", self.ID, self.chromosome, "start:",
                         str(self.start_pos), "end:", str(self.end_pos), "pool",
                         str(self.pool), "length", str(self.length), "gene symbol" ,self.gene_symbol, ">>>"])

    def ampl_intersection(self, other_ampl):
        start_ampl1 = self.start_pos
        end_ampl1 = self.end_pos
        start_ampl2 = other_ampl.start_pos
        end_ampl2 = other_ampl.end_pos
        if (end_ampl1 <= start_ampl2 or start_ampl1 >= end_ampl2):
            return 0
        return max(0, max(start_ampl1, start_ampl2) - min(end_ampl1, end_ampl2))


    def intersection(self, tup, perc):
        """
        finds intersection between read and amplicon (in bp or in percents)
        params: tuple (start of read - end of read), percentage mode (True or False)
        trivial algo O(n^2), n - len( max (set_of_reads, set_of_amplicons) )
        """
        start_read = tup[1]
        end_read = tup[2]
        start_ampl = self.start_pos
        end_ampl = self.end_pos

        if (end_read <= start_ampl or start_read >= end_ampl):
            return 0

        first_diff = start_ampl - start_read
        second_diff = end_ampl - end_read
        common_part = 0 # relatively to amplicon

        if first_diff >= 0 and second_diff >= 0:
            common_part = end_read - start_ampl

        elif first_diff >= 0 and second_diff <= 0:
            common_part = end_ampl - start_ampl

        elif first_diff <= 0 and second_diff >= 0:
            common_part = end_read - start_read

        elif first_diff <= 0 and second_diff <= 0:
            common_part = end_ampl - start_read

        if perc:
            return 100 * float(common_part) / (end_ampl - start_ampl)
        else:
            return common_part


def parse_bed_file(directory, filename):
    """

    :param directory: path to directory (absolute or relative)
    :param filename: name of file with .bed file
    :return: panel of amplicons and counter of .bed files in folder with .bam files
    """
    bed_file = ""
    path_to_bed = ""
    if filename == "":
        counter_of_beds = 0
        directory = os.path.abspath(directory)
        for filename in os.listdir(directory):
            if filename.endswith((".bed",".BED")):
                path_to_bed = os.path.join(directory, filename)
                counter_of_beds += 1
                if counter_of_beds > 1:
                    logger.warn("Too much BED files. chimeric_solver call is ambiguous.")
                    return None, 0
        if counter_of_beds == 0:
            logger.warn("There is no BED file in directory.")
            return None, 0
    else:
        if not filename.endswith((".bed",".BED")):
            logger.warn(bed_file)
            logger.warn("Your .bed file should have .bed file extenstion.")
            return None, 0
        else:
            counter_of_beds = 1
            path_to_bed = os.path.abspath(filename)

    panel_of_amplicons = defaultdict(list)
    with open(path_to_bed) as bed_f:
        bed_f.readline()
        for line in bed_f:
            tmp_amplicon = Amplicon(line)
            panel_of_amplicons[tmp_amplicon.chromosome].append(tmp_amplicon)

    return panel_of_amplicons, counter_of_beds



def additional_info(panel_of_amplicons, add_file):
    """
    :param panel_of_amplicons: dict of amplicons from .bed file
    :param add_file: file with the information about pools
    :return: panel of amplicons with additional fields added to each amplicon, type of pool
    """
    with open(add_file) as f:
        array_of_string = f.read().split("\r")
        for tmp_line in array_of_string:
                tmp_list = tmp_line.split()
                if len(tmp_list) <= 6:
                    continue
                else:
                    elem.pool = tmp_list[0][-3]
                    elem.fwd_primer = tmp_list[4]
                    elem.rev_primer = tmp_list[5]
                    elem.reference = ""
                    if tmp_list[8] in ("exon","intron"):
                        """
                        two types of markings in pools list
                        """
                        chromosome = tmp_list[11]
                        for elem in panel_of_amplicons[chromosome]:
                            if elem.ID == tmp_list[3]:
                                elem.length = int(tmp_list[15]) - int(tmp_list[12])
                                elem.gene_symbol = tmp_list[7] + tmp_list[8] + tmp_list[9]
                    else:
                        chromosome = tmp_list[9]
                        for elem in panel_of_amplicons[chromosome]:
                            if elem.ID == tmp_list[3]:
                                elem.length = int(tmp_list[13]) - int(tmp_list[10])
                                elem.gene_symbol = tmp_list[7]

    return panel_of_amplicons


def bam_to_sam(directory, outputdir):
    """
    :param directory: path to directory with .bam files
    :param outputdir: output directory for .sam files
    :return:
    """
    if not (os.path.exists(outputdir)):
        logger.warn("Making directory... " + outputdir)
        os.makedirs(outputdir)
    else:
        for the_file in os.listdir(outputdir):
            file_path = os.path.join(outputdir, the_file)
            try:
                if os.path.isfile(file_path):
                    os.unlink(file_path)
            except Exception as e:
                logger.warn(str(e))
    number_of_bam_files = 0
    for filename in os.listdir(directory):
        if filename.endswith((".bam",".BAM")) and not filename.startswith("nomatch"):
            number_of_bam_files += 1

    if number_of_bam_files == 0:
        logger.warn("There is no bam file in directory. Stopped.")
        sys.exit()

    counter_of_progress = 0
    for filename in os.listdir(directory):
        tmp_string = ""
        if filename.endswith((".bam",".BAM"))  and not filename.startswith("nomatch"):
            out_file_name = filename[:-4] + ".sam"
            tmp_string = ("samtools view -h " + os.path.join(directory,
                          filename) + " > " +
                          os.path.join(outputdir, out_file_name))
            os.system(tmp_string)
            counter_of_progress += 1
            logger.info(" ".join(["Progress in converting to SAM:", str(counter_of_progress),
                   "files of", str(number_of_bam_files), "total"]))


def process_cigar(cigar):
    """
    :param cigar: CIGAR-string
    :return: list with [numer of bps, CIGAR-symbol, ..., number of bps, CIGAR-symbol]
    """
    cigar_list = []
    prev_num = ""
    option = ""
    for ch in cigar:
        if (ch.isdigit()):
            if option:
                cigar_list.append(option)
            option = ""
            prev_num += ch
        if (ch.isalpha()):
            if prev_num:
                cigar_list.append(int(prev_num))
            prev_num = ""
            option += ch
    cigar_list.append(option)
    return cigar_list

counter_glob = 0

def extract_info_of_read(tmp_list, clip_cutoff, canonical_strings_for_each_amplicon, score_mq):
    """
    this function takes splitted read (just splitted read string)
    and return tuples of values:
    (CHR, START, END, %GC, [CLIPPED_START_PART, CLIPPED_END_PART])
    """
    flag = tmp_list[1]
    chromosome = tmp_list[2]
    start_pos = int(tmp_list[3])
    cigar_string = process_cigar(tmp_list[5])
    string = tmp_list[9]
    clipped_parts = []

    true_length = 0
    """
    types of clipping (hard and soft)
    """
    clipping_letters = ("S", "H")

    first_number = cigar_string[0]
    last_number = cigar_string[-2]

    first_option = cigar_string[1]
    last_option = cigar_string[-1]

    number_of_clipping_letters = 0

    number_of_matches = 0

    for i in xrange(len(cigar_string) - 1):
        if type(cigar_string[i]) == int:
            option = cigar_string[i + 1]
            if option in ("M","D","X","="):
                true_length += cigar_string[i]
            elif option in clipping_letters:
                number_of_clipping_letters += cigar_string[i]
                true_length += cigar_string[i]
            if option == "M":
                number_of_matches += cigar_string[i]

    not_chimera_flag = True

    if first_option in clipping_letters:
        start_pos += first_number
        true_length -= first_number
        if first_number > clip_cutoff:
            clipped_parts.append(string[:first_number])
            not_chimera_flag = False
    if last_option in clipping_letters:
        true_length -= last_number
        if last_number > clip_cutoff:
            clipped_parts.append((string[len(string) - last_number - 1:]))
            not_chimera_flag = False

    if not_chimera_flag and int(tmp_list[4]) < score_mq:
        return []

    if first_option in clipping_letters:
        processed_string = string[first_number:first_number + true_length]
    else:
        processed_string = string[:true_length]

    GC = 0.0
    tot_len = len(tmp_list[9])
    for i in tmp_list[9]:
        if i in ("G","C","g","c"):
            GC += 1
    GC /= tot_len

    end_pos = start_pos + true_length
    result_tuple = (chromosome, start_pos, end_pos, GC, true_length, clipped_parts)

    # canonical string, length of intersection, number of Ms, abs(start_read - start_ampl) + abs(end_read - end_ampl)
    alphabet = ("A","C","G","T","-","I")
    if chromosome in canonical_strings_for_each_amplicon and number_of_clipping_letters == 0 and int(tmp_list[4]) >= 60:
        for ampl in canonical_strings_for_each_amplicon[chromosome]:
            intersection_length = ampl.intersection(result_tuple, False)
            if intersection_length >= ampl.length and true_length <= ampl.length + 60:
                left_primer_length = ampl.start_pos - start_pos
                position = 60 - left_primer_length
                counter_of_steps = 0
                counter_of_seq_steps = 0
                for i in xrange(len(cigar_string) - 1):
                    if type(cigar_string[i]) == int:
                        counter_of_interval = 0
                        length_of_interval = cigar_string[i]
                        option = cigar_string[i + 1]
                        if option in ("M", "X", "="):
                            while counter_of_interval < length_of_interval:
                                canonical_strings_for_each_amplicon[chromosome][ampl][string[counter_of_seq_steps]][position + counter_of_steps] += 1
                                counter_of_steps += 1
                                counter_of_seq_steps += 1
                                counter_of_interval += 1
                        elif option == "D":
                            while counter_of_interval < length_of_interval:
                                canonical_strings_for_each_amplicon[chromosome][ampl]["-"][position + counter_of_steps] += 1
                                counter_of_steps += 1
                                counter_of_interval += 1
                        elif option == "I":
                            while counter_of_interval < length_of_interval:
                                canonical_strings_for_each_amplicon[chromosome][ampl]["I"][position + counter_of_steps] += 1
                                counter_of_interval += 1
                                counter_of_seq_steps += 1
                        else:
                            logger.warn(" ".joind["Unknown option in CIGAR:", option])
    return result_tuple


def if_len_of_list_intersected_equal_two(list_of_potentially_intersected, intersection, coverage_of_amplicons):
    """
    :param list_of_potentially_intersected: list of amplicons that intersects
    :param intersection: dict, {amplicon name : intersection length}
    :param coverage_of_amplicons: dict, {name of ampl : coverage}
    :return:
    """
    first = list_of_potentially_intersected[0]
    second = list_of_potentially_intersected[1]
    first_intersection = intersection[first]
    second_intersection = intersection[second]
    if first_intersection > second_intersection:
        coverage_of_amplicons[first] += 1
    else:
        coverage_of_amplicons[second] += 1

def get_coverage(reads_after_final_processing, cutoff,
                 panel_of_amplicons, percentage_mode_on):
    """
    :param reads_after_final_processing: reads withou chimeras
    :param cutoff: length for intersection
    :param panel_of_amplicons: panel of amplicons from .bed file
    :param percentage_mode_on: intersection in percents of in base pairs
    :return: calculated coverages of amplicons
    """
    coverage_of_amplicons = defaultdict(int)
    total_num_of_chimeras_of_first_type = 0
    for key, item in reads_after_final_processing.iteritems():
        if key in panel_of_amplicons:
            for read in item:
                list_of_potentially_intersected = []
                intersection = {}
                for amplicon in panel_of_amplicons[key]:
                   if not read[0] == amplicon.chromosome:
                       continue
                   inter = amplicon.intersection(read, percentage_mode_on)
                   if (inter > cutoff):
                       list_of_potentially_intersected.append(amplicon)
                       intersection[amplicon] = inter
                if len(list_of_potentially_intersected) == 1:
                    amplicon = list_of_potentially_intersected[0]
                    coverage_of_amplicons[amplicon] += 1
                elif len(list_of_potentially_intersected) == 2:
                    if_len_of_list_intersected_equal_two(list_of_potentially_intersected, intersection, coverage_of_amplicons)
                elif len(list_of_potentially_intersected) == 3:
                            left_coord = float('Inf')
                            right_coord = 0
                            right_ampl = None
                            left_ampl = None
                            left_and_right = {}
                            center_amplicon = None
                            for elem in list_of_potentially_intersected:
                                if elem.start_pos < left_coord:
                                    left_coord = elem.start_pos
                                    left_ampl = elem
                                    left_and_right["left"] = elem
                                if elem.end_pos > right_coord:
                                    right_coord = elem.end_pos
                                    right_ampl = elem
                                    left_and_right["right"] = elem
                            for elem in list_of_potentially_intersected:
                                if elem == left_and_right["left"] or elem == left_and_right["right"]:
                                    continue
                                else:
                                    center_amplicon = elem
                            left_interstection = center_amplicon.ampl_intersection(left_ampl)
                            right_interstection = center_amplicon.ampl_intersection(right_ampl)

                            if left_interstection + 30 > intersection[left_ampl]:
                                list_of_potentially_intersected.remove(left_ampl)
                            if right_interstection + 30 > intersection[right_ampl]:
                                list_of_potentially_intersected.remove(right_ampl)
                                if len(list_of_potentially_intersected) == 1:
                                    coverage_of_amplicons[list_of_potentially_intersected[0]] += 1
                                else:
                                    if_len_of_list_intersected_equal_two(list_of_potentially_intersected, intersection, coverage_of_amplicons)
                            if len(list_of_potentially_intersected) == 3:
                                coverage_of_amplicons[left_ampl] += 1
                                coverage_of_amplicons[right_ampl] += 1
                                total_num_of_chimeras_of_first_type += 1
                elif len(list_of_potentially_intersected) == 4:
                    logger.warn(" ".joind("We can not decide the right mapping of read ", str(list_of_potentially_intersected)))
                    logger.ward(read)
    logger.info(" ".join(["TOTAL AMOUNT OF I TYPE CHIMERAS:", str(total_num_of_chimeras_of_first_type)]))
    return coverage_of_amplicons


def rev_compl(string):
    """
    :param string: DNA containing 4 nucleotides.
    :return: reverse compimentary DNA
    """
    res_string = ""
    for i in xrange(len(string) - 1, -1, -1):
        letter = string[i]
        if letter == "G":
            res_string += "C"
        elif letter == "C":
            res_string += "G"
        elif letter == "A":
            res_string += "T"
        elif letter == "T":
            res_string += "A"
        else:
            res_string += letter
    return res_string

def remove_homopolymers(string, max_homopolymer_length):
    """
    :param string: DNA sequence
    :param max_homopolymer_length: maximum allowed length of homopolymer
    :return: sequence without homopolymers
    """
    result_string = ""
    homo_symbol = "$"
    homo_len = 1
    for i in xrange(len(string)):
        if homo_symbol == string[i] and homo_len < max_homopolymer_length:
            homo_len += 1
            result_string += string[i]
        elif not homo_symbol == string[i]:
            homo_symbol = string[i]
            result_string += string[i]
            homo_len = 1
    return result_string

def calculate_corrected_reads(outputdir, panel_of_amplicons, min_length, mq,
                              percentage_mode_on, clip_cutoff, cutoff, num_of_reads,
                              res_file):
    """
    We calculate the coverages of each amplicon with our metrics of quality.
    Percentage mode / BP mode - different metrics of intersection.
    Output: file with samples and their coverages in each amplicon.
    """
    string_to_cmd = "javac ./parseq/chimeric_solver/Main.java"
    os.system(string_to_cmd)
    if mq > 0:
        score_of_mq = floor(-10 * log10(mq)) # prob -> score
    else:
        score_of_mq = 0

    dict_to_know_what_amplicons_are_more_chimeric = defaultdict(int)

    chromosomes_in_bed = list(panel_of_amplicons.iterkeys())

    string_to_output = {}
    top_string = "Gene\tTarget\t"
    for key in panel_of_amplicons:
        for ampl in panel_of_amplicons[key]:
            string_to_output[ampl] = "N/A\t"
            string_to_output[ampl] += ampl.ID + "\t"


    if not os.path.exists(outputdir):
        os.makedirs(outputdir)


    for filename in os.listdir(outputdir):
        logger.info(filename)
        max_length = 0
        starts_end_ends = defaultdict(list)

        if not filename.endswith((".sam",".SAM")):
            continue

        total_reads_in_file = 0
        after_processing = 0

        reads_after_final_processing = defaultdict(list)

        with open(os.path.join(outputdir, filename)) as sam_to_correct:
            alphabet = ("A","C","G","T","-","I")
            nucleotides_alphabet = ("A","C","G","T")
            canonical_strings_for_each_amplicon = {}
            canonical_strings_for_each_amplicon_with_primers = {}
            for key in panel_of_amplicons:
                canonical_strings_for_each_amplicon[key] = {}
                canonical_strings_for_each_amplicon_with_primers[key] = {}
                for ampl in panel_of_amplicons[key]:
                    len_of_consensus = ampl.length + 120
                    # canonical string, length of intersection, number of Ms, abs(start_read - start_ampl) + abs(end_read - end_ampl)
                    canonical_strings_for_each_amplicon[key][ampl] = {}
                    canonical_strings_for_each_amplicon_with_primers[key][ampl] = {}
                    for letter in alphabet:
                        canonical_strings_for_each_amplicon[key][ampl][letter] = [0 for i in xrange(len_of_consensus)]

            for tmp_line in sam_to_correct:
                clipping = 0

                if not tmp_line.startswith("@"):
                    clipping = 0
                    tmp_list = tmp_line.split()
                    total_reads_in_file += 1
                    if (int(tmp_list[4]) > 0 and
                        len(tmp_list[9]) >= min_length): # and
                        # tmp_list[2] in chromosomes_in_bed): # if throw away, then reads
                        processed_read = extract_info_of_read(tmp_list, int(clip_cutoff), canonical_strings_for_each_amplicon, score_of_mq)
                        if len(processed_read) > 1:
                            reads_after_final_processing[processed_read[0]].append(processed_read)
                            after_processing += 1
            for key in panel_of_amplicons:
                for ampl in panel_of_amplicons[key]:
                    heterozigouthy_state = False
                    canonical_read_with_primers = ""
                    summa = []
                    for i in xrange(len(canonical_strings_for_each_amplicon[key][ampl][letter])):
                        summa_in = 0
                        for letter in nucleotides_alphabet:
                            summa_in += canonical_strings_for_each_amplicon[key][ampl][letter][i]
                        summa.append(summa_in)
                    for i in xrange(len(summa)):
                        if summa[i] > 0:
                            max_amount = 0
                            canonical_letter = ""
                            lower_threshold_for_amigous_case = 0.4 * summa[i]
                            upper_threshold_for_amigous_case = 0.6 * summa[i]
                            for letter in alphabet:
                                if (canonical_strings_for_each_amplicon[key][ampl][letter][i] > lower_threshold_for_amigous_case and
                                    canonical_strings_for_each_amplicon[key][ampl][letter][i] < upper_threshold_for_amigous_case):
                                    canonical_letter = "*"
                                    break
                                if (canonical_strings_for_each_amplicon[key][ampl][letter][i] > max_amount and
                                not letter in ("-","I")):
                                    canonical_letter = letter
                                    max_amount = canonical_strings_for_each_amplicon[key][ampl][letter][i]
                            canonical_read_with_primers += canonical_letter
                    canonical_strings_for_each_amplicon_with_primers[key][ampl] = canonical_read_with_primers

            chimeras_list = []
            for chr in reads_after_final_processing:
                for read in reads_after_final_processing[chr]:
                    if not read[5] == []:
                        for elem in read[5]:
                            without_homopolymers = remove_homopolymers(elem, 4)
                            chimeras_list.append(without_homopolymers)


        if total_reads_in_file < num_of_reads:
            logger.warn('\n'.join(["WARNING",
                             "In file " + str(filename), "only " + str(total_reads_in_file) + " reads! And "
                             + str(num_of_reads) + " is minimum required. You can change" ,
                             "this option with the key -n.", "Stop analysis of this file.", ""]))
            continue

        size_of_list = len(chimeras_list)
        for i in xrange(size_of_list):
            chimera = chimeras_list[i]
            chimeras_list.append(rev_compl(chimera))


        # total reads
        coverage_of_amplicons = get_coverage(reads_after_final_processing, cutoff,
                                             panel_of_amplicons, percentage_mode_on)
        counter = 0
        list_of_ampls_to_search_chimera_in = []

        with open("tmp_chimeras.txt", "wb") as f:
            for chimera in chimeras_list:
                f.write(chimera + "\n")
        with open("tmp_references.txt", "wb") as f:
            for chr in canonical_strings_for_each_amplicon_with_primers:
                for ampl in canonical_strings_for_each_amplicon_with_primers[chr]:
                    if not canonical_strings_for_each_amplicon_with_primers[chr][ampl] == "":
                        reference = canonical_strings_for_each_amplicon_with_primers[chr][ampl]
                        f.write(remove_homopolymers(reference, 4) + "\n")
                        list_of_ampls_to_search_chimera_in.append(ampl)


        string_to_cmd = ("").join(["java parseq/chimeric_solver/Main" ])
        os.system(string_to_cmd)

        with open("tmp_output.txt", "r") as f:
            for ampl in list_of_ampls_to_search_chimera_in:
                chimeric_increase = 0
                try:
                    chimeric_increase = int(f.readline())
                except:
                    chimeric_increase = 0
                coverage_of_amplicons[ampl] += chimeric_increase
                dict_to_know_what_amplicons_are_more_chimeric[ampl] += chimeric_increase
                counter += chimeric_increase

        logger.info(" ".join(["TOTAL AMOUNT OF CHIMERAS OF II TYPE: ", str(len(chimeras_list) / 2)]))
        logger.info(" ".join(["TOTAL AMOUNT OF CHIMERAS OF II TYPE RE-ALIGNED: ", str(counter)]))


        top_string += (filename + "\t")
        for key in panel_of_amplicons:
            for ampl in panel_of_amplicons[key]:
                string_to_output[ampl] += str(coverage_of_amplicons[ampl])+ "\t"

    for key in dict_to_know_what_amplicons_are_more_chimeric:
        logger.info(" ".join(["AMPL", str(key.ID), str(key.chromosome), "HAS THIS AMOUNT OF CHIMERAS II TYPE", str(dict_to_know_what_amplicons_are_more_chimeric[key])]))

    if not (os.path.exists("result")):
            logger.warn("Making directory result")
            os.makedirs("result")
    newpath = os.path.join(os.path.abspath("result"), res_file)
    if "/" in res_file:
        newpath = os.path.abspath(res_file)
    with open(newpath, "wb") as output_file:
            output_file.write(top_string + "\n")
            for key in panel_of_amplicons.iterkeys():
                for amplicon in panel_of_amplicons[key]:
                    output_file.write(string_to_output[amplicon] + "\n")


def main():
    parser = argparse.ArgumentParser(description="""Convert everything from BAM to SAM.
Then count corrected_reads.
Needs input directory with bam files and file with coordinates of amplicons.""")
    parser.add_argument('--dir','-d', action="store", dest = "directory",
                        required=True,
                        help = 'Directory with 1 BED file and BAM files.')
    parser.add_argument('--bed','-b', action="store", dest = "bed",
                        default = "",
                        help = 'BED file (if not in directory with dataset)')
    parser.add_argument('--out','-o', action="store", dest = "outputdir",
                        default = "output_with_SAMs",
                        help = 'Output directory for all SAM files')
    parser.add_argument('--len','-l', action="store", dest="len_threshold",
                        type=int,
                        default="1")
    # in form "probability mapping position is wrong"
    parser.add_argument('--mq','-m', action="store", dest="MQ", type=float, default="0.9")
    parser.add_argument('--mode','-md', action="store_true", dest="mode", default=False)
    parser.add_argument('--cutoff','-c', action="store", dest="overlap", type=int, default=20)
    parser.add_argument('--convert','-conv',action="store_true",dest="converter", default=False)
    parser.add_argument('--clipping','-clip',action="store", dest="clip_amount", type=int, default=30)
    parser.add_argument('--numOfReads','-n',action="store", dest="min_number_of_reads", type=int, default=15000)
    parser.add_argument('--addInfo','-i',action="store", dest="additional_file", type=str, default=None)

    parser.add_argument('--resFile','-r',action="store", dest="result_file", type=str, default="result.xls")

    args = parser.parse_args()
    # .bam files which names starts with "nomatch" will not be considered!
    directory = args.directory
    outputdir = args.outputdir

    len_threshold = (args.len_threshold)
    mq = (args.MQ)
    cutoff = (args.overlap)
    clip_cutoff = (args.clip_amount)
    num_of_reads = (args.min_number_of_reads)
    bed_file = (args.bed)

    panel_of_amplicons, counter_of_bed = parse_bed_file(directory, bed_file)

    if not counter_of_bed == 1:
        logger.warn("The chimeric_solver can not be started because the .bed file was not found. Quit.")
        sys.exit()

    pool_info = ""

    if not args.additional_file == None:
        additional_info(panel_of_amplicons, args.additional_file)

    if args.converter:
        bam_to_sam(directory, outputdir)

    # two modes - percentage (coverage defined by the % of amplicon
    # covered with the read or
    # - base mode (number of bp for coverage)

    percentage_mode_on = args.mode

    res_file = args.result_file

    calculate_corrected_reads(outputdir, panel_of_amplicons,
                              len_threshold, mq, percentage_mode_on,
                              clip_cutoff, cutoff, num_of_reads,
                              res_file)


if __name__ == "__main__":
	main()
