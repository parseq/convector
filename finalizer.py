#!/usr/bin/env python2

import argparse
from collections import defaultdict
import os



def parser_results(file_with_results):
    """
    :param file_with_results: file with matrix of CNVs, tab delimited
    :return: tuple, list of names of amplicons, their order, list of samples that contain DEL marks, names of samples,
             names of exons.
    """
    amplicons = {}
    order_of_amplicons = defaultdict(list)
    sample_amplicon_coverage = defaultdict(list)
    dirty_samples = defaultdict(dict)
    exons = []
    with open(file_with_results) as f:
        content = f.read()
        # uncomment for usage with other format
        # lines = content.split("\r")
        lines = content.split("\n")
        samples = []
        for line in lines:
            if line.startswith("Chr"):
                tmp_line = line.replace("\t","*")
                tmp_line = tmp_line.split("*")
                samples = tmp_line[5:]
            if line.startswith("chr"):
                tmp_line = line.replace("\t","*")
                tmp_line = tmp_line.split("*")
                if not tmp_line[5] == "EX" and not tmp_line[5] == "NE":
                    amplicons[tmp_line[3]] = (tmp_line[0], tmp_line[1],
                                              tmp_line[2], tmp_line[4])
                    if tmp_line[4] not in exons:
                        exons.append(tmp_line[4])
                    order_of_amplicons[tmp_line[0]].append(tmp_line[3])
                    for i in range(5, len(samples) + 5):
                        sample_amplicon_coverage[samples[i - 5]].append(tmp_line[i])
                        if tmp_line[i] == "DEL" or tmp_line[i] == "HOMO DEL":# or tmp_line[i] == "^":
                            if tmp_line[4] not in dirty_samples[samples[i - 5]]:
                                dirty_samples[samples[i - 5]][tmp_line[4]] = [tmp_line[3]]
                            else:
                                dirty_samples[samples[i - 5]][tmp_line[4]].append(tmp_line[3])
    return amplicons, order_of_amplicons, dirty_samples, samples, exons
                
def determine_neighbors(exons):
    """
    :param exons: list of exons
    :return: dictionary that contains information about left and right neighboor of exon
    """
    neighboors = defaultdict(list)
    for i in range(len(exons)):
        current_exon = exons[i].split()
        left_neighboor = ["LEFT"]
        right_neighboor = ["RIGHT"]
        if i > 0:
            left_neighboor = exons[i-1].split()
        if i < len(exons) - 1:
            right_neighboor = exons[i+1].split()
        if left_neighboor[0] == current_exon[0]:
            neighboors[exons[i]].append(exons[i - 1])
        if right_neighboor[0] == current_exon[0]:
            neighboors[exons[i]].append(exons[i + 1])
    return neighboors

            
def generate_final_output(file_with_results, file_with_deletions):
    """
    :param file_with_results: name of file with result report
    :param file_with_deletions: input matrix with CNVs
    :return:
    """
    amplicons, order_of_amplicons, dirty_samples, samples, order_of_exons = parser_results(file_with_results)
    exons_before_length_marking = defaultdict(list)
    exons = defaultdict(list)
    for amplicon, lst in amplicons.iteritems():
        exons[lst[3]].append(amplicon)

        
    neighboors = determine_neighbors(order_of_exons)
    
    final_output = defaultdict(list)

    for suspicious, dirty_exons in dirty_samples.iteritems():
        for dirty_exon in dirty_exons:
            exon = exons[dirty_exon]
            if len(exon) == 1:
                if len(dirty_exons[dirty_exon]) == len(exon):
                    final_output[suspicious].append(str(dirty_exon))
            elif len(exon) == 2:
                if len(dirty_exons[dirty_exon]) == len(exon):
                    final_output[suspicious].append(dirty_exon)
            else:
                if len(dirty_exons[dirty_exon]) >= ((len(exon) + 1)// 2):
                    final_output[suspicious].append(dirty_exon)
    changes = True

    while changes:
        changes = False
        for suspicious, dirty_exons in final_output.iteritems():
            for dirty_exon in dirty_exons:
                    for elem in neighboors[dirty_exon]:
                        if elem not in final_output[suspicious] and elem in dirty_samples[suspicious]:
                            if len(dirty_samples[suspicious][elem]) >= ((len(exons[elem]) + 1)// 2):
                                final_output[suspicious].append(elem)
                                changes = True
                            
    # determine neighboring with less strict criteria
    with open(os.path.join(file_with_deletions), "wb") as f:
        for elem in sorted(final_output.iterkeys()):
            tmp_string = str(elem) + "\t"
            for ex in order_of_exons:
                if ex in final_output[elem]:
                        if len(exons[ex]) == 1:
                            tmp_string += ex + " : SH\t"
                        else:
                            tmp_string += ex + "\t"
            tmp_string += "\n"
            f.write(tmp_string)
               




def main():
    parser = argparse.ArgumentParser(description="""Output result of detection of deletions.""")
    parser.add_argument('--input','-i',action="store", dest = "input_file",
                        help=".xls file with results of CNVdetector",
                        required=True)
    parser.add_argument('--output','-o',action="store", dest = "output_file",
                        help="File to output results",
                        required=True)

    parser.add_argument('--folder','-f',action="store", dest = "output_folder",
                        help="Folder to output results",
                        required=False, default=None)

    args = parser.parse_args()

    output_folder = args.output_folder
    generate_final_output(args.input_file, args.output_file)

main()
