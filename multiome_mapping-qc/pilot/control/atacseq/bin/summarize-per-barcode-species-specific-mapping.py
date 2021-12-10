#!/usr/bin/env python
# coding: utf-8

import os
import sys
import bam_utils
import ChunkedFileReader
import argparse
import logging

parser = argparse.ArgumentParser()
parser.add_argument('file_in_1', help = 'For genome 1, list of read, barcode, flag')
parser.add_argument('file_in_2', help = 'For genome 2, list of read, barcode, flag')
parser.add_argument('genome_1', help = 'Name of genome 1')
parser.add_argument('genome_2', help = 'Name of genome 2')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

FILE_1 = args.file_in_1
FILE_2 = args.file_in_2

mapping_outcomes = dict() # barcode --> [total, mapped_to_species_1, mapped_to_species_2, mapped_only_to_species_1, mapped_only_to_species_2]
reads_processed = 0


def is_primary_alignment_of_read_1(flag):
    x = bam_utils.interpret_flag(flag)
    return x['is_first'] and x['is_paired'] and not x['is_secondary'] and not x['is_supplementary']


with open(FILE_1, 'r') as fh_species_1:
    with open(FILE_2, 'r') as fh_species_2:
        s1_file_reader = ChunkedFileReader.ChunkedFileReader(fh_species_1, 0)
        s2_file_reader = ChunkedFileReader.ChunkedFileReader(fh_species_2, 0)
        for species_1 in s1_file_reader:
            species_2 = next(s2_file_reader)
            reads_processed += 1
            if reads_processed % 1000000 == 0:
                logging.info(f'Processed {reads_processed} reads')
            # verify the read name and barcodes are the same
            species_1_alignments = [x.split('\t') for x in species_1]
            species_2_alignments = [x.split('\t') for x in species_2]
            assert(species_1_alignments[0][0] == species_2_alignments[0][0])
            assert(species_1_alignments[0][1] == species_2_alignments[0][1])
            barcode = species_1_alignments[0][1]
            # get primary alignment of read 1 for species 1
            species_1_flags = [int(i[2]) for i in species_1_alignments]
            species_1_primary_alignment = [i for i in species_1_flags if is_primary_alignment_of_read_1(i)]
            assert(len(species_1_primary_alignment) == 1)
            species_1_primary_alignment = species_1_primary_alignment[0]
            # get primary alignment of read 1 for species 2
            species_2_flags = [int(i[2]) for i in species_2_alignments]
            species_2_primary_alignment = [i for i in species_2_flags if is_primary_alignment_of_read_1(i)]
            assert(len(species_2_primary_alignment) == 1)
            species_2_primary_alignment = species_2_primary_alignment[0]
            # compare outcomes
            species_1_primary_alignment_interpretation = bam_utils.interpret_flag(species_1_primary_alignment)
            species_2_primary_alignment_interpretation = bam_utils.interpret_flag(species_2_primary_alignment)
            species_1_is_mapped_in_proper_pair = species_1_primary_alignment_interpretation['is_proper_pair'] and not species_1_primary_alignment_interpretation['is_unmapped'] and not species_1_primary_alignment_interpretation['mate_unmapped']
            species_2_is_mapped_in_proper_pair = species_2_primary_alignment_interpretation['is_proper_pair'] and not species_2_primary_alignment_interpretation['is_unmapped'] and not species_2_primary_alignment_interpretation['mate_unmapped']
            if barcode not in mapping_outcomes:
                mapping_outcomes[barcode] = [0, 0, 0, 0, 0]
            mapping_outcomes[barcode][0] += 1
            if species_1_is_mapped_in_proper_pair:
                mapping_outcomes[barcode][1] += 1
            if species_2_is_mapped_in_proper_pair:
                mapping_outcomes[barcode][2] += 1
            if species_1_is_mapped_in_proper_pair and not species_2_is_mapped_in_proper_pair:
                mapping_outcomes[barcode][3] += 1
            if species_2_is_mapped_in_proper_pair and not species_1_is_mapped_in_proper_pair:
                mapping_outcomes[barcode][4] += 1


print('\t'.join(['barcode', 'total_read_pairs', f'number_mapped_to_{args.genome_1}', f'number_mapped_to_{args.genome_2}', f'number_mapped_only_to_{args.genome_1}', f'number_mapped_only_to_{args.genome_2}']))

for barcode, counts in mapping_outcomes.items():
    print('\t'.join([barcode] + [str(i) for i in counts]))

