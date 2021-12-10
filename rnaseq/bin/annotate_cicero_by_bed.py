#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pybedtools as bt
import pandas as pd
import argparse

parser = argparse.ArgumentParser(
    description='Reduce a cicero co-accessibility file to a file containing only interactions in which at least one TSS is involved', add_help=True)
parser.add_argument('cicero', help='Cicero co-accessibility file')
parser.add_argument(
    'bed', help='BED file (TSS for example) (chrom, start, end, name, score, strand)')
parser.add_argument('chrom_sizes', help='Chrom size file')
parser.add_argument('-r', '--right', dest='right', type=int, default=0,
                    help='Extend intervals by distance (bp) towards right. Default: 0')

parser.add_argument('--min-score', dest='min_score', type=float, default=0,
                    help='Filter connections for the minimum score. Default: 0')

command_group = parser.add_mutually_exclusive_group()
command_group.add_argument('-b', '--both', dest='both', type=int, default=0,
                           help='Extend intervals by distance in bp in both directions. Default: 0')
command_group.add_argument('-l', '--left', dest='left', type=int, default=0,
                           help='Extend intervals by distance (bp) towards left. Default: 0')

args = parser.parse_args()

CICERO_IN = args.cicero
BED_FILE = args.bed
CHROM_SIZES = args.chrom_sizes

# first, determine all peaks present
KEEP_PEAKS = set()
ALL_PEAKS = set()

with open(CICERO_IN, 'r') as f:
    for line in f:
        first, second, score = line.rstrip().split()
        # if float(score) < args.min_score:
        #     continue
        ALL_PEAKS.add(first)
        ALL_PEAKS.add(second)

peaks = pd.DataFrame([i.split('_') for i in ALL_PEAKS],
                     columns=['chrom', 'start', 'end'])
peaks.start = peaks.start.astype(int)
peaks.end = peaks.end.astype(int)
peaks_bt = bt.BedTool().from_dataframe(peaks).sort()

bed_bt = bt.BedTool(BED_FILE).sort()
peak2annot = dict()  # chr_start_end --> [tss_1, tss_2, ...]

if args.both:
    print("INFO: --both specified; using that.", file=sys.stderr)
    intersection = peaks_bt.intersect(bed_bt.slop(
        b=args.both, g=CHROM_SIZES), wa=True, wb=True)
else:
    print("INFO: --left or --right specified; using that.", file=sys.stderr)
    intersection = peaks_bt.intersect(bed_bt.slop(
        l=args.left, r=args.right, s=True, g=CHROM_SIZES), wa=True, wb=True)

for i in intersection:
    peak_chrom, peak_start, peak_end, tss_chrom, tss_start, tss_end, gene, score, tss_strand = i
    peak = f'{peak_chrom}_{peak_start}_{peak_end}'
    if peak not in peak2annot:
        peak2annot[peak] = set()
    peak2annot[peak].add(gene)

print('\t'.join(['peak_1', 'peak_2', 'score', 'peak_1_name', 'peak_2_name']))
with open(CICERO_IN, 'r') as f:
    for line in f:
        peak_1, peak_2, score = line.strip().split()
        peak_1_tss = 'None' if peak_1 not in peak2annot else ','.join(
            set(peak2annot[peak_1]))
        peak_2_tss = 'None' if peak_2 not in peak2annot else ','.join(
            set(peak2annot[peak_2]))
        print('\t'.join([peak_1, peak_2, score, peak_1_tss, peak_2_tss]))
