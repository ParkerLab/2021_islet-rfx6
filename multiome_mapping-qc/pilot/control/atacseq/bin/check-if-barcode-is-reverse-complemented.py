#!/usr/bin/python

import os
import sys
import gzip
import logging
import argparse

parser = argparse.ArgumentParser(description='Compare nuclear barcodes to the barcode whitelist to infer whether they need to be reverse complemented', add_help=True)
parser.add_argument('--check-first', default=1000000, type=int, help='Check no more than this number of barcodes (default: 1000000)')
parser.add_argument('fastq', help='Fastq file containing nuclear barcodes (gzipped)')
parser.add_argument('whitelist', help='Nuclear barcode whitelist')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

COMPLEMENTS = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'N': 'N'
}


def reverse_complement(s):
    return ''.join([COMPLEMENTS[x.upper()] for x in s][::-1])

whitelist_barcodes = set()
whitelist_barcodes_rc = set()

logging.info('Reading whitelist barcodes and determining their reverse complements')

with open(args.whitelist, 'r') as f:
    for line in f:
        line = line.rstrip()
        #whitelist_barcodes.add(line)
        whitelist_barcodes.add(reverse_complement(reverse_complement(line)))
        whitelist_barcodes_rc.add(reverse_complement(line))

logging.info('Finished reading barcodes')
logging.info('Processing fastq file {}'.format(args.fastq))

line_count = 0
record_count = 0
in_whitelist_count = 0
in_whitelist_rc_count = 0
with gzip.open(args.fastq, 'rt') as f:
    for line in f:
        line_count += 1
        if (line_count - 2) % 4 != 0: # skip everything but the sequence itself
            continue
        record_count += 1
        line = line.rstrip()
        if line in whitelist_barcodes:
            in_whitelist_count += 1
        if line in whitelist_barcodes_rc:
            in_whitelist_rc_count += 1
        if record_count == args.check_first:
            break
        if record_count % 1000000 == 0:
            percent_in_whitelist = round(100*float(in_whitelist_count) / record_count, 2)
            percent_in_whitelist_rc = round(100*float(in_whitelist_rc_count) / record_count, 2)
            logging.info('Processed {} records so far. {} ({}%) in whitelist, {} ({}%) in whitelist reverse complemented'.format(record_count, in_whitelist_count,  percent_in_whitelist, in_whitelist_rc_count, percent_in_whitelist_rc))

logging.info('Finished processing barcodes.')
logging.info('{} of {} barcodes ({}%) were in the whitelist'.format(in_whitelist_count, record_count, round(100*in_whitelist_count/record_count, 2)))
logging.info('{} of {} barcodes ({}%) were in the reverse-complemented whitelist'.format(in_whitelist_rc_count, record_count, round(100*in_whitelist_rc_count/record_count, 2)))
decision = 'Do not reverse complement' if in_whitelist_count >= in_whitelist_rc_count else 'Reverse complement'
print(decision)
