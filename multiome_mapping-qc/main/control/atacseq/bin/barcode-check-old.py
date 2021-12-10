#!/usr/bin/env python
# coding: utf-8

import os
import sys
import gzip
import logging
import argparse
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd

parser = argparse.ArgumentParser(description='Compare nuclear barcodes to the barcode whitelist to infer whether they need to be reverse complemented', add_help=True)
parser.add_argument('--check-first', default=1000000, type=int, help='Check no more than this number of barcodes (default: 1000000)')
parser.add_argument('--prefix', default='barcode-check.', help='Prefix for output files (default: barcode-check.)')
parser.add_argument('fastq', help='Fastq file containing nuclear barcodes (gzipped)')
parser.add_argument('whitelist', help='Nuclear barcode whitelist')
#args = parser.parse_args(['/lab/work/porchard/Nova-303/data/fastq/1846_ATAC_b.2.fastq.gz', '/home/porchard/github/snATACseq-NextFlow/737K-arc-v1.txt'])
args = parser.parse_args()

PREFIX = args.prefix

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

COMPLEMENTS = {
    'A': 'T',
    'C': 'G',
    'G': 'C',
    'T': 'A',
    'N': 'N'
}


def reverse_complement(s):
    return ''.join([COMPLEMENTS[x.upper()] for x in s][::-1])


def complement(s):
    return ''.join([COMPLEMENTS[x.upper()] for x in s])


def reverse(s):
    return ''.join(x.upper() for x in s[::-1])


whitelist_barcodes = set()
whitelist_barcodes_rc = set()
whitelist_barcodes_r = set()
whitelist_barcodes_c = set()
logging.info('Reading whitelist barcodes and determining their reverse complements')

#count = 0
with open(args.whitelist, 'r') as f:
    for line in f:
        line = line.rstrip()
#        count += 1
#        if count == 1:
#            print(line)
#            print(reverse_complement(line))
#            print(complement(line))
#            print(reverse(line))
        whitelist_barcodes.add(line)
        whitelist_barcodes_rc.add(reverse_complement(line))
        whitelist_barcodes_c.add(complement(line))
        whitelist_barcodes_r.add(reverse(line))

logging.info('Finished reading barcodes')


# In[35]:


# determine barcode length
barcode_lengths = set()
for i in whitelist_barcodes:
    barcode_lengths.add(len(i))
if not len(barcode_lengths) == 1:
    raise ValueError('Barcodes are not all of uniform length')
barcode_length = list(barcode_lengths)[0]


# In[36]:


logging.info('Processing fastq file {}'.format(args.fastq))

match_counts = dict() # offset --> [number_in_whitelist, number_in_whitelist_rc]
line_count = 0
record_count = 0
barcode_counts = dict()
with gzip.open(args.fastq, 'rt') as f:
    for line in f:
        line_count += 1
        if (line_count - 2) % 4 != 0: # skip everything but the sequence itself
            continue
        record_count += 1
        line = line.rstrip()
        if line not in barcode_counts:
            barcode_counts[line] = 0
        barcode_counts[line] += 1
        if len(line) < barcode_length:
            raise ValueError('Barcode reads are shorted than the barcodes in the whitelist')
        for offset in range(0, len(line) - barcode_length + 1):
            if offset not in match_counts:
                match_counts[offset] = [0, 0, 0, 0]
            potential_barcode = line[offset:(offset+barcode_length)]
            if potential_barcode in whitelist_barcodes:
                match_counts[offset][0] += 1
            if potential_barcode in whitelist_barcodes_rc:
                match_counts[offset][1] += 1
            if potential_barcode in whitelist_barcodes_r:
                match_counts[offset][2] += 1
            if potential_barcode in whitelist_barcodes_c:
                match_counts[offset][3] += 1
        if record_count == args.check_first:
            break
        if record_count % 1000000 == 0:
            for offset in match_counts:
                in_whitelist_count = match_counts[offset][0]
                in_whitelist_rc_count = match_counts[offset][1]
                percent_in_whitelist = round(100*float(in_whitelist_count) / record_count, 2)
                percent_in_whitelist_rc = round(100*float(in_whitelist_rc_count) / record_count, 2)
                logging.info('Processed {} records so far. For offset {}, {} ({}%) in whitelist, {} ({}%) in whitelist reverse complemented'.format(record_count, offset, in_whitelist_count,  percent_in_whitelist, in_whitelist_rc_count, percent_in_whitelist_rc))


# In[39]:


df = [[i] + j for i, j in match_counts.items()]
df = pd.DataFrame(df, columns=['start_index', 'in_whitelist', 'in_whitelist_rc', 'in_whitelist_r', 'in_whitelist_c'])
df.in_whitelist = 100 * df.in_whitelist / record_count
df.in_whitelist_rc = 100 * df.in_whitelist_rc / record_count
df.in_whitelist_r = 100 * df.in_whitelist_r / record_count
df.in_whitelist_c = 100 * df.in_whitelist_c / record_count
fig, ax = plt.subplots()
if len(df) > 1:
    ax.plot('start_index', 'in_whitelist', data=df, color='black', label='Using original whitelist')
    ax.plot('start_index', 'in_whitelist_rc', data=df, color='blue', label='Using reverse complemented whitelist')
    ax.plot('start_index', 'in_whitelist_r', data=df, color='green', label='Using reverse whitelist')
    ax.plot('start_index', 'in_whitelist_c', data=df, color='red', label='Using complemented whitelist')
else:
    ax.scatter('start_index', 'in_whitelist', data=df, color='black', label='Using original whitelist', alpha=0.5)
    ax.scatter('start_index', 'in_whitelist_rc', data=df, color='blue', label='Using reverse complemented whitelist', alpha=0.5)
    ax.scatter('start_index', 'in_whitelist_r', data=df, color='green', label='Using reverse whitelist', alpha=0.5)
    ax.scatter('start_index', 'in_whitelist_c', data=df, color='red', label='Using complemented whitelist', alpha=0.5)
ax.set_xlabel('Offset X from barcode read start')
ax.set_ylabel('% observed barcodes in whitelist\n(using {}M reads)'.format(record_count/1e6))
ax.legend()
fig.tight_layout()
fig.savefig(f'{PREFIX}percent-barcodes-in-whitelist.png')
fig.clf()
