#!/usr/bin/env python
# coding: utf-8


#!/usr/bin/python

import os
import sys
import gzip
import logging
import argparse
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd

parser = argparse.ArgumentParser(description='Compare nuclear barcodes to the barcode whitelist to infer whether they need to be reverse complemented', add_help=True)
parser.add_argument('--check-first', default=1000000, type=int, help='Check no more than this number of barcodes (default: 1000000)')
parser.add_argument('--barcode-length', dest='barcode_length', default=16, type=int, help='Barcode length (default: 16)')
parser.add_argument('--prefix', default='barcode-check.', help='Prefix for output files (default: barcode-check.)')
parser.add_argument('fastq', help='Fastq file containing nuclear barcodes (gzipped)')
#args = parser.parse_args(['/lab/work/porchard/Nova-303/data/fastq/1846_ATAC_b.2.fastq.gz'])
args = parser.parse_args()

PREFIX = args.prefix

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')


# In[7]:


logging.info('Processing fastq file {}'.format(args.fastq))

line_count = 0
record_count = 0
barcode_counts = dict() # offset --> barcode --> count
with gzip.open(args.fastq, 'rt') as f:
    for line in f:
        line_count += 1
        if (line_count - 2) % 4 != 0: # skip everything but the sequence itself
            continue
        record_count += 1
        line = line.rstrip()
        if len(line) < args.barcode_length:
            raise ValueError('Barcode reads are shorted than the barcodes in the whitelist')
        for offset in range(0, len(line) - args.barcode_length + 1):
            if offset not in barcode_counts:
                barcode_counts[offset] = dict()
            potential_barcode = line[offset:(offset+args.barcode_length)]
            if potential_barcode not in barcode_counts[offset]:
                barcode_counts[offset][potential_barcode] = 0
            barcode_counts[offset][potential_barcode] += 1
        if record_count == args.check_first:
            break
        if record_count % 1000000 == 0:
            logging.info('Processed {} records so far'.format(record_count))


# In[11]:


for offset in list(barcode_counts.keys()):
    df = pd.DataFrame([[barcode, count] for barcode, count in barcode_counts[offset].items()], columns=['barcode', 'count'])
    df = df.sort_values('count', ascending=False)
    df['rank'] = list(range(1, len(df) + 1))
    fig, ax = plt.subplots()
    ax.plot('rank', 'count', data=df)
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_xlabel('Barcode')
    ax.set_ylabel('Barcode count')
    fig.tight_layout()
    fig.savefig(f'{PREFIX}offset_{offset}.elbow-plot.png')
    fig.clf()
