#!/usr/bin/env python

import sys
import json
import gzip
import numpy
import argparse
import re
import pysam
import logging

parser = argparse.ArgumentParser(description='Get read name, barcode, and flag for all reads with barcode in the whitelist', add_help=True)
parser.add_argument('whitelist', help = '10X whitelist')
parser.add_argument('bam', help = 'Bam file')
args = parser.parse_args()

logging.basicConfig(level=logging.DEBUG, format='%(asctime)s - %(levelname)s: %(message)s')

whitelist = set()
with open(args.whitelist, 'r') as f:
    for line in f:
        line = line.rstrip()
        whitelist.add(line)

with pysam.AlignmentFile(args.bam, 'rb') as f:
    for read in f.fetch(until_eof=True):
        if not read.has_tag('CB'):
            continue
        read_name = read.query_name
        barcode = read.get_tag('CB')
        flag = read.flag
        if barcode not in whitelist:
            continue
        print(f'{read_name}\t{barcode}\t{flag}')
