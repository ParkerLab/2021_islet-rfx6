#!/usr/bin/env python

import os
import sys
import argparse
import gzip
from Bio import SeqIO
from Bio.SeqIO.QualityIO import FastqGeneralIterator
import re


parser = argparse.ArgumentParser(description="Append UMIs from one fastq file to the read names of another fastq file")
parser.add_argument('umi_fastq', action='store', type=str)
parser.add_argument('read_fastq', action='store', type=str)

args = parser.parse_args()

READ_NAME_RE = re.compile('(.*)\s.*')

with gzip.open(args.umi_fastq, 'rt') as umi_fastq:
    with gzip.open(args.read_fastq, 'rt') as read_fastq:
        x = FastqGeneralIterator(read_fastq)
        for umi_read_name, barcode, umi_phred in FastqGeneralIterator(umi_fastq):
            read_name, read, phred = next(x)
            read_name = re.match(READ_NAME_RE, read_name).group(1)
            print('@' + read_name + '_' + barcode)
            print(read)
            print('+')
            print(phred)
