#!/usr/bin/env python
# coding: utf-8

import pandas
import math
import sys
import os
import glob
from itertools import cycle
import argparse

def getOpts():
    parser = argparse.ArgumentParser(description='Split droplet barcodes into set of similarly sized subset bams ')
    parser.add_argument('--qc', required=True, help="""qc file""")
    parser.add_argument('--umi-threshold', type=int, required=True, help ="""UMI threshold""")
    parser.add_argument('--mitochondrial-fraction', type=float, required=True, help ="""Max mitochondrial fraction""")
    parser.add_argument('--nbarcodes', type=int, help ="""number of barcodes per subset file""")
    parser.add_argument('--out-prefix', type=str, help ="""Output prefix""")
    args = parser.parse_args()
    return args


if __name__ == '__main__':
    
    args = getOpts()

    q = pandas.read_csv(args.qc, sep='\t',
            usecols=['barcode', 'total_reads', 'umis', 'fraction_mitochondrial'],
            na_values="None"
        )
    q = q[(q['umis']>=args.umi_threshold) & (q['fraction_mitochondrial']<args.mitochondrial_fraction)]
    q.sort_values('total_reads', ascending=False)
    total_barcodes = len(q.index)
    
    nfiles = round(total_barcodes/args.nbarcodes)
    cycles = cycle(list(range(nfiles)))

    files = [next(cycles) for i in list(range(total_barcodes))]
    assert len(files) == total_barcodes
    
    q['subset_file'] = files

    for index, grp in q.groupby('subset_file'):
        filename = f"{args.out_prefix}_{index}"
        print(f"outputting {len(grp.index)} barcodes in file {index}, total reads = {grp['total_reads'].sum()}")
        grp[['barcode']].to_csv(filename, header=False, index=False)
