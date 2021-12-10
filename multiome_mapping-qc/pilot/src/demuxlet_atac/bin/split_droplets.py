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
    parser.add_argument('--hqaa-threshold', type=int, required=True, help ="""Min HQAA per barcode""")
    parser.add_argument('--fraction-mitochondrial-threshold', type=float, required=True, help ="""Max mitochondrial fraction""")
    parser.add_argument('--tss-threshold', type=float, required=True, help ="""Min TSS enrichment""")
    parser.add_argument('--nbarcodes', type=int, help ="""number of barcodes per subset file""")
    parser.add_argument('--out-prefix', type=str, help ="""Output prefix""")
    args = parser.parse_args()
    return args

def getlib(x):
    s = x.split('-')
    batch = f"{s[0]}-{s[1]}"
    library = s[2]
    return [batch, library]

if __name__ == '__main__':
    
    args = getOpts()

    q = pandas.read_csv(args.qc, sep='\t', header=None, names=['library', 'barcode', 'metric', 'value'])
    metrics = ['hqaa', 'tss_enrichment', 'percent_mitochondrial', 'total_autosomal_reads']
    q = q[q['metric'].isin(metrics)]
    q = q[q['value']!='None']
    q['value'] = q['value'].astype(float)
    q['batch'], q['library'] = zip(*q['library'].map(lambda x: getlib(x.replace("-hg19.ataqv.json.gz", ""))))
    d = q.pivot_table(index=['batch','library','barcode'], columns='metric', values='value').reset_index()
    d['fraction_mitochondrial'] = d['percent_mitochondrial'].map(lambda x: x/100)
    d = d[(d['hqaa']>=args.hqaa_threshold) & (d['fraction_mitochondrial']<args.fraction_mitochondrial_threshold) & (d['tss_enrichment']>=args.tss_threshold)]
    
    d.sort_values('total_autosomal_reads', ascending=False)
    total_barcodes = len(d.index)
    print(f"total selected barcodes = {total_barcodes}")
    nfiles = round(total_barcodes/args.nbarcodes)
    cycles = cycle(list(range(nfiles)))

    files = [next(cycles) for i in list(range(total_barcodes))]
    assert len(files) == total_barcodes
    
    d['subset_file'] = files

    for index, grp in d.groupby('subset_file'):
        filename = f"{args.out_prefix}_{index}"
        print(f"outputting {len(grp.index)} barcodes in file {index}, total reads = {grp['total_autosomal_reads'].sum()}")
        grp[['barcode']].to_csv(filename, header=False, index=False)
