#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse

import pandas as pd
from scipy.io import mmread

# convert RNA STARsolo counts to HDF5 format

#RNA_MTX = '/lab/work/porchard/Nova-315/work/rnaseq/results/starsolo/1846_RNA-hg19/Solo.out/GeneFull/raw/matrix.mtx'
#RNA_BARCODES = '/lab/work/porchard/Nova-315/work/rnaseq/results/starsolo/1846_RNA-hg19/Solo.out/GeneFull/raw/barcodes.tsv'
#RNA_FEATURES = '/lab/work/porchard/Nova-315/work/rnaseq/results/starsolo/1846_RNA-hg19/Solo.out/GeneFull/raw/features.tsv'
#OUT = 'rna.hdf5'

parser = argparse.ArgumentParser()
parser.add_argument(
    "--solo-dir",
    required=True,
    help="Paths to the counts file"
)
parser.add_argument(
    "--barcodes",
    required=True,
    help="Path to the text file containing QC-pass barcodes"
)
parser.add_argument(
    "--library",
    required=True,
    help="Library name/label to use"
)
parser.add_argument(
    "--out",
    required=True,
    help="Path to the output HDF5 file"
)

args = parser.parse_args()

include_barcodes = pd.read_csv(args.barcodes, header=None).iloc[:, 0].to_list()

print("Reading data..", file=sys.stderr)

rna_counts = mmread(f"{args.solo_dir}/GeneFull/raw/matrix.mtx").tocsc()
rna_barcodes = pd.read_csv(f"{args.solo_dir}/GeneFull/raw/barcodes.tsv", header=None).iloc[:,0]
rna_features = pd.read_csv(f"{args.solo_dir}/GeneFull/raw/features.tsv", sep='\t', header=None).iloc[:,-1]

KEEP_BARCODES = rna_barcodes.isin(include_barcodes)
KEEP_BARCODES = KEEP_BARCODES[KEEP_BARCODES].index.to_list()

rna_counts = pd.DataFrame(rna_counts[:, KEEP_BARCODES].todense())
rna_counts.index = rna_features
rna_counts.columns = rna_barcodes[KEEP_BARCODES]
rna_counts = rna_counts.T

rna_counts.index = ['{}-{}'.format(args.library, i) for i in rna_counts.index]
rna_counts.index = rna_counts.index.rename('barcode')
rna_counts.columns = rna_counts.columns.rename('gene')

# check for and sum duplicate features
feature_counts = rna_features.value_counts()
duplicate_features = feature_counts[feature_counts>1].index.to_list()
if len(duplicate_features) > 0:
    nonduplicate_features = feature_counts[feature_counts==1].index.to_list()
    deduplicated_features = pd.concat([rna_counts[[i]].sum(axis=1).to_frame(name=i) for i in duplicate_features], axis=1)
    rna_counts = pd.concat([rna_counts[nonduplicate_features], deduplicated_features], axis=1)
    assert(len(rna_counts.columns) == len(feature_counts))

print("Writing output..", file=sys.stderr)
rna_counts.to_hdf(args.out, key=os.path.basename(args.out).replace('.hdf5', '').replace('-', '_'))
