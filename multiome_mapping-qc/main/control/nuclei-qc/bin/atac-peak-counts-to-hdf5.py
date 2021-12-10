#!/usr/bin/env python
# coding: utf-8

import os
import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "--counts",
    required=True,
    help="Paths to the counts file"
)
parser.add_argument(
    "--barcodes",
    required=True,
    help="Path to the text file containing QC-pass barcodes"
)
parser.add_argument(
    "--out",
    required=True,
    help="Path to the output HDF5 file"
)

args = parser.parse_args()

include_barcodes = pd.read_csv(args.barcodes, header=None).iloc[:,0]

atac_counts = pd.read_csv(args.counts, header=None, sep='\t', names = ["library", "barcode", "feature", "fragments"])
atac_counts = atac_counts[atac_counts.barcode.isin(include_barcodes)]
atac_counts = atac_counts.pivot(index='barcode', columns='feature', values='fragments').fillna(0)
atac_counts = atac_counts.astype(int)
atac_counts.to_hdf(args.out, key=os.path.basename(args.out).replace('.hdf5', ''))

