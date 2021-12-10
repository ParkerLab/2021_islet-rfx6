#!/usr/bin/env python
# coding: utf-8

import os
import sys
import argparse
import pandas


def getOpts():
    parser = argparse.ArgumentParser(
        description='Format SNP position and annotation overlap file')
    parser.add_argument('--link', required=True,
                        help="""Annotation link file to fetch the order of annotations""")
    parser.add_argument('--annots', nargs='+', required=True,
                        help="""Annotation overlap bed files""")
    parser.add_argument('--snp', required=True, help="""SNP pos file""")
    parser.add_argument('--output', required=True, help="""output file name""")
    args = parser.parse_args()
    return args


def read_annot(fname):
    a = pandas.read_csv(fname, sep='\t', header=None, index_col=[0])
    return a


if __name__ == '__main__':
    args = getOpts()

    print("Reading SNP file..", file=sys.stderr)
    d = pandas.read_csv(args.snp, sep=' ', header=None,
                        usecols=[0]).set_index(0)

    link = pandas.read_csv(args.link, sep='\t', usecols=['Index', 'Path'], dtype={
                           'Index': int}).sort_values('Index')
    annotation_order = link['Path'].map(os.path.basename).tolist()

    annot_file_basenames = [os.path.basename(f) for f in args.annots]
    assert sorted(annotation_order) == sorted(annot_file_basenames)

    annotation_files = [
        f for f in args.annots for af in annotation_order if os.path.basename(f) == af]

    print("Creating annotation data frame..", file=sys.stderr)
    for f in annotation_files:
        print(f"Reading {f}..")
        d = pandas.concat([d] + [read_annot(f)], axis=1, join="inner")

    d = d.astype(str)
    print("Writing output csv file..", file=sys.stderr)
    d.apply(lambda x: "".join(x), axis=1).to_csv(
        args.output, header=False, sep=" ")
