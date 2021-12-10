import os
import sys
import glob
import re

import argparse

parser = argparse.ArgumentParser()
parser.add_argument(
    "-fq",
    "--fastq-dir",
    required=True,
    help="""Path to the parent directory containing per-sample FASTQ directories"""
)
parser.add_argument(
    "-o",
    "--out-dir",
    required=True,
    help="Output directory where to symlink files"
)

args = parser.parse_args()

FASTQ_FILES = glob.glob(f'{args.fastq_dir}/**/*')

def parse_file_name(f):
    RE = f'{args.fastq_dir}/Sample_(\d+)-CV-([0-9])-ATAC_([a-z])/.*_S\d+_(R\d+)_.*.fastq.gz'
    sample, well, rg, read = re.match(RE, f).groups()
    return {'sample': sample, 'well': f"CV{well}", 'rg': rg, 'read': read}


for f in FASTQ_FILES:
    info = parse_file_name(f)
    info['library'] = '{}_ATAC'.format(info['sample'])
    new_file = os.path.join(args.out_dir, '{library}_{well}_{rg}.{read}.fastq.gz'.format(**info))
    try:
        os.symlink(f, new_file)
    except:
        print("fatal: target already exists")
        exit(1)
