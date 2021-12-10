import os
import sys
import glob
import re

FASTQ_FILES = glob.glob('/lab/data/seqcore/3528-CV/ATAC/fastqs_3528-CV/Sample_3528-CV-*/*')

TARGET_DIR = sys.argv[1]

def parse_file_name(f):
    RE = '/lab/data/seqcore/3528-CV/ATAC/fastqs_3528-CV/Sample_(\d+)-CV-([0-9])-ATAC_([a-z])/.*_S\d+_(R\d+)_.*.fastq.gz'
    sample, well, rg, read = re.match(RE, f).groups()
    return {'sample': sample, 'well': f"CV{well}", 'rg': rg, 'read': read}


for f in FASTQ_FILES:
    info = parse_file_name(f)
    info['library'] = '{}_ATAC'.format(info['sample'])
    new_file = os.path.join(TARGET_DIR, '{library}_{well}_{rg}.{read}.fastq.gz'.format(**info))
    os.symlink(f, new_file)
