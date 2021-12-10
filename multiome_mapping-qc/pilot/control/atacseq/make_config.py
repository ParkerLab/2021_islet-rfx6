#!/usr/bin/env python

import os
import sys
import re
import glob
import yaml

ROOT = sys.argv[1]

DATA = os.path.join(ROOT, 'data')
FASTQ_DIR = os.path.join(DATA, 'fastq', 'atac')
FASTQ_FILES = glob.glob(os.path.join(FASTQ_DIR, '*.fastq.gz'))

def parse_fastq_name(f):
    RE = '^(.*)_([a-z]).R(\d+).fastq.gz'
    library, readgroup, read = re.match(RE, os.path.basename(f)).groups()
    return {'library': library, 'readgroup': readgroup, 'read': read}


libraries = set([parse_fastq_name(f)['library'] for f in FASTQ_FILES])
libraries = {library: {'genome': ['hg19-mCherry-mKate2'], 'readgroups': {}} for library in libraries}

for f in FASTQ_FILES:
    info = parse_fastq_name(f)
    library = info['library']
    readgroup = '{library}_{readgroup}'.format(**info)
    if readgroup not in libraries[library]['readgroups']:
        libraries[library]['readgroups'][readgroup] = dict()
    read = info['read']
    if read == '1':
        read = '1'
    elif read == '2':
        read = 'index'
    elif read == '3':
        read = '2'
    libraries[library]['readgroups'][readgroup][read] = f

CONFIG = {'libraries': libraries}
print(yaml.dump(CONFIG, sort_keys = True, indent = 4))
