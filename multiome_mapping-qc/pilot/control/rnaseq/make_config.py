#!/usr/bin/env python

import os
import sys
import re
import glob
import yaml

ROOT = sys.argv[1]

DATA = os.path.join(ROOT, 'data')
FASTQ_DIR = os.path.join(DATA, 'fastq', 'rna')
FASTQ_FILES = glob.glob(os.path.join(FASTQ_DIR, '*.fastq.gz'))

def parse_fastq_name(f):
    RE = '(.*)_R(\d+).fastq.gz'
    library, read = re.match(RE, os.path.basename(f)).groups()
    readgroup = 'rg1'
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
    libraries[library]['readgroups'][readgroup][read] = f

CONFIG = {
    'libraries': libraries
}
print(yaml.dump(CONFIG, sort_keys = True, indent = 4))
