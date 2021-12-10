#! /usr/bin/env python3

import sys
import re
import yaml
import pathlib
import glob

ROOT = "/home/vivekrai/analyses/2020-01_vanderbilt_rna"
BAM_RE = "^(.*)_(.*)-(.*)-(\w+).*"

LIBRARIES = {}


def parse_bam(f):
    """ Return FASTQ "library" name and whether it is "first" or "second" set
    of reads (paired-end)."""

    m = re.search(BAM_RE, pathlib.Path(f).name)
    sample_id = m.group(1)
    sample_type = m.group(2)
    status = m.group(3)
    age = m.group(4)
    return [pathlib.Path(f).name, sample_type, status]


def create_library_item():
    files = glob.glob(f"{ROOT}/data/bams/*.bam")

    for bam in files:
        lib, s, status = parse_bam(bam)

        if s == 'i':
            s = 'Whole.islet'
        elif s == 'a':
            s = "Alpha.cell"
        elif s == 'b':
            s = "Beta.cell"

        if lib not in LIBRARIES:
            LIBRARIES[lib] = {}
            LIBRARIES[lib]['genome'] = 'hg38'
            LIBRARIES[lib]['type'] = s
            LIBRARIES[lib]['status'] = status
            LIBRARIES[lib]['bam'] = bam


if __name__ == "__main__":
    create_library_item()
    print(yaml.dump({"libraries": LIBRARIES}, indent=4, default_flow_style=False))
