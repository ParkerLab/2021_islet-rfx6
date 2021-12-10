#!/usr/bin/env python

import sys
import json
import gzip
import numpy
import argparse
import Ataqv

parser = argparse.ArgumentParser()
parser.add_argument('--files', nargs = '*', help = 'Names of ataqv metric files from which to extract metrics. These files will be within the data/ directory of an ataqv app.')
args = parser.parse_args()


for ataqv_data_file in args.files:
    d = Ataqv.load_ataqv_json(ataqv_data_file)[0] # TODO: this assumes that only one dict is in the json...which won't always be the case
    for bp in d.tss_coverage:
        position = int(bp[0]) - 1001
        coverage = bp[1]
        print('{}\t{}\t{}'.format(ataqv_data_file, position, coverage))
