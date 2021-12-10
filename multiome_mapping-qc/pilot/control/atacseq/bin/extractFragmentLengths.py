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
    for fragment_length, count in sorted(d.fragment_length_counts.items()):
        print('{}\t{}\t{}'.format(ataqv_data_file, fragment_length, count))
