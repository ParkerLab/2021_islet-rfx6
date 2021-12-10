#!/usr/bin/env python

import sys
import json
import gzip
import numpy
import argparse
import Ataqv

parser = argparse.ArgumentParser()
parser.add_argument('--metrics', nargs = '*', help = 'Names of metrics to extract (acceptable metrics listed below).')
parser.add_argument('--files', nargs = '*', help = 'Names of ataqv metric files from which to extract metrics. These files will be within the data/ directory of an ataqv app.')
args = parser.parse_args()


for ataqv_data_file in args.files:
    d = Ataqv.load_ataqv_json(ataqv_data_file)[0]
    for metric in args.metrics:
        if not hasattr(d, metric):
            raise ValueError('Metric {} not present'.format(metric))
        print('{}\t{}\t{}'.format(ataqv_data_file, metric, getattr(d, metric)))
