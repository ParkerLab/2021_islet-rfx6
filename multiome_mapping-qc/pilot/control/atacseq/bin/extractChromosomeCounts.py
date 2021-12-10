#!/usr/bin/env python

import sys
import argparse
import Ataqv

parser = argparse.ArgumentParser()
parser.add_argument('--files', nargs = '*', help = 'Names of ataqv metric files from which to extract metrics. These files will be within the data/ directory of an ataqv app.')
args = parser.parse_args()

for ataqv_data_file in args.files:
    print("Loading Ataqv json file..", file=sys.stderr)
    d = Ataqv.load_ataqv_json(ataqv_data_file)
    print("done. writing data now..", file=sys.stderr)
    for item in d:
        for chrom, count in item.chromosome_counts:
            print('{}\t{}\t{}\t{}'.format(ataqv_data_file, item.name, chrom, count))
