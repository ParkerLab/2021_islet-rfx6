#!/usr/bin/env python

import os
import sys
import glob
import pandas as pd
import numpy as np
import json
import gzip

JSON, OUT = sys.argv[1:]

qc = None
with gzip.open(JSON, 'rt') as f:
    qc = json.load(f)

for i in range(len(qc)):
    number_mitochondrial = qc[i]['chromosome_read_counts']['chrM'] if 'chrM' in qc[i]['chromosome_read_counts'] else 0
    number_total = sum(qc[i]['chromosome_read_counts'].values())
    qc[i]['fraction_mitochondrial'] = number_mitochondrial / number_total if number_total != 0 else np.nan


# what if we use marker gene specificity as well? if plot the % of reads assigned to each
cols = [k for k, v in qc[0].items() if isinstance(v, int) or isinstance(v, float) or isinstance(v, str)]
df = pd.DataFrame([[x[i] for i in cols] for x in qc])
df.columns = list(cols)
df = df.drop(columns=['reads_no_umi', 'reads_passing_filters'])
df['percent_mapped'] = df.uniquely_mapped_reads / df.total_reads
df['percent_no_gene'] = df.reads_no_gene / df.total_reads
df.to_csv(OUT, index=False, sep = '\t')
