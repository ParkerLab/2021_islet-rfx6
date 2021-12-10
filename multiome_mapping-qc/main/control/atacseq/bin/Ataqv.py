#!/usr/bin/env python
# coding: utf-8

import sys
import json
import gzip
import numpy
import logging
import dict_utils

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s: %(message)s')

class Ataqv():
    def __init__(self, metrics_file = None, d = None, peaks = False):
        """Load the metrics_file or json_string into an ataqv object. If peaks = False (default), do not load the peak list. If json_string is not None then it is assumed to be a json string representing a previously-created Ataqv object, and it will be loaded"""
        self.reads_with_mate_mapped_to_different_reference = None
        self.duplicate_fraction_not_in_peaks = None
        self.second_reads = None
        self.total_peaks = None
        self.hqaa = None
        self.tss_coverage = None
        self.ff_reads = None
        self.library = None
        self.first_reads = None
        self.fragment_length_counts = None
        self.reverse_reads = None
        self.maximum_proper_pair_fragment_size = None
        self.properly_paired_and_mapped_reads = None
        self.ppm_not_in_peaks = None
        self.fragment_length_counts_fields = None
        self.paired_reads = None
        self.duplicate_reads = None
        self.peaks = None
        self.reads_with_mate_too_distant = None
        self.duplicate_mitochondrial_reads = None
        self.hqaa_in_peaks = None
        self.hqaa_tf_count = None
        self.reads_mapped_and_paired_but_improperly = None
        self.metrics_url = None
        self.total_peak_territory = None
        self.forward_reads = None
        self.peak_duplicate_ratio = None
        self.rr_reads = None
        self.mapq_counts_fields = None
        self.fr_reads = None
        self.peak_percentiles = None
        self.unmapped_reads = None
        self.forward_mate_reads = None
        self.duplicate_autosomal_reads = None
        self.duplicates_not_in_peaks = None
        self.description = None
        self.rf_reads = None
        self.median_mapq = None
        self.secondary_reads = None
        self.total_mitochondrial_reads = None
        self.reverse_mate_reads = None
        self.unmapped_mate_reads = None
        self.unclassified_reads = None
        self.ppm_in_peaks = None
        self.mapq_counts = None
        self.unpaired_reads = None
        self.total_reads = None
        self.mean_mapq = None
        self.hqaa_mononucleosomal_count = None
        self.peaks_fields = None
        self.duplicates_in_peaks = None
        self.supplementary_reads = None
        self.name = None
        self.reads_mapped_with_zero_quality = None
        self.url = None
        self.organism = None
        self.hqaa_overlapping_peaks_percent = None
        self.total_autosomal_reads = None
        self.tss_enrichment = None
        self.short_mononucleosomal_ratio = None
        self.fragment_length_distance = None
        self.qcfailed_reads = None
        self.duplicate_fraction_in_peaks = None
        self.median_fragment_length = None
        self.fragment_length_distance = None
        self.percent_mitochondrial = None
        self.percent_autosomal_duplicate = None
        self.percent_mitochondrial_duplicate = None
        self.fragment_length_distribution = dict()
        self.chromosome_counts = dict()
        self.max_fraction_reads_from_single_autosome = None
        self.percent_hqaa = None
        self.percent_properly_paired_and_mapped = None
        self.percent_secondary = None
        self.percent_supplementary = None
        self.percent_duplicate = None
        self.percent_unmapped = None
        self.percent_unmapped_mate = None
        self.percent_qcfailed = None
        self.percent_unpaired = None
        self.percent_mapq_0 = None
        self.percent_rf = None
        self.percent_ff = None
        self.percent_rr = None
        self.percent_autosomal = None
        self.percent_mate_separate_chromosome = None
        self.percent_mate_too_distant = None
        self.percent_improperly_paired = None

        if metrics_file is not None:
            load_ataqv_json(f)
            with gzip.open(metrics_file, 'rb') as f:
                line = f.read()
                if isinstance(line, bytes):
                    line = line.decode('utf-8')
                d = json.loads(line)[0]
                if 'metrics' in d:
                    d = d['metrics']

                for metric, value in d.items():
                    if metric == 'peaks' and not peaks:
                        continue
                    setattr(self, metric, value)

                fragment_length_counts = dict()
                if type(self.fragment_length_counts) == type(dict()):
                    for i, vals in self.fragment_length_counts.items():
                        fragment_length_counts[int(i)] = vals[0]
                else:
                    for i, count, proportion in self.fragment_length_counts:
                        fragment_length_counts[int(i)] = int(count)
                self.fragment_length_counts = fragment_length_counts
                self.library = self.library['library']

                # add some additional metrics
                self.percent_mitochondrial = 100.0 * float(self.total_mitochondrial_reads) / self.total_reads if self.total_reads > 0 else None
                self.percent_mitochondrial_duplicate = 100.0 * float(self.duplicate_mitochondrial_reads) / self.total_mitochondrial_reads if self.total_mitochondrial_reads > 0 else None
                self.percent_autosomal_duplicate = 100.0 * float(self.duplicate_autosomal_reads) / self.total_autosomal_reads if self.total_autosomal_reads > 0 else None
                self.median_fragment_length = dict_utils.median_from_dict(self.fragment_length_counts)
                self.percent_hqaa = 100.0 * float(self.hqaa) / self.total_reads
                self.percent_properly_paired_and_mapped = 100.0 * float(self.properly_paired_and_mapped_reads) / self.total_reads
                self.percent_secondary = 100.0 * float(self.secondary_reads) / self.total_reads
                self.percent_supplementary = 100.0 * float(self.supplementary_reads) / self.total_reads
                self.percent_duplicate = 100.0 * float(self.duplicate_reads) / self.total_reads
                self.percent_unmapped = 100.0 * float(self.unmapped_reads) / self.total_reads
                self.percent_unmapped_mate = 100.0 * float(self.unmapped_mate_reads) / self.total_reads
                self.percent_qcfailed = 100.0 * float(self.qcfailed_reads) / self.total_reads
                self.percent_unpaired = 100.0 * float(self.unpaired_reads) / self.total_reads
                self.percent_mapq_0 = 100.0 * float(self.reads_mapped_with_zero_quality) / self.total_reads
                self.percent_rf = 100.0 * float(self.rf_reads) / self.total_reads
                self.percent_ff = 100.0 * float(self.ff_reads) / self.total_reads
                self.percent_rr = 100.0 * float(self.rr_reads) / self.total_reads
                self.percent_autosomal = 100.0 * float(self.total_autosomal_reads) / self.total_reads if self.total_reads > 0 else None
                self.percent_mate_separate_chromosome = 100.0 * float(self.reads_with_mate_mapped_to_different_reference) / self.total_reads
                self.percent_mate_too_distant = 100.0 * float(self.reads_with_mate_too_distant) / self.total_reads
                self.percent_improperly_paired = 100.0 * float(self.reads_mapped_and_paired_but_improperly) / self.total_reads

        elif d is not None:
            for key, val in d.items():
                setattr(self, key, val)
        else:
            logging.error('Need to pass a metrics file or JSON string to Ataqv()')

    def as_dict(self):
        attributes = [x for x in dir(self) if not x.startswith('__') and not callable(getattr(self, x))]
        tmp = {a: getattr(self, a) for a in attributes}
        return tmp


def _process_dict(d, peaks = False):
    assert(type(d) == type(dict()))
    if not peaks:
        if 'peaks' in d:
            del d['peaks']
    fragment_length_counts = dict()
    if type(d['fragment_length_counts']) == type(dict()):
        fragment_length_counts = {int(i): vals[0] for i, vals in d['fragment_length_counts'].items()}
    else:
        fragment_length_counts = {int(i): int(count) for i, count, proportion in d['fragment_length_counts']}
    d['fragment_length_counts'] = fragment_length_counts
    d['library'] = d['library']['library']
    d['percent_mitochondrial'] = 100 * float(d['total_mitochondrial_reads']) / d['total_reads'] if d['total_reads'] > 0 else None
    d['percent_mitochondrial_duplicate'] = 100 * float(d['duplicate_mitochondrial_reads']) / d['total_mitochondrial_reads'] if d['total_mitochondrial_reads'] > 0 else None
    d['percent_autosomal_duplicate'] = 100 * float(d['duplicate_autosomal_reads']) / d['total_autosomal_reads'] if d['total_autosomal_reads'] > 0 else None
    d['median_fragment_length'] = dict_utils.median_from_dict(d['fragment_length_counts'])
    d['percent_hqaa'] = 100.0 * float(d['hqaa']) / d['total_reads']
    d['percent_properly_paired_and_mapped'] = 100.0 * float(d['properly_paired_and_mapped_reads']) / d['total_reads']
    d['percent_secondary'] = 100.0 * float(d['secondary_reads']) / d['total_reads']
    d['percent_supplementary'] = 100.0 * float(d['supplementary_reads']) / d['total_reads']
    d['percent_duplicate'] = 100.0 * float(d['duplicate_reads']) / d['total_reads']
    d['percent_unmapped'] = 100.0 * float(d['unmapped_reads']) / d['total_reads']
    d['percent_unmapped_mate'] = 100.0 * float(d['unmapped_mate_reads']) / d['total_reads']
    d['percent_qcfailed'] = 100.0 * float(d['qcfailed_reads']) / d['total_reads']
    d['percent_unpaired'] = 100.0 * float(d['unpaired_reads']) / d['total_reads']
    d['percent_mapq_0'] = 100.0 * float(d['reads_mapped_with_zero_quality']) / d['total_reads']
    d['percent_rf'] = 100.0 * float(d['rf_reads']) / d['total_reads']
    d['percent_ff'] = 100.0 * float(d['ff_reads']) / d['total_reads']
    d['percent_rr'] = 100.0 * float(d['rr_reads']) / d['total_reads']
    d['percent_autosomal'] = 100.0 * float(d['total_autosomal_reads']) / d['total_reads'] if d['total_reads'] > 0 else None
    d['percent_mate_separate_chromosome'] = 100.0 * float(d['reads_with_mate_mapped_to_different_reference']) / d['total_reads']
    d['percent_mate_too_distant'] = 100.0 * float(d['reads_with_mate_too_distant']) / d['total_reads']
    d['percent_improperly_paired'] = 100.0 * float(d['reads_mapped_and_paired_but_improperly']) / d['total_reads']
    return d




def load_ataqv_json(metrics_file, peaks = False):
    with gzip.open(metrics_file, 'rb') as f:
        line = f.read()
        if isinstance(line, bytes):
            line = line.decode('utf-8')
        d = json.loads(line)
        if 'metrics' in d[0]:
            return [Ataqv(d=_process_dict(i['metrics'], peaks = peaks)) for i in d]
        else:
            return [Ataqv(d=_process_dict(i, peaks=peaks)) for i in d]
