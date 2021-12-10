#!/usr/bin/env python
# coding: utf-8

import os
import sys
import pandas as pd

flagstat = dict()

def _interpret_flag(d):
    binary = bin(d).replace('0b', '').rjust(12, '0')
    x = {interpretation: bit == '1' for interpretation, bit in zip(['is_paired', 'is_proper_pair', 'is_unmapped', 'mate_unmapped', 'is_reverse', 'mate_reverse', 'is_first', 'is_second', 'is_secondary', 'fail_qc', 'is_duplicate', 'is_supplementary'], binary[::-1])}
    return x


def interpret_flag(d):
    assert(isinstance(d, int))
    if d not in flagstat:
        flagstat[d] = _interpret_flag(d)
    return flagstat[d]
