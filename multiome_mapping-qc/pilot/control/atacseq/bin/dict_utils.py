#!/usr/bin/env python

import numpy

def find_index(l, i):
    """Input: [3, 2, 6], 4
    Output: 1
    Logic:
    [3, 2, 6] --> [0, 0, 0, 1, 1, 2, 2, 2, 2, 2, 2]
    4th value in that list is 1
    """

    assert(sum(l) >= i)

    running_sum = 0
    current_index = 0

    for j in l:
        if running_sum + j < i:
            running_sum += j
            current_index += 1
        else:
            return current_index


def median_from_dict(d):
    """Given a dictionary representing count data (d[value] --> count), return the median value"""

    keys = sorted(d.keys())
    values = [d[i] for i in keys]

    assert(sum(values) > 0)

    # remove 0 values...
    while 0 in values:
        i = values.index(0)
        del keys[i]
        del values[i]

    count = sum(values)

    middle = None

    if count % 2 == 0:
        middle = [int(count / 2) - 1, int(count / 2)]  # we want these indices, when the values are rep()'ed and ranked
    else:
        middle = [int(count) / 2]

    indices = [find_index(values, i + 1) for i in middle]
    corresponding_keys = [keys[i] for i in indices]

    return numpy.mean(corresponding_keys)
