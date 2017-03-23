"""
Utils for Range Minimum Queries
===============================

This module implements the algorithm by Fischer and Heun [FISCHER2006]_ to
answer range minimum/maximum queries in constant time, with linear time
preprocessing.

.. [FISCHER2006] Fischer, J., & Heun, V. (2006). Theoretical and Practical
                 Improvements on the RMQ-Problem, with Applications to LCA
                 and LCE. Proc. Cpm. Volume 4009 of Lncs, 36–48.
                 http://doi.org/10.1007/11780441_5
"""

import math
from typing import List

import numpy


def _array_power2_idx(i: int, j: int):
    return i, i + 2**(j+1)


def _precompute_rmq_table(array: numpy.array, operator=numpy.argmin):
    size = len(array)
    table = -numpy.ones((size, int(numpy.log2(size))), dtype=numpy.int32)

    # Base case, initialise dynamic programming table
    for i in range(size):
        s, e = _array_power2_idx(i, 0)
        table[i, 0] = i + operator(array[s:e])

    k = numpy.zeros((2,), dtype=numpy.int32)
    values = numpy.zeros((2,), dtype=numpy.int32)
    for j in range(1, int(numpy.log2(size))):
        for i in range(size):
            k[0] = table[i, j-1]
            row = min(i + 2**j, size-1)
            k[1] = table[row, j-1]

            values[0] = array[k[0]]
            values[1] = array[k[1]]

            table[i, j] = k[operator(values)]

    return table


class BaseRangeQueryDP:
    """Implements a O(n log n) preprocessing algorithm to answer range minimum
    or maximum queries in constant time.

    For an actual application you'll probably want to use
    :class:`RangeMinimumQuery` because that implementation only requires linear
    preprocessing time. This class exists because :class:`RangeMinimumQuery`
    uses this implementation on a smaller sized array to answer queries.

    .. seealso:: RangeMinimumQuery
    """

    def __init__(self, array: List[int], func):
        self.array = numpy.array(array)
        self.func = func
        self.table = _precompute_rmq_table(self.array, func)

    def query(self, i: int, j: int):
        if j <= i:
            return -1

        if j - i == 1:
            return i

        j = min(j, len(self.array))

        k = numpy.zeros((2,), dtype=numpy.int32)
        values = numpy.zeros((2,), dtype=numpy.int32)

        l = int(numpy.log2(j - i))
        col = max(0, l-1)
        k[0] = self.table[i, col]
        k[1] = self.table[j - 2**l, col]

        values[0] = self.array[k[0]]
        values[1] = self.array[k[1]]

        return k[self.func(values)]


class RangeMinimumQueryDP(BaseRangeQueryDP):
    def __init__(self, array: List[int]):
        super().__init__(array, numpy.argmin)


class RangeMaximumQueryDP(BaseRangeQueryDP):
    def __init__(self, array: List[int]):
        super().__init__(array, numpy.argmax)


class RangeMinimumQuery:
    def __init__(self, array: List[int]):
        self.array = array
        self.block_size = math.log(len(self.array))/4
