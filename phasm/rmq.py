"""
Utils for Range Minimum/Maximum Queries
=======================================

This module implements the algorithm by Fischer and Heun [FISCHER2006]_ to
answer range minimum/maximum queries in constant time, with linear time
preprocessing.

This module is a bit of a mess in terms of functions that could be part of
class, but this is to leave the possibility open for Numba JIT compilation.
Numba JIT complilation does not work very well with class methods, because it
does not know how to optimise the first `self` parameter of methods.

.. [FISCHER2006] Fischer, J., & Heun, V. (2006). Theoretical and Practical
                 Improvements on the RMQ-Problem, with Applications to LCA
                 and LCE. Proc. Cpm. Volume 4009 of Lncs, 36–48.
                 http://doi.org/10.1007/11780441_5
"""

import math
from typing import List, Tuple, Callable

import numpy


def _array_power2_idx(i: int, j: int) -> Tuple[int, int]:
    return i, i + 2**(j+1)


def _precompute_rmq_table(array: numpy.array,
                          operator: Callable=numpy.argmin) -> numpy.array:
    """This function fills the dynamic programming table with the indices
    of the computed minimum or maximum in a given range."""

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

    This is a base class: it needs to know which kind of operation you want to
    perform on the ranges, for example minimum or maximum. For this there are
    two convencience classes available: `RangeMinimumQueryDP` and
    `RangeMaximumQueryDP`.

    For an actual application you'll probably want to use
    `RangeMinimumQuery` or `RangeMaximumQuery` because those implementations
    only require linear preprocessing time. This class exists because
    `BaseRangeQuery` uses this implementation on a smaller sized array to
    answer queries.

    .. seealso:: RangeMinimumQueryDP, RangeMaximumQueryDP, BaseRangeQuery,
                 RangeMaximumQuery, RangeMinimumQuery
    """

    def __init__(self, array: List[int], func: Callable):
        self.array = numpy.array(array)
        self.func = func
        self.table = _precompute_rmq_table(self.array, func)

    def query(self, i: int, j: int) -> int:
        if j <= i:
            return -1

        if j - i == 1:
            return i

        j = min(j, len(self.array))

        k = numpy.zeros((2,), dtype=numpy.int32)
        values = numpy.zeros((2,))

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


def _calculate_block_values(array: numpy.array, block_size: int,
                            func: Callable) -> numpy.array:
    num_blocks = len(array) // block_size
    values = numpy.zeros((num_blocks,))
    for i in range(num_blocks):
        bs = i * block_size
        be = min(len(array), (i+1) * block_size)

        values[i] = func(array[bs:be])

    return values


class BaseRangeQuery:
    """This is an implementation for the range minimum/maximum query problem
    with O(n) preprocessing time, and answers each query in O(log n) time.
    It is possible to construct an algorithm which has a theoretical constant
    runtime for answering queries, in practice this was not faster than the
    current implementation [FISCHER2006]_.

    The linear preprocessing time is obtained by dividing the original array
    into multiple blocks, and calculate the minimum/maximum value for each
    block. These minima/maxima are stored in a new array, and this new array is
    then preprocessed using the dynamic programming implementation for RMQ.

    Now when we need to query a range, we only need to check the boundary
    blocks, and if there more blocks in between, we can quickly query which
    block is the minimum/maximum.

    This is the base class, so you should provide the "operator" function, like
    argmin or argmax. There are two convenience classes: `RangeMinimumQuery`
    and `RangeMaximumQuery`.

    .. seealso:: RangeMinimumQuery, RangeMaximumQuery, RangeMinimumQueryDP,
                 RangeMaximumQueryDP, BaseRangeQueryDP
    """

    def __init__(self, array: List[int], func: Callable, argfunc: Callable):
        self.array = numpy.array(array)
        self.block_size = int(math.ceil(math.log(len(self.array))/4))

        self.func = func
        self.argfunc = argfunc
        self.block_values = _calculate_block_values(self.array,
                                                    self.block_size,
                                                    self.func)
        self.rmq_block = BaseRangeQueryDP(self.block_values, self.argfunc)

    def query(self, i: int, j: int) -> int:
        if j <= i:
            return -1

        if j - i == 1:
            return i

        # Calculate the most left and most right block of the given range
        left_boundary = i // self.block_size
        right_boundary = (j-1) // self.block_size

        if (left_boundary == right_boundary or
                right_boundary - left_boundary == 1):
            # Range is within the same block, or we have two consecutive
            # blocks, calculate using naive method
            return i + self.argfunc(self.array[i:j])
        else:
            k = numpy.zeros((3,), dtype=numpy.int32)
            values = numpy.zeros((3,))

            le = (left_boundary+1) * self.block_size
            k[0] = i + self.argfunc(self.array[i:le])
            values[0] = self.array[k[0]]

            # Value of blocks in between, can be easily queried using
            # our preprocessed block values
            block_winner = self.rmq_block.query(left_boundary+1,
                                                right_boundary)
            bs = block_winner * self.block_size
            be = (block_winner+1) * self.block_size
            k[1] = bs + self.argfunc(self.array[bs:be])
            values[1] = self.array[k[1]]

            rs = right_boundary * self.block_size
            k[2] = rs + self.argfunc(self.array[rs:j])
            values[2] = self.array[k[2]]

            return k[self.argfunc(values)]


class RangeMinimumQuery(BaseRangeQuery):
    def __init__(self, array: List[int]):
        super().__init__(array, numpy.amin, numpy.argmin)


class RangeMaximumQuery(BaseRangeQuery):
    def __init__(self, array: List[int]):
        super().__init__(array, numpy.amax, numpy.argmax)
