import pytest  # noqa

from phasm.rmq import (RangeMinimumQueryDP, RangeMaximumQueryDP,
                       RangeMinimumQuery, RangeMaximumQuery)


def test_rmq_min_dp():
    a = [2, 4, 6, 8, 1]

    rmq = RangeMinimumQueryDP(a)
    k = 2

    min_values = [rmq.query(i, i+k) for i in range(len(a) - k + 1)]
    assert min_values == [0, 1, 2, 4]

    a = [8, 6, 4, 6, 10, 2, 3, 5, 0]
    rmq = RangeMinimumQueryDP(a)
    k = 3
    min_values = [rmq.query(i, i+k) for i in range(len(a) - k + 1)]
    assert min_values == [2, 2, 2, 5, 5, 5, 8]
    assert rmq.query(0, len(a)) == 8

    assert rmq.query(4, 2) == -1
    assert rmq.query(1, 2) == 1


def test_rmq_max_dp():
    a = [2, 4, 6, 8, 1]

    rmq = RangeMaximumQueryDP(a)
    k = 2

    max_values = [rmq.query(i, i+k) for i in range(len(a) - k + 1)]
    assert max_values == [1, 2, 3, 3]

    a = [8, 6, 4, 6, 10, 2, 3, 5, 0]
    rmq = RangeMaximumQueryDP(a)
    k = 3
    min_values = [rmq.query(i, i+k) for i in range(len(a) - k + 1)]
    assert min_values == [0, 1, 4, 4, 4, 7, 7]
    assert rmq.query(0, len(a)) == 4

    assert rmq.query(4, 2) == -1
    assert rmq.query(1, 2) == 1


def test_rmq_min():
    a = [2, 4, 6, 8, 1]

    rmq = RangeMinimumQuery(a)
    k = 2

    min_values = [rmq.query(i, i+k) for i in range(len(a) - k + 1)]
    assert min_values == [0, 1, 2, 4]

    a = [8, 6, 4, 6, 10, 2, 3, 5, 0]
    rmq = RangeMinimumQueryDP(a)
    k = 3
    min_values = [rmq.query(i, i+k) for i in range(len(a) - k + 1)]
    assert min_values == [2, 2, 2, 5, 5, 5, 8]
    assert rmq.query(0, len(a)) == 8

    assert rmq.query(4, 2) == -1
    assert rmq.query(1, 2) == 1


def test_rmq_max():
    a = [2, 4, 6, 8, 1]

    rmq = RangeMaximumQuery(a)
    k = 2

    max_values = [rmq.query(i, i+k) for i in range(len(a) - k + 1)]
    assert max_values == [1, 2, 3, 3]

    a = [8, 6, 4, 6, 10, 2, 3, 5, 0]
    rmq = RangeMaximumQueryDP(a)
    k = 3
    min_values = [rmq.query(i, i+k) for i in range(len(a) - k + 1)]
    assert min_values == [0, 1, 4, 4, 4, 7, 7]
    assert rmq.query(0, len(a)) == 4

    assert rmq.query(4, 2) == -1
    assert rmq.query(1, 2) == 1
