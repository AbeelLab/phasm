import random
import string


def round_up(number: int, multiple: int) -> int:
    """Round a number up to a multiple of another number. Only works for
    positive numbers.

    Example:

    >>> round_up(57, 100)
    100
    >>> round_up(102, 100)
    200
    >>> round_up(43, 50)
    50
    >>> round_up(77, 50)
    100
    """
    assert multiple != 0

    return int((number + multiple - 1) / multiple) * multiple


def random_string(n, alphabet=string.ascii_letters):
    return ''.join(
        random.choice(alphabet) for _ in range(n)
    )
