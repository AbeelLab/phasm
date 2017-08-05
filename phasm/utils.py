import json
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


class DebugDataLogger:
    def __init__(self, f=None):
        self.outfile = f

    def new_bubble(self, entrance, exit, start_of_block, rel_read_info):
        if not self.outfile:
            return

        json.dump({
            'type': 'new_bubble',
            'entrance': str(entrance),
            'exit': str(exit),
            'start_of_block': start_of_block,
            'rel_read_info': {str(r): {
                'overlap_max': i.overlap_max,
                'alignments': [a.as_tuple() for a in i.alignments]
            } for r, i in rel_read_info.items()}
        }, self.outfile)
        self.outfile.write("\n")

    def candidate_set(self, haplotype_set, extension, read_probs, p_sr, prior):
        if not self.outfile:
            return

        json.dump({
            'type': 'candidate_set',
            'haplotype_set': [
                [str(n) for n in h] for h in haplotype_set.haplotypes
            ],
            'extension': [
                [str(n) for n in h] for h in extension
            ],
            'read_probs': read_probs,
            'p_sr': p_sr,
            'prior': prior
        }, self.outfile)
        self.outfile.write("\n")

    def haploblock(self, hs, name):
        if not self.outfile:
            return

        json.dump({
            'type': 'haploblock',
            'name': name,
            'from_large_bubble': hs.from_large_bubble,
            'haplotypes': [
                [str(n) for n in h] for h in hs.haplotypes
            ]
        }, self.outfile)
        self.outfile.write("\n")
