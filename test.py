from itertools import islice
import os
import random
import unittest

import sais


# Functions for generating test sequences
def generate_sequence(alphabet='ACTG', length=None):
    def generate():
        while True:
            yield random.choice(alphabet)

    if length:
        return islice(generate(), length)
    else:
        return generate()


def create_sequence_file(path, sequences=1, sequence_length=100000):
    f = open(path, 'wb')
    with open(path, 'wb') as f:
        for _ in range(sequences):
            for c in generate_sequence(length=sequence_length):
                f.write(c)
            f.write(os.linesep)


class TestSuffixArrays(unittest.TestCase):

    def testSuffixArray(self):
        seq = sais.Sequence('banana$')
        sa = seq.suffix_array
        self.assertEqual([6, 5, 3, 1, 0, 4, 2], sa[:len(seq)])

    def testLongestCommonPrefix(self):
        seq = sais.Sequence('banana$')
        lcp = seq.lcp_array
        self.assertEqual([-1, 0, 1, 3, 0, 0, 2], lcp[:len(seq)])


if __name__ == '__main__':
    unittest.main()
