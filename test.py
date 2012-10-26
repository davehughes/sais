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



class TestSequenceMatching(unittest.TestCase):

    def testSimpleMatch(self):
        s1 = sais.Sequence(''.join(generate_sequence(length=1000)))
        s2 = sais.Sequence(''.join(generate_sequence(length=1000)))
        match_tuples = sais.find_best_subsequence_matches(s1, s2)

        # scan to make sure lcp matches
        for idx, (match, lcp) in enumerate(match_tuples):
            match_idx1 = s1.suffix_array[match]
            match_idx2 = s2.suffix_array[idx]
            self.assertEqual(s1.seq[match_idx1: match_idx1 + lcp],
                             s2.seq[match_idx2: match_idx2 + lcp])


if __name__ == '__main__':
    unittest.main()
