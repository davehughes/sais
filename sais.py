import os
from ctypes import *

here = os.path.dirname(__file__)
SAIS = cdll.LoadLibrary(os.path.join(here, 'lib/libsais.so'))

LIBC = cdll.LoadLibrary('libc.so.6')
LIBC.malloc.argtypes = [c_int]


def malloc(size, type=None):
    '''
    Helper function for calling libc's malloc() and appropriately casting the
    result.  If `type` is specified, the return value is cast as a pointer to
    that type and the size is multiplied by the sizeof() that type.
    '''
    if type:
        LIBC.malloc.restype = POINTER(type)
        size = size * sizeof(type)
    else:
        LIBC.malloc.restype = c_void_p
    return LIBC.malloc(size)


class Sequence(object):

    def __init__(self, seq, terminator='$'):
        self.seq = seq

    @property
    def suffix_array(self):
        if not getattr(self, '_suffix_array', None):
            suffix_array = malloc(len(self.seq), c_int)
            status = SAIS.sais(c_char_p(self.seq), suffix_array, len(self.seq))
            if status:
                raise Exception("Error code %s encountered while constructing suffix array", status)
            self._suffix_array = suffix_array
        return self._suffix_array

    @property
    def lcp_array(self):
        if not getattr(self, '_lcp_array', None):
            lcp_array = malloc(len(self), c_int)
            status = SAIS.compute_lcp(c_char_p(self.seq),
                                      self.suffix_array,
                                      lcp_array,
                                      len(self))
            if status:
                raise Exception("Error code %s encountered while constructing lcp array", status)
            self._lcp_array = lcp_array
        return self._lcp_array

    def suffix(self, start=0, length=-1):
        if length < 0:
            length = len(self) - start
        return self.seq[start:start+length]

    def __len__(self):
        return len(self.seq)

def find_best_subsequence_matches(seq1, seq2):
    matches = malloc(len(seq2), c_int)
    lcps = malloc(len(seq2), c_int)
    SAIS.find_best_subsequence_matches(seq1.seq, seq1.suffix_array, len(seq1),
                                       seq2.seq, seq2.suffix_array, len(seq2),
                                       matches, lcps)
    return [(matches[i], lcps[i]) for i in range(len(seq2))]

