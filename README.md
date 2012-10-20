This is a simple Python ctypes-based wrapper around [this speedy suffix array construction implementation in C.](https://sites.google.com/site/yuta256/sais)  It also includes an implementation of Kasai's method for constructing the longest common prefix array from a suffix array in linear time.

Building:
---------
```
cd sais/sais-lite
make dll
```

Basic Usage:
------------
```python
import sais

seq = sais.Sequence('banana$')

print seq.suffix_array[:len(myseq)]   # -> [6, 5, 3, 1, 0, 4, 2]
print seq.lcp_array[:len(myseq)]      # -> [-1, 0, 1, 3, 0, 0, 2]

def longest_repeated_substring(seq):
    idx, lcp = max(seq.lcp_array[:len(myseq)], key=lambda (_, v): v)
    suffix_start = seq.suffix_array[idx]
    suffix_end = suffix_start + lcp
    return seq.seq[suffix_start:suffix_end]

print longest_repeated_substring(seq)     # -> 'ana'
```

References:
-----------

* [Suffix arrays](http://en.wikipedia.org/wiki/Suffix_array)
* [Longest common prefix](http://en.wikipedia.org/wiki/LCP_array)