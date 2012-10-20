import sys
import time
import sais


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print "Usage: %s [SEQFILE]" % sys.argv[0]
        sys.exit()

    seq = sais.Sequence(open(sys.argv[1]).read())
    startTime = time.time()
    print seq.suffix_array[0] + seq.lcp_array[0]
    endTime = time.time()
    print 'Construction took %s ms' % ((endTime - startTime) * 1000.0)
