/*
 * sais.h for sais-lite
 * Copyright (c) 2008-2010 Yuta Mori All Rights Reserved.
 *
 * Permission is hereby granted, free of charge, to any person
 * obtaining a copy of this software and associated documentation
 * files (the "Software"), to deal in the Software without
 * restriction, including without limitation the rights to use,
 * copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the
 * Software is furnished to do so, subject to the following
 * conditions:
 *
 * The above copyright notice and this permission notice shall be
 * included in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
 * EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
 * OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
 * NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
 * HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
 * WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
 * FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
 * OTHER DEALINGS IN THE SOFTWARE.
 */

#ifndef _SAIS_H
#define _SAIS_H 1

#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

#include <sys/queue.h>

/* find the suffix array SA of T[0..n-1]
   use a working space (excluding T and SA) of at most 2n+O(lg n) */
int
sais(const unsigned char *T, int *SA, int n);

/* find the suffix array SA of T[0..n-1] in {0..k-1}^n
   use a working space (excluding T and SA) of at most MAX(4k,2n) */
int
sais_int(const int *T, int *SA, int n, int k);

/* burrows-wheeler transform */
int
sais_bwt(const unsigned char *T, unsigned char *U, int *A, int n);
int
sais_int_bwt(const int *T, int *U, int *A, int n, int k);

/* longest common prefix generation */
int
compute_lcp(const unsigned char * T, const int *SA, int *LCP, int n);

//int compute_lcp_int(const int *SA, int *LCP, int start, int lengt


typedef struct Run {
    int start;
    int length;
    STAILQ_ENTRY(Run) list_entry;
} Run;

typedef STAILQ_HEAD(RunList, Run) RunList;

RunList *
findCommonSeqsByLength(const int *LCP, int n);

void testCommonSeqs();

void
find_best_subsequence_matches(const unsigned char *S1, const int *SA1, int n1, 
                              const unsigned char *S2, const int *SA2, int n2,
                              int *bestMatchIndices, int *bestMatchLCPs);

#ifdef __cplusplus
} /* extern "C" */
#endif /* __cplusplus */

#endif /* _SAIS_H */
