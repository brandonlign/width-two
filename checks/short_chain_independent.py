#!/usr/bin/env python3
"""Short-chain theorem: independent certificate.

Recomputes h_j by listing every linear extension and counting its descents
under the earliest-extension labeling, then repeats the binomial inversion.
Usage: python3 short_chain_independent.py 4 --multiply
"""
from itertools import combinations, combinations_with_replacement, product
from functools import lru_cache
from math import comb, factorial, prod
import sys

k = int(sys.argv[1])
multiplied = '--multiply' in sys.argv
dmax = k + multiplied
for L in combinations(range(2*k), k):
    U = sorted(set(range(2*k)) - set(L))
    if any(l > u for l, u in zip(L,U)):
        continue

    @lru_cache(None)
    def value(s):
        cuts = [sum(s[:i+1]) for i in range(2*k)]
        m = sum(s)
        total = 0
        for xs in combinations_with_replacement(range(m+1),k):
            if not all(cuts[l] <= x <= cuts[u] for l,x,u in zip(L,xs,U)):
                continue
            # Sort the full extension by the execution order, then read its
            # natural labels as keys supplied by the earliest extension.
            word = [(x,1,i,(cuts[L[i]],1,i)) for i,x in enumerate(xs)]
            word += [(j,0,0,(j,0,0)) for j in range(1,m+1)]
            labels = [row[3] for row in sorted(word)]
            j = sum(a > b for a,b in zip(labels,labels[1:]))
            total += (-1)**j * factorial(j) * prod(range(m,m+k-j))
        return total*(m+1) if multiplied else total

    coefficients = []
    by_degree = []
    for d in range(dmax+1):
        block = []
        for r in combinations_with_replacement(range(2*k+1), d):
            s = [r.count(i) for i in range(2*k+1)]
            entry = sum(
                (-1)**(d-sum(v))*prod(comb(a,b) for a,b in zip(s,v))*value(v)
                for v in product(*(range(a+1) for a in s)))
            block.append(entry)
        by_degree.append(min(block))
        coefficients += block
    print(L,len(coefficients),min(coefficients),max(coefficients),coefficients.count(0),by_degree)
    assert min(coefficients) >= 0
