#!/usr/bin/env python3
"""Short-chain theorem: binomial-basis certificate for k = 1, 2, 3, 4.

For each ordering (chamber) of the 2k interval endpoints, the quantity
(m+1) N (N-1)_k [t]Omega is a polynomial in the 2k+1 gap sizes. Its coefficients
in the multivariate binomial basis are recovered by exact binomial inversion and
checked to be nonnegative. Usage: python3 short_chain_certificate.py k
"""
from itertools import combinations,combinations_with_replacement
from math import comb,factorial,prod
from functools import lru_cache
import sys
k=int(sys.argv[1])
def solve():
 for L in combinations(range(2*k),k):
  U=[i for i in range(2*k) if i not in L]
  if any(l>u for l,u in zip(L,U)):continue
  @lru_cache(None)
  def f(s):
   cuts=[sum(s[:i+1]) for i in range(2*k)]
   m=sum(s); n=m+k
   dp={(0,0):1}
   for l,u in zip(L,U):
    nex=dp.copy()
    for (j,x),v in dp.items():
     for y in range(max(x,cuts[l])+1,cuts[u]+1):
      nex[j+1,y]=nex.get((j+1,y),0)+v
    dp=nex
   return (m+1)*sum((-1)**j*factorial(j)*v*prod(range(n-k,n-j)) for (j,x),v in dp.items())
  coeff={}
  for d in range(k+2):
   for t in combinations_with_replacement(range(2*k+1),d):
    s=tuple(t.count(i) for i in range(2*k+1))
    val=f(s)
    from itertools import product
    for r in product(*(range(i+1) for i in s)):
     if r!=s: val-=coeff[r]*prod(comb(a,b) for a,b in zip(s,r))
    coeff[s]=val
    if val<0: raise AssertionError((L,t,val))
 print('verified',k,flush=True)
solve()
