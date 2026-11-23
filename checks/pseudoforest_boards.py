#!/usr/bin/env python3
"""Pseudoforest theorem: exact check of (N-3)h_2 - 3h_3 >= C(N-2,2) and of the
strictness criterion on every board of size <= 10 whose components are unicyclic."""
from fractions import Fraction
from itertools import combinations, combinations_with_replacement
from math import comb


def all_unicyclic(m, k, cells):
    adj = [set() for _ in range(m+k)]
    for i, j in cells:
        adj[j-1].add(m+i)
        adj[m+i].add(j-1)
    unseen = set(range(m+k))
    while unseen:
        seed = unseen.pop()
        stack = [seed]
        vertices = degree_sum = 0
        while stack:
            v = stack.pop()
            vertices += 1
            degree_sum += len(adj[v])
            new = adj[v] & unseen
            stack.extend(new)
            unseen.difference_update(new)
        if degree_sum != 2*vertices:
            return False
    return True


def increasing(s):
    s = sorted(s)
    return all(i0 < i1 and j0 < j1
               for (i0,j0),(i1,j1) in zip(s,s[1:]))


def main():
    count = 0
    for n in range(4,11):
        for k in range(2,n//2+1):
            m = n-k
            endpoints = list(combinations_with_replacement(range(m+1), k))
            for lower in endpoints:
                for upper in endpoints:
                    if any(l>u for l,u in zip(lower,upper)):
                        continue
                    if sum(u-l for l,u in zip(lower,upper)) != n:
                        continue
                    cells = [(i,j) for i,(l,u) in enumerate(zip(lower,upper))
                             for j in range(l+1,u+1)]
                    if not all_unicyclic(m,k,cells):
                        continue
                    h = [sum(increasing(s) for s in combinations(cells,j))
                         for j in range(k+1)]
                    excess = (n-3)*h[2] -3*(h[3] if k>=3 else 0) - comb(n-2,2)
                    scaled_linear = sum(Fraction((-1)**j*h[j],comb(n-1,j))
                                        for j in range(k+1))
                    assert excess >= 0
                    assert (scaled_linear > 0) == (n > 4)
                    count += 1
    print(f"Checked {count} exact unicyclic-union boards of total size <= 10.")


if __name__ == "__main__":
    main()
