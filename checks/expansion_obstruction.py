#!/usr/bin/env python3
"""Final remark: the fixed-weight rectangle expansion fails for a width-two poset.

Two chains of length three with l = (0,0,0), u = (2,3,3). Computes W_P by
counting increasing matchings in the board and solves for its coordinates in the
basis (1-z)^r W_{3-r,3-r}(z), r = 0..3, in exact arithmetic.
"""
from fractions import Fraction
from itertools import combinations
from math import comb


def poly_mul(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, a in enumerate(p):
        for j, b in enumerate(q):
            out[i + j] += a * b
    return out


def solve(A, y):
    n = len(y)
    M = [[Fraction(A[i][j]) for j in range(n)] + [Fraction(y[i])] for i in range(n)]
    for c in range(n):
        p = next(i for i in range(c, n) if M[i][c] != 0)
        M[c], M[p] = M[p], M[c]
        for i in range(n):
            if i != c and M[i][c] != 0:
                f = M[i][c] / M[c][c]
                M[i] = [a - f * b for a, b in zip(M[i], M[c])]
    return [M[i][n] / M[i][i] for i in range(n)]


def main():
    lower, upper = (0, 0, 0), (2, 3, 3)
    cells = [(i, x) for i in range(3) for x in range(lower[i] + 1, upper[i] + 1)]
    W = [sum(1 for s in combinations(cells, j)
             if all(a[0] < b[0] and a[1] < b[1] for a, b in zip(s, s[1:])))
         for j in range(4)]
    assert W == [1, 8, 9, 1]
    basis = []
    for r in range(4):
        rect = [comb(3 - r, j) ** 2 for j in range(4 - r)]
        one_minus_z = [1]
        for _ in range(r):
            one_minus_z = poly_mul(one_minus_z, [1, -1])
        b = poly_mul(one_minus_z, rect) + [0] * 4
        basis.append(b[:4])
    coords = solve([[basis[r][j] for r in range(4)] for j in range(4)], W)
    assert coords == [Fraction(57, 60), Fraction(-5, 60), Fraction(3, 60), Fraction(5, 60)]
    print("W_P = 1 + 8z + 9z^2 + z^3; coordinates (57, -5, 3, 5)/60")


if __name__ == "__main__":
    main()
