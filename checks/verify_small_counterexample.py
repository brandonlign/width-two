#!/usr/bin/env python3
"""Exact checks for the 43-element counterexample Q of Section 4.

Q has chains A_1<...<A_22 and B_1<...<B_21 and the relations A_1<B_17 and
A_13<B_18.  The linear coefficient is computed three independent ways, all in
exact arithmetic:

  1. formula  - the descent numbers h_d of Lemma 4.2 and equation (2.1);
  2. ideals   - counting multichains of order ideals, then Newton interpolation;
  3. descents - a dynamic program over the shuffles of the two chains.
"""
from fractions import Fraction
from math import comb

from verify_counterexample import all_coefficients, linear_coefficient

M, K = 22, 21
RELS = [(1, 17), (13, 18)]                       # A_i < B_j
EXPECTED = Fraction(-2098, 3179576907)
TABLE = [1, 409, 42000, 1785800, 38729530, 481892958, 3701653032, 18420331776,
         61340152896, 139684812696, 220506645788, 243008774736, 187186137658,
         100328434948, 37022208048, 9239548416, 1517780145, 157608990, 9728950,
         323400, 4851, 22]


def C(n, k):
    return comb(n, k) if 0 <= k <= n else 0


def h_formula(d):
    total = sum(C(13, a) * C(9, d - a) * C(17, b) * C(4, d - b)
                for a in range(d + 1) for b in range(a, d + 1))
    return total - (C(9, d - 1) * C(4, d - 1) if d >= 1 else 0)


def lin_from_h(h):
    N = M + K
    return sum(Fraction((-1) ** d * v, N * comb(N - 1, d)) for d, v in enumerate(h))


def ideal_values(m, k, rels, tmax):
    """Omega(t), t = 0..tmax, as multichains of ideals (a, b) = prefix sizes."""
    ok = [[all(b < j or a >= i for i, j in rels) for b in range(k + 1)]
          for a in range(m + 1)]
    f = [[0] * (k + 1) for _ in range(m + 1)]
    f[0][0] = 1
    values = [f[m][k]]
    for _ in range(tmax):
        g = [[0] * (k + 1) for _ in range(m + 1)]
        for a in range(m + 1):
            for b in range(k + 1):
                g[a][b] = f[a][b]
                if a:
                    g[a][b] += g[a - 1][b]
                if b:
                    g[a][b] += g[a][b - 1]
                if a and b:
                    g[a][b] -= g[a - 1][b - 1]
        f = [[g[a][b] if ok[a][b] else 0 for b in range(k + 1)]
             for a in range(m + 1)]
        values.append(f[m][k])
    return values


def descent_counts(m, k, rels):
    """h_d by a DP over shuffles; a descent is a B immediately followed by an A."""
    need = [0] * (k + 1)                         # A's required before B_j
    for i, j in rels:
        for jj in range(j, k + 1):
            need[jj] = max(need[jj], i)
    state = {(0, 0, 0, 0): 1}                    # (#A, #B, last, descents)
    for _ in range(m + k):
        nxt = {}
        for (a, b, last, d), v in state.items():
            if a < m:
                key = (a + 1, b, 1, d + (last == 2))
                nxt[key] = nxt.get(key, 0) + v
            if b < k and a >= need[b + 1]:
                key = (a, b + 1, 2, d)
                nxt[key] = nxt.get(key, 0) + v
        state = nxt
    h = [0] * (m + 1)
    for (_, _, _, d), v in state.items():
        h[d] += v
    while h and h[-1] == 0:
        h.pop()
    return h


def shift_one(poly):
    """Coefficients of f(t+1) from those of f(t)."""
    out = [Fraction(0)] * len(poly)
    for i, c in enumerate(poly):
        for j in range(i + 1):
            out[j] += c * comb(i, j)
    return out


def main():
    N = M + K
    h = [h_formula(d) for d in range(K + 1)]
    assert h == TABLE == descent_counts(M, K, RELS)
    assert h_formula(K + 1) == 0
    from_formula = lin_from_h(h)

    values = ideal_values(M, K, RELS, N)
    coeffs = all_coefficients(values)
    assert from_formula == linear_coefficient(values) == coeffs[1] == EXPECTED
    assert [i for i, c in enumerate(coeffs) if c < 0] == [1]
    assert all(c >= 0 for c in shift_one(coeffs))
    print(f"Q: |Q| = {N}, [t]Omega = {EXPECTED}")
    print(f"  Omega(0..3) = {values[:4]}; formula, ideal count and descent DP agree")
    print("  t is the only negative coefficient; every coefficient of Omega(t+1) is >= 0")

    for rels in ([RELS[1]], [RELS[0]]):
        c1 = linear_coefficient(ideal_values(M, K, rels, N))
        assert c1 > 0
        print(f"  with only {rels}: [t]Omega = {c1}")


if __name__ == "__main__":
    main()
