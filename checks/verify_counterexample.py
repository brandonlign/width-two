#!/usr/bin/env python3
"""Exact checks for the counterexample, the infinite family and Proposition 3.1.

P(m,k,r,h) is the poset with chains A_1<...<A_m and B_1<...<B_k plus the
relations A_1 < B_{k-h+1} and A_r < B_k.  The linear coefficient of its order
polynomial is computed in four independent ways, all in exact arithmetic:

  1. formula   - the descent polynomial of Proposition 3.1 and the functional L_N;
  2. maps      - counting weakly order-preserving maps P -> [t] by summing over
                 the values of A_1 and A_r, then Newton interpolation;
  3. ideals    - counting multichains of order ideals, then Newton interpolation;
  4. descents  - a dynamic program over linear extensions that counts descents.
"""
from fractions import Fraction
from math import comb, factorial


# ---------------------------------------------------------------- helpers

def linear_coefficient(values):
    """[t] of the polynomial with values f(0), f(1), ..., f(D)."""
    diffs, out, k = list(values), Fraction(0), 0
    while diffs:
        if k:
            out += Fraction(diffs[0] * (-1) ** (k - 1), k)
        diffs = [b - a for a, b in zip(diffs, diffs[1:])]
        k += 1
    return out


def all_coefficients(values):
    """Monomial coefficients of the polynomial with values f(0), ..., f(D)."""
    diffs, newton = list(values), []
    while diffs:
        newton.append(diffs[0])
        diffs = [b - a for a, b in zip(diffs, diffs[1:])]
    poly = [Fraction(0)] * len(newton)
    falling = [Fraction(1)]                      # t(t-1)...(t-k+1)
    for k, d in enumerate(newton):
        for i, c in enumerate(falling):
            poly[i] += Fraction(d, factorial(k)) * c
        nxt = [Fraction(0)] * (len(falling) + 1)
        for i, c in enumerate(falling):
            nxt[i + 1] += c
            nxt[i] -= k * c
        falling = nxt
    return poly


# ---------------------------------------------------------------- method 1

def by_formula(m, k, r, h):
    N = m + k

    def L(a, b, shift):          # L_N applied to z^shift * W_{a,b}
        return sum(Fraction((-1) ** (j + shift) * comb(a, j) * comb(b, j),
                            N * comb(N - 1, j + shift))
                   for j in range(min(a, b) + 1))

    return (-sum(L(k - 1, q, 1) for q in range(r))
            - sum(L(m - 1, q, 1) for q in range(h))
            + L(0, 0, 1)
            + sum(L(i, j, 2) for i in range(h - 1) for j in range(r - 1)))


# ---------------------------------------------------------------- method 2

def count_maps(m, k, r, h, t):
    if t == 0:
        return 0
    p = k - h + 1                                # A_1 < B_p
    total = 0
    for a1 in range(1, t + 1):
        for ar in range(a1, t + 1):
            chain_a = comb(ar - a1 + r - 2, r - 2) * comb(t - ar + m - r, m - r)
            chain_b = 0
            for v in range(a1, t + 1):           # v = f(B_p)
                below = comb(v + p - 2, p - 1)
                n = k - p
                if n == 0:
                    above = 1 if v >= ar else 0
                else:
                    above = comb(t - v + n, n)
                    if ar - 1 >= v:              # remove f(B_k) < f(A_r)
                        above -= comb(ar - 1 - v + n, n)
                chain_b += below * above
            total += chain_a * chain_b
    return total


# ---------------------------------------------------------------- method 3

def count_by_ideals(m, k, r, h, tmax):
    """Omega(t) for t = 0..tmax as multichains of ideals (a, b) = prefix sizes."""
    ok = [[(b < k - h + 1 or a >= 1) and (b < k or a >= r)
           for b in range(k + 1)] for a in range(m + 1)]
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


# ---------------------------------------------------------------- method 4

def by_descents(m, k, r, h):
    """Labels list every A before every B; a descent is a B followed by an A."""
    state = {(0, 0, 0, 0): 1}                     # (#A, #B, last, descents)
    for _ in range(m + k):
        nxt = {}
        for (a, b, last, d), v in state.items():
            if a < m:
                key = (a + 1, b, 1, d + (last == 2))
                nxt[key] = nxt.get(key, 0) + v
            need = max(1 if b + 1 >= k - h + 1 else 0, r if b + 1 == k else 0)
            if b < k and a >= need:
                key = (a, b + 1, 2, d)
                nxt[key] = nxt.get(key, 0) + v
        state = nxt
    N = m + k
    return sum(Fraction((-1) ** d * v, N * comb(N - 1, d))
               for (_, _, _, d), v in state.items())


# ---------------------------------------------------------------- main

def family_closed_form(m):
    num = 10 * comb(m - 1, 4) - comb(2 * m - 6, 4) + 4 * comb(2 * m - 6, 3)
    return Fraction(num, comb(2 * m - 2, 4))


def main():
    m, k, r, h = 28, 33, 5, 7
    N = m + k
    expected = Fraction(-233, 872541560)

    maps = [count_maps(m, k, r, h, t) for t in range(N + 1)]
    ideals = count_by_ideals(m, k, r, h, N)
    assert maps == ideals
    coeffs = all_coefficients(maps)
    assert by_formula(m, k, r, h) == linear_coefficient(maps) == by_descents(m, k, r, h) == expected
    assert coeffs[1] == expected
    assert [i for i, c in enumerate(coeffs) if c < 0] == [1]
    print(f"Counterexample: P(28,33,5,7) has [t]Omega = {expected}")
    print(f"  first values Omega(0..3) = {maps[:4]}")
    print(f"  formula, map count, ideal count and descent DP agree;"
          f" t is the only negative coefficient")

    for mm in range(5, 46):
        NN = 2 * mm
        direct = NN * (NN - 1) * linear_coefficient(count_by_ideals(mm, mm, 5, 5, NN))
        assert direct == family_closed_form(mm) == NN * (NN - 1) * by_formula(mm, mm, 5, 5)
        assert (direct < 0) == (mm >= 35)
    for mm in range(46, 2001):
        assert family_closed_form(mm) < 0
    for x in range(6):                           # degree-4 numerator identity
        mm = 35 + x
        num = 10 * comb(mm - 1, 4) - comb(2 * mm - 6, 4) + 4 * comb(2 * mm - 6, 3)
        assert num == (Fraction(-x ** 4, 4) - Fraction(143 * x ** 3, 6)
                       - Fraction(3049 * x ** 2, 4) - Fraction(25001 * x, 3) - 4960)
    lim = 2 * sum(Fraction((-1) ** j * (j + 1) * comb(q, j), 2 ** j)
                  for q in range(5) for j in range(q + 1)) - 1
    assert lim == Fraction(-3, 8)
    print("Infinite family: closed form matches direct counts for 5 <= m <= 45;"
          " negative exactly for m >= 35 (checked to m = 2000); limit -3/8")

    for (mm, kk, rr, hh) in [(30, 31, 6, 6), (35, 35, 5, 5), (50, 50, 5, 5)]:
        NN = mm + kk
        assert by_formula(mm, kk, rr, hh) == linear_coefficient(
            count_by_ideals(mm, kk, rr, hh, NN)) == by_descents(mm, kk, rr, hh)
        print(f"  P{(mm, kk, rr, hh)}: [t]Omega = {by_formula(mm, kk, rr, hh)}")


if __name__ == "__main__":
    main()
