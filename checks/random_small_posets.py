#!/usr/bin/env python3
"""Brute-force tests of the positive theorems on random small width-two posets.

Each poset is two chains plus random cross relations, transitively closed. Its
order polynomial is computed exactly from multichains of order ideals, and every
coefficient is checked in the three positive classes:
  series-parallel (comparability graph has no induced P4),
  short chain (one chain has at most four elements),
  pseudoforest (each incomparability component has at most one cycle).
The strictness criterion (linear coefficient > 0 iff the comparability graph is
connected) is checked in the series-parallel and pseudoforest classes. Also
checks the order polynomial of the complement zigzag Z_4 and nonnegativity for
complement zigzags with N <= 40.
"""
import random
from fractions import Fraction
from math import factorial


def coefficients(values):
    diffs, newton = list(values), []
    while diffs:
        newton.append(diffs[0])
        diffs = [b - a for a, b in zip(diffs, diffs[1:])]
    poly = [Fraction(0)] * len(newton)
    falling = [Fraction(1)]
    for k, d in enumerate(newton):
        for i, c in enumerate(falling):
            poly[i] += Fraction(d, factorial(k)) * c
        nxt = [Fraction(0)] * (len(falling) + 1)
        for i, c in enumerate(falling):
            nxt[i + 1] += c
            nxt[i] -= k * c
        falling = nxt
    return poly


def order_polynomial(n, lt, A, B):
    """Exact coefficients of Omega for a poset covered by chains A and B."""
    ideals = []
    for i in range(len(A) + 1):
        for j in range(len(B) + 1):
            members = set(A[:i] + B[:j])
            if all(y in members for x in members for y in range(n) if lt[y][x]):
                ideals.append(frozenset(members))
    supersets = {I: [J for J in ideals if I <= J] for I in ideals}
    full = frozenset(range(n))
    vec, values = {frozenset(): 1}, [1 if n == 0 else 0]
    for _ in range(n):
        new = {}
        for I, c in vec.items():
            for J in supersets[I]:
                new[J] = new.get(J, 0) + c
        vec = new
        values.append(vec.get(full, 0))
    return coefficients(values)


def closure(n, rel):
    lt = [[False] * n for _ in range(n)]
    for a, b in rel:
        lt[a][b] = True
    for x in range(n):
        for a in range(n):
            if lt[a][x]:
                for b in range(n):
                    if lt[x][b]:
                        lt[a][b] = True
    return None if any(lt[i][i] for i in range(n)) else lt


def random_poset(rng, m, k, p):
    n = m + k
    A, B = list(range(m)), list(range(m, n))
    rel = [(A[i], A[i + 1]) for i in range(m - 1)] + [(B[i], B[i + 1]) for i in range(k - 1)]
    for a in A:
        for b in B:
            x = rng.random()
            if x < p:
                rel.append((a, b))
            elif x < 2 * p:
                rel.append((b, a))
    lt = closure(n, rel)
    return (n, lt, A, B) if lt else None


def connected(n, adj):
    seen, stack = {0}, [0]
    while stack:
        v = stack.pop()
        for w in range(n):
            if adj(v, w) and w not in seen:
                seen.add(w)
                stack.append(w)
    return len(seen) == n


def series_parallel(n, comp):
    for a in range(n):
        for b in range(n):
            if a == b or not comp[a][b]:
                continue
            for c in range(n):
                if c in (a, b) or not comp[b][c] or comp[a][c]:
                    continue
                for d in range(n):
                    if d not in (a, b, c) and comp[c][d] and not comp[a][d] and not comp[b][d]:
                        return False
    return True


def pseudoforest(n, edges):
    parent = list(range(n))

    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for a, b in edges:
        parent[find(a)] = find(b)
    vertices, count = {}, {}
    for v in range(n):
        vertices[find(v)] = vertices.get(find(v), 0) + 1
    for a, b in edges:
        count[find(a)] = count.get(find(a), 0) + 1
    return all(count.get(c, 0) <= vertices[c] for c in vertices)


def main():
    rng = random.Random(1)
    tested = {"series-parallel": 0, "short chain": 0, "pseudoforest": 0, "none": 0}
    for _ in range(4000):
        m, k = rng.randint(1, 9), rng.randint(1, 9)
        if m + k > 14:
            continue
        got = random_poset(rng, m, k, rng.choice([0.02, 0.05, 0.1, 0.2]))
        if not got:
            continue
        n, lt, A, B = got
        coeffs = order_polynomial(n, lt, A, B)
        comp = [[lt[i][j] or lt[j][i] for j in range(n)] for i in range(n)]
        inc = [(i, j) for i in range(n) for j in range(i + 1, n) if not comp[i][j]]
        conn = connected(n, lambda v, w: comp[v][w])
        classes = []
        if series_parallel(n, comp):
            classes.append("series-parallel")
        if min(m, k) <= 4:
            classes.append("short chain")
        if pseudoforest(n, inc):
            classes.append("pseudoforest")
        for c in classes or ["none"]:
            tested[c] += 1
        if classes:
            assert min(coeffs) >= 0, (classes, m, k, coeffs)
        if "series-parallel" in classes or "pseudoforest" in classes:
            assert (coeffs[1] > 0) == conn
    print("random width-two posets with N <= 14, all coefficients nonnegative in each class:")
    for c, v in tested.items():
        print(f"  {c}: {v}")

    def zigzag(N):
        lt = [[j - i >= 2 for j in range(N)] for i in range(N)]
        A, B = list(range(0, N, 2)), list(range(1, N, 2))
        return order_polynomial(N, lt, A, B)

    assert zigzag(4) == [0, Fraction(2, 24), Fraction(7, 24), Fraction(10, 24), Fraction(5, 24)]
    for N in range(4, 41):
        c = zigzag(N)
        assert all(x > 0 for x in c[1:]), N
    print("complement zigzags: Omega(Z_4) = (5t^4+10t^3+7t^2+2t)/24;"
          " coefficients in degrees 1..N positive for 4 <= N <= 40")


if __name__ == "__main__":
    main()
