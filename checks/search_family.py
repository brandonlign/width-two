#!/usr/bin/env python3
"""Remark 3.2: smallest size of a negative member of the family P(m,k,r,h).

Uses the exact formula of Proposition 3.1, after checking it against direct
counts of order-preserving maps on random small parameters.
"""
import random
from verify_counterexample import by_formula, count_by_ideals, linear_coefficient


def main():
    rng = random.Random(3)
    for _ in range(200):
        m, k = rng.randint(2, 12), rng.randint(2, 12)
        r, h = rng.randint(2, m), rng.randint(2, k)
        N = m + k
        assert by_formula(m, k, r, h) == linear_coefficient(count_by_ideals(m, k, r, h, N))
    print("Proposition 3.1 formula matches direct counts on 200 random parameter sets")

    for N in range(4, 62):
        found = [(m, N - m, r, h, by_formula(m, N - m, r, h))
                 for m in range(2, N - 1)
                 for r in range(2, m + 1)
                 for h in range(2, N - m + 1)
                 if by_formula(m, N - m, r, h) < 0]
        if found:
            print(f"smallest size with a negative linear coefficient: {N}")
            for f in found:
                print(f"  P{f[:4]}: {f[4]}")
            assert N == 61
            return
    raise AssertionError("no negative member found up to 61")


if __name__ == "__main__":
    main()
