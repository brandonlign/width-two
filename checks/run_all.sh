#!/usr/bin/env sh
set -e
cd "$(dirname "$0")"
python3 verify_counterexample.py
python3 verify_small_counterexample.py
python3 search_family.py
for k in 1 2 3 4; do python3 short_chain_certificate.py "$k"; done
python3 short_chain_independent.py 4 --multiply | tail -n 1
python3 pseudoforest_boards.py
python3 random_small_posets.py
python3 expansion_obstruction.py
