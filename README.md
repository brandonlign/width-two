# Order polynomials of width-two posets can have negative coefficients

Brandon Li

Paper and verification code disproving Conjecture 7.1 of Ferroni, Morales and
Panova ([arXiv:2503.16403](https://arxiv.org/abs/2503.16403)), which asserts
that every poset of width two has an order polynomial with nonnegative
coefficients.

Let $P$ consist of two chains $A_1<\dots<A_{28}$ and $B_1<\dots<B_{33}$ with
the two extra relations $A_1<B_{27}$ and $A_5<B_{33}$. Then $P$ has width two,
61 elements, and

$$
[t]\,\Omega_P(t)=-\frac{233}{872541560}.
$$

The linear coefficient is the only negative one. In the skew-shape language of
the conjecture, $P$ is $R_{\lambda/\mu}$ with $\lambda=(33^{28})$ and
$\mu=(7,1,1,1,1)$. Removing the hook $(5,1,1,1,1)$ from an $m\times m$ square
gives a counterexample for every $m\ge35$, and within this two-relation family
61 elements is the minimum.

The paper also proves nonnegativity for three classes of width-two posets:
series-parallel posets, posets covered by two chains one of which has at most
four elements, and posets whose incomparability graph has at most one cycle in
each component. For the first and third classes, the linear coefficient is
positive exactly when the comparability graph is connected.

- [Paper (PDF)](paper/main.pdf)
- [LaTeX source](paper/)
- [Verification code](checks/)

## Reproduce

```sh
sh checks/run_all.sh
```

The scripts use only the Python standard library and exact rational
arithmetic. `verify_counterexample.py` computes the linear coefficient of the
counterexample in four independent ways and checks the infinite family.
`search_family.py` searches the two-relation family for the smallest negative
member. The remaining scripts check the short-chain certificate, the
pseudoforest boards, random small posets in each positive class, and the
expansion obstruction in the final remark. Their output is in
[`checks/output.txt`](checks/output.txt). The PDF builds with Tectonic
(`cd paper && tectonic main.tex`) or `pdflatex`.

## License

Code: MIT. Paper: CC BY 4.0.
