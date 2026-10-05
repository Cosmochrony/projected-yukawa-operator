# Projected Yukawa Operator and Chiral Polar Factor

This repository contains the source of the Cosmochrony fermionic-matter paper
*Projected Yukawa Operator and Chiral Polar Factor — The Squared Observable
$H_\Pi = (Y^{\mathrm{gen}}_\Pi)^{\dagger}Y^{\mathrm{gen}}_\Pi$ and the Undetermined Chiral Polar Map $U_\Pi$*.

This note takes the second step of the mass-sector frontier of the fermionic-matter
sub-programme, and begins with a negative result. The companion line note types the Yukawa coupling line
as the determinant line $L_Y = \wedge^2(E_{\mathrm{weak}})$ of a distinct weak factor, conditionally on the hypothesis
[H-Weak] of Q14 (the fermion sectors need [H-Spin] as well), and reads generation levels on the model operator
$\operatorname{diag}(1, \tfrac12+u, \tfrac12-u)$ of $\mathbb{C}^3_{\mathrm{gen}}$.

Version 2.1.

## Core Result

The paper asks whether the squared projective residue determines the Yukawa morphism itself, and
shows that **it does not**. Four premises are named, none supplied by a source:
[H-Res], the identification of the model operator on the generation factor with the restriction
$E_\Pi^2|_{\mathbb{C}^3_{\mathrm{gen}}}$;
[H-Sq], the dictionary $H_\Pi = (Y^{\mathrm{gen}}_\Pi)^{\dagger}Y^{\mathrm{gen}}_\Pi
= \lambda_Y^2 E_\Pi^2|_{\mathbb{C}^3_{\mathrm{gen}}}$;
[H-Fac], made of two separate premises, the existence of a weak linking carrier and the multiplicity at most one of the
invariant couplings, under which the full coupling $Y_\Pi$ factorises through its generation block;
and [H-Grad], the $J_3$ grading of each generation carrier.
The results are stated for the generation block
$Y^{\mathrm{gen}}_\Pi : \mathbb{C}^3_{\mathrm{gen}}{}^{(L)} \to \mathbb{C}^3_{\mathrm{gen}}{}^{(R)}$
and transfer to $Y_\Pi$ only under [H-Fac].
The polar decomposition $Y^{\mathrm{gen}}_\Pi = U_\Pi H_\Pi^{1/2}$ is linear algebra, without any of the
premises; under [H-Res] and [H-Sq], $H_\Pi$ has the squared Yukawa levels
$\lambda_Y^2\{1, \tfrac12+u, \tfrac12-u\}$ and the **chiral polar factor $U_\Pi$ is undetermined**.
The undetermined object carries the chiral orientation that later inputs must supply; the squared observable
alone cannot resolve it.
Under [H-Grad] the polar factor is a rephasing class for the unitaries commuting with the Cartan generator $J_3$ of
each carrier, non-trivial at first order exactly when its polar generator $\Omega$ has a transverse part.
For a smooth real-parameter family $Y(\gamma) = U(\gamma)P(\gamma)$ of invertible morphisms between fixed orthonormal
frames, $P > 0$, the polar generator is $\Omega = U^\dagger \partial_\gamma U$ and it solves the Sylvester equation
$\Omega P + P\Omega = U^\dagger \partial_\gamma Y - (\partial_\gamma Y)^\dagger U$, uniquely since $P > 0$
(singular supports and moving frames need a separate treatment).
The projection $A^{\mathrm{odd}}_\Pi$ of the $\mathfrak{sl}_2$ lift onto the $J_\Pi$-odd anti-hermitian part vanishes
identically for every complex coefficient, because with the internal antilinear parity $J_\Pi$ of the generation copy
(Q14 Section 6) the $J_\Pi$-odd part of the lift is its Hermitian part; this is a statement about
$A^{\mathrm{odd}}_\Pi$ as defined, and $A^{\mathrm{odd}}_\Pi$ does not determine the polar factor.
It is not the anomaly density of Q14, and it is not the polar generator: at the identity, for
$Y(0) = 1$ and $\partial_\gamma Y(0) = L(M)$, $\Omega(0) = \tfrac12(L(M) - L(M)^\dagger)$,
the $J_\Pi$-even anti-Hermitian part of the lift, non-zero in general.
For the real cascade step $g = \exp(tE)\exp(sF)$, the unitary polar factor of $\mathrm{Sym}^2 g$ (not a Yukawa
morphism; the step is not at the identity) has the external entry $(e_+, e_-)$ equal to
$(t-s)^2/((2+ts)^2+(t-s)^2)$, that is $1/5$ at $(1, 0)$ and $1/17$ at $(2, 1)$, and the Jarlskog invariant of the real
factor is zero.
The absence of an external entry in an infinitesimal generator does not exclude one in the finite factor.
For fixed $D > 0$, $Y_0 = D$ and $Y_1 = VD$ have the same square for every unitary $V$: the square does not determine
the polar data.
The corpus supplies no family of Yukawa morphisms whose polar factor is $U_\Pi$: whether such a family exists, and which
polar class it has, is open (a missing element, not a refutation).
No source justifies the selection of the $J_\Pi$-odd part.
The $J_\Pi$-even anti-Hermitian part of the lift has a non-zero internal entry $\tfrac{\sqrt2}{2}(q - \bar p)$, non-zero
for real $p \neq q$, so the image of $\mathfrak{sl}_2$ does reach the internal block $e_0 \leftrightarrow e_\pm$
through that part, while the external block $R_{\mathrm{mix}}$, along which PRS expands $E_\Pi^2$, has no matrix
element in the image of $\mathfrak{sl}_2$ (Q14 Remark 6.4, a statement on matrix elements of the lift);
Q14 Proposition 6.3 (i)-(ii) excludes the internal block only by hypothesis, for the oriented odd lift it considers.
For the cascade step $g = \exp(tE)\exp(sF)$ with real $t$, $s$, the $J_\Pi$-even internal entry is
$\tfrac{\sqrt2}{2}\,c\,(s - t)$ with $c = \theta/\sinh\theta$: it vanishes at $t = s$, where the whole
anti-Hermitian part of the lift is zero.
No physical exclusion of the internal block is claimed; the vanishing of $A^{\mathrm{odd}}_\Pi$ is
representation-theoretic (the $J_\Pi$-odd anti-Hermitian operators are spin 0 and spin 2, the $\mathfrak{sl}_2$
image is spin 1).
On the orientation-compatible branch of Q14 (an input), in the framework of PRS a non-zero diagonal split requires a locking operator that
breaks $J_\Pi$-equivariance (in the framework of PRS, under its chiral splitting and [H-Gen](i), a $J_\Pi$-commuting locking gives $E[P] = 0$ on
that branch, which concerns Q14's $E_\Pi$ only when $\Pi_S D^2 \Pi_S^* = \mathrm{Lich}$); the loss of
$J_\Pi$-compatibility alone neither builds a split nor its spectral reading, and no mechanism for it is derived.
The branch is not the $V{-}A$ structure of [H-Weak].

## Keywords

projected Yukawa operator; hermitian square; polar decomposition; chiral polar factor; polar class;
rephasing invariants; Jarlskog phase; generation mixing; generation levels; partial no-go; Schur residue;
CP-real branch; non-injective projection.

## Repository Contents

```
projected-yukawa-operator/
├── code/        # Diagnostic scripts (+ requirements.txt)
├── tex/         # LaTeX sources (main + cosmochrony-bibliography.bib)
├── out/         # Compiled paper PDF (ProjectedYukawaOperator.pdf)
├── zenodo.json  # Zenodo deposition metadata
└── README.md
```

## Reproduction

```
pip install -r code/requirements.txt
python code/front3b_yukawa_operator.py
python code/front_polar_class_rephasing.py
python code/front_reality_structure_sl2_lift.py
python code/front_polar_generator.py
python code/front_hfac_multiplicity.py
python code/front_level_assignment_nogo.py
python code/front_depths_ng_reconciliation.py
```

The scripts combine exact arithmetic with delimited numerical checks.
The symbolic checks use SymPy; the level-assignment script uses exact rational levels (`Fraction`), floating-point
square roots, and a `1e-12` tolerance for its root-ratio check.
The multiplicity script uses numerical Weyl integration, with comparison thresholds `1e-9` and `1e-12`;
the rephasing script runs random-point
rank tests in 60-digit mpmath arithmetic with a fixed seed.
Step 2 of `front_depths_ng_reconciliation.py` includes floating-point ratios, logarithms and depth estimates.
`front3b_yukawa_operator.py` uses floating-point sorting keys and compares the sorted eigenvalues exactly.
Their printed outputs are committed beside them in `code/`.
The depth-bracketing audits recorded in the paper are not reproduced in this repository.

## Links

- 🔗 DOI: [10.5281/zenodo.20767498](https://doi.org/10.5281/zenodo.20767498)
- 🌐 Website: https://cosmochrony.org/science/fermionic-matter/projected-yukawa-operator/

## Citation

> J. Beau, *Projected Yukawa Operator and Chiral Polar Factor*, Zenodo, 2026.
> DOI: 10.5281/zenodo.20767498.

## Acknowledgements

Portions of the editorial refinement benefited from iterative interactions with large
language models, used as analytical assistants. All claims and final formulations remain
the sole responsibility of the author.
