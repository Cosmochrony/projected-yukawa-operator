# Projected Yukawa Operator and Chiral Polar Factor

This repository contains the source of the Cosmochrony fermionic-matter paper
*Projected Yukawa Operator and Chiral Polar Factor — The Squared Observable
$H_\Pi = (Y^{\mathrm{gen}}_\Pi)^{\dagger}Y^{\mathrm{gen}}_\Pi$ and the Undetermined Chiral Polar Map $U_\Pi$*.

This note takes the second step of the mass-sector frontier of the fermionic-matter
sub-programme, and begins with a negative result. The companion line note types the Yukawa coupling line
as the determinant line $L_Y = \wedge^2(E_{\mathrm{weak}})$ of a distinct weak factor, conditionally on the hypotheses
[H-Spin] and [H-Weak] of Q14, and reads generation levels on the model operator
$\operatorname{diag}(1, \tfrac12+u, \tfrac12-u)$ of $\mathbb{C}^3_{\mathrm{gen}}$.

Version 2.0.

## Core Result

The paper asks whether the squared projective residue determines the Yukawa morphism itself, and
shows that **it does not**. Four premises are named, none supplied by a source:
[H-Res], the identification of the model operator on the generation factor with the restriction
$E_\Pi^2|_{\mathbb{C}^3_{\mathrm{gen}}}$;
[H-Sq], the dictionary $H_\Pi = (Y^{\mathrm{gen}}_\Pi)^{\dagger}Y^{\mathrm{gen}}_\Pi
= \lambda_Y^2 E_\Pi^2|_{\mathbb{C}^3_{\mathrm{gen}}}$;
[H-Fac], made of two separate premises, the existence of a weak linking carrier and the multiplicity one of the
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
each carrier, non-trivial exactly when its generator has a transverse part.
For the $J_\Pi$-odd anti-hermitian part $A_\Pi$ of the $\mathfrak{sl}_2$ lift, defined in the paper, the transverse
part vanishes identically for every complex coefficient: neither the internal block $e_0 \leftrightarrow e_\pm$ nor
the external block $R_{\mathrm{mix}}$ is sourced by the lift, and a source needs a generator in the spin-two sector.
The real diagonal split is carried, on the orientation-compatible branch of Q14 (an input), by the
Lorentz-chirality imbalance, not by the $V{-}A$ structure of [H-Weak].

## Keywords

Projected Yukawa operator, chiral polar factor, squared observable, Schur residue, mass sector,
generation levels, hermitian square.

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
python code/front_hfac_multiplicity.py
python code/front_level_assignment_nogo.py
python code/front_depths_ng_reconciliation.py
```

The scripts are exact (sympy); the rephasing script also computes one rank in 60-digit arithmetic with a fixed seed,
and the multiplicity script uses numerical Weyl integration.
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
