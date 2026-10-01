# Projected Yukawa Operator and Chiral Polar Factor

This repository contains the source of the Cosmochrony fermionic-matter paper
*Projected Yukawa Operator and Chiral Polar Factor — The Squared Observable
$H_\Pi = Y_\Pi^{\dagger}Y_\Pi$ and the Undetermined Chiral Polar Map $U_\Pi$*.

This note takes the second step of the mass-sector frontier of the fermionic-matter
sub-programme, and begins with a negative result. The companion line note types the Yukawa coupling line
as the determinant line $L_Y = \wedge^2(E_{\mathrm{weak}})$ of a distinct weak factor, conditionally on the hypotheses
[H-Spin] and [H-Weak] of Q14, and reads generation levels on the model operator
$\operatorname{diag}(1, \tfrac12+u, \tfrac12-u)$ of $\mathbb{C}^3_{\mathrm{gen}}$.

Version 2.0.

## Core Result

The paper asks whether the squared projective residue determines the Yukawa morphism itself, and
shows that **it does not**. Three premises are named, none supplied by a source: [H-Res], the identification of the
model operator on the generation factor with the restriction $E_\Pi^2|_{\mathbb{C}^3_{\mathrm{gen}}}$; [H-Sq], the
dictionary $H_\Pi = Y_\Pi^{\dagger}Y_\Pi = \lambda_Y^2 E_\Pi^2|_{\mathbb{C}^3_{\mathrm{gen}}}$; and [H-Fac], the
factorisation of the full Yukawa coupling through its generation block. The results are stated for the generation block
$Y_\Pi : \mathbb{C}^3_{\mathrm{gen}}{}^{(L)} \to \mathbb{C}^3_{\mathrm{gen}}{}^{(R)}$ and transfer to the full operator
only under [H-Fac]. The polar decomposition $Y_\Pi = U_\Pi H_\Pi^{1/2}$ is linear algebra, without any of the three
premises; under [H-Res] and [H-Sq], $H_\Pi$ has the squared Yukawa levels
$\lambda_Y^2\{1, \tfrac12+u, \tfrac12-u\}$ and the **chiral polar factor $U_\Pi$ is undetermined**. The undetermined
object carries the chiral orientation that later inputs must supply; the squared observable alone cannot resolve it.
The polar factor is a rephasing class for the unitaries commuting with the Cartan generator $J_3$ of each carrier. On the
real cascade of the model the class is trivial and no mixing is produced; the real diagonal split is carried, on the
orientation-compatible branch of Q14 (an input), by the Lorentz-chirality imbalance, not by the $V{-}A$ structure of
[H-Weak]. Q14 excludes a complex metaplectic phase as a source of the external block $R_{\mathrm{mix}}$ within its
model; the internal block $e_0 \leftrightarrow e_\pm$ is neither excluded by Q14 nor constructed here.

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
```

Both scripts are exact (sympy); the second also computes one rank in 60-digit arithmetic with a fixed seed.
Their printed outputs are committed beside them in `code/`.

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
