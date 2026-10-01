# Projected Yukawa Operator and Chiral Polar Factor

This repository contains the source of the Cosmochrony fermionic-matter paper
*Projected Yukawa Operator and Chiral Polar Factor — The Squared Observable
$H_\Pi = Y_\Pi^{\dagger}Y_\Pi$ and the Undetermined Chiral Polar Map $U_\Pi$*.

This note takes the second step of the mass-sector frontier of the fermionic-matter
sub-programme, and draws a **sharp negative conclusion** at its entry. The companion line note
types the Yukawa coupling line as the determinant line $L_Y = \wedge^2(E_{\mathrm{weak}})$ of a distinct weak
doublet, conditionally on the hypotheses [H-Spin] and [H-Weak] of Q14, and reads the generation-level assignment
$E_\Pi^2|_{\mathrm{gen}} = \operatorname{diag}(1, \tfrac12+u, \tfrac12-u)$.

Version 2.0.

## Core Result

The paper asks whether the squared projective residue determines the Yukawa morphism itself, and
shows that **it does not**. The generation-level observable produced by the Schur-residue sector
is, under a named identification (the operators on the generation factor are restrictions of
$E_\Pi^2$, not supplied), the positive hermitian square $H_\Pi = Y_\Pi^{\dagger}Y_\Pi$, which fixes the
singular values of the Yukawa map but leaves the **chiral polar factor $U_\Pi$ undetermined**. The undetermined
object carries the chiral orientation that later inputs must supply; the squared observable alone
cannot resolve it. The polar-factor statements are algebraic consequences using only the generation factor and
positive hermitian squares; the spinor, weak and line factors are spectators. The absence of mixing on the real
cascade reads the diagonal split as the Lorentz-chirality (left-admissibility) imbalance on the
orientation-compatible branch, an input, and not as $V{-}A$ selection.

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
