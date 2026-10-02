"""Representation theory behind the premise [H-Fac] (weak linking carrier K, multiplicity at most one).

Setting. G = Spin(3,1) x U(2) = SL(2,C) x U(2). Left sector F_L = P_L S (x) (E (x) det^k), right sector F_R = P_R S (x)
det^m, E the weak doublet of U(2), L_Y = det. G acts trivially on the generation factor. The coupling is a G-invariant
sesquilinear pairing, antilinear in F_L, linear in F_R, with a Lorentz-trivial U(2)-module K (the weak linking carrier):
I_K = (conj(F_L) (x) F_R (x) K)^G.  Claim of the paper: for K Lorentz-trivial,
    dim I_K = multiplicity of E (x) det^(k-m) in K.
Checks.
  1  Lorentz: (0,1/2) (x) (0,1/2) = (0,0) + (0,1) (SU(2) character integration for the compact form): the scalar occurs
     once, so the Lorentz factor of the sesquilinear pairing is not an obstruction.
  2  (1/2,0) and (0,1/2) are inequivalent SL(2,C) modules: their characters tr g and conj(tr g) differ at an
     explicit g = diag(1+i, 1/(1+i)); hence the LINEAR Hom_G(P_L S, P_R S) is zero before any weak factor.
  3  Weak side, Weyl integration on U(2) (grid, exact to rounding for these Laurent polynomials): for K = E (x) det^a
     and k, m, a in {-2..2}, dim I_K (weak part) = 1 iff a = k - m, else 0 (125 cases).
  4  K a U(2)-character det^a: the invariant space is zero for all k, m, a in {-2..2} (125 cases): a doublet component
     in K is required, in the sesquilinear convention.
  5  Linear convention (control): Hom_{U(2)}(E (x) det^k (x) K, det^m) is non-zero for K = E (x) det^a iff
     a = m - k - 1, so a doublet in K is required in both conventions (the twist is convention dependent).
  6  Multiplicity: for K = (E (x) det^(k-m)) + (E (x) det^(k-m)) the dimension is 2 (F1 holds, F2 fails); for K = 0 it is 0
     (F2 holds, F1 fails): the two premises are independent.
Exit status 0 iff every check passes.
"""

from itertools import product

import numpy as np

N = 64
TH = 2 * np.pi * np.arange(N) / N
T1, T2 = np.meshgrid(TH, TH, indexing="ij")
Z1, Z2 = np.exp(1j * T1), np.exp(1j * T2)
W = np.abs(Z1 - Z2) ** 2 / 2            # Weyl density of U(2) in eigenphase coordinates
CHI_E = Z1 + Z2
DET = Z1 * Z2


def integrate(ch):
    return (ch * W).mean()


checks = {}

# 1 Lorentz scalar multiplicity in (0,1/2)(x)(0,1/2): SU(2) characters on the compact form
tt = np.linspace(0.0, 2 * np.pi, 4096, endpoint=False)
chi_half = 2 * np.cos(tt / 2)
haar = (1 / np.pi) * np.sin(tt / 2) ** 2   # SU(2) Haar density in the class angle t = 2*theta in [0, 2pi)
dt = tt[1] - tt[0]
mult_scalar = np.sum(chi_half * chi_half * haar) * dt
mult_vector = np.sum(chi_half * chi_half * (1 + 2 * np.cos(tt)) * haar) * dt   # chi_1 = 1 + 2 cos t
checks["1_lorentz_scalar_multiplicity_one"] = abs(mult_scalar - 1) < 1e-9
checks["1_lorentz_spin_one_multiplicity_one"] = abs(mult_vector - 1) < 1e-9

# 2 inequivalence of (1/2,0) and (0,1/2) for SL(2,C)
a = 1 + 1j
g = np.diag([a, 1 / a])
checks["2_characters_differ"] = abs(np.trace(g) - np.conj(np.trace(g))) > 1e-9

rng = range(-2, 3)

# 3 sesquilinear convention, K = E (x) det^a
ok3 = True
for k, m, aa in product(rng, repeat=3):
    d = integrate(np.conj(CHI_E * DET**k) * DET**m * (CHI_E * DET**aa))
    ok3 = ok3 and abs(d - (1 if aa == k - m else 0)) < 1e-9
checks["3_doublet_K_multiplicity_formula_125_cases"] = ok3

# 4 characters: zero
ok4 = all(abs(integrate(np.conj(CHI_E * DET**k) * DET**m * DET**aa)) < 1e-9 for k, m, aa in product(rng, repeat=3))
checks["4_character_K_gives_zero_125_cases"] = ok4

# 5 linear convention control: Hom(E (x) det^k (x) K, det^m) = invariants of conj(E det^k K) (x) det^m ... computed
# as the multiplicity of det^m in E (x) det^k (x) K via the inner product with conj(det^m)
ok5 = True
for k, m, aa in product(rng, repeat=3):
    d = integrate(np.conj(DET**m) * (CHI_E * DET**k) * (CHI_E * DET**aa))
    ok5 = ok5 and abs(d - (1 if aa == m - k - 1 else 0)) < 1e-9
checks["5_linear_convention_twist_is_m_minus_k_minus_1"] = ok5

# 6 independence of F1 and F2
k0, m0 = 1, -1
chi_K2 = 2 * CHI_E * DET ** (k0 - m0)
d2 = integrate(np.conj(CHI_E * DET**k0) * DET**m0 * chi_K2)
d0 = integrate(np.conj(CHI_E * DET**k0) * DET**m0 * 0 * DET)
checks["6_two_copies_dimension_two"] = abs(d2 - 2) < 1e-9
checks["6_zero_carrier_dimension_zero"] = abs(d0) < 1e-12

print("[H-Fac]: representation theory of the weak linking carrier (Weyl integration on U(2), SU(2) characters)")
print("=" * 96)
allok = True
for k, v in checks.items():
    ok = bool(v)
    allok = allok and ok
    print(f"  [{'PASS' if ok else 'FAIL'}]  {k}")
print(f"  {sum(bool(v) for v in checks.values())}/{len(checks)} checks pass")
print("=" * 96)
print("RESULT: the Lorentz factor of the sesquilinear pairing carries one invariant; the weak factor requires a doublet")
print("        component E (x) det^(k-m) in the linking carrier K (a character gives zero); (F1) existence and (F2)")
print("        multiplicity at most one are independent premises. No source supplies K.")
print("ALL CHECKS PASS" if allok else "SOME CHECKS FAILED")
raise SystemExit(0 if allok else 1)
