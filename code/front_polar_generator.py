"""Polar generator of a morphism family versus the vanishing projection A_Pi^odd (exact sympy, no numerics).

Objects. L(M) is the lift of M in sl_2(C) to Sym^2(V_gen) in the orthonormal weight basis (e_0, e_+, e_-);
S = Sym^2(eps) and J_Pi z = S conj(z) is the internal antilinear parity of the generation copy (Q14 Section 6).
A_Pi^odd := antiherm(Pi_odd(L(M))) is the projection of the paper (it vanishes identically).
For a smooth real-parameter family Y(gamma) = U(gamma) P(gamma) of invertible morphisms between FIXED orthonormal
frames, P > 0, the polar generator is Omega = U^dagger dU/dgamma (anti-hermitian).
The checks apply only to a real parameter, fixed orthonormal frames and P > 0: singular supports and moving frames
(frame transport, covariant derivative) are not covered.

Checks.
  M1  the Sylvester equation Omega P + P Omega = U^dagger dY - dY^dagger U holds for the family Y = Sym^2 g(gamma),
      g = exp(2 gamma E) exp(gamma F), at gamma = 1 (exact); its linear operator is invertible (P > 0, leading
      minors), and the solution equals the eigenbasis formula Omega_ij = (C - C^dagger)_ij / (p_i + p_j).
  M2  at the identity Omega(0) = antiherm(L(M)) for dY = L(M) with symbolic complex p, q, r; it is J_Pi-even and
      non-zero, whereas A_Pi^odd = antiherm(Pi_odd(L(M))) = 0: the two objects differ.
  M3  exact real counterexample g = exp(tE) exp(sF): the unitary polar factor is
      [[2 + ts, t - s], [s - t, 2 + ts]] / sqrt((2 + ts)^2 + (t - s)^2); Sym^2 of it is the polar unitary factor of
      Sym^2 g; its (e_+, e_-) entry is (t - s)^2 / ((2 + ts)^2 + (t - s)^2): 1/5 at (1, 0), 1/17 at (2, 1), zero at
      t = s; the tangent at 0 is L((E - F)/2) (times t0 - s0), non-zero and J_Pi-even, while A_Pi^odd(L(E)) = 0;
      the Jarlskog invariant of the real factor is zero.
  M4  non-identifiability: Y0 = D and Y1 = V D with V = Sym^2 u have the same square D^2, and different polar factors.
  M6  first order versus finite: the (e_+, e_-) tangent component of Omega(0) is zero, the finite entry is
      O(gamma^2); a family with Omega(0) = 0 (U = exp(gamma^2 K)) has an off-diagonal entry of order gamma^2.
Exit status 0 iff every check passes.
"""

import sympy as sp

checks = {}


def zero(M):
    return all(sp.simplify(sp.expand(e)) == 0 for e in M)


B = sp.Matrix.hstack(sp.Matrix([0, 1, 1, 0]) / sp.sqrt(2), sp.Matrix([1, 0, 0, 0]), sp.Matrix([0, 0, 0, 1]))


def rho(g):
    return sp.simplify(B.H * sp.kronecker_product(g, g) * B)


def lift(M):
    I2 = sp.eye(2)
    return sp.simplify(B.H * (sp.kronecker_product(M, I2) + sp.kronecker_product(I2, M)) * B)


def antiherm(X):
    return (X - X.H) / 2


E = sp.Matrix([[0, 1], [0, 0]])
F = sp.Matrix([[0, 0], [1, 0]])
S = rho(sp.Matrix([[0, 1], [-1, 0]]))


def T(X):
    return S * X.conjugate() * S.inv()


def polar_u(t, s):
    n = sp.sqrt((2 + t * s) ** 2 + (t - s) ** 2)
    return sp.Matrix([[2 + t * s, t - s], [s - t, 2 + t * s]]) / n


p, q, r = sp.symbols("p q r", complex=True)
Mgen = sp.Matrix([[r, p], [q, -r]])
L = lift(Mgen)
checks["0_S_is_the_paper_matrix"] = S == sp.Matrix([[-1, 0, 0], [0, 0, 1], [0, 1, 0]])
checks["0_lift_has_paper_entries"] = zero(L - sp.Matrix([[0, sp.sqrt(2) * q, sp.sqrt(2) * p],
                                                          [sp.sqrt(2) * p, 2 * r, 0],
                                                          [sp.sqrt(2) * q, 0, -2 * r]]))

# M3: real polar factor, general (t, s)
t, s = sp.symbols("t s", real=True)
g = sp.Matrix([[1 + t * s, t], [s, 1]])
checks["3_g_equals_exp_tE_exp_sF"] = g == (sp.eye(2) + t * E) * (sp.eye(2) + s * F)
u = polar_u(t, s)
Pm = sp.simplify(u.T * g)
checks["3_u_orthogonal"] = zero(u.T * u - sp.eye(2))
checks["3_P_symmetric_det1"] = zero(Pm - Pm.T) and sp.simplify(Pm.det() - 1) == 0
for (t0, s0) in ((1, 0), (2, 1)):
    P0 = Pm.subs({t: t0, s: s0})
    checks[f"3_P_positive_definite_at_t{t0}_s{s0}"] = sp.simplify(P0[0, 0]) > 0 and sp.simplify(P0.det() - 1) == 0
Gs, Us, Ps = rho(g), rho(u), rho(Pm)
checks["3_sym2_polar_product"] = zero(Us * Ps - Gs)
checks["3_sym2_u_unitary"] = zero(Us.H * Us - sp.eye(3))
checks["3_sym2_P_hermitian"] = zero(Ps - Ps.H)
checks["3_external_entry_formula"] = sp.simplify(
    Us[1, 2] - (t - s) ** 2 / ((2 + t * s) ** 2 + (t - s) ** 2)) == 0
checks["3_external_entry_1_5_at_t1_s0"] = sp.simplify(Us[1, 2].subs({t: 1, s: 0})) == sp.Rational(1, 5)
checks["3_external_entry_1_17_at_t2_s1"] = sp.simplify(Us[1, 2].subs({t: 2, s: 1})) == sp.Rational(1, 17)
checks["3_zero_rotation_at_t_equal_s"] = zero((Us - sp.eye(3)).subs(s, t))
g21 = g.subs({t: 2, s: 1})
checks["3_t2_s1_cascade_matrix"] = g21 == sp.Matrix([[3, 2], [1, 1]])
u21 = polar_u(2, 1)
checks["3_t2_s1_polar_factor"] = zero(u21 - sp.Matrix([[4, 1], [-1, 4]]) / sp.sqrt(17))
P21 = sp.simplify(Ps.subs({t: 2, s: 1}))
checks["3_t2_s1_sym2_P_positive_definite"] = (P21[0, 0] > 0 and P21[:2, :2].det() > 0 and P21.det() > 0
                                              and zero(P21 - P21.H))
U21 = sp.simplify(Us.subs({t: 2, s: 1}))
checks["3_t2_s1_internal_entries_pm_4sqrt2_over_17"] = (
    sp.simplify(U21[0, 1] + 4 * sp.sqrt(2) / 17) == 0 and sp.simplify(U21[0, 2] - 4 * sp.sqrt(2) / 17) == 0
    and sp.simplify(U21[1, 0] - 4 * sp.sqrt(2) / 17) == 0 and sp.simplify(U21[2, 0] + 4 * sp.sqrt(2) / 17) == 0)
checks["3_jarlskog_zero_for_real_factor"] = sp.simplify(
    (U21[0, 0] * U21[1, 1] * sp.conjugate(U21[0, 1]) * sp.conjugate(U21[1, 0])).as_real_imag()[1]) == 0
checks["3_polar_factor_J_even"] = zero(T(Us) - Us)

# M1: Sylvester equation for the family g(gamma), (t, s) = (2 gamma, gamma)
gam = sp.symbols("gamma", real=True)
Yfam = rho(g.subs({t: 2 * gam, s: gam}))
Ufam = rho(u.subs({t: 2 * gam, s: gam}))
Pfam = sp.simplify(Ufam.H * Yfam)
Omega = sp.simplify(Ufam.H * Ufam.diff(gam))
dY = Yfam.diff(gam)
checks["1_Y_equals_U_P_for_all_gamma"] = zero(Ufam * Pfam - Yfam)
checks["1_Omega_antihermitian"] = zero(Omega + Omega.H)
C = Ufam.H * dY
checks["1_Sylvester_identity_all_gamma"] = zero(Omega * Pfam + Pfam * Omega - (C - C.H))
Y1 = sp.simplify(Yfam.subs(gam, 1))
P1 = sp.simplify(Pfam.subs(gam, 1))
O1 = sp.simplify(Omega.subs(gam, 1))
C1 = sp.simplify(C.subs(gam, 1))
checks["1_P_at_gamma1_positive_definite"] = (P1[0, 0] > 0 and P1[:2, :2].det() > 0 and P1.det() > 0
                                             and zero(P1 - P1.H))
# Sylvester linear operator X -> X P + P X as a 9x9 matrix: invertible, hence the solution is unique.
sylv = sp.kronecker_product(sp.eye(3), P1.T) + sp.kronecker_product(P1, sp.eye(3))
checks["1_Sylvester_operator_invertible"] = sp.simplify(sylv.det()) != 0
# eigenbasis formula, exactly: P1 = W diag(w) W^T with W orthogonal (rational P1)
ev = P1.eigenvects()
vals, vecs = [], []
for val, mult, vs in ev:
    for v in vs:
        vals.append(val)
        vecs.append(v / sp.sqrt((v.H * v)[0]))
W = sp.Matrix.hstack(*vecs)
checks["1_P1_has_three_positive_eigenvalues"] = len(vals) == 3 and all(sp.simplify(v) > 0 for v in vals)
Cw = sp.simplify(W.H * (C1 - C1.H) * W)
Ow = sp.Matrix(3, 3, lambda i, j: sp.simplify(Cw[i, j] / (vals[i] + vals[j])))
checks["1_eigenbasis_formula_matches_definition"] = zero(W * Ow * W.H - O1)
# Without Sylvester: one can solve X P + P X = C - C^dagger and recover Omega
X = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"x{i}{j}"))
sol = sp.solve(list(X * P1 + P1 * X - (C1 - C1.H)), list(X), dict=True)
checks["1_Sylvester_solution_unique_and_equal_Omega"] = len(sol) == 1 and zero(X.subs(sol[0]) - O1)

# M2: at the identity
gam_s = sp.symbols("gamma_s", real=True)
Omega0 = (L - L.H) / 2
checks["2_Omega0_equals_antiherm_L"] = zero(Omega0 - antiherm(L))
checks["2_Omega0_is_J_even"] = zero(T(Omega0) - Omega0)
odd = (L - T(L)) / 2
checks["2_odd_projection_is_hermitian_part"] = zero(odd - (L + L.H) / 2)
checks["2_Aodd_vanishes_for_all_complex_pqr"] = zero(antiherm(odd))
checks["2_Omega0_nonzero_example"] = (Omega0.subs({p: 1, q: 0, r: 0}) != sp.zeros(3, 3))
# Omega(0) obtained from the Sylvester formula applied to Y = I + gamma L(M) (P = U = I at gamma = 0)
Yl = sp.eye(3) + gam_s * L
Ul0, Pl0 = sp.eye(3), sp.eye(3)
Cl = Ul0.H * Yl.diff(gam_s)
Om_syl = sp.Matrix(3, 3, lambda i, j: (Cl - Cl.H)[i, j] / 2)
checks["2_Sylvester_at_identity_gives_antiherm"] = zero(Om_syl - antiherm(L))
# polar tangent of the cascade at 0: dU/dgamma|0 = lift((E - F)/2) for (t, s) = (2 gamma, gamma)
tang = sp.simplify(Ufam.diff(gam).subs(gam, 0))
checks["3_tangent_is_lift_of_E_minus_F_half"] = zero(tang - lift((E - F) / 2))
checks["3_tangent_nonzero_J_even"] = tang != sp.zeros(3, 3) and zero(T(tang) - tang)
checks["3_Aodd_of_L_E_is_zero"] = zero(antiherm((lift(E) - T(lift(E))) / 2))
checks["3_external_tangent_component_zero"] = sp.simplify(tang[1, 2]) == 0 and sp.simplify(tang[2, 1]) == 0
checks["3_internal_tangent_component_nonzero"] = sp.simplify(tang[0, 1]) != 0

# M4: same square, different polar factors
D = sp.diag(1, sp.sqrt(sp.Rational(3, 5)), sp.sqrt(sp.Rational(2, 5)))
Y0 = D
Y1m = U21 * D
checks["4_same_square"] = zero(Y1m.H * Y1m - D * D) and zero(Y0.H * Y0 - D * D)
checks["4_different_morphisms"] = not zero(Y1m - Y0)
checks["4_polar_factors_differ"] = zero(Y0 * D.inv() - sp.eye(3)) and not zero(Y1m * D.inv() - sp.eye(3))
checks["4_V_unitary"] = zero(U21.H * U21 - sp.eye(3))
V2 = rho(sp.Matrix([[0, 1], [-1, 0]]))
checks["4_second_countermodel_any_unitary"] = zero((V2 * D).H * (V2 * D) - D * D)

# M6: first order versus finite
checks["6_external_entry_O_gamma2"] = (sp.simplify(sp.series(Ufam[1, 2], gam, 0, 2).removeO()) == 0
                                       and sp.simplify(sp.series(Ufam[1, 2], gam, 0, 3).removeO()
                                                       - gam ** 2 / 4) == 0)
Uk = sp.Matrix([[sp.cos(gam ** 2), sp.sin(gam ** 2)],
                                                                  [-sp.sin(gam ** 2), sp.cos(gam ** 2)]])
Omk = sp.simplify(Uk.H * Uk.diff(gam))
checks["6_vanishing_tangent_at_0"] = zero(Omk.subs(gam, 0))
checks["6_finite_exit_at_second_order"] = sp.simplify(sp.series(Uk[0, 1], gam, 0, 3).removeO() - gam ** 2) == 0

print("Polar generator of a morphism family versus the vanishing projection A_Pi^odd (exact symbolic)")
print("=" * 92)
allok = True
for k, val in checks.items():
    ok = bool(val)
    allok = allok and ok
    print(f"  [{'PASS' if ok else 'FAIL'}]  {k}")
print(f"  {sum(bool(v) for v in checks.values())}/{len(checks)} checks pass")
print("=" * 92)
print("RESULT: A_Pi^odd, the projection of the sl_2 lift onto the J_Pi-odd anti-hermitian part, vanishes and does not")
print("        determine the polar factor. The polar generator Omega of a morphism family solves the Sylvester")
print("        equation Omega P + P Omega = U^dagger dY - dY^dagger U (real parameter, fixed orthonormal frames,")
print("        P > 0); at the identity it is antiherm(L(M)), J_Pi-even and non-zero. For the Sym^2 lift of the real")
print("        cascade the finite polar factor has the non-zero external entry (t - s)^2 / ((2 + ts)^2 + (t - s)^2)")
print("        (1/5 at (1, 0), 1/17 at (2, 1)); the Jarlskog invariant of the real factor is zero. The square")
print("        D^2 does not determine the polar factor. Whether the physical Yukawa morphism family exists and which")
print("        polar class it has is open: the corpus supplies no such family.")
print("ALL CHECKS PASS" if allok else "SOME CHECKS FAILED")
raise SystemExit(0 if allok else 1)
