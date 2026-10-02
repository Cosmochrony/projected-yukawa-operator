"""Reality structure of the sl_2 lift on C^3_gen = Sym^2(V_gen) under the ANTILINEAR internal parity J_Pi.

Objects. M = pE + qF + rH in sl_2(C), p, q, r complex symbols. L(M) is its lift to Sym^2(V_gen) in the orthonormal
weight basis e_0 = (v+v- + v-v+)/sqrt2, e_+ = v+v+, e_- = v-v-. J_Pi z = S conj(z), S = Sym^2(eps), eps = [[0,1],[-1,0]]
(the internal antilinear parity of Q14 Section 6, J_Pi^2 = -1 on V_gen, so S real, S^2 = +1 on Sym^2).
J_Pi-odd projection: Pi_odd(X) = (X - S conj(X) S^-1)/2. The chiral generator of the paper is
A_Pi := antiherm(Pi_odd(L(M))), a DEFINITION of the paper (no source selects the J_Pi-odd part for the generator of U_Pi).

Checks (exact, sympy).
  1  derived lift: L(H) = 2 J_3 and L(M) has the printed entries; S = Sym^2(eps) is real with S^2 = 1.
  2  Proposition (reality structure) (i)  : S conj(L(M)) S^-1 = -L(M)^dagger for symbolic complex p, q, r.
  3  (ii) : Pi_odd(L) = (L + L^dagger)/2 and Pi_even(L) = (L - L^dagger)/2.
  4  (iii): A_Pi = antiherm(Pi_odd(L)) = 0 identically, for all complex p, q, r.
  5  the J_Pi-even anti-hermitian part antiherm(L) has the internal entry sqrt2 (q - conj p)/2, non-zero for real p != q:
     it is NOT part of A_Pi (the question of the paper's definition, flagged in the paper).
  6  (iv): the J_Pi-odd anti-hermitian subspace of u(3) has real dimension 6 and is Hilbert-Schmidt orthogonal to
     L(E), L(F), L(H), hence to the complex span L(sl_2(C)) (the spin-one component of End(Sym^2 V) = 1 + 3 + 5).
  7  NEGATIVE CONTROL: the LINEAR involution X -> S X S^-1 (not the antilinear parity of the paper) gives
     antiherm(Pi_lin_odd(L)) with internal entries (sqrt2/4)(p + q - conj p - conj q) = i (sqrt2/2) Im(p + q): non-zero
     iff Im(p + q) != 0 (so not "Im p or Im q": Im p = -Im q gives zero); the two parities agree on real L.
  8  detectors of the polar non-triviality proposition: for U(g) = exp(g A), A anti-hermitian with generic entries,
     |U_12|^2 = g^2 |A_12|^2 + O(g^3) and the quartet phase is -g^3 Im(A_12 A_23 A_31) + O(g^4) (exact series, to order 3).
Exit status 0 iff every check passes.
"""

import sympy as sp

p, q, r = sp.symbols("p q r", complex=True)
checks = {}


def zero(M):
    return all(sp.simplify(sp.expand(e)) == 0 for e in M)


def sym2_group(g):
    B = sp.Matrix.hstack(sp.Matrix([0, 1, 1, 0]) / sp.sqrt(2), sp.Matrix([1, 0, 0, 0]), sp.Matrix([0, 0, 0, 1]))
    return sp.simplify(B.H * sp.kronecker_product(g, g) * B)


def lift(M):
    I2 = sp.eye(2)
    D = sp.kronecker_product(M, I2) + sp.kronecker_product(I2, M)
    B = sp.Matrix.hstack(sp.Matrix([0, 1, 1, 0]) / sp.sqrt(2), sp.Matrix([1, 0, 0, 0]), sp.Matrix([0, 0, 0, 1]))
    return sp.simplify(B.H * D * B)


def antiherm(X):
    return (X - X.H) / 2


E = sp.Matrix([[0, 1], [0, 0]])
F = sp.Matrix([[0, 0], [1, 0]])
H = sp.Matrix([[1, 0], [0, -1]])
eps = sp.Matrix([[0, 1], [-1, 0]])
S = sym2_group(eps)
L = lift(p * E + q * F + r * H)
Sc = S.inv()

# 1
checks["1_L_H_equals_2J3"] = lift(H) == sp.diag(0, 2, -2)
checks["1_L_entries"] = L == sp.Matrix([[0, sp.sqrt(2) * q, sp.sqrt(2) * p],
                                        [sp.sqrt(2) * p, 2 * r, 0],
                                        [sp.sqrt(2) * q, 0, -2 * r]])
checks["1_S_real_involution"] = S == sp.Matrix([[-1, 0, 0], [0, 0, 1], [0, 1, 0]]) and S * S == sp.eye(3)

# 2-4
conjugated = sp.simplify(S * L.conjugate() * Sc)
checks["2_i_S_conjL_Sinv_equals_minus_Ldagger"] = zero(conjugated + L.H)
odd = (L - conjugated) / 2
even = (L + conjugated) / 2
checks["3_ii_odd_is_hermitian_part"] = zero(odd - (L + L.H) / 2)
checks["3_ii_even_is_antihermitian_part"] = zero(even - (L - L.H) / 2)
A = sp.simplify(antiherm(odd))
checks["4_iii_A_Pi_vanishes_for_all_complex_pqr"] = (A == sp.zeros(3, 3))

# 5
internal_even = sp.simplify(antiherm(even)[0, 1])
checks["5_even_internal_entry_formula"] = sp.simplify(internal_even - sp.sqrt(2) * (q - sp.conjugate(p)) / 2) == 0
real_pq = {p: sp.Rational(1), q: sp.Rational(3)}
checks["5_even_internal_entry_nonzero_for_real_p_neq_q"] = sp.simplify(internal_even.subs(real_pq)) != 0

# 6
xs = sp.symbols("x0:9", real=True)
d, a, b = xs[0:3], xs[3:6], xs[6:9]
X = sp.zeros(3, 3)
for k, (i, j) in enumerate([(0, 1), (0, 2), (1, 2)]):
    X[i, j] = a[k] + sp.I * b[k]
    X[j, i] = -(a[k] - sp.I * b[k])
for k in range(3):
    X[k, k] = sp.I * d[k]
cond = X + S * X.conjugate() * Sc
eqs = []
for e in cond:
    eqs += [sp.re(sp.expand(e)), sp.im(sp.expand(e))]
eqs = [sp.simplify(e) for e in eqs if sp.simplify(e) != 0]
sol = sp.solve(eqs, xs, dict=True)[0]
free = [x for x in xs if x not in sol]
checks["6_odd_antiherm_subspace_dim_6"] = len(free) == 6
Xs = X.subs(sol)
checks["6_odd_antiherm_family_matches_paper"] = (
    Xs[1, 1] == Xs[2, 2] and sp.simplify(Xs[0, 2] - sp.conjugate(Xs[0, 1])) == 0
    and sp.simplify(Xs[2, 0] + Xs[0, 1]) == 0 and sp.simplify(Xs[1, 0] + sp.conjugate(Xs[0, 1])) == 0
    and sp.simplify(Xs[2, 1] + sp.conjugate(Xs[1, 2])) == 0)
checks["6_orthogonal_to_L_sl2"] = all(sp.simplify(sp.expand((lift(G).H * Xs).trace())) == 0 for G in (E, F, H))

# 7 negative control: linear involution
lin = sp.simplify((L - S * L * Sc) / 2)
A_lin = sp.simplify(antiherm(lin))
expected_internal = sp.sqrt(2) * (p + q - sp.conjugate(p) - sp.conjugate(q)) / 4
checks["7_linear_control_internal_is_Im_p_plus_q"] = (sp.simplify(A_lin[0, 1] - expected_internal) == 0)
pq_anti = {p: sp.I, q: -sp.I}
checks["7_linear_control_vanishes_for_Im_p_eq_minus_Im_q"] = sp.simplify(A_lin[0, 1].subs(pq_anti)) == 0
checks["7_linear_control_nonzero_for_Im_p_only"] = sp.simplify(A_lin[0, 1].subs({p: sp.I, q: 0})) != 0
checks["7_linear_antilinear_agree_on_real_L"] = zero((lin - odd).subs({p: sp.Rational(2), q: sp.Rational(5), r: sp.Rational(7)}))

# 8 detectors of the polar non-triviality proposition (exact, to order 3 in g)
g = sp.symbols("g", real=True)
a12, a13, a23 = sp.symbols("a12 a13 a23", complex=True)
d1, d2, d3 = sp.symbols("d1 d2 d3", real=True)
Agen = sp.Matrix([[sp.I * d1, a12, a13], [-sp.conjugate(a12), sp.I * d2, a23],
                  [-sp.conjugate(a13), -sp.conjugate(a23), sp.I * d3]])
U = sp.eye(3) + g * Agen + g**2 * Agen**2 / 2 + g**3 * Agen**3 / 6
U12, U21 = U[0, 1], U[1, 0]
mod2 = sp.expand(U12 * sp.conjugate(U12))
c2 = sp.simplify(sp.expand(mod2).coeff(g, 2))
checks["8_modulus_leading_order"] = sp.simplify(c2 - a12 * sp.conjugate(a12)) == 0
quartet = sp.expand(U[0, 0] * U[1, 1] * sp.conjugate(U12) * sp.conjugate(U21))
Jq = sp.simplify((quartet - sp.conjugate(quartet)) / (2 * sp.I))
Jg3 = sp.simplify(sp.expand(Jq).coeff(g, 3))
a31 = -sp.conjugate(a13)
target = -sp.im(a12 * a23 * a31)
checks["8_jarlskog_orders_0_1_2_vanish"] = all(sp.simplify(sp.expand(Jq).coeff(g, k)) == 0 for k in (0, 1, 2))
checks["8_jarlskog_cubic_coefficient"] = sp.simplify(sp.expand_complex(Jg3 - target)) == 0

print("Reality structure of the sl_2 lift under the antilinear parity J_Pi (exact symbolic)")
print("=" * 92)
print("  L(M) = [[0, r2 q, r2 p], [r2 p, 2r, 0], [r2 q, 0, -2r]],  S = [[-1,0,0],[0,0,1],[0,1,0]],  J_Pi z = S conj(z)")
print("  A_Pi := antiherm(Pi_odd(L(M))), Pi_odd(X) = (X - S conj(X) S^-1)/2  (a definition of the paper)")
print("-" * 92)
allok = True
for k, val in checks.items():
    ok = bool(val)
    allok = allok and ok
    print(f"  [{'PASS' if ok else 'FAIL'}]  {k}")
print(f"  {sum(bool(v) for v in checks.values())}/{len(checks)} checks pass")
print("=" * 92)
print("RESULT: A_Pi, the anti-hermitian part of the J_Pi-odd part of the sl_2 lift of the step generator, vanishes")
print("        identically for all complex (p, q, r), because with the antilinear J_Pi the J_Pi-odd part of the lift")
print("        is its hermitian part. This is a statement about A_Pi as defined. The J_Pi-even anti-hermitian part of the")
print("        lift has the non-zero internal entry (sqrt2/2)(q - conj p), non-zero for real p != q, so the image of")
print("        sl_2 does reach the internal block e_0 <-> e_pm through that part; the external block R_mix is not")
print("        reached by the image of sl_2 (Q14 Rem. 6.4, Prop. 6.3). No physical exclusion of the internal block is")
print("        claimed: whether A_Pi, rather than the J_Pi-even part, is the right object is a modelling choice of the")
print("        paper that no source justifies. The linear involution is a different operation (negative control).")
print("ALL CHECKS PASS" if allok else "SOME CHECKS FAILED")
raise SystemExit(0 if allok else 1)
