"""Polar class of the generation block: rephasing groups defined by the J3 grading, not by the eigenbasis of Y Y^dag.

Purpose. Proposition (Polar class) of the paper quotients the polar factor U of the generation block
Y = U H^{1/2} by the rephasing groups G_L, G_R, each the group of unitaries commuting with the Cartan generator
J3 = diag(0, 1, -1) of its own carrier (a copy of C^3_gen = Sym^2(V_gen); the grading of each carrier is the premise
[H-Grad] of the paper). The groups are defined by the J3 grading of the carriers and not by the eigenbasis of Y Y^dag.
Here Y is the generation block Ygen = Y_Pi^gen of the full morphism Y_Pi. This script checks, with exact arithmetic over Q(i) (sympy) and one
high-precision numerical rank computation (mpmath, 60 digits), every step of the argument.

Setup. H = lambda^2 diag(1, 1/2 + u, 1/2 - u) is a positive operator commuting with J3 (the model operator on C^3_gen).
No identification with E_Pi^2 and no Yukawa dictionary is used here: the checks concern the algebra of Y = U H^{1/2}.

Checks.
  1  G_L = G_R = diagonal unitaries: the commutant of J3 (distinct eigenvalues 0, 1, -1) is the diagonal matrices.
  2  Generic distinct H, exact example (u = 7/50, so that H^{1/2} = diag(1, 4/5, 3/5) is rational):
     a  U is a non-monomial exact unitary and Y = U H^{1/2} has Y^dag Y = H, Y H^{-1/2} = U;
     b  for diagonal unitary V_L, V_R: Y' = V_R Y V_L^{-1} has Y'^dag Y' = H and U' = V_R U V_L^{-1};
     c  the moduli |U_ij|^2 and the quartet J are invariant, exactly;
     d  Y Y^dag = U H U^dag is NOT diagonal in the right J3 basis, and (P_j)_ii = |U_ij|^2 for the spectral
        projectors P_j of Y Y^dag: the moduli measure the misalignment of the right level eigenbasis with the right
        J3 grading, with the grading fixed independently of Y Y^dag;
     e  the circular choice V_R := U (right basis pinned by the eigenbasis of Y Y^dag) removes U for every Y, which is
        why the right group is not defined that way.
  3  Degenerate H, u = 0: H = diag(1, 1/2, 1/2). The J3 groups still act and H and the moduli are invariant, but the
     commutant of H is strictly larger than the diagonal group, and a unitary of that commutant changes the moduli of U
     while leaving H invariant: the J3 class is then finer than the class under Stab(H).
  4  Singular H, u = 1/2: H = diag(1, 1, 0). The unitary polar factor is unique only up to the phase of the column
     paired with the zero singular value; that ambiguity is a J3 rephasing of the left carrier, so the class is well
     defined (the spectrum is however degenerate there, as at u = 0).
  5  Local rank: at random points of U(3) the differential of (|U_ij|^2, J) has rank 4 = 9 - 5 (60-digit arithmetic),
     so the class space has dimension four and (moduli, J) are local coordinates on it at generic points; the
     rephasing orbit has dimension 5 (stabiliser the common phase), computed from the differential of the action,
     and orbit plus class dimensions add to nine.
Deterministic: fixed seed for the random points of check 5. Exit status 0 iff every check passes.
"""

import random

import mpmath as mp
import sympy as sp

I = sp.I


def dag(M):
    return M.conjugate().T


def zero(M):
    return all(sp.simplify(sp.expand(e)) == 0 for e in M)


def phase_matrix(a, b, c):
    return sp.diag(a, b, c)


def quartet(U):
    return sp.im(sp.expand(U[0, 0] * U[1, 1] * sp.conjugate(U[0, 1]) * sp.conjugate(U[1, 0])))


def moduli2(U):
    return sp.Matrix(3, 3, lambda i, j: sp.simplify(sp.expand(U[i, j] * sp.conjugate(U[i, j]))))


def check1():
    W = sp.Matrix(3, 3, lambda i, j: sp.Symbol(f"w{i}{j}"))
    J3 = sp.diag(0, 1, -1)
    eqs = list(W * J3 - J3 * W)
    sol = sp.solve(eqs, list(W), dict=True)[0]
    Wc = W.subs(sol)
    off = [Wc[i, j] for i in range(3) for j in range(3) if i != j]
    return all(sp.simplify(e) == 0 for e in off)


def exact_unitary():
    K = sp.Matrix([[0, 1 + I, 2], [1 - I, 1, -I], [2, I, -1]])
    assert zero(K - dag(K))
    U = (sp.eye(3) - I * K) * (sp.eye(3) + I * K).inv()
    U = U.applyfunc(lambda e: sp.simplify(sp.expand_complex(e)))
    return U


def check2():
    out = {}
    U = exact_unitary()
    out["2a_unitary"] = zero(dag(U) * U - sp.eye(3))
    u = sp.Rational(7, 50)
    Hh = sp.diag(1, sp.Rational(4, 5), sp.Rational(3, 5))
    H = Hh * Hh
    out["2a_levels"] = H == sp.diag(1, sp.Rational(1, 2) + u, sp.Rational(1, 2) - u)
    Y = U * Hh
    out["2a_Y_dag_Y_is_H"] = zero(dag(Y) * Y - H)
    out["2a_U_recovered"] = zero(Y * Hh.inv() - U)
    off = [U[i, j] for i in range(3) for j in range(3) if i != j]
    out["2a_U_not_diagonal"] = any(sp.simplify(e) != 0 for e in off)
    zeros_per_row = [sum(1 for j in range(3) if sp.simplify(U[i, j]) == 0) for i in range(3)]
    out["2a_U_not_monomial"] = max(zeros_per_row) < 2

    VL = phase_matrix((3 + 4 * I) / 5, (5 + 12 * I) / 13, (8 - 15 * I) / 17)
    VR = phase_matrix((20 + 21 * I) / 29, (7 - 24 * I) / 25, (-12 + 5 * I) / 13)
    out["2b_phases_unit"] = zero(dag(VL) * VL - sp.eye(3)) and zero(dag(VR) * VR - sp.eye(3))
    Yp = VR * Y * VL.inv()
    out["2b_Hprime_equals_H"] = zero(dag(Yp) * Yp - H)
    out["2b_H_commutes_with_VL"] = zero(VL * H - H * VL)
    Up = Yp * Hh.inv()
    out["2b_U_transforms"] = zero(Up - VR * U * VL.inv())
    out["2b_right_square_transforms"] = zero(Yp * dag(Yp) - VR * (Y * dag(Y)) * dag(VR))
    out["2c_moduli_invariant"] = zero(moduli2(Up) - moduli2(U))
    out["2c_quartet_invariant"] = sp.simplify(quartet(Up) - quartet(U)) == 0
    out["2c_quartet_nonzero"] = sp.simplify(quartet(U)) != 0
    out["2c_phases_actually_change_U"] = not zero(Up - U)

    YYd = Y * dag(Y)
    out["2d_right_square_is_UHUdag"] = zero(YYd - U * H * dag(U))
    off = [YYd[i, j] for i in range(3) for j in range(3) if i != j]
    out["2d_right_square_not_J3_diagonal"] = any(sp.simplify(e) != 0 for e in off)
    h = [H[k, k] for k in range(3)]
    ok = True
    for j in range(3):
        P = sp.eye(3)
        for k in range(3):
            if k != j:
                P = P * (YYd - h[k] * sp.eye(3)) / (h[j] - h[k])
        P = P.applyfunc(lambda e: sp.simplify(sp.expand_complex(e)))
        for i in range(3):
            ok = ok and sp.simplify(P[i, i] - U[i, j] * sp.conjugate(U[i, j])) == 0
    out["2d_projector_diagonals_are_moduli"] = ok

    Ypin = dag(U) * Y
    out["2e_circular_pinning_trivialises_U"] = zero(Ypin - Hh)
    return out


def check3():
    out = {}
    U = exact_unitary()
    Hh = sp.diag(1, 1 / sp.sqrt(2), 1 / sp.sqrt(2))
    H = Hh * Hh
    Y = U * Hh
    out["3_H_is_u0_level_operator"] = zero(H - sp.diag(1, sp.Rational(1, 2), sp.Rational(1, 2)))
    VL = phase_matrix((3 + 4 * I) / 5, (5 + 12 * I) / 13, (8 - 15 * I) / 17)
    VR = phase_matrix((20 + 21 * I) / 29, (7 - 24 * I) / 25, (-12 + 5 * I) / 13)
    Yp = VR * Y * VL.inv()
    Up = Yp * Hh.inv()
    out["3_J3_groups_still_act"] = zero(dag(Yp) * Yp - H) and zero(moduli2(Up) - moduli2(U))
    c, s = sp.Rational(3, 5), sp.Rational(4, 5)
    V = sp.Matrix([[1, 0, 0], [0, c, -s], [0, s, c]])
    out["3_V_in_commutant_of_H"] = zero(V * H - H * V)
    out["3_V_not_diagonal"] = V[1, 2] != 0
    Y2 = Y * V.inv()
    out["3_H_invariant_under_V"] = zero(dag(Y2) * Y2 - H)
    U2 = Y2 * Hh.inv()
    out["3_moduli_not_invariant_under_V"] = not zero(moduli2(U2) - moduli2(U))
    return out


def check4():
    out = {}
    U = exact_unitary()
    Hh = sp.diag(1, 1, 0)
    H = Hh * Hh
    Y = U * Hh
    out["4_H_is_u_half_operator"] = zero(H - sp.diag(1, sp.Rational(1, 2) + sp.Rational(1, 2),
                                                     sp.Rational(1, 2) - sp.Rational(1, 2)))
    phi = (8 - 15 * I) / 17
    VL = phase_matrix(1, 1, phi)
    Uprime = U * VL
    out["4_same_Y_two_unitaries"] = zero(Uprime * Hh - Y) and zero(dag(Uprime) * Uprime - sp.eye(3))
    out["4_ambiguity_is_a_J3_rephasing"] = zero(Uprime - U * VL) and zero(VL * dag(VL) - sp.eye(3))
    out["4_class_unchanged"] = zero(moduli2(Uprime) - moduli2(U)) and sp.simplify(quartet(Uprime) - quartet(U)) == 0
    return out


def check5(npoints=6, seed=20261001):
    mp.mp.dps = 60
    rnd = random.Random(seed)
    sv_counts = []
    orbit_dims = []
    for _ in range(npoints):
        Z = mp.matrix(3, 3)
        for i in range(3):
            for j in range(3):
                Z[i, j] = mp.mpc(rnd.gauss(0, 1), rnd.gauss(0, 1))
        Q, _R = mp.qr(Z)
        U = Q

        def f(Um):
            vals = [abs(Um[i, j]) ** 2 for i in range(3) for j in range(3)]
            q = Um[0, 0] * Um[1, 1] * mp.conj(Um[0, 1]) * mp.conj(Um[1, 0])
            return vals + [mp.im(q)]

        basis = []
        for i in range(3):
            A = mp.matrix(3, 3)
            A[i, i] = 1j
            basis.append(A)
        for i in range(3):
            for j in range(i + 1, 3):
                A = mp.matrix(3, 3)
                A[i, j], A[j, i] = 1, -1
                basis.append(A)
                B = mp.matrix(3, 3)
                B[i, j], B[j, i] = 1j, 1j
                basis.append(B)
        h = mp.mpf(10) ** -25
        D = mp.matrix(10, 9)
        for k, A in enumerate(basis):
            Up = U * mp.expm(h * A)
            Um = U * mp.expm(-h * A)
            fp, fm = f(Up), f(Um)
            for r in range(10):
                D[r, k] = (fp[r] - fm[r]) / (2 * h)
        S = mp.svd_r(D, compute_uv=False)
        sv = [S[i] for i in range(len(S))]
        sv_counts.append(sum(1 for s in sv if s > mp.mpf(10) ** -12))
        # orbit dimension: rank of (a, b) -> d/dt [exp(i t a) U exp(-i t b)] over the 6 real phase parameters
        T = mp.matrix(18, 6)
        for k in range(6):
            a = [0, 0, 0]
            b = [0, 0, 0]
            if k < 3:
                a[k] = 1
            else:
                b[k - 3] = 1
            dU = mp.matrix(3, 3)
            for i in range(3):
                for j in range(3):
                    dU[i, j] = 1j * (a[i] - b[j]) * U[i, j]
            for i in range(3):
                for j in range(3):
                    T[2 * (3 * i + j), k] = mp.re(dU[i, j])
                    T[2 * (3 * i + j) + 1, k] = mp.im(dU[i, j])
        So = mp.svd_r(T, compute_uv=False)
        orbit_dims.append(sum(1 for t in range(len(So)) if So[t] > mp.mpf(10) ** -12))
    return {"5_rank_is_four_at_all_points": all(c == 4 for c in sv_counts),
            "5_orbit_dimension_is_five_at_all_points": all(c == 5 for c in orbit_dims),
            "5_orbit_plus_class_dimension_is_nine": all(c + d == 9 for c, d in zip(sv_counts, orbit_dims)),
            "_ranks": sv_counts, "_orbits": orbit_dims}


def main():
    results = {"1_commutant_of_J3_is_diagonal": check1()}
    results.update(check2())
    results.update(check3())
    results.update(check4())
    r5 = check5()
    ranks = r5.pop("_ranks")
    orbits = r5.pop("_orbits")
    results.update(r5)
    print("Polar class of the generation block: exact checks over Q(i) and a 60-digit rank computation")
    print("=" * 96)
    allok = True
    for k, v in results.items():
        ok = bool(v)
        allok = allok and ok
        print(f"  [{'PASS' if ok else 'FAIL'}]  {k}")
    print(f"  ranks of d(moduli, J) at the random points: {ranks}")
    print(f"  dimensions of the rephasing orbit at the random points: {orbits}")
    print("=" * 96)
    print("The rephasing groups are the J3-commutants of the two carriers, defined independently of Y Y^dag.")
    print("Restriction flagged: distinct spectrum of H (u not in {0, +-1/2}) for the J3 class to coincide with the")
    print("class under the commutant of H; the statement is made for 0 < u < 1/2.")
    print("ALL CHECKS PASS" if allok else "SOME CHECKS FAILED")
    return allok


if __name__ == "__main__":
    raise SystemExit(0 if main() else 1)
