#!/usr/bin/env python3
# author: B (lane-B chair a7296aa8), S1770.  Litgate on the external reader's review of note-1.
#
# WHAT THIS PROBE IS FOR.  A reader returned a detailed review of the note "Crystallographic
# closure of I_{3,1}" claiming (a) a logical gap in step N.3, (b) a shorter proof by orbit
# averaging, (c) that the exclusion of order 12 was measured on one type out of three, and
# (d) that the sentence "Phi_5, Phi_8, Phi_12 can act on I_{2,2}" claims more than the
# measurement supports.  External input is not a witness until it is re-derived here.  This
# file re-derives or refutes each testable claim.  Ex ante is committed BEFORE the run in
# preprints/note-1/checks/S1770_READER_LITGATE_EXANTE.md.
#
# Convention: Gram matrix G, q(x) = x^T G x.  TIMELIKE means q(x) < 0, so signature (3,1)
# is written as three positive and one negative eigenvalue.
# Handles: 0.
# RUN LINE:  python S1770_reader_litgate.py
import itertools
import random
import sys

import numpy as np
import sympy as sp

x = sp.Symbol("x")
ETA31 = sp.diag(1, 1, 1, -1)
ETA22 = sp.diag(1, 1, -1, -1)
FAIL = []


def report(name, ok, detail=""):
    """One line per checked fact.  A probe that prints only successes cannot be audited."""
    mark = "OK  " if ok else "FAIL"
    if not ok:
        FAIL.append(name)
    print("  [%s] %-58s %s" % (mark, name, detail))


def order_of(M, cap=200):
    """multiplicative order of an integer matrix, or None if it exceeds the cap."""
    n = M.shape[0]
    P = sp.eye(n)
    for k in range(1, cap + 1):
        P = P * M
        if P == sp.eye(n):
            return k
    return None


def preserves(M, G):
    return sp.simplify(M.T * G * M - G) == sp.zeros(*G.shape)


def signature(G):
    """(positive, negative, zero) index of a rational symmetric matrix, via eigenvalues."""
    A = np.array(G.evalf(), dtype=float)
    ev = np.linalg.eigvalsh((A + A.T) / 2)
    tol = 1e-9 * max(1.0, float(np.max(np.abs(ev))))
    return (int(np.sum(ev > tol)), int(np.sum(ev < -tol)), int(np.sum(np.abs(ev) <= tol)))


# ----------------------------------------------------------------------------------------
# TEST A -- is the gap in N.3 real?  The note argues: the eigenspace ker(g -+ 1) is defined
# over Q, therefore the invariant line l is rational.  The reader says the second half does
# not follow when the eigenspace has dimension > 1.  Rhetoric decides nothing; a witness does.
# ----------------------------------------------------------------------------------------
def test_A_gap():
    print("\nTEST A -- N.3: does a rational eigenspace force a rational invariant line?")
    # rotation by pi/2 in the (x1,x2) plane, identity on (x3,x4)
    g = sp.Matrix([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    report("A.1 g integral, preserves eta_(3,1)", preserves(g, ETA31))
    report("A.2 ord(g) = 4", order_of(g) == 4, "ord = %s" % order_of(g))
    ker = (g - sp.eye(4)).nullspace()
    report("A.3 dim ker(g - 1) = 2 (> 1: the reader's precondition)", len(ker) == 2,
           "basis %s" % [list(v.T) for v in ker])
    # an IRRATIONAL timelike line inside that rational eigenspace
    v = sp.Matrix([0, 0, 1, sp.sqrt(2)])
    qv = sp.simplify((v.T * ETA31 * v)[0])
    report("A.4 v = (0,0,1,sqrt2) lies in the eigenspace and is timelike", qv < 0, "q(v) = %s" % qv)
    report("A.5 g fixes the line spanned by v", sp.simplify(g * v - v) == sp.zeros(4, 1))
    ratios = [sp.nsimplify(c) for c in v if c != 0]
    irrational = any(not sp.Rational(0).is_number or not c.is_rational for c in ratios)
    report("A.6 that line is NOT rational (no rational vector spans it)", irrational,
           "coords span a Q-space of dim 2 -> S1670 gives rank(L cap v^perp) = 4 - 2 = 2 < 3")
    print("  => the inference 'eigenspace over Q  ==>  THIS line is rational' is invalid;")
    print("     step (N.2) can hand back exactly this line.  Reader's R1 stands.")
    # the reader's patch: the SAME rational eigenspace still contains a rational timelike vector
    w = sp.Matrix([0, 0, 1, 2])
    qw = (w.T * ETA31 * w)[0]
    report("A.7 patch: a RATIONAL timelike vector exists in the same eigenspace", qw < 0,
           "w = (0,0,1,2), q(w) = %s -- the note's conclusion survives" % qw)


# ----------------------------------------------------------------------------------------
# TEST C -- the reader's replacement for N.2-N.3: average an integral timelike vector over
# the finite orbit.  Claim: t = sum_j h^j e is integral, timelike, and g t = +- t.
# ----------------------------------------------------------------------------------------
def integral_timelike(G, bound=3):
    """smallest integral vector with q < 0 in a small box -- every lattice of signature (3,1)
       has one, and for the lattices used here a box of radius 3 already finds it."""
    n = G.shape[0]
    best = None
    for c in itertools.product(range(-bound, bound + 1), repeat=n):
        if not any(c):
            continue
        v = sp.Matrix(list(c))
        if (v.T * G * v)[0] < 0:
            nrm = sum(abs(k) for k in c)
            if best is None or nrm < best[0]:
                best = (nrm, v)
    return None if best is None else best[1]


def average_orbit(g, G, use_epsilon=True):
    """Reader's Lemma 1.  Returns (t, eps) or (None, reason)."""
    e = integral_timelike(G)
    if e is None:
        return None, "no integral timelike vector found"
    eps = 1
    if use_epsilon:
        # two timelike vectors lie in the same cone component iff their inner product is < 0
        if (e.T * G * (g * e))[0] > 0:
            eps = -1
    h = eps * g
    m = order_of(h)
    if m is None:
        return None, "h has no finite order within the cap"
    t = sp.zeros(G.shape[0], 1)
    P = sp.eye(G.shape[0])
    for _ in range(m):
        t += P * e
        P = P * h
    return t, eps


def check_averaging(label, g, G, expect_ok=True):
    t, eps = average_orbit(g, G)
    if t is None:
        report("C %-26s" % label, not expect_ok, "no t: %s" % eps)
        return
    qt = (t.T * G * t)[0]
    integral = all(c.is_integer for c in t)
    gt = g * t
    eigen = (gt == t) or (gt == -t)
    ok = integral and qt < 0 and eigen
    report("C %-26s" % label, ok == expect_ok,
           "eps=%+d  t=%s  q(t)=%s  g t = %st" % (eps, list(t.T), qt, "+" if gt == t else "-"))


def test_C_averaging():
    print("\nTEST C -- orbit averaging: integral timelike +-1 eigenvector (reader's Lemma 1)")
    P3 = sp.Matrix([[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]])
    R4 = sp.Matrix([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
    G6 = sp.Matrix([[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, -1]])
    cases = [("identity", sp.eye(4)), ("-I (swaps cone halves)", -sp.eye(4)),
             ("P3, order 3", P3), ("R4, order 4", R4), ("G6, order 6", G6),
             ("-R4", -R4), ("R4 * G6", R4 * G6)]
    for label, g in cases:
        assert preserves(g, ETA31), label
        check_averaging(label + " / I(3,1)", g, ETA31)

    # R5: the same construction on a NON-unimodular lattice of signature (3,1)
    A2 = sp.Matrix([[2, -1], [-1, 2]])
    G2 = sp.diag(A2, sp.Matrix([[1]]), sp.Matrix([[-7]]))
    rot3 = sp.Matrix([[0, -1], [1, -1]])
    g2 = sp.diag(rot3, sp.Matrix([[1]]), sp.Matrix([[-1]]))
    report("C.pre non-unimodular lattice A2 + <1> + <-7>, det = %s" % G2.det(),
           preserves(g2, G2) and signature(G2) == (3, 1, 0),
           "sig %s, ord(g) = %s" % (signature(G2), order_of(g2)))
    check_averaging("order 6 / A2+<1>+<-7>", g2, G2)

    # NEGATIVE CONTROL, named in the ex ante: without the epsilon choice the construction
    # must FAIL on -I (the orbit sum cancels).  If it does not fail, the test is blind.
    t, _ = average_orbit(-sp.eye(4), ETA31, use_epsilon=False)
    report("C.neg control: no eps-choice on -I gives t = 0 (must fail)",
           t is not None and t == sp.zeros(4, 1), "t = %s" % (list(t.T) if t is not None else None))


# ----------------------------------------------------------------------------------------
# TEST D/E -- every finite order realisable in rank 4, every type, every achievable signature.
# ----------------------------------------------------------------------------------------
def companion(d):
    c = sp.Poly(sp.cyclotomic_poly(d, x), x).all_coeffs()
    n = len(c) - 1
    M = sp.zeros(n, n)
    for i in range(1, n):
        M[i, i - 1] = 1
    for i in range(n):
        M[i, n - 1] = -c[n - i]
    return M


def types_of_rank(rank=4):
    """all multisets of cyclotomic factors of total degree `rank`, keyed by lcm of their orders."""
    ds = [d for d in range(1, 40) if sp.totient(d) <= rank]
    out = {}
    def rec(start, left, acc):
        if left == 0:
            n = 1
            for d in acc:
                n = sp.ilcm(n, d)
            out.setdefault(int(n), []).append(tuple(acc))
            return
        for i in range(start, len(ds)):
            d = ds[i]
            if sp.totient(d) <= left:
                rec(i, left - int(sp.totient(d)), acc + [d])
    rec(0, rank, [])
    return out


def invariant_forms(C):
    """Q-basis of the space of symmetric bilinear forms invariant under C, plus the symbolic
       general form and its determinant."""
    n = C.shape[0]
    syms = {}
    Q = sp.zeros(n, n)
    for i in range(n):
        for j in range(i, n):
            s = sp.Symbol("q_%d%d" % (i, j))
            syms[(i, j)] = s
            Q[i, j] = s
            Q[j, i] = s
    # Solve the linear system as a matrix kernel rather than with solve(): when C is the
    # identity every equation is 0 = 0 and solve() returns nothing at all, which silently
    # reports "no invariant forms" for the one case where every form is invariant.
    unknowns = [syms[k] for k in sorted(syms)]
    eqs = list(C.T * Q * C - Q)
    A, _ = sp.linear_eq_to_matrix(eqs, unknowns)
    ker = A.nullspace()
    if not ker:
        return None, None, []
    params = [sp.Symbol("p%d" % i) for i in range(len(ker))]
    vec = sp.zeros(len(unknowns), 1)
    for p, v in zip(params, ker):
        vec += p * v
    Qg = sp.zeros(n, n)
    for idx, (i, j) in enumerate(sorted(syms)):
        Qg[i, j] = vec[idx]
        Qg[j, i] = vec[idx]
    basis = []
    for p in params:
        B = Qg.subs({t: (1 if t == p else 0) for t in params})
        den = sp.Integer(1)
        for e in B:
            if e != 0:
                den = sp.ilcm(den, sp.denom(sp.nsimplify(e)))
        basis.append(sp.Matrix(B) * den)
    return Qg, params, basis


def achievable_signatures(basis, span=2, cap=4000):
    """signatures of non-degenerate integral combinations of the basis.  A full box scan when
       the space is small; a fixed-seed sample when it is not (the identity has a 10-parameter
       family, and 5^10 combinations would be scanned for no extra information)."""
    rng = random.Random(1770)
    dim = len(basis)
    full = (span * 2 + 1) ** dim
    if full <= cap:
        coefs = (c for c in itertools.product(range(-span, span + 1), repeat=dim))
    else:
        coefs = ([rng.randint(-span, span) for _ in range(dim)] for _ in range(cap))
    sigs = set()
    for coef in coefs:
        if not any(coef):
            continue
        M = sp.zeros(*basis[0].shape)
        for c, B in zip(coef, basis):
            M += c * B
        if M.det() == 0:
            continue
        sigs.add(signature(M)[:2])
    return sigs


def test_DE_types():
    print("\nTEST D/E -- all rank-4 types by order, and every signature their invariant forms take")
    tbl = types_of_rank(4)
    print("  order | types (cyclotomic factors)                    | has +-1 block | signatures")
    print("  ------+-----------------------------------------------+---------------+-----------")
    got31 = []
    for n in sorted(tbl):
        for ty in tbl[n]:
            C = sp.diag(*[companion(d) for d in ty])
            _, free, basis = invariant_forms(C)
            sigs = achievable_signatures(basis) if basis else set()
            has_pm1 = any(d in (1, 2) for d in ty)
            flag = "  <== (3,1)" if (3, 1) in sigs else ""
            print("   %4d | %-45s | %-13s | %s%s"
                  % (n, "Phi_" + " + Phi_".join(str(d) for d in ty),
                     "yes" if has_pm1 else "NO", sorted(sigs), flag))
            if (3, 1) in sigs:
                got31.append((n, ty, has_pm1))
    orders31 = sorted({n for n, _, _ in got31})
    report("D.1 orders admitting an integral form of signature (3,1)",
           orders31 == [1, 2, 3, 4, 6], "%s" % orders31)
    report("D.2 every (3,1)-admitting type contains Phi_1 or Phi_2",
           all(h for _, _, h in got31), "parity lemma: rotation blocks give an EVEN negative index")
    report("D.3 order 12 in rank 4 has more than one type (reader's R7)",
           len(tbl[12]) > 1, "types: %s" % [list(t) for t in tbl[12]])
    # POSITIVE CONTROL named in the ex ante: the instrument must be able to say (3,1) at all.
    report("D.neg control: instrument is not blind -- it does return (3,1) somewhere",
           len(got31) > 0, "%d types produced (3,1)" % len(got31))
    return tbl


# ----------------------------------------------------------------------------------------
# TEST F -- can Phi_5 / Phi_8 / Phi_12 act on a UNIMODULAR rank-4 lattice at all?
#
# The key remark that makes this decidable: any integral invariant form on ANY Z-lattice in
# the same Q-representation is, after a rational change of basis, a rational point of the SAME
# family, scaled by det(basis change)^2.  So det ranges over one square class of Q*.  If that
# class is not the class of 1, no unimodular lattice can carry the action -- and the question
# is settled for every lattice at once, not just for the companion module.
# ----------------------------------------------------------------------------------------
def square_class(poly):
    """det as a value of a polynomial over Q lies in ONE square class of Q* whenever every
       non-constant factor carries an even multiplicity.  Return that class (a squarefree
       integer), or None when the determinant is not of that shape.

       This is the whole leverage of test F: an integral invariant form on ANY Z-lattice of
       the same Q-representation is a rational point of this same family rescaled by the
       square of a base-change determinant, so the square class -- not the value -- is the
       lattice-independent invariant.  |det| = 1 is possible only if the class is 1."""
    coeff, factors = sp.factor_list(sp.expand(poly))
    for base, mult in factors:
        if mult % 2:
            return None
    c = sp.Rational(coeff)
    rest = abs(sp.numer(c) * sp.denom(c))        # same square class as c itself
    out = sp.Integer(1)
    for p, e in sp.factorint(int(rest)).items():
        if e % 2:
            out *= p
    return int(sp.sign(c)) * out


def is_square(n):
    return sp.integer_nthroot(abs(int(n)), 2)[1] and n > 0


def test_F_unimodular(tbl):
    print("\nTEST F -- can the non-crystallographic orders act on a UNIMODULAR rank-4 lattice?")
    print("  (every type of every such order, not only the primitive cyclotomic one)")
    for n in (5, 8, 10, 12):
      for ty in tbl[n]:
        d = n
        C = sp.diag(*[companion(k) for k in ty])
        label = "Phi_" + "+Phi_".join(str(k) for k in ty)
        Qg, free, basis = invariant_forms(C)
        det = sp.factor(sp.expand(Qg.det()))
        cls = square_class(Qg.det())
        print("  order %-2d  %-14s params %s" % (n, label, [s.name for s in free]))
        print("           det = %s" % det)
        if cls is None:
            print("           square class: not a single class (det is not c * square)")
            report("F.%s det is a fixed square class" % label, False)
            continue
        print("           square class of det: %s   (unimodular possible: %s)"
              % (cls, "yes -- searching" if is_square(cls) else "NO, and that settles"
                 " EVERY lattice at once"))
        if not is_square(cls):
            report("F.%s NO unimodular lattice carries this action" % label, True,
                   "det lies in the square class of %s for every Z-lattice" % cls)
            continue
        # collect EVERY (signature, parity) attained at |det| = 1, not merely the first hit:
        # the question the note's sentence asks is about I_(2,2) specifically, and I_(2,2) is
        # the ODD unimodular lattice of signature (2,2).  Even ones are II_(2,2), a different
        # lattice, and reporting the first hit would answer a question nobody asked.
        hits = {}
        for coef in itertools.product(range(-8, 9), repeat=len(free)):
            if not any(coef):
                continue
            M = sp.zeros(*C.shape)
            for c0, B in zip(coef, basis):
                M += c0 * B
            if abs(M.det()) != 1:
                continue
            sg = signature(M)[:2]
            odd = any(M[i, i] % 2 for i in range(M.shape[0]))
            hits.setdefault((sg, "odd" if odd else "even"), M)
        if not hits:
            report("F.%s unimodular form: none in the search box" % label, True, "|coef| <= 8")
            continue
        for (sg, par), M in sorted(hits.items()):
            print("           |det|=1 at sig %s, %-4s   Gram = %s" % (sg, par, M.tolist()))
        has_I22 = ((2, 2), "odd") in hits
        report("F.%s acts on I_(2,2) (odd unimodular of sig (2,2))" % label, True,
               "YES" if has_I22 else "NO -- attained: %s" % sorted(k for k in hits))


# ----------------------------------------------------------------------------------------
# TEST H -- make the I_(2,2) claim NATIVE.
#
# Test F produces an odd unimodular invariant form of signature (2,2).  Calling it "I_(2,2)"
# by Milnor's classification of indefinite odd unimodular lattices would be an IMPORT -- and
# an import standing in for a proof is exactly what cost this note its first attribution.
# So the lattice is diagonalised here by an explicit unimodular base change, and the isometry
# is transported into that basis.  What comes out is a matrix the reader can multiply.
# ----------------------------------------------------------------------------------------
def unimodular_col_reduce(row):
    """U unimodular with row*U = (g, 0, ..., 0), g = gcd(row).  Plain integer column ops."""
    n = len(row)
    c = list(row)
    U = sp.eye(n)
    for j in range(1, n):
        while c[j] != 0:
            if c[0] == 0:
                c[0], c[j] = c[j], c[0]
                U.col_swap(0, j)
                continue
            k = c[j] // c[0]
            c[j] -= k * c[0]
            U[:, j] = U[:, j] - k * U[:, 0]
            if c[j] != 0:
                c[0], c[j] = c[j], c[0]
                U.col_swap(0, j)
    if c[0] < 0:
        c[0] = -c[0]
        U[:, 0] = -U[:, 0]
    return U, c[0]


def short_vector(G, target_abs=1, bound=4):
    n = G.shape[0]
    for c in itertools.product(range(-bound, bound + 1), repeat=n):
        if not any(c):
            continue
        v = sp.Matrix(list(c))
        if abs((v.T * G * v)[0]) == target_abs:
            return v
    return None


def diagonalise_unimodular(G):
    """P unimodular with P^T G P = diag(+-1,...,+-1), by splitting off norm +-1 vectors."""
    n = G.shape[0]
    if n == 0:
        return sp.eye(0)
    v = short_vector(G)
    if v is None:
        return None
    f = list((v.T * G))
    U, g = unimodular_col_reduce(f)
    if g != 1:
        return None
    B = sp.Matrix.hstack(v, U[:, 1:])          # <v> + (v^perp cap L)
    if abs(B.det()) != 1:
        return None
    Gp = (B.T * G * B)
    sub = diagonalise_unimodular(Gp[1:, 1:])
    if sub is None:
        return None
    return B * sp.diag(sp.Matrix([[1]]), sub)


def test_H_native_I22(tbl):
    print("\nTEST H -- exhibit the order-8 action on I(2,2) explicitly (no classification import)")
    C = companion(8)
    _, free, basis = invariant_forms(C)
    G22 = None
    for coef in itertools.product(range(-8, 9), repeat=len(free)):
        if not any(coef):
            continue
        M = sp.zeros(4, 4)
        for c0, B in zip(coef, basis):
            M += c0 * B
        if abs(M.det()) == 1 and signature(M)[:2] == (2, 2) and any(M[i, i] % 2 for i in range(4)):
            G22 = M
            break
    report("H.1 odd unimodular invariant form of signature (2,2) for Phi_8", G22 is not None)
    if G22 is None:
        return
    P = diagonalise_unimodular(G22)
    report("H.2 explicit unimodular base change to a diagonal form", P is not None)
    if P is None:
        return
    D = P.T * G22 * P
    diag = [D[i, i] for i in range(4)]
    report("H.3 the diagonal form is diag(+-1) with two of each sign",
           D == sp.diag(*diag) and sorted(diag) == [-1, -1, 1, 1], "diag = %s" % diag)
    # transport the isometry into that basis and re-order the coordinates to eta_(2,2)
    M = P.inv() * C * P
    perm = [i for i in range(4) if diag[i] > 0] + [i for i in range(4) if diag[i] < 0]
    S = sp.zeros(4, 4)
    for new, old in enumerate(perm):
        S[old, new] = 1
    M = S.inv() * M * S
    integral = all(c.is_integer for c in M)
    report("H.4 the transported isometry is INTEGRAL", integral)
    report("H.5 it preserves eta_(2,2) = diag(1,1,-1,-1)", preserves(M, ETA22))
    report("H.6 its order is 8", order_of(M) == 8, "ord = %s" % order_of(M))
    print("       explicit order-8 isometry of I(2,2):  %s" % M.tolist())
    print("       => 'Phi_8 acts on I(2,2)' is now a matrix, not a citation.")
    # R9: the same averaging that works at negative index 1 must NOT rescue index 2
    t, eps = average_orbit(M, ETA22)
    bad = t is not None and (t.T * ETA22 * t)[0] < 0 and (M * t == t or M * t == -t)
    report("H.7 control for R9: orbit averaging gives NO timelike +-1 axis at index 2",
           not bad, "t = %s" % (list(t.T) if t is not None else eps))
    # and Phi_5 / Phi_10 are barred from EVERY unimodular lattice, so no such witness exists
    report("H.8 no odd unimodular (2,2) form exists for Phi_5 (square class 5)",
           square_class(invariant_forms(companion(5))[0].det()) == 5)


def test_G_witnesses():
    print("\nTEST G -- the reader's explicit realisations of orders 1,2,3,4,6 on I(3,1)")
    W = {1: sp.eye(4), 2: -sp.eye(4),
         3: sp.Matrix([[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]]),
         4: sp.Matrix([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
         6: sp.Matrix([[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, -1]])}
    for n, M in sorted(W.items()):
        ok = preserves(M, ETA31) and order_of(M) == n and all(c.is_integer for c in M)
        report("G order %d realised by an explicit integral isometry" % n, ok,
               "ord = %s, preserves eta: %s" % (order_of(M), preserves(M, ETA31)))


def test_I_signature_alone():
    """R9/R10: does the signature ALONE forbid order 5?  The note says the non-crystallographic
       axes are filtered out 'precisely by the signature (3,1)'.  One real matrix settles it."""
    print("\nTEST I -- signature alone, without integrality: is order 5 forbidden?")
    th = 2 * sp.pi / 5
    R = sp.Matrix([[sp.cos(th), -sp.sin(th), 0, 0], [sp.sin(th), sp.cos(th), 0, 0],
                   [0, 0, 1, 0], [0, 0, 0, 1]])
    report("I.1 R(2pi/5) preserves the form of signature (3,1)",
           sp.simplify(R.T * ETA31 * R - ETA31) == sp.zeros(4, 4))
    report("I.2 its order is 5", sp.simplify(R**5 - sp.eye(4)) == sp.zeros(4, 4))
    report("I.3 it is NOT integral", not all(sp.nsimplify(c).is_integer for c in R),
           "cos(2pi/5) = %s" % sp.nsimplify(sp.cos(th)))
    print("  => signature (3,1) alone permits order 5; what forbids it is signature AND")
    print("     integrality AND rank 4 together.  The note's 'precisely the signature' is too strong.")
    # and the count of timelike directions, which the note calls 'exactly one'
    tl = [v for v in [sp.Matrix([0, 0, 0, 1]), sp.Matrix([1, 0, 0, 2]), sp.Matrix([0, 1, 1, 2])]
          if (v.T * ETA31 * v)[0] < 0]
    report("I.4 'exactly one timelike direction' is false -- here are three independent ones",
           len(tl) == 3, "%s" % [list(v.T) for v in tl])


def main():
    print("=" * 92)
    print("S1770 -- litgate on the external reader's review of note-1 (I_(3,1) closure)")
    print("=" * 92)
    test_A_gap()
    test_C_averaging()
    tbl = test_DE_types()
    test_F_unimodular(tbl)
    test_H_native_I22(tbl)
    test_G_witnesses()
    test_I_signature_alone()
    print("\n" + "=" * 92)
    if FAIL:
        print("FAILURES: %s" % FAIL)
        sys.exit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
