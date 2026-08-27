#!/usr/bin/env python3
# author: B (lane-B chair a7296aa8), S1771.  Finite orders of integral isometries, by rank and
# signature -- the instrument behind the note's scope table.
#
# WHAT IT COMPUTES.  For a given rank n: every type (multiset of cyclotomic factors of total
# degree n, i.e. every rational representation a finite-order integral matrix can carry), the
# order it realises, and EVERY signature an integral invariant form on it can take.  The answer
# to "which orders live on a lattice of signature (p,q)" is then read off the table.
#
# WHY STRUCTURALLY AND NOT BY SEARCH.  S1770 sampled signatures inside a box.  A sample can
# witness "achievable"; it can never establish "not achievable", and reading silence as a
# verdict is exactly the VOID-BLIND failure of S1756.  Here the signature set is DERIVED:
#
#   * isotypic components of distinct real irreducible constituents are orthogonal for any
#     invariant form, so the signature is the sum of the components' signatures;
#   * the Phi_1 component (eigenvalue +1) of multiplicity m carries an arbitrary symmetric
#     form -- signature (a, m-a) for every a;  same for Phi_2 (eigenvalue -1);
#   * a rotation component (d >= 3, angle != 0, pi) of multiplicity m has invariant symmetric
#     forms in bijection with HERMITIAN forms on C^m, so its signature is (2a, 2(m-a)) --
#     always EVEN on both sides;
#   * Phi_d for d >= 3 splits into phi(d)/2 components, one per rotation angle, and distinct d
#     never share an angle.
#
# The parity lemma follows: an odd negative index forces a +-1 block.  The whole content of
# "signature (n-1,1) is crystallographic" sits in that one word, odd.
#
# CONTROL.  The structural count is cross-checked, type by type, against the linear-algebra
# measurement of S1770 (solve C^T Q C = Q, then sample).  Disagreement on any type voids the
# table.  See --check.
# Handles: 0.
# RUN LINE:  python S1771_lattice_orders.py              (table for ranks 2..8 + controls)
#            python S1771_lattice_orders.py --md FILE    (also write the table as markdown)
import argparse
import itertools
import random
import sys

import numpy as np
import sympy as sp

x = sp.Symbol("x")
FAIL = []


def report(name, ok, detail=""):
    mark = "OK  " if ok else "FAIL"
    if not ok:
        FAIL.append(name)
    print("  [%s] %-64s %s" % (mark, name, detail))


# ----------------------------------------------------------------------------------------
# types
# ----------------------------------------------------------------------------------------
def types_of_rank(rank):
    """every multiset of cyclotomic factors of total degree exactly `rank`, keyed by the order
       (the lcm) it realises."""
    ds = [d for d in range(1, 200) if sp.totient(d) <= rank]
    out = {}

    def rec(start, left, acc):
        if left == 0:
            n = 1
            for d in acc:
                n = int(sp.ilcm(n, d))
            out.setdefault(n, []).append(tuple(acc))
            return
        for i in range(start, len(ds)):
            d = ds[i]
            if sp.totient(d) <= left:
                rec(i, left - int(sp.totient(d)), acc + [d])

    rec(0, rank, [])
    return out


def components(ty):
    """the real-irreducible components of a type, as a list of (kind, multiplicity).
       kind 'pm1' contributes any signature (a, m-a); kind 'rot' contributes (2a, 2(m-a))."""
    mult = {}
    for d in ty:
        mult[d] = mult.get(d, 0) + 1
    comps = []
    for d, m in sorted(mult.items()):
        if d in (1, 2):
            comps.append(("pm1", m))
        else:
            for _ in range(int(sp.totient(d)) // 2):   # one component per rotation angle
                comps.append(("rot", m))
    return comps


def signatures_of(ty):
    """EVERY signature an integral invariant form on this type can have.  Derived, not sampled."""
    sigs = {(0, 0)}
    for kind, m in components(ty):
        step = 1 if kind == "pm1" else 2
        span = m if kind == "pm1" else m
        add = {(step * a, step * (span - a)) for a in range(span + 1)}
        sigs = {(p + dp, q + dq) for (p, q) in sigs for (dp, dq) in add}
    return sigs


# ----------------------------------------------------------------------------------------
# the control: the same question answered by linear algebra, the S1770 way
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


def invariant_basis(C):
    """Z-scaled basis of the space of C-invariant symmetric bilinear forms."""
    n = C.shape[0]
    keys = [(i, j) for i in range(n) for j in range(i, n)]
    syms = {k: sp.Symbol("q_%d_%d" % k) for k in keys}
    Q = sp.zeros(n, n)
    for (i, j), s in syms.items():
        Q[i, j] = s
        Q[j, i] = s
    A, _ = sp.linear_eq_to_matrix(list(C.T * Q * C - Q), [syms[k] for k in keys])
    basis = []
    for v in A.nullspace():
        B = sp.zeros(n, n)
        for idx, (i, j) in enumerate(keys):
            B[i, j] = v[idx]
            B[j, i] = v[idx]
        den = sp.Integer(1)
        for e in B:
            if e != 0:
                den = sp.ilcm(den, sp.denom(sp.nsimplify(e)))
        basis.append(np.array((B * den).tolist(), dtype=float))
    return basis


def sampled_signatures(basis, span=2, cap=6000, seed=1771):
    """the S1770 method: integral combinations in a box, signature by eigenvalues."""
    rng = random.Random(seed)
    dim = len(basis)
    full = (2 * span + 1) ** dim
    if full <= cap:
        coefs = itertools.product(range(-span, span + 1), repeat=dim)
    else:
        coefs = ([rng.randint(-span, span) for _ in range(dim)] for _ in range(cap))
    out = set()
    for coef in coefs:
        if not any(coef):
            continue
        M = sum(c * B for c, B in zip(coef, basis))
        ev = np.linalg.eigvalsh((M + M.T) / 2)
        tol = 1e-9 * max(1.0, float(np.max(np.abs(ev))))
        if np.any(np.abs(ev) <= tol):
            continue
        out.add((int(np.sum(ev > tol)), int(np.sum(ev < -tol))))
    return out


def control_rank4():
    """structural vs linear-algebraic, type by type, at rank 4 -- the rank S1770 measured."""
    print("\nCONTROL -- structural signature set vs the linear-algebra measurement (rank 4)")
    tbl = types_of_rank(4)
    bad = []
    for n in sorted(tbl):
        for ty in tbl[n]:
            derived = signatures_of(ty)
            measured = sampled_signatures(invariant_basis(sp.diag(*[companion(d) for d in ty])))
            if not measured <= derived:
                bad.append((ty, sorted(measured - derived)))
    report("C.1 every measured signature is predicted by the structural count", not bad,
           "%d types checked%s" % (sum(len(v) for v in tbl.values()),
                                   "" if not bad else "; offenders %s" % bad))
    # and the other direction, where sampling can reach it: the small types must be attained
    miss = []
    for n in sorted(tbl):
        for ty in tbl[n]:
            if len(ty) > 2:
                continue
            derived = signatures_of(ty)
            measured = sampled_signatures(invariant_basis(sp.diag(*[companion(d) for d in ty])))
            if measured != derived:
                miss.append((ty, sorted(derived - measured)))
    report("C.2 on the small types the measurement attains the whole predicted set", not miss,
           "" if not miss else "unattained: %s" % miss)


# ----------------------------------------------------------------------------------------
# the table
# ----------------------------------------------------------------------------------------
def orders_with_signature(rank, q):
    """orders of finite-order integral isometries realisable on a lattice of signature
       (rank - q, q), with the types that realise them."""
    tbl = types_of_rank(rank)
    out = {}
    target = (rank - q, q)
    for n in sorted(tbl):
        for ty in tbl[n]:
            if target in signatures_of(ty):
                out.setdefault(n, []).append(ty)
    return out


def fmt_type(ty):
    return "+".join("P%d" % d for d in ty)


def build_table(ranks):
    rows = []
    for rank in ranks:
        lor = orders_with_signature(rank, 1)
        defi = orders_with_signature(rank, 0)
        rows.append((rank, sorted(defi), sorted(lor), lor))
    return rows


def order_of(M, cap=200):
    n = M.shape[0]
    P = sp.eye(n)
    for k in range(1, cap + 1):
        P = P * M
        if P == sp.eye(n):
            return k
    return None


def control_inclusion():
    """The '<=' half of Theorem 1: a Lorentzian order is ALREADY realised on the definite
       complement.  The only case needing an argument is eps = -1 with k odd, where the
       realiser is -g|_L; the step it rests on is  ord(-A) = 2k  for  ord(A) = k  odd.
       Stated in the text as two lines -- checked here, because a freshly added leg of a
       proof is exactly the kind of thing that gets asserted and never run."""
    print("\nCONTROL -- Theorem 1, the inclusion that needs -g|_L  (ord(-A) = 2k for k odd)")
    cases = [("identity, rank 3", sp.eye(3)),
             ("3-cycle, rank 3", sp.Matrix([[0, 0, 1], [1, 0, 0], [0, 1, 0]])),
             ("companion(Phi_5), rank 4", companion(5)),
             ("companion(Phi_3), rank 2", companion(3))]
    bad = []
    for label, A in cases:
        k = order_of(A)
        km = order_of(-A)
        want = 2 * k if k % 2 else None      # even k is not the case the text claims
        if k % 2 and km != want:
            bad.append((label, k, km))
        print("       %-26s ord(A) = %-3s  ord(-A) = %-3s  %s"
              % (label, k, km, "(k odd -> expect %d)" % (2 * k) if k % 2 else "(k even)"))
    report("I.1 ord(-A) = 2k whenever ord(A) = k is odd", not bad, "" if not bad else "%s" % bad)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--md", default=None, help="also write the table as markdown")
    ap.add_argument("--max-rank", type=int, default=8)
    args = ap.parse_args()

    print("=" * 94)
    print("S1771 -- finite orders of integral isometries by rank and signature")
    print("=" * 94)

    control_rank4()
    control_inclusion()

    ranks = list(range(2, args.max_rank + 1))
    rows = build_table(ranks)

    print("\nTABLE -- which orders a lattice of the given signature admits")
    print("  rank | definite (n,0)                     | LORENTZIAN (n-1,1)")
    print("  -----+------------------------------------+-----------------------------------")
    for rank, defi, lor, _ in rows:
        print("  %4d | %-34s | %s" % (rank, defi, lor))

    # THE IDENTITY.  The Lorentzian question at rank n is the DEFINITE question at rank n-1 --
    # which is the general theorem itself, read off the table: an odd negative index forces a
    # +-1 axis, its complement is definite of one rank less, and lcm(2, .) adds nothing new.
    by_rank = {r: (defi, lor) for r, defi, lor, _ in rows}
    ident = [r for r in ranks if r - 1 in by_rank and by_rank[r][1] != by_rank[r - 1][0]]
    report("T.0 orders on (n-1,1) = orders on the definite lattice of rank n-1", not ident,
           "checked ranks %s%s" % ([r for r in ranks if r - 1 in by_rank],
                                   "" if not ident else "; mismatch at %s" % ident))

    cry = [1, 2, 3, 4, 6]
    report("T.rank 2 (signature (1,1)) is {1,2}: the complement has rank 1",
           by_rank[2][1] == [1, 2], "%s" % by_rank[2][1])
    for rank, defi, lor, _ in rows:
        if 3 <= rank <= 4:
            report("T.rank %d Lorentzian set is crystallographic" % rank, lor == cry, "%s" % lor)
        if rank == 5:
            report("T.rank 5 BREAKS the crystallographic set", lor != cry, "%s" % lor)
            report("T.rank 5 gains exactly 5, 8, 10, 12",
                   sorted(set(lor) - set(cry)) == [5, 8, 10, 12], "%s" % sorted(set(lor) - set(cry)))
    if len(rows) >= 5:
        r5 = dict((r, l) for r, _, l, _ in rows)
        if 6 in r5:
            report("T.rank 6 is a plateau (same set as rank 5)", r5[6] == r5[5], "%s" % r5[6])
        if 7 in r5:
            report("T.rank 7 breaks again (order 7 appears)", 7 in r5[7], "%s" % r5[7])

    # the parity lemma, checked on every type of every rank in range
    off = []
    for rank in ranks:
        for n, tys in types_of_rank(rank).items():
            for ty in tys:
                odd_q = any(q % 2 for (_, q) in signatures_of(ty))
                has_pm1 = any(d in (1, 2) for d in ty)
                if odd_q != has_pm1:
                    off.append((rank, ty))
    report("P.1 odd negative index  <=>  the type contains Phi_1 or Phi_2", not off,
           "all types, ranks %d..%d%s" % (ranks[0], ranks[-1],
                                          "" if not off else "; offenders %s" % off[:4]))

    if args.md:
        with open(args.md, "w", encoding="utf-8") as fh:
            fh.write("# Finite orders of integral isometries, by rank and signature\n\n")
            fh.write("Generated by `src/S1771_lattice_orders.py` -- do not edit by hand.\n\n")
            fh.write("An order appears in a row when some integral lattice of that rank and\n")
            fh.write("signature carries an isometry of that order.  The negative index, not the\n")
            fh.write("rank, is what decides: an odd negative index forces a `+-1` eigenvector,\n")
            fh.write("and the orthogonal complement is a definite lattice of one rank less.\n\n")
            fh.write("| rank | definite `(n,0)` | Lorentzian `(n-1,1)` |\n|---:|:--|:--|\n")
            for rank, defi, lor, _ in rows:
                fh.write("| %d | `%s` | `%s` |\n" % (rank, defi, lor))
            fh.write("\n## Which types realise each Lorentzian order\n\n")
            for rank, _, lor, det in rows:
                fh.write("\n**rank %d, signature (%d,1)**\n\n" % (rank, rank - 1))
                fh.write("| order | types |\n|---:|:--|\n")
                for n in sorted(det):
                    fh.write("| %d | %s |\n" % (n, ", ".join("`%s`" % fmt_type(t)
                                                             for t in det[n])))
        print("\nwritten: %s" % args.md)

    print("\n" + "=" * 94)
    if FAIL:
        print("FAILURES: %s" % FAIL)
        sys.exit(1)
    print("all checks passed")


if __name__ == "__main__":
    main()
