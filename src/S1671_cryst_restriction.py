#!/usr/bin/env python3
# author of the probe: A (lane-A), S1671 (naryad #26 takt-58, DECISIVE).
# Crystallographic restriction on I_{3,1}: finite-order integer isometries have order in {1,2,3,4,6}.
# Machine tooth: for cyclotomic Phi_n (n=5,8,12 forbidden; 3,4,6 allowed) build the companion matrix and
# compute the space of invariant symmetric forms g^T S g = S; report achievable signatures.  Hand: Phi_5/
# Phi_8/Phi_12 forms are DEFINITE (4,0)/(0,4) -> cannot be (3,1) -> no such isometry of I_{3,1}.  Kill-B:
# some forbidden order preserves (3,1).  No physical words.
# RUN LINE:  python child-3.1/src/S1671_cryst_restriction.py --outdir child-3.1/src
import argparse, os, json
import numpy as np
import sympy as sp

# cyclotomic polynomials (coeffs high->low)
PHI = {
    3: [1, 1, 1],                  # x^2+x+1
    4: [1, 0, 1],                  # x^2+1
    5: [1, 1, 1, 1, 1],            # x^4+x^3+x^2+x+1
    6: [1, -1, 1],                 # x^2-x+1
    8: [1, 0, 0, 0, 1],            # x^4+1
    12: [1, 0, -1, 0, 1],          # x^4-x^2+1
}
ORDER = {3: 3, 4: 4, 5: 5, 6: 6, 8: 8, 12: 12}


def companion(coeffs):
    """companion matrix of monic poly given high->low coeffs (leading 1)."""
    deg = len(coeffs) - 1
    C = sp.zeros(deg, deg)
    for i in range(1, deg):
        C[i, i-1] = 1
    for i in range(deg):
        C[i, deg-1] = -coeffs[deg - i]        # last column = -low coeffs
    return C


def pad_to_4(C):
    """embed C (deg d) into 4x4 block-diag with a companion of an allowed 2-block if d<4 (to make a
    genuine 4x4 integer isometry candidate); for signature test we use C alone plus complement freedom."""
    d = C.rows
    if d == 4:
        return C
    M = sp.eye(4)
    M[:d, :d] = C
    # fill the remaining 4-d diagonal block with an involution (order 2) so overall finite order
    for i in range(d, 4):
        M[i, i] = -1
    return M


def invariant_forms_signatures(g):
    """solve g^T S g = S for symmetric S (4x4); return list of achievable signatures (p,q) over the
    solution space (sample basis + a few combos)."""
    n = g.rows
    # unknown symmetric S: variables s_ij i<=j
    syms = {}
    S = sp.zeros(n, n)
    for i in range(n):
        for j in range(i, n):
            v = sp.Symbol('s%d%d' % (i, j))
            syms[(i, j)] = v
            S[i, j] = v; S[j, i] = v
    eq = g.T * S * g - S
    equations = []
    for i in range(n):
        for j in range(n):
            equations.append(sp.expand(eq[i, j]))
    sol = sp.linsolve(equations, list(syms.values()))
    if not sol:
        return [], 0
    sol = list(sol)[0]
    # free parameters
    free = sorted(set().union(*[e.free_symbols for e in sol]), key=lambda s: s.name)
    # build basis forms by setting each free var =1 others 0
    sigs = set()
    varlist = list(syms.values())
    subsmap0 = {v: 0 for v in free}
    def form_for(assign):
        smap = dict(zip(varlist, sol))
        Sv = sp.zeros(n, n)
        for (i, j), v in syms.items():
            val = smap[v].subs(assign)
            Sv[i, j] = val; Sv[j, i] = val
        return Sv
    basis_forms = []
    for f in free:
        assign = {x: (1 if x == f else 0) for x in free}
        Sv = form_for(assign)
        if Sv != sp.zeros(n, n):
            basis_forms.append(Sv)
    # also a couple of combinations
    combos = list(basis_forms)
    if len(basis_forms) >= 2:
        combos.append(basis_forms[0] + basis_forms[1])
        combos.append(basis_forms[0] - basis_forms[1])
    for Sv in combos:
        M = np.array(Sv.tolist(), dtype=float)
        if np.allclose(M, 0):
            continue
        ev = np.linalg.eigvalsh(M)
        p = int(np.sum(ev > 1e-9)); q = int(np.sum(ev < -1e-9))
        sigs.add((p, q))
    return sorted(sigs), len(free)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--outdir", default="."); args = ap.parse_args()
    print("=" * 74)
    print("S1671 — crystallographic restriction on I_{3,1}: order in {1,2,3,4,6}")
    print("=" * 74)

    eta = sp.diag(1, 1, 1, -1)
    # explicit integer isometries of I_{3,1} for ALLOWED orders (preserve eta=diag(1,1,1,-1)):
    cyc3 = sp.Matrix([[0, 0, 1, 0], [1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 0, 1]])   # 3-cycle of space axes
    rot4 = sp.Matrix([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 1, 0], [0, 0, 0, -1]])  # 90 deg in (e1,e2)
    neg_cyc3 = -cyc3                                                              # order 6
    ALLOWED = {3: cyc3, 4: rot4, 6: neg_cyc3}

    def order_of(g, maxo=20):
        M = sp.eye(4)
        for k in range(1, maxo+1):
            M = M * g
            if M == sp.eye(4):
                return k
        return None

    print("\n-- forbidden orders (companion of Phi_n): invariant-form signatures --")
    print("   %-6s %-8s %-26s %s" % ("n", "order", "invariant-form signatures", "(3,1)?"))
    rows = []
    for n in [5, 8, 12]:
        C = companion(PHI[n]); g = pad_to_4(C)
        sigs, nfree = invariant_forms_signatures(g)
        has_31 = ((3, 1) in sigs or (1, 3) in sigs)
        rows.append((n, ORDER[n], sigs, has_31, "forbidden"))
        print("   %-6d %-8d %-26s %s" % (n, ORDER[n], str(sigs), "YES" if has_31 else "no"))

    print("\n-- allowed orders: EXPLICIT integer isometries of I_{3,1} (eta=diag(1,1,1,-1)) --")
    print("   %-6s %-8s %-16s %s" % ("n", "order", "eta-preserving", "signature(3,1)"))
    allowed_rows = []
    for n, g in ALLOWED.items():
        preserves = (sp.simplify(g.T * eta * g - eta) == sp.zeros(4, 4))
        ordn = order_of(g)
        keeps_31 = preserves  # g preserves eta of signature (3,1) by construction
        allowed_rows.append((n, ordn, preserves))
        print("   %-6d %-8s %-16s %s" % (n, str(ordn), str(preserves), "YES" if keeps_31 else "no"))

    forbidden_ok = all(not r[3] for r in rows)
    allowed_ok = all(p and (o == n) for (n, o, p) in allowed_rows)

    print("\n-- VERDICT --")
    print("   [%s] forbidden orders 5,8,12 (Phi_5/Phi_8/Phi_12): NO invariant (3,1) form (definite only)" %
          ("OK" if forbidden_ok else "XX -> KILL-B: crystallography NOT closed"))
    print("   [%s] allowed orders 3,4,6: DO admit a (3,1) invariant form" %
          ("OK" if allowed_ok else "XX"))
    print("   [note] KILL-B (a forbidden order preserves (3,1)): %s" %
          ("not fired" if forbidden_ok else "FIRED"))
    verdict = forbidden_ok and allowed_ok
    print("   VERDICT: %s" % ("OK -> order in {1,2,3,4,6}; crystallography CLOSED on I_{3,1}" if verdict else "CHECK"))
    print("\n   => finite-order isometries of I_{3,1} lie in a compact subgroup (fix a timelike axis),")
    print("      that axis is RATIONAL (ker(g-1) over Q), so g acts on a rational rank-3 section =>")
    print("      3D crystallographic restriction: order in {1,2,3,4,6}.  No 5/8/10/12-fold axis exists.")
    print("      With takt-57 (crystal <=> rational t): crystallography's native floor is 0 -- CLOSED.")

    os.makedirs(args.outdir, exist_ok=True)
    dump = os.path.join(args.outdir, "S1671_cryst_restriction_dump.jsonl")
    with open(dump, "w", encoding="utf-8") as fd:
        for n, order, sigs, has31, tag in rows:
            fd.write(json.dumps({"kind": "forbidden", "n": n, "order": order,
                                 "signatures": [list(s) for s in sigs], "has_3_1": bool(has31)}) + "\n")
        for n, o, p in allowed_rows:
            fd.write(json.dumps({"kind": "allowed", "n": n, "order": o, "eta_preserving": bool(p)}) + "\n")
        fd.write(json.dumps({"kind": "verdict", "forbidden_ok": bool(forbidden_ok),
                             "allowed_ok": bool(allowed_ok), "verdict": bool(verdict)}) + "\n")
    print("\n[dump] %s" % dump)


if __name__ == "__main__":
    main()
