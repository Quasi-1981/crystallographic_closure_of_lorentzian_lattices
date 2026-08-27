#!/usr/bin/env python3
# author of the probe: A (lane-A), S1670 (naryad #26 takt-57).
# Theorem: rank(Lambda cap t^perp) = 4 - r, r = dim_Q<t_1..t_4>.  Irrational t (r>=2) => no rank-3
# rational section => no 3D crystal without a window.  Crystallography's native floor is 0 (rational).
# Machine: symbolic rank of the integer kernel of x -> eta(x,t)=0 over Q(sqrt2,sqrt3) for controls.
# No physical words.  RUN LINE:  python child-3.1/src/S1670_rank_drop.py --outdir child-3.1/src
import argparse, os, json
import sympy as sp

s2, s3, s6 = sp.sqrt(2), sp.sqrt(3), sp.sqrt(6)
eta = sp.diag(1, 1, 1, -1)


def qdim(coords):
    """dim_Q of the Q-span of the coordinates, using a fixed Q-basis {1,sqrt2,sqrt3,sqrt6}."""
    basis = [sp.Integer(1), s2, s3, s6]
    # express each coord in basis; build rational coefficient matrix, rank over Q
    rows = []
    for c in coords:
        c = sp.expand(c)
        # coefficients wrt basis
        coeffs = []
        for bvec in basis:
            if bvec == 1:
                # rational part = c with all radicals removed
                coeffs.append(c.as_independent(s2, s3, s6, as_Add=True)[0])
            else:
                coeffs.append(c.coeff(bvec))
        rows.append([sp.nsimplify(x) for x in coeffs])
    M = sp.Matrix(rows)
    return M.rank()


def section_rank(t):
    """rank over Z (=over Q) of {x in Q^4 : eta(x,t)=0}: 4 - (number of independent rational constraints)."""
    t = sp.Matrix(t)
    c = eta * t                                  # coefficient vector of eta(x,t) = c . x
    # c . x = 0 ; split c_i into Q-basis components -> rational constraint rows
    basis = [sp.Integer(1), s2, s3, s6]
    constraint_rows = []
    for bvec in basis:
        row = []
        for i in range(4):
            ci = sp.expand(c[i])
            if bvec == 1:
                row.append(ci.as_independent(s2, s3, s6, as_Add=True)[0])
            else:
                row.append(ci.coeff(bvec))
        if any(x != 0 for x in row):
            constraint_rows.append([sp.nsimplify(x) for x in row])
    if not constraint_rows:
        return 4
    C = sp.Matrix(constraint_rows)
    return 4 - C.rank()


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--outdir", default="."); args = ap.parse_args()
    print("=" * 74)
    print("S1670 — rank(Lambda cap t^perp) = 4 - r  (irrational observer has no 3D crystal)")
    print("=" * 74)

    controls = [
        ("(1,1,1,2)      rational",        [1, 1, 1, 2],            1, 3),
        ("(1,sqrt2,0,0)  r=2",             [1, s2, 0, 0],           2, 2),
        ("(1,sqrt2,sqrt3,sqrt6)  r=4",     [1, s2, s3, s6],         4, 0),
        ("(2,3,-1,5)     rational (extra)",[2, 3, -1, 5],           1, 3),
        ("(1,sqrt2,1,0)  r=2 (extra)",     [1, s2, 1, 0],           2, 2),
    ]
    print("\n   %-32s %-6s %-8s %-8s %s" % ("t", "r(hand)", "rank(4-r)", "rank(meas)", "match"))
    ok = True
    rows = []
    for name, t, r_hand, rank_hand in controls:
        r = qdim(t)
        rk = section_rank(t)
        match = (int(r) == r_hand and int(rk) == rank_hand and int(rk) == 4 - int(r))
        ok = ok and match
        rows.append((name, int(r), int(rk), match))
        print("   %-32s %-6d %-8d %-8d %s" % (name, r, 4-r, rk, "OK" if match else "XX"))

    print("\n-- VERDICT --")
    print("   [%s] rank(Lambda cap t^perp) = 4 - r on all controls" % ("OK" if ok else "XX"))
    print("   [%s] 3D crystal (rank 3) <=> r=1 <=> t rational (up to scalar)" % ("OK" if ok else "XX"))
    print("   [note] KILL-A (rank != 4-r): %s" % ("not fired" if ok else "FIRED"))
    print("   VERDICT: %s" % ("OK -> irrational t (r>=2) has NO 3D crystal; needs a window" if ok else "CHECK"))
    print("\n   => crystallography's NATIVE FLOOR is 0 (rational/field level): a rank-3 crystal section")
    print("      exists iff the timelike direction t is RATIONAL (r=1).  Irrational (algebraic) t gives")
    print("      rank 4-r < 3 => no 3D lattice; a quasicrystal needs an external WINDOW (cut-and-project).")

    os.makedirs(args.outdir, exist_ok=True)
    dump = os.path.join(args.outdir, "S1670_rank_drop_dump.jsonl")
    with open(dump, "w", encoding="utf-8") as fd:
        for name, r, rk, match in rows:
            fd.write(json.dumps({"kind": "control", "t": name, "r": r, "rank": rk,
                                 "match": bool(match)}) + "\n")
        fd.write(json.dumps({"kind": "verdict", "all_match": bool(ok)}) + "\n")
    print("\n[dump] %s" % dump)


if __name__ == "__main__":
    main()
