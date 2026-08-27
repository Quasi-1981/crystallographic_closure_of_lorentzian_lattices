# Crystallographic Closure of Integral Lorentzian Lattices and the Rank at Which It Ends

**A finite-order integral isometry of any lattice of signature `(3,1)` has order in
`{1, 2, 3, 4, 6}`; for `I_{3,1}` all five are realized. On signature `(4,1)` this is no longer true.**

**Vladimir Sobol** · ORCID 0009-0006-9829-7931 · Independent Researcher

---

## 1 · Statement

> **Lemma 1 (integral timelike axis).** Let `Λ` be an integral lattice of signature `(n,1)`, `g ∈ O(Λ)` a
> finite-order isometry. Then there exists an **integral** vector `t ∈ Λ` with `q(t) < 0` and
> `ε ∈ {±1}` such that `g·t = ε·t`.

> **Theorem 1 (reduction).** In the notation of Lemma 1, the lattice `L := t^⊥ ∩ Λ` is **positive definite**,
> has rank `n`, is invariant under `g`, and
> `ord(g) = lcm( ord(g|_L), ord(ε) )`.
> Hence the set of orders on signature `(n,1)` **equals** the set of orders on positive definite
> lattices of rank `n`.

> **Theorem 2 (`n = 3`).** For **any** integral lattice of signature `(3,1)`:
> `ord(g) ∈ {1, 2, 3, 4, 6}`.

> **Corollary (sharpness).** For `I_{3,1}` **all five** are realized:
> `{ ord(g) : g ∈ O(I_{3,1}), ord(g) < ∞ } = {1, 2, 3, 4, 6}`.

**Unimodularity does no work anywhere.** It is needed only for the corollary — to name a specific
lattice on which all five orders are present. The upper bound rests on **integrality** and
**negative index 1**, and on nothing else.

This is the same restriction as in three-dimensional crystallography — but here it is proved for
a **Lorentzian** lattice, and the reason it survives is the content of §4. The reason it
**ends** at the next rank is the content of §5.

---

## 2 · Proof

### Lemma 1 — averaging over the finite orbit

**(L.1.1) `Λ` contains an integral timelike vector.** `Λ` spans `Λ⊗ℝ`, the condition `q < 0` is **open** and
nonempty (signature `(n,1)`), and `Λ⊗ℚ` is dense in `Λ⊗ℝ` ⟹ there exists a rational vector with `q < 0`;
multiplying by a common denominator gives an integral one. Denote it `e`.

**(L.1.2) Choice of sign.** The set `{q < 0}` has exactly **two** connected components (a cone and
its opposite). An isometry either preserves both or swaps them; for two timelike vectors
`u, v` they lie in the same component if and only if `⟨u,v⟩ < 0`. Hence there exists
`ε ∈ {±1}` such that `h := ε·g` **preserves** the component of `e`. Together with `g`, the matrix `h`
has finite order; let `m := ord(h)`.

**(L.1.3) Averaging.** Set

```
   t := Σ_{j=0}^{m−1} h^j · e
```

- `t` is **integral**: `h` is integral, `e` is integral;
- `t` is **timelike and nonzero**: all `h^j·e` lie in **one** component, and it is a
  **convex** cone, so the sum stays in it ⟹ `q(t) < 0`, in particular `t ≠ 0`;
- `h·t = t`, because `h^m = 1` only permutes the summands. Since `h = ε·g` and `ε² = 1`, we get
  `g·t = ε·t`. ∎

**Why averaging specifically.** The seemingly more natural route — from the compactness of the generated group to
an invariant **real** line, and then to its rationality — has a gap: rationality of a proper
*subspace* does not imply rationality of a specific *line* within it. A witness — rotation by
`π/2` in the plane `(x₁,x₂)` on `I_{3,1}`: the kernel `ker(g−1) = ⟨e₃,e₄⟩` is rational, but contains
an invariant timelike line `(0,0,1,√2)`, which is not rational. Along with it, the rank-3
section would fall too. Averaging gives an axis that is **integral** from the start, constructively, and without any theory
of compact groups.

### Theorem 1 — reduction

`t` is timelike ⟹ `t^⊥` is positive definite, so `L = t^⊥ ∩ Λ` is a positive definite integral lattice.
`t` is integral ⟹ `t^⊥` is rational ⟹ `rank L = n`. From `g·t = ±t` it follows that `g(t^⊥) = t^⊥`, hence
`g·L = L`.

The sum `L ⊕ ℤt` has **finite index** in `Λ`, and an isometry that is the identity on a finite-index
sublattice is the identity on all of `Λ` (both span the same ℚ-space). Therefore
`ord(g) = lcm( ord(g|_L), ord(g|_{ℤt}) )`, and `g|_{ℤt} = ε` has order 1 or 2. ∎

**Equality of sets — both inclusions.**

**(⊆)** Let `k := ord(g|_L)`. If `ε = +1` **or** `k` is even, then `ord(g) = lcm(k, ord ε) = k`,
and this order is already carried by `g|_L` itself on the positive definite `L`. If `ε = −1` **and** `k` is odd,
then `ord(g) = 2k`, and the isometry needed is given by `−g|_L`: it is integral, preserves the same `L`, and from
`(−g|_L)^m = (−1)^m · (g|_L)^m`, with `k` odd, the identity is first reached at `m = 2k`, i.e.
`ord(−g|_L) = 2k = ord(g)`.

**(⊇)** If `L₀` is a positive definite integral lattice of rank `n` with an isometry of order `k`, then
`L₀ ⊕ ⟨−1⟩` has signature `(n,1)` and carries the same isometry `⊕ 1`. ∎

### Lemma 2 — definite rank 3: a trace argument

> An integral isometry `h` of finite order of a **positive definite** lattice of rank 3 has
> `ord(h) ∈ {1,2,3,4,6}`.

Over `ℝ` such `h` lies in `O(3)`, hence is conjugate to `R_θ ⊕ (δ)`, `δ = det(h) ∈ {±1}`, where `R_θ` is
a rotation of the plane. The trace does not depend on the basis, and in the lattice basis the matrix is **integral**:

```
   tr(h) = δ + 2·cos θ ∈ ℤ   ⟹   2·cos θ ∈ ℤ   ⟹   cos θ ∈ {0, ±½, ±1}
```

These are exactly the possible orders of the rotation `r := ord(R_θ) ∈ {1,2,3,4,6}`. Next `ord(h) = lcm(r, ord(δ))`, and

```
     r        :  1   2   3   4   6
     lcm(2, r):  2   2   6   4   6      — all in {1,2,3,4,6}
```

so the set remains closed in the case `δ = −1` as well. ∎

**Three different quantities, not one.** `r` is the order of the planar **rotation**; `k = ord(h)` is the order of
the **whole** rank-3 isometry; `ord(g)` is the order of the original Lorentzian isometry. `lcm(2,·)`
is applied twice for different reasons: inside Lemma 2 (the `δ` block) and in Theorem 1 (the axis `ε`).

**Zero imports.** Here a reference to the classical crystallographic restriction could have stood
(Scherrer 1946, Coxeter 1969). It does not: the argument uses only the **integrality of the trace**, and this is
three lines. The rest of the classical facts are also derived on the spot: a finite-order integral matrix is semisimple with
roots of unity (the minimal polynomial divides `x^m − 1`, which is separable over `ℚ`);
a finite subgroup of `GL_n(ℤ)` preserves a positive definite form (averaging). Scherrer and Coxeter
remain a **historical reference, not a leg of the proof**.

### Theorem 2 and the corollary

Theorem 1 at `n = 3` together with Lemma 2 gives `ord(g) ∈ {1,2,3,4,6}`. Sharpness for `I_{3,1}` — by
exhibition; `η = diag(1,1,1,−1)`:

| order | integral isometry | why this order |
|---:|:--|:--|
| 1 | `I₄` | — |
| 2 | `−I₄` | — |
| 3 | cyclic permutation `x₁→x₂→x₃→x₁` | order of the cycle |
| 4 | rotation by `π/2` in the plane `(x₁,x₂)` | `lcm(4,1)` |
| 6 | the same permutation `⊕ (−1)` on `x₄` | `lcm(3,2) = 6` |

All five are integral, all preserve `η`. ∎

---

## 3 · Exact scope — in two tags

**The text must not claim more than the probe measured.** Hence:

| statement | tag | what supports it |
|:--|:--|:--|
| Lemma 1, Theorem 1, Lemma 2 | **derived** | §2, in full; zero imports |
| orders `1,2,3,4,6` are realized on `I_{3,1}` | **derived** | the five explicit matrices above — multiplied by hand |
| orders `5, 8, 10, 12` do not occur on signature `(3,1)` | **derived** | Theorem 2 (Lemma 1 + Lemma 2). The parity lemma (§4) gives only a **structural reason** — the mandatory `±1` block; the prohibition itself is completed by the fact that a **definite rank 3** remains in the orthogonal complement |
| the same, by an **independent route** | **measured** | enumeration of **all 24 types** of rank 4 and all reachable signatures |
| `Φ₅`, `Φ₁₀` do not occur on **any unimodular** lattice of rank 4 | **measured** | square class of `det` = 5 (§4) |
| `Φ₈` acts on `I_{2,2}` | **measured** | explicit matrix of order 8 (§4) |
| table across ranks 2–8 | **measured** | structural count, checked against linear algebra type by type |

**Machine check against a known answer:** enumerating integral isometries gave `ℤ³` — 48 elements,
orders `{1,2,3,4,6}`; `A₃` — 48, the same; hexagonal — 24, `{1,2,3,6}`. The union is
`{1,2,3,4,6}`, **no order 5 at all**. This is a **check**, not a leg: the verdict is carried by the trace argument.

---

## 4 · What actually filters — and what does not

### The mechanism — three parts, not one

> **negative index 1** + **integrality** + **rank 4**

The roles are separated:

1. **negative index 1** makes the negative block one-dimensional ⟹ a finite-order action on it is
   only `±1` ⟹ **the parity lemma**: the isotypic components of different real irreducible constituents
   are orthogonal for any invariant form, while on a rotation block (`θ ≠ 0, π`) invariant
   symmetric forms correspond to **Hermitian** forms, and hence give an **even** negative index.
   An odd index ⟹ a `±1` block **must** be present;
2. **integrality** turns this axis into an **arithmetic** one (Lemma 1) and leaves an integral
   lattice in the orthogonal complement;
3. **rank 4** turns this lattice into rank 3, where the trace quantizes `2cos θ`.

**Signature alone filters nothing.** In the real group `O(3,1)` an element of order 5 does exist:
rotation by `2π/5` in a spatial plane preserves `diag(1,1,1,−1)` and has order exactly 5. It is
simply **not integral** (`cos(2π/5) = (√5−1)/4`).

**"Exactly one timelike direction" is false.** There are infinitely many timelike directions
(`(0,0,0,1)`, `(1,0,0,2)`, `(0,1,1,2)` — already three independent ones). The correct formulation is:
**the maximal negative definite subspace is one-dimensional**, i.e. `negative index = 1`.

### Non-crystallographic axes: three different verdicts, not one

The tempting phrase "`Φ₅`, `Φ₈`, `Φ₁₂` can act on `I_{2,2}`, but not on `I_{3,1}`" merges
**three statements of different truth-status**. One invariant distinguishes them.

**Square class.** An integral invariant form on **any** ℤ-lattice within the same
ℚ-reduction is a rational point of the same family, scaled by `det(P)²` of the change of basis. Hence the class of
`det` in `ℚ*/(ℚ*)²` **does not depend on the lattice**, and the question "does a unimodular one occur" is settled for
all lattices at once — by a single number, without enumeration:

| type | `det` of the family of invariant forms | class | unimodular lattice |
|:--|:--|:--:|:--|
| `Φ₅` | `5(4p₀² + 2p₀p₁ − p₁²)² / 16` | **5** | **impossible, none** |
| `Φ₁₀` | `5(4p₀² − 2p₀p₁ − p₁²)² / 16` | **5** | **impossible, none** |
| `Φ₃+Φ₄`, `Φ₄+Φ₆` | `3p₀²p₁² / 4` | **3** | **impossible, none** |
| `Φ₈` | `(2p₀² − p₁²)²` | 1 | **exists**, and odd at `(2,2)` |
| `Φ₁₂` | `(4p₀² − 3p₁²)² / 16` | 1 | exists, but found to be **even** |

Hence:

- **order 5** (and 10) does not act on `I_{2,2}` — and does not act on any unimodular lattice of rank 4;
- **order 8** does act, and precisely on `I_{2,2}`. To avoid relying on the classification of indefinite
  odd unimodular lattices, the form was diagonalized by an **explicit unimodular change of basis**
  to `diag(1,1,−1,−1)`, and the isometry was carried into this basis:

```
   M = [ −3  −3  −1  −4 ]
       [  3  −3  −4   1 ]      MᵀηM = η for η = diag(1,1,−1,−1),  ord(M) = 8
       [ −4   1   3  −3 ]
       [  1   4   3   3 ]
```

- **order 12**: a unimodular invariant form of signature `(2,2)` exists, but the one found is
  **even**, i.e. it is `II_{2,2}`, not `I_{2,2}`. Whether an odd one exists is **open**: parity
  also depends on the lattice, so a search within a bounded box proves nothing here.

**Why this is not a bookkeeping detail.** The "safe" formulation ("they admit integral nondegenerate
invariant forms of signature (2,2), but not (3,1)") is true — and it blurs precisely this distinction,
leaving the reader with the impression that all three cases are the same.

**Answer to the natural question "what about quasiperiodic axes?"**: in this arithmetic there are none, and
the reason is named explicitly — not "not found", but **negative index 1 together with integrality in
rank 4 leaves no room for them**.

---

## 5 · Where the closure ends

Theorem 1 reduces the Lorentzian question of rank `n+1` to the definite question of rank `n`. Measured for
ranks 2–8; at rank 4 the structural count was checked against linear algebra type by type:

| rank | signature | orders |
|---:|:--|:--|
| 2 | `(1,1)` | `1, 2` |
| 3 | `(2,1)` | `1, 2, 3, 4, 6` |
| **4** | **`(3,1)`** | **`1, 2, 3, 4, 6`** |
| 5 | `(4,1)` | `1, 2, 3, 4, 5, 6, 8, 10, 12` |
| 6 | `(5,1)` | `1, 2, 3, 4, 5, 6, 8, 10, 12` |
| 7 | `(6,1)` | `1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 12, 14, 15, 18, 20, 24, 30` |
| 8 | `(7,1)` | same as at rank 7 |

**The closure breaks exactly at signature `(4,1)`**, and it breaks substantively: `Φ₅ ⊕ Φ₁` has
degree 5 and signature `(4,0) + (0,1) = (4,1)` — the fivefold axis goes through as soon as the definite
complement has enough rank to hold it. Next comes a **plateau** (rank 6 adds nothing), and the next
break is at `(6,1)`, where order 7 appears, because `φ(7) = 6`.

**An identity, read off from the table and verified for ranks 3–8:**

```
   { orders on signature (n−1,1) }  =  { orders on a definite lattice of rank n−1 }
```

This is Theorem 1 seen from the other side: **for the spectrum of finite orders**, the Lorentzian problem
`(n,1)` adds nothing beyond the definite question of rank `n`. **The formulation is deliberately narrow:** it is
about the **set of orders**, not about the classification of Lorentzian lattices in general. Within these
limits, "the crystallographicity of signature `(3,1)`" is not a property of Lorentzian-ness, but a property of **the number 3**:
rank 3 turns out to be the **last definite rank** where the spectrum still equals the classical
`{1,2,3,4,6}`.

---

## 6 · Limits

- **Multiplicity 1.** This statement and the three-dimensional crystallographic restriction are **not two witnesses**:
  Theorem 1 **reduces** the former to the latter, a common ancestor ⟹ multiplicity 1. The same holds for the two routes
  to the prohibition of order 5 — the signature one and the trace one: both come out of **one and the same** integral
  arithmetic. This is **corroboration**, not independent confirmation.
- **Window (quasiperiodicity) is a separate input, not a consequence.** If a structure admits a window, this is
  a condition introduced **from outside**; what is measured here is what integral Lorentzian arithmetic allows.
- **No physical reading whatsoever.** The term "timelike" is used **only** in the sense of the sign of the
  quadratic form; no interpretation of a time coordinate is introduced.
- **Open:** whether order 12 acts on `I_{2,2}` (§4). A general description of the square class of `det`
  for an arbitrary cyclotomic type was not sought here.
- **Zero imports in the chain.** What is telling is **how** the last debt was closed: the first
  attribution attached to the step of Lemma 2 ("Hermann–Mauguin") turned out to be **the name of a notation mistaken for
  a theorem**, and was removed; instead of looking for a replacement, the step was derived. It proved cheaper to **prove**
  than to cite correctly. The same thing happened a second time in §4: the classification of unimodular
  lattices was replaced by explicit diagonalization.

---

## 7 · Relation to preprint-7

| what | where |
|:--|:--|
| full statement, proof, scope | **here** |
| a line about the scope, **without a reference here** | P7 §11(d) |
| a reference **from here to P7** — legitimate | P7 is published, concept DOI `10.5281/zenodo.22125370` |
| a reference **from P7 to here** — not yet | waits for this work's address |

**A forward-reference is not placed:** a reference to a nonexistent artifact is broken from day one and
is fixed only by a new version. The asymmetry is deliberate and temporary: when this work receives an address,
the line at P7 §11(d) will receive a reference **via a new version of P7**, not by editing what is published.

**Preprint-7:** *Section Criteria for Integral Lorentzian Lattices and the Cost of Realization*,
Zenodo [10.5281/zenodo.22125370](https://doi.org/10.5281/zenodo.22125370).

---

## Acknowledgement

This text went through **two rounds of external review**. The first pointed to a gap in the proof of the
rationality of the axis, to an incomplete measurement for order 12, and to an overstated formulation about `I_{2,2}`;
the second — to the fact that Theorem 1 claimed equality of sets while only one inclusion was proved, and to two
overstated formulations in §3 and §5. All remarks were accepted: the proof route was rebuilt, §4
was rewritten, the inclusion `⊆` was added. The errors that remain are the author's.
