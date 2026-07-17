# Known-sire / full-sib support — core design (branch `m6-known-sire-core`)

*Development design note for independent review. This describes the CORE oracle-phase
extension only. It is not a validated biological method and does not change any
published claim. Items marked "deferred" belong to a later branch.*

## 0. Goal and non-goals

Add support for **known-sire full-sib crosses** alongside the existing
**open-pollinated (OP)** half-sib model, sharing a single maternal recombination map
and (for full-sib crosses) estimating a **separate paternal recombination map**.

A known heterozygous sire transmits *linked* paternal haplotypes, so a correct model
needs a **paternal hidden inheritance path**, not an independent per-marker paternal
allele frequency. Fixing paternal `q ∈ {0, 0.5, 1}` per marker is explicitly **not** a
correct full-sib model and is not used here.

**In scope:** cross-aware data model; sire retention/validation; full-sib + mixed
simulator; oracle-phase four-state full-sib HMM; joint engine combining OP and
known-sire crosses; tests; docs; this design.

**Deferred (next branch):** full-sib two-point estimation; automatic sire-phase
inference; full-sib linkage-group construction; real-data analysis; publication
simulations; parent-specific *maps per parent instance* beyond one shared maternal +
one shared paternal map; polyploidy.

## 1. Supported family types

| `family_type` | father | father genotype | model used |
|---|---|---|---|
| `open_pollinated`   | unknown           | —          | 2-state maternal HMM + dam-specific paternal `q` (current model) |
| `known_sire_genotyped` | known ID       | available  | **4-state full-sib HMM** (this design) |
| `known_sire_untyped`   | known ID       | missing    | falls back to the OP model **only on explicit request**, with a status/warning; never silently treated as full-sib |

`known_sire_untyped` is never silently promoted to full-sib. The reader records the
type; the fitting API refuses to treat an untyped sire as genotyped unless the caller
opts into the OP fallback (`untyped_sire = "open_pollinated"`), which is reported.

## 2. Cross vs. parent identity

A population keyed only by mother ID is insufficient when a mother has offspring from
several fathers. We introduce an explicit **cross**:

```
cross_id = mother_id  ×  father_id
```

with a defined unknown-sire token for OP crosses (`father_id = "<unknown>"`,
`cross_id = "<mother>__x__<unknown>"`). Parent genotypes are stored **once per parent
ID**; crosses reference parents by ID. Consequences enforced by construction:

- the **same mother** across crosses shares one genotype, one phase vector, and one
  maternal recombination map;
- the **same sire** across crosses shares one genotype, one phase vector, and one
  paternal recombination map.

## 3. Phase representation (applied ONCE)

The current OP engine uses the **allele-state** parameterization: the hidden state is
the transmitted maternal *allele*, the transition carries the phase (coupling preserves
/ repulsion switches), and the emission is phase-free. This cannot factor into a
maternal⊗paternal Kronecker transition.

The full-sib engine therefore uses the **homolog-state** parameterization (the other
standard, mathematically equivalent HMM form):

- hidden state = transmitted **homolog index** (per parent), so the transition is
  **phase-free** and factors as a Kronecker product (Section 5);
- **phase enters only the emission**, through the parent's *phased haplotype allele
  labels*.

Phase is thus applied exactly once (in the emission), never twice. A parent's oracle
phase is represented as a `2 × z` **allele matrix** `H` (`H[j,k] ∈ {0,1}` = allele on
homolog `j` at marker `k`), obtained from the parent genotype `g` and an adjacent phase
vector `ψ` (`1` = coupling, `0` = repulsion) by the package's standard anchored
construction:

```
orient[1] = 0;  orient[k] = orient[k-1] XOR (ψ[k-1] == repulsion)
g[k]==2 (AA): H[,k] = (1,1)
g[k]==0 (aa): H[,k] = (0,0)
g[k]==1 (Aa): H[1,k] = (orient[k]==0 ? 1 : 0), H[2,k] = 1 - H[1,k]
```

A helper `phase_to_haplotypes(g, ψ)` implements this; the simulator emits the true `H`
directly; a test verifies the two agree. Because a physical maternal crossover is a
homolog-index switch here and an allele-switch-under-coupling in the allele-state OP
engine, **the two forms count the same maternal recombination events**, so maternal
counts pool validly across OP and full-sib crosses (Section 6).

## 4. Hidden states, initial probabilities, emissions

For a genotyped mother (haplotypes `Hm`, `2×z`) and genotyped sire (`Hp`, `2×z`), the
per-offspring hidden state at marker `k` is

```
Z_ik = (h_m, h_p),  h_m,h_p ∈ {1,2}   → 4 states (1,1),(1,2),(2,1),(2,2)
```

encoded `s = 2*(h_m-1) + (h_p-1) ∈ {0,1,2,3}`.

**Initial:** `P(Z_i1) = 1/4` for each state (uninformative; documented default).

**Emission.** The expected offspring genotype (dosage) in state `(h_m,h_p)` at marker
`k` is the sum of the transmitted alleles `d = Hm[h_m,k] + Hp[h_p,k] ∈ {0,1,2}`. With a
symmetric error `ε`,

```
b_k(y | s) = 1                         if y is missing
           = 1 - ε                     if y == d
           = ε/2                       otherwise (each of the two other classes)
```

No paternal `q` appears in the full-sib emission.

## 5. Transitions (Kronecker)

With maternal and paternal per-interval recombination fractions `r_m,k`, `r_p,k`,

```
T_m,k = [[1-r_m,k, r_m,k],[r_m,k, 1-r_m,k]]
T_p,k = [[1-r_p,k, r_p,k],[r_p,k, 1-r_p,k]]
T_k   = T_m,k ⊗ T_p,k                                  (4×4)
T_k( (h_m',h_p') | (h_m,h_p) ) = T_m,k(h_m'|h_m) · T_p,k(h_p'|h_p)
```

Maternal and paternal `r` are **separate**; they are never assumed equal.

## 6. EM sufficient statistics and parameter sharing

Per offspring, from the scaled forward–backward posteriors `ξ_k(s,s')`:

```
maternal switch  m_switch(k) = Σ_i Σ_{s,s': h_m≠h_m'} ξ_k(s,s')
paternal switch  p_switch(k) = Σ_i Σ_{s,s': h_p≠h_p'} ξ_k(s,s')
total(k)         = Σ_i Σ_{s,s'} ξ_k(s,s')   (= #offspring at observed loci)
```

**M-step (consensus maps, shared across meioses):**

```
r_m,k = ( Σ_{OP crosses} diff_OP(k) + Σ_{full-sib crosses} m_switch(k) )
        / ( Σ_{OP} (same_OP+diff_OP)(k) + Σ_{full-sib} total(k) )
r_p,k = ( Σ_{full-sib crosses} p_switch(k) ) / ( Σ_{full-sib} total(k) )
```

where `same_OP/diff_OP` are the existing OP engine's expected non-recombinant/
recombinant maternal counts (a maternal crossover = allele-switch-under-coupling =
homolog-index switch, so they are the same quantity). OP crosses contribute **only** to
`r_m` and continue to update their dam-specific paternal `q`. Full-sib crosses
contribute to **both** `r_m` and `r_p`. Bounds: `r ∈ [ε_r, 0.5]`, numerically capped
just below 0.5, consistent with the OP engine.

**Joint likelihood** = sum of the correct family-specific observed-data log-likelihoods
(4-state for full-sib, 2-state for OP). **Active objective** = observed-data
log-likelihood (`λ=0`) or observed-data LL + the OP `q` pseudocount penalty (`λ>0`);
convergence requires the relative active-objective change `< tol` **and**
`max|Δr_m|, max|Δr_p| (, max|Δq|) < tol`, matching the OP engine's contract. The final
returned log-likelihood is evaluated at the final parameters.

## 7. Missing-data behavior

- **Missing offspring genotype** at a marker → emission `= 1` (neutral); the marker is
  informative through neighbors via the transition. (Same as the OP engine.)
- **Missing required sire (or mother) genotype** at a fitted marker of a full-sib
  cross → the phased allele label is unknown, so paternal (or maternal) transmission
  cannot be scored. The core does **not** invent an allele and does **not** substitute
  `q=0.5`. Documented safe behavior: **stop with an informative validation error**
  naming the cross and markers (`on_missing_parent = "error"`, the default). A future
  branch may add per-cross blockwise exclusion; that is deferred.

## 8. Compatibility with the current open-pollinated model

- A dataset with only `open_pollinated` crosses reproduces the current model exactly:
  the mixed API **dispatches to the existing `hmm_map`/`hmm_map_joint`** engine, so
  OP-only results are numerically unchanged (test 11).
- The reader remains backward compatible: existing pedigree/genotype CSVs (mother +
  offspring, no father) produce an `HSMap.data` with the same `G_list`/`M_list`, plus
  new cross-aware fields (`cross_table`, `crosses`, `parent_genotypes`, `F_list`) that
  default to the OP interpretation. A migration note documents old-object handling.
- A homozygous sire at all markers makes the paternal emission independent of `h_p`, so
  the 4-state model collapses to the maternal-only 2-state model with fixed paternal
  alleles (test 5); symmetrically for a homozygous mother (test 6).

## 9. Public API (Milestone 4)

- `hmm_map_fullsib(x, phased_m, phased_p, ...)` — known-sire crosses only.
- `hmm_map_mixed(x, phased_m, phased_p = NULL, ...)` — OP + known-sire; dispatches to
  the OP engine when no full-sib cross is present.

Output (both): maternal `r_m`, paternal `r_p`, maternal & paternal map distances,
observed log-likelihood, active objective, convergence status/reason, iteration count,
objective and parameter traces, contributing crosses/mothers/sires, per-cross family
type used, parent phase supplied, and (on request) posterior inheritance probabilities.
The paternal map is a genuine sire recombination map and is **not** labelled a
pollen-pool map.

## 10. Complexity

Per full-sib offspring: forward–backward over 4 states and `z` markers is `O(z·4²)`
time and `O(z·4)` memory (streamed per offspring, no `n×z×16` dense object). A cross of
`n` offspring is `O(n·z)`; the Kronecker transition is applied factored (two 2-state
steps) to avoid materializing 4×4 per interval where possible. Benchmark vs. the
2-state OP HMM is reported in Milestone 6.
