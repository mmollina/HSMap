# Full-sib engine — performance note (m6-known-sire-core)

*Implementation-check benchmarks, not a publication result. Single-threaded, serial C++
(`src/hmm_fullsib.cpp`); correctness was prioritized over optimization in this first
implementation.*

## Complexity

For a cross of `n` offspring and `z` markers, per EM iteration:

| engine | states | time / iter | working memory |
|---|---|---|---|
| open-pollinated (2-state) | 2 | `O(n · z · 2²)` | `O(z · 2)` per offspring (streamed) + `O(z)` reductions |
| full-sib (4-state) | 4 | `O(n · z · 4²)` | `O(z · 4)` per offspring (streamed) + `O(z)` reductions |

The four-state forward–backward is applied per offspring and reduced immediately into
per-interval switch/total accumulators; no `n × z × 16` dense object is materialized.
Memory is therefore `O(z)` resident plus the per-offspring `O(4z)` scratch, independent
of `n`. Posterior inheritance probabilities are returned only on request (`n × 4z`).

## Benchmark vs. the two-state E-step (1 thread)

Per single E-step call (seconds), and the 4-state / 2-state ratio:

| z | n | fs_estep (4-state) | op_estep (2-state) | ratio |
|---:|---:|---:|---:|---:|
| 50  | 500  | 0.0075 | 0.0034 | 2.2 |
| 100 | 500  | 0.0152 | 0.0067 | 2.3 |
| 200 | 500  | 0.0300 | 0.0136 | 2.2 |
| 100 | 2000 | 0.0589 | 0.0266 | 2.2 |
| 500 | 1000 | 0.1525 | 0.0908 | 1.7 |

The four-state engine costs ~2.2× the two-state engine — well under the naive 16/4 = 4×,
because emission construction, scaling, and reductions are shared. Timing is **linear in
both `z` and `n`** (4× markers → 4× time at fixed `n`; 4× offspring → 4× time at fixed `z`),
confirming `O(n·z)`.

## Full-fit runtime (`hmm_map_fullsib`, 1 thread)

| z | n | wall time | iterations | ms/iter |
|---:|---:|---:|---:|---:|
| 100 | 500  | 0.68 s | 46 | 14.8 |
| 200 | 500  | 1.50 s | 50 | 30.0 |
| 100 | 2000 | 2.13 s | 36 | 59.2 |

A full-sib map over ~100–200 markers and hundreds–thousands of offspring fits in ~1–2 s.
Cost scales linearly with markers, offspring, crosses, and iterations.

## Future optimization (deferred)

Not done here (correctness first): (i) parallelize the per-offspring FB across cores with
RcppParallel, as the two-state OP engine already does — the offspring loop is embarrassingly
parallel; (ii) apply the Kronecker transition factored as two 2-state steps (8 vs 16 mults
per interval); (iii) cache emission columns shared by offspring with identical genotype rows.
These are performance-only and must preserve the brute-force-verified numerics.
