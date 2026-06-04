# The Tuple-Size Law at m=4: A Prime-Quadruplet Confirmation (Part XVIII)

Part XVII established that in the proper-depth coordinate `ell = ln ln X` a prime
`m`-tuple obeys `d ln rho_m / d ell -> -m`, and staked a falsifiable prediction:
the prime quadruplet must give **slope -4** and **ratio 2:1** against the twin.

**This is that test, and it passes.** Recomputing the quadruplet
`(6N-1, 6N+1, 6N+5, 6N+7)` from sieve primitives gives **152,141** centres on
S10 (reproduced on a second machine), a deep-shell slope `-4.175`, and a ratio to
the twin of **2.0145** — inside the bound `2 +/- 0.05`.

| pattern | m | deep slope (S9->S10) | intercept (count-aware) | ratio to twin |
|---|---:|---:|---:|---:|
| twin | 2 | -2.072 | -1.96 | 1.0000 |
| triplet | 3 | -3.111 | -3.01 | 1.5011 |
| **quadruplet** | **4** | **-4.175** | **-3.94** | **2.0145** |

The three points are **collinear through the origin**: `|d ln rho_m / d ell| =
1.041 m`, a single shared logarithmic tail `kappa = 1.041`, so the slope ratios
are the bare integers `1 : 1.5011 : 2.0145` against `1 : 3/2 : 2`, to 0.1%. Three
points on a line through the origin is qualitatively stronger than two: the
topological size `m` of the constellation, and nothing else, sets the decay rate.

> **The sparse-regime instrument.** Reaching m=4 forced the estimator to evolve.
> The quadruplet is ~160x rarer than the twin (only 26 members on S5, 128 on S6),
> so the unweighted deepest-5 intercept of Part XVII — fine in the dense m=2,3
> regime — is dominated by small-count noise and returns a spurious -5.12. A
> count-aware (Poisson-weighted) fit, `w = dell^2 / (1/c_a + 1/c_b)`, restores it,
> returning -1.96, -3.01, -3.94 for m=2,3,4 (each better than the unweighted
> value). This is **not** a correction to Part XVII: in the dense regime the two
> estimators agree. It is the new physics of rare constellations that demanded the
> sharper instrument. The robust, sparsity-independent statement is the **ratio**,
> which the absolute slopes inherit a tail and the ratios do not.
>
> **Scope (as Part XVII).** First moments only; the omega-variance / Erdos-Kac law
> on twin centres is open (Part XII). The `(ln X)^-m` law is Hardy-Littlewood; rho
> is the Part I scale. No infinitude is claimed.

**The next gauntlet (m=5).** The prime quintuplet `(6N-1, 6N+1, 6N+5, 6N+7,
6N+11)` must give slope **-5** and ratio **5:2** against the twin (`= 2.5`), 5:3
against the triplet, 5:4 against the quadruplet. ~10x rarer than the quadruplet,
it will lean entirely on the count-aware estimator and may need shells beyond S10.
A twin ratio outside `2.5 +/- 0.06` would falsify the tuple-size law.

Part I: doi:10.5281/zenodo.20470367 · XV: doi:10.5281/zenodo.20534965 ·
XVII: doi:10.5281/zenodo.20541575

---

## Layout

```
.
├── README.md
├── LICENSE                 (MIT)
├── CITATION.cff
├── .zenodo.json
├── data/
│   ├── tuple_flow.csv          shell, ell, and (count, rho, ln_rho) for twin/triplet/quad
│   └── quadruplet_counts.csv   shell, quad_count, centres, rho_quad, source
├── code/
│   ├── quadruplet_geodynamics.py   three-pattern verdict + figure; count-aware estimator
│   └── quadruplet_deep.py          from-scratch quadruplet recompute at S9, S10
├── figures/                fig_paper18_tuplelaw.{pdf,png}
└── paper/                  Chen_6N_Paper18.{tex,pdf} + figure
```

## Reproducing

Requirements: Python 3.8+, `numpy`, `matplotlib`.

```bash
pip install numpy matplotlib

# Three-pattern verdict + the two-panel figure (uses the embedded counts).
python code/quadruplet_geodynamics.py

# From-scratch quadruplet recompute at the deep shells (the bedrock check).
python code/quadruplet_deep.py                 # S9 and S10 (S10 a few minutes)
SHELLS=10 python code/quadruplet_deep.py        # S10 only
```

On Windows PowerShell set env vars separately, e.g.
`$env:SHELLS="10"; python code/quadruplet_deep.py`.

### Conventions (same as Parts I-XVII)

- Skeleton 6N±1; N the centre; rho_m = #{m-tuple centres} / #{centres N} per shell.
- ell = ln ln X (proper depth). Shells S_k : 6N in [10^(k-1), 10^k);
  centres per shell = 1.5e(k-1).
- quadruplet (6N-1, 6N+1, 6N+5, 6N+7) = two consecutive twin centres N, N+1.
- Engine: segmented value sieve; every count recomputed from primitives.

## License

MIT — see `LICENSE`.
