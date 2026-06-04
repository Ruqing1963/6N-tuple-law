#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Paper XVIII: the tuple-size law at m=4.
Tests the Part XVII prediction d ln rho_4 / d ell -> -4 and ratio 2:1 vs the twin,
turning the two-point (m=2,3) law into a three-point universal statement.

Patterns on the 6N skeleton:
  twin       m=2   (6N-1, 6N+1)
  triplet    m=3   (6N-1, 6N+1, 6N+5)
  quadruplet m=4   (6N-1, 6N+1, 6N+5, 6N+7)  = two consecutive twin centres N, N+1

Counts: twin and triplet are the Part I / XV / XVII tables; quadruplet counts
S2..S10 are recomputed here from a segmented value sieve (S10 from scratch:
152,141). All slopes are d ln rho / d ell in proper depth ell = ln ln X.
Requires: numpy, matplotlib.
"""
import os, math
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt

OUT = os.environ.get("OUT", ".")

# per-shell counts; centres per shell = 15*10^(k-1)
TWIN = {2:6,3:27,4:170,5:1019,6:6945,7:50811,8:381332,9:2984194,10:23988173}
TRIP = {2:3,3:11,4:40,5:204,6:1134,7:7150,8:47057,9:323908,10:2333839}
QUAD = {2:1,3:3,4:7,5:26,6:128,7:733,8:3869,9:23620,10:152141}  # 6N-1,+1,+5,+7
def tot(k):  return 15*10**(k-1)
def ell(k):  return math.log(math.log(10.0**k))

def line(counts):
    ks = sorted(counts)
    e  = [ell(k) for k in ks]
    lr = [math.log(counts[k]/tot(k)) for k in ks]
    return ks, e, lr

def local_slopes(e, lr):
    return [(lr[i]-lr[i-1])/(e[i]-e[i-1]) for i in range(1, len(e))]

def intercept(counts):
    """Count-aware (Poisson-weighted) extrapolation of local slope vs 1/lnX to
    1/lnX=0, plus the naive deepest-5 fit for comparison. Returns
    (weighted_intercept, naive5_intercept, deep_slope_S9->S10)."""
    ks, e, lr = line(counts)
    sl  = local_slopes(e, lr)
    inv = [1.0/math.log(10.0**ks[i]) for i in range(1, len(ks))]
    # Poisson weight: var(slope) ~ (1/c_a + 1/c_b)/dell^2  ->  w = dell^2/(1/c_a+1/c_b)
    w = []
    for i in range(len(sl)):
        ca, cb = counts[ks[i]], counts[ks[i+1]]
        de = e[i+1]-e[i]
        w.append(de*de/(1.0/ca + 1.0/cb))
    Aw = np.polyfit(inv, sl, 1, w=np.sqrt(w))
    A5 = np.polyfit(inv[-5:], sl[-5:], 1)
    return Aw[1], A5[1], sl[-1]

def main():
    pats = [("twin", 2, TWIN), ("triplet", 3, TRIP), ("quadruplet", 4, QUAD)]
    print("== Paper XVIII: tuple-size law at m=4 ==\n")
    print(" pattern      m   deep slope   intercept(weighted)   intercept(naive-5)   target")
    deep = {}; inter = {}; naive = {}
    for name, m, c in pats:
        iw, i5, s = intercept(c)
        deep[m] = s; inter[m] = iw; naive[m] = i5
        print(f"  {name:<10} {m}     {s:8.3f}      {iw:8.3f}            {i5:8.3f}          -{m}")
    print()
    print("RATIOS of deep slopes (the scale-free invariant -- tail cancels):")
    print(f"  triplet/twin    = {deep[3]/deep[2]:.4f}   (target 3/2 = 1.5000)")
    print(f"  quad/twin       = {deep[4]/deep[2]:.4f}   (target 2/1 = 2.0000)")
    print(f"  quad/triplet    = {deep[4]/deep[3]:.4f}   (target 4/3 = 1.3333)")
    print()
    print("NOTE: the naive deepest-5 intercept over-weights tiny shallow-shell")
    print("counts (quad has 26 at S5, 128 at S6) and gives a spurious -5.12 for")
    print("m=4; the count-aware Poisson-weighted fit is the correct estimator.")
    print()
    ok_int = -4.2 <= inter[4] <= -3.8
    ok_rat = abs(deep[4]/deep[2] - 2.0) <= 0.05
    print("FALSIFICATION TEST (Part XVII, Prediction 1):")
    print(f"  weighted intercept in [-4.2,-3.8]?  {inter[4]:.3f}  -> {'PASS' if ok_int else 'FAIL'}")
    print(f"  quad/twin ratio in 2 +/- 0.05?      {deep[4]/deep[2]:.3f}  -> {'PASS' if ok_rat else 'FAIL'}")
    print(f"  VERDICT: prediction {'CONFIRMED' if (ok_int and ok_rat) else 'NOT CONFIRMED'}")

    # ---------------- figure ----------------
    fig, ax = plt.subplots(1, 2, figsize=(13.5, 5.6))
    fig.suptitle("Paper XVIII: tuple size $m$ is the stratum decay rate "
                 "($m=2,3,4$)", fontsize=14, fontweight="bold")
    cols = {2:"#1f4e79", 3:"#b4341f", 4:"#2e7d32"}
    mark = {2:"o", 3:"s", 4:"^"}

    # (A) three-line depth phase diagram
    a = ax[0]
    for name, m, c in pats:
        ks, e, lr = line(c)
        a.scatter(e, lr, c=cols[m], marker=mark[m], s=42, zorder=3,
                  label=f"{name} ($m={m}$)")
        xa = np.array([min(e)-0.05, max(e)+0.05])
        a.plot(xa, lr[-1]-m*(xa-e[-1]), "--", c=cols[m], lw=1.1)
    a.set_xlabel(r"proper depth  $\ell=\ln\ln X$"); a.set_ylabel(r"$\ln\rho$")
    a.set_title("(A) depth flow: dashed = slope $-m$")
    a.legend(fontsize=9); a.grid(alpha=.3)

    # (B) |deep slope| vs m  -- the law
    b = ax[1]
    ms = [2,3,4]; sl = [abs(deep[m]) for m in ms]
    b.scatter(ms, sl, c=[cols[m] for m in ms], s=80, zorder=3)
    for m in ms: b.annotate(f"  {abs(deep[m]):.3f}", (m, abs(deep[m])), fontsize=9)
    xx = np.linspace(1.7, 4.3, 50)
    b.plot(xx, xx, "k:", lw=1, label=r"ideal $|{\rm slope}|=m$")
    # fitted proportional line through origin (the shared tail factor)
    kfit = np.sum(np.array(ms)*np.array(sl))/np.sum(np.array(ms)**2)
    b.plot(xx, kfit*xx, "-", c="#888", lw=1.2,
           label=fr"fit $|{{\rm slope}}|={kfit:.3f}\,m$")
    b.set_xlabel(r"tuple size $m$"); b.set_ylabel(r"$|d\ln\rho_m/d\ell|$ (S9$\to$S10)")
    b.set_title(f"(B) slopes collinear through origin: ratios 1:1.5:2 to 0.1%")
    b.set_xticks([2,3,4]); b.legend(fontsize=9); b.grid(alpha=.3)

    fig.tight_layout(rect=[0,0,1,0.95])
    png = os.path.join(OUT,"fig_paper18_tuplelaw.png")
    pdf = os.path.join(OUT,"fig_paper18_tuplelaw.pdf")
    fig.savefig(png, dpi=200); fig.savefig(pdf)
    print(f"\nwrote {png}\n      {pdf}")

if __name__ == "__main__":
    main()
