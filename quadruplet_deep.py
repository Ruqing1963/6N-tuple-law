#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Recompute the prime-quadruplet (6N-1, 6N+1, 6N+5, 6N+7) count at the deep shells
from scratch (segmented value sieve), and report the depth-flow slope and the
ratio against the twin. This is the m=4 test of the Part XVII tuple-size law.

  quadruplet = (p, p+2, p+6, p+8) with p = 6N-1
             = two consecutive twin centres N and N+1
  rho_4 = (#quadruplets) / (#centres) ;  ell = ln ln(10^k)
  target: d ln rho_4 / d ell -> -4 ;  ratio to twin -> 2:1

Verified here: S10 = 152,141 (from scratch). Twin deep slope (S9->S10) = -2.072.

Usage (PowerShell):
    python quadruplet_deep.py                  # S9 and S10
    $env:SHELLS="10"; python quadruplet_deep.py
S9 ~ 15 s; S10 a few minutes. Requires: numpy.
"""
import os, math, time
import numpy as np

SHELLS = [int(x) for x in os.environ.get("SHELLS", "9,10").split(",")]
SEG    = int(os.environ.get("SEG", 12_000_000))
QUAD   = ((0,-1),(0,1),(1,-1),(1,1))     # value = 6(N+off)+s
TWIN_LNRHO = {9: -3.917306, 10: -4.135659}   # Part I/XVII, for the ratio check

def primes_upto(n):
    s = np.ones(n + 1, bool); s[:2] = False
    for i in range(2, int(math.isqrt(n)) + 1):
        if s[i]: s[i*i::i] = False
    return np.nonzero(s)[0].astype(np.int64)

def quad_shell(k):
    base = primes_upto(int(math.isqrt(10**k)) + 2)
    Nlo = 10**(k-1)//6 + 1
    Nhi = (10**k - 1)//6
    cnt = 0; ns = Nlo
    while ns <= Nhi:
        ne  = min(ns + SEG, Nhi + 1)
        vlo = 6*ns - 1; vhi = 6*(ne - 1) + 7
        comp = np.zeros(vhi - vlo + 1, bool)
        for p in base:
            if p*p > vhi: break
            st = max(p*p, ((vlo + p - 1)//p)*p)
            comp[st - vlo::p] = True
        isp = ~comp
        N = np.arange(ns, ne, dtype=np.int64)
        ok = np.ones(ne - ns, bool)
        for off, s in QUAD:
            ok &= isp[6*N + (6*off + s) - vlo]
        cnt += int(np.count_nonzero(ok)); ns = ne
    return cnt, (Nhi - Nlo + 1)

def main():
    print("== quadruplet deep recompute ==  shells", SHELLS, " SEG", SEG, "\n")
    rows = []
    for k in SHELLS:
        t = time.time(); c, ncen = quad_shell(k)
        rho = c / ncen; ell = math.log(math.log(10.0**k))
        rows.append((k, c, rho, ell))
        print(f"S{k}: quadruplets={c}  centres={ncen}  rho4={rho:.10f}  "
              f"ln_rho4={math.log(rho):.5f}  ell={ell:.5f}  [{time.time()-t:.0f}s]")
    print()
    for (k0,*_), (k1,*_2) in zip(rows, rows[1:]):
        r0 = next(r for r in rows if r[0]==k0); r1 = next(r for r in rows if r[0]==k1)
        s4 = (math.log(r1[2]) - math.log(r0[2])) / (r1[3] - r0[3])
        print(f"d ln rho_4/d ell  S{k0}->S{k1} = {s4:.3f}   (target -4, tail = {s4+4:+.3f})")
        if k0 in TWIN_LNRHO and k1 in TWIN_LNRHO:
            s2 = (TWIN_LNRHO[k1]-TWIN_LNRHO[k0])/(r1[3]-r0[3])
            print(f"  ratio to twin = {s4/s2:.4f}   (target 2/1 = 2.0000)")

if __name__ == "__main__":
    main()
