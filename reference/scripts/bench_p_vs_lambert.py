#!/usr/bin/env python3
"""Micro-benchmark: Izzo Lambert solve vs closed-form p-parametrised velocities.

Usage (from the repo root):   python3 scripts/bench_p_vs_lambert.py [N] [REPEATS]
Reports the best (lowest) mean time per call over REPEATS repeats, which is the standard way to
reduce noise from other processes.  What is measured: ONE function call per evaluation, same
geometry, Python + numpy overhead included.  It does NOT measure a full minimum-delta-v search.
"""
import os, platform, sys, time
import numpy as np
sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from planechange.lambert_izzo import lambert
from planechange.pparam import velocities_from_p, time_of_flight

N = int(sys.argv[1]) if len(sys.argv) > 1 else 2000
REPEATS = int(sys.argv[2]) if len(sys.argv) > 2 else 5
MU = 3.986004418e5
r1 = np.array([7000.0, 0, 0])
r2 = np.array([-3000.0, 12000, 4000])
r1n, r2n = np.linalg.norm(r1), np.linalg.norm(r2)


def best_mean(fn):
    fn(0)  # warm-up
    runs = []
    for _ in range(REPEATS):
        t = time.perf_counter()
        for i in range(N):
            fn(i)
        runs.append((time.perf_counter() - t) / N)
    return min(runs)


def f_lambert(i):
    lambert(r1, r2, 3000.0 + i, MU)          # tof varies a little so nothing is cached


def f_p(i):
    velocities_from_p(r1, r2, 9000.0 + i, MU)


def f_p_plus_tof(i):
    v1, _, th = velocities_from_p(r1, r2, 9000.0 + i, MU)
    time_of_flight(r1n, r2n, th, 9000.0 + i, v1, r1, MU)


tl, tp, tpt = best_mean(f_lambert), best_mean(f_p), best_mean(f_p_plus_tof)
print(f"python {platform.python_version()}  numpy {np.__version__}  {platform.machine()}  N={N} repeats={REPEATS}")
print(f"lambert()                         {tl*1e6:8.1f} us/call")
print(f"velocities_from_p()               {tp*1e6:8.1f} us/call   ratio {tl/tp:5.1f}x")
print(f"velocities_from_p()+time_of_flight {tpt*1e6:7.1f} us/call   ratio {tl/tpt:5.1f}x")
