#!/usr/bin/env python3
"""
Energy exchange of a cold beam with a diffusely reflecting moving piston
(in.micro_energy): momentum, energy, work and heat per collision compared
with the kernel theory

usage: python3 verify_micro_energy.py OUTDIR

OUTDIR holds log.micro_energy.yes and log.micro_energy.no.  Pure Python.
Exit code 1 if a check fails.

compute surf tallies every piston collision with the lab-frame velocities
before (v) and after (v') the bounce.  Summed over all collisions, gas ->
piston:
  Gamma = number of collisions per time               nflux_incident
  Fx    = sum m (vx - vx') per time                   fx
  E     = sum 1/2 m (|v|^2 - |v'|^2) per time         etot
The piston moves at U along its normal +x, so per collision
  W = U Fx      work done by the gas on the piston
  Q = E - W     energy lost by the gas in the piston frame (heat)
because 1/2 m |v - U|^2 = 1/2 m |v|^2 - U m v + 1/2 m U^2 (Galilean).

Theory, per collision (w = v - U, wall frame, normal +x into the gas):
  incident: beam with drift -V (V = 1000 m/s) and temperature Tb along x,
    s = -vx ~ N(V, sb^2), sb^2 = k Tb / m.  The emitted beam is flux
    weighted, so its density is ~ f(s) and the piston, closing in at U,
    is hit at a rate ~ f(s) (U + s):  <g> = <g (U+s)>_f / <U+s>_f
      <U+s>      = U + (U V + V^2 + sb^2) / (U + V)
      <(U+s)^2>  = U^2 + 2U<s> + (U (V^2 + sb^2) + V^3 + 3 V sb^2) / (U + V)
    tangential: 1/2 m <vy^2 + vz^2> = k Tb
  outgoing: half-range flux Maxwellian at Tw in the wall frame
    <w'_n> = sqrt(pi)/2 c,  1/2 m <|w'|^2> = 2 k Tw,  c = sqrt(2 k Tw / m)
  Fx/Gamma = -m (<U+s> + sqrt(pi)/2 c)
  Q/Gamma  = 1/2 m <(U+s)^2> + k Tb - 2 k Tw
  W/Gamma  = U Fx/Gamma
  E/Gamma  = Q/Gamma + W/Gamma
Pass: |measured - theory| < 4 standard errors (block averages of the 20-step
tallies) + 0.1% of the magnitude (beam tail and finite-sample margin).
wallvel no is reported for comparison only (INFO).
"""

import math
import os
import sys

KB = 1.380649e-23
HERE = os.path.dirname(os.path.abspath(__file__))

TWALL = 300.0
TBEAM = 30.0
VBEAM = 1000.0
U = 1.0 / (1000 * 1.0e-5)           # trans 1.0 in 1000 steps of 1e-5 s
STEP0 = 500                         # tallies of steps 501-1000 are used

results = []


def report(status, name, msg):
    results.append(status)
    print(f"{status:4s}  {name:34s} {msg}")


def species_mass(name):
    for line in open(os.path.join(HERE, "air.species")):
        s = line.split()
        if s and not line.startswith("#") and s[0] == name:
            return float(s[2])
    raise ValueError(name)


def stats_rows(logfile):
    """rows of all stats blocks in a SPARTA log as dicts keyed by column"""
    rows, hdr = [], None
    for line in open(logfile):
        s = line.split()
        if s and s[0] == "Step":
            hdr = s
            continue
        if hdr and s and s[0] == "Loop":
            hdr = None
            continue
        if hdr and len(s) == len(hdr):
            try:
                rows.append(dict(zip(hdr, map(float, s))))
            except ValueError:
                continue
    return rows


def theory(m):
    sb2 = KB * TBEAM / m
    c = math.sqrt(2.0 * KB * TWALL / m)
    es = (U * VBEAM + VBEAM**2 + sb2) / (U + VBEAM)
    es2 = (U * (VBEAM**2 + sb2) + VBEAM**3 + 3.0 * VBEAM * sb2) / (U + VBEAM)
    rel = U + es
    rel2 = U * U + 2.0 * U * es + es2
    fx = -m * (rel + 0.5 * math.sqrt(math.pi) * c)
    q = 0.5 * m * rel2 + KB * TBEAM - 2.0 * KB * TWALL
    w = U * fx
    return {"Fx": fx, "W": w, "Q": q, "E": q + w}


def measured(logfile):
    """per-collision Fx, E, W, Q: total and block values"""
    rows = [r for r in stats_rows(logfile) if r["Step"] > STEP0]
    if not rows:
        return None
    blocks = []
    for r in rows:
        g, fx, e = r["c_tot[1]"], r["c_tot[2]"], r["c_tot[3]"]
        blocks.append((g, fx, e))
    gsum = sum(b[0] for b in blocks)
    tot = {"Fx": sum(b[1] for b in blocks) / gsum,
           "E": sum(b[2] for b in blocks) / gsum}
    tot["W"] = U * tot["Fx"]
    tot["Q"] = tot["E"] - tot["W"]

    # standard error of the ratio estimate from the block scatter
    se = {}
    n = len(blocks)
    for key in ("Fx", "E", "W", "Q"):
        vals = []
        for g, fx, e in blocks:
            v = {"Fx": fx, "E": e, "W": U * fx, "Q": e - U * fx}[key]
            vals.append(v - tot[key] * g)
        gmean = gsum / n
        var = sum(v * v for v in vals) / (n - 1)
        se[key] = math.sqrt(var / n) / gmean
    deleted = max(r["f_10"] for r in rows)
    return tot, se, n, gsum / n, deleted


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    outdir = sys.argv[1]
    m = species_mass("N")
    th = theory(m)
    unit = m * 1.0e4        # print energies per collision in m * 1e4 m^2/s^2

    print("micro energy balance per piston collision "
          f"(U = {U:.0f} m/s, Twall = {TWALL:.0f} K, beam {VBEAM:.0f} m/s "
          f"{TBEAM:.0f} K)")
    print(f"      energies in units of m x 1e4 m^2/s^2 = {unit:.4g} J, "
          f"momentum in m x 1 m/s")

    for wv in ("yes", "no"):
        log = os.path.join(outdir, f"log.micro_energy.{wv}")
        if not os.path.exists(log):
            report("FAIL", f"micro energy wallvel {wv}", f"missing {log}")
            continue
        res = measured(log)
        if res is None:
            report("FAIL", f"micro energy wallvel {wv}", "no stats rows")
            continue
        tot, se, n, grate, deleted = res
        print(f"      wallvel {wv}: {n} blocks of 20 steps, "
              f"{grate:.4g} collisions/s, deleted {deleted:.0f}")
        strict = wv == "yes"
        for key in ("Fx", "W", "Q", "E"):
            scale = m if key == "Fx" else unit
            dev = tot[key] - th[key]
            tol = 4.0 * se[key] + 0.001 * abs(th[key])
            ok = abs(dev) < tol
            if strict:
                status = "PASS" if ok else "FAIL"
            else:
                status = "INFO"
            flag = "" if ok else "  <- differs"
            report(status, f"micro wallvel {wv} {key}/collision",
                   f"{tot[key]/scale:10.3f}  theory {th[key]/scale:10.3f}  "
                   f"(tol {tol/scale:.3f}){flag}")
        if strict:
            report("PASS" if deleted == 0 else "FAIL",
                   f"micro wallvel {wv} deleted particles", f"{deleted:.0f}")
        # balance sheet: lab energy = work + heat, with signs
        print(f"      balance: E = W + Q  ->  {tot['E']/unit:.3f} = "
              f"{tot['W']/unit:.3f} + {tot['Q']/unit:.3f}   "
              f"(wall does work {-tot['W']/unit:.3f}, "
              f"gas loses heat {tot['Q']/unit:.3f} in the wall frame)")

    nfail = results.count("FAIL")
    npass = results.count("PASS")
    print()
    print(f"{npass} passed, {nfail} failed")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
