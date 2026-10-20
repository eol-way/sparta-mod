#!/usr/bin/env python3
"""
Long-time energy budget of a gas in a rotating square container
(in.macro_heating): first law for the gas, temperature plateau, drift,
work and heat at steady state

usage: python3 verify_macro_heating.py OUTDIR [PLOTDIR]

OUTDIR holds log.macro_heating.yes / .no and grid.macro_heating.yes / .no
.txt.  Checks are pure Python.  If matplotlib is available, the time
histories are plotted to PLOTDIR/macro_heating.png (default OUTDIR).
Exit code 1 if a check fails.

Also reported (INFO): the angular momentum not exchanged with the walls,
  tau = dLgas/dt + tz, Lgas = fnum sum m ((x-5) vy - (y-5) vx).
  DSMC pair collisions do not conserve Lz exactly: a pair at r1, r2 with
  momentum exchange dp changes Lz by ((r1 - r2) x dp)_z, a drag in rigid
  rotation, smaller with nearcp yes.  Particles pushed out of the walls
  (engulf project) change r without a wall tally and add to tau as well.

Quantities (gas -> walls, lab frame, every wall collision, compute surf):
  E = etot              energy flow from the gas into the walls
  T = tz                torque of the gas on the walls about the axis
  W = omega T           work done by the gas on the rotating walls
  Q = E - W             heat lost by the gas in the wall frame
  Egas = fnum sum 1/2 m v^2   (monatomic gas, no internal energy)
Checks for wallvel yes (wallvel no is reported as INFO, for comparison):
  run         the run reached its last step (no abort)
  particles   particle count constant, no particle deleted
  first law   Egas(t) - Egas(t - 20 dt) = -E dt * 20 for every 20-step
              window, to round-off (1e-6 of the total energy change)
  plateau     mean T over the second half = Twall within 4 standard
              errors (1000-step batch means) + 0.3 K
  drift       linear trend of T over the second half, change over the
              half < 4 standard errors + 0.3 K
  W, Q        means over the second half = 0 within 4 standard errors +
              0.1 % of the largest |E| during the transient
              (in any closed steady state Egas is constant, so E = W + Q
              averages to 0 by itself; W = 0 in addition needs zero net
              torque, i.e. no artificial torque (friction) on the gas)
  rotation    least-squares u_theta = Omega r (r < 1.75) of the cell
              velocities averaged over the second half, |Omega/omega - 1|
              < 3 %
"""

import math
import os
import sys

KB = 1.380649e-23
HERE = os.path.dirname(os.path.abspath(__file__))

TWALL = 300.0
DT = 1.0e-4
OMEGA = 2.0 * math.pi / (1000 * DT)     # rotate 360 in 1000 steps
SIDE = 3.9
BATCH = 50                              # rows (20 steps each) per batch

results = []


def report(status, name, msg):
    results.append(status)
    print(f"{status:4s}  {name:32s} {msg}")


def species_mass(name):
    for line in open(os.path.join(HERE, "air.species")):
        s = line.split()
        if s and not line.startswith("#") and s[0] == name:
            return float(s[2])
    raise ValueError(name)


def stats_rows(logfile):
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


def run_info(logfile):
    text = open(logfile).read()
    nrun = None
    for line in text.splitlines():
        s = line.split()
        if len(s) >= 4 and s[0] == "variable" and s[1] == "nrun":
            nrun = int(float(s[3]))
    return nrun, text.count("Loop time") >= 2


def mean_se(vals):
    """mean and standard error from batch means of BATCH rows (1000 steps)"""
    n = len(vals)
    mean = sum(vals) / n
    nb = n // BATCH
    if nb < 2:
        return mean, float("nan")
    bm = [sum(vals[i * BATCH:(i + 1) * BATCH]) / BATCH for i in range(nb)]
    bmean = sum(bm) / nb
    var = sum((b - bmean) ** 2 for b in bm) / (nb - 1)
    return mean, math.sqrt(var / nb)


def trend(steps, vals):
    """least-squares slope and its standard error (from batch means)"""
    nb = len(vals) // BATCH
    xs = [sum(steps[i * BATCH:(i + 1) * BATCH]) / BATCH for i in range(nb)]
    ys = [sum(vals[i * BATCH:(i + 1) * BATCH]) / BATCH for i in range(nb)]
    xm, ym = sum(xs) / nb, sum(ys) / nb
    sxx = sum((x - xm) ** 2 for x in xs)
    slope = sum((x - xm) * (y - ym) for x, y in zip(xs, ys)) / sxx
    res = [y - ym - slope * (x - xm) for x, y in zip(xs, ys)]
    se = math.sqrt(sum(r * r for r in res) / (nb - 2) / sxx)
    return slope, se


def box_rotation(gridfile, rmax=1.75):
    """least-squares slope of u_theta(r) in cells within rmax of center"""
    lines = open(gridfile).read().splitlines()
    i = [k for k, l in enumerate(lines) if l.startswith("ITEM: CELLS")][0]
    num = den = 0.0
    for l in lines[i + 1:]:
        if not l.strip():
            continue
        _, x, y, n, u, v = (float(t) for t in l.split())
        dx, dy = x - 5.0, y - 5.0
        r = math.hypot(dx, dy)
        if r > rmax or r == 0.0 or n <= 0.0:
            continue
        num += (-dy * u + dx * v) / r * r
        den += r * r
    return num / den


def analyse(outdir, wv, m, strict):
    tag = f"wallvel {wv}"
    log = os.path.join(outdir, f"log.macro_heating.{wv}")
    if not os.path.exists(log):
        report("FAIL" if strict else "INFO", tag, f"missing {log}")
        return None
    rows = stats_rows(log)
    nrun, complete = run_info(log)
    last = int(rows[-1]["Step"])
    ok = complete and last == nrun
    report(("PASS" if ok else "FAIL") if strict else "INFO", f"{tag} run",
           f"reached step {last} of {nrun}"
           + ("" if ok else " (aborted, checks below use the steps done)"))
    st = (lambda c: ("PASS" if c else "FAIL")) if strict else \
        (lambda c: "INFO")

    steps = [r["Step"] for r in rows]
    egas = [r["v_egas"] for r in rows]
    eflow = [r["c_wall[1]"] for r in rows]
    work = [OMEGA * r["c_wall[2]"] for r in rows]
    heat = [e - w for e, w in zip(eflow, work)]
    tgas = [r["v_tgas"] for r in rows]
    lgas = [r["v_lgas"] for r in rows]

    # particles

    npart = {r["Np"] for r in rows}
    ndel = sum(r["f_10"] for r in rows)
    report(st(len(npart) == 1 and ndel == 0), f"{tag} particles",
           f"{int(rows[0]['Np'])} -> {int(rows[-1]['Np'])}, "
           f"deleted {int(ndel)}")

    # transition regime: mean free path from the collision rate

    tail = rows[len(rows) // 2:]
    ncol = sum(r["Ncoll"] for r in tail) / len(tail)
    nu = 2.0 * ncol / (rows[0]["Np"] * DT)
    cbar = math.sqrt(8.0 * KB * TWALL / (math.pi * m))
    lam = cbar / nu
    print(f"      {tag}: mean free path {lam:.3f} m, Kn = lambda/side "
          f"= {lam / SIDE:.2f}")

    # first law, window by window: dEgas = -E dt

    dsum = 0.0
    maxres = 0.0
    total = egas[-1] - egas[0]
    cum_pred = [0.0]
    for k in range(1, len(rows)):
        nstep = steps[k] - steps[k - 1]
        pred = -eflow[k] * nstep * DT
        dsum += pred
        cum_pred.append(dsum)
        maxres = max(maxres, abs((egas[k] - egas[k - 1]) - pred))
    rel = maxres / abs(total)
    report(st(rel < 1.0e-6), f"{tag} first law per window",
           f"max |dEgas + E dt| = {maxres:.3e} J, {rel:.1e} of "
           f"Egas change {total:.4e} J")
    report(st(abs((egas[-1] - egas[0]) - dsum) < 1.0e-6 * abs(total)),
           f"{tag} first law cumulative",
           f"Egas(end)-Egas(0) = {total:.6e} J, -int E dt = {dsum:.6e} J")

    # plateau, drift, W and Q over the second half

    half = len(rows) // 2
    s2 = steps[half:]
    tau = (lgas[-1] - lgas[half]) / ((steps[-1] - steps[half]) * DT) + \
        sum(r["c_wall[2]"] for r in rows[half + 1:]) / (len(rows) - half - 1)
    print(f"      {tag}: angular momentum not exchanged with the walls "
          f"{tau:+.3e} N m (equivalent work omega tau {OMEGA * tau:+.3e} W)")
    tmean, tse = mean_se(tgas[half:])
    tol = 4.0 * tse + 0.3
    report(st(abs(tmean - TWALL) < tol), f"{tag} temperature plateau",
           f"{tmean:8.3f} K (Twall {TWALL:.0f} K, tol {tol:.2f} K), "
           f"steps {int(s2[0])}-{int(s2[-1])}")
    slope, sse = trend(s2, tgas[half:])
    span = s2[-1] - s2[0]
    drift, dse = slope * span, sse * span
    tol = 4.0 * dse + 0.3
    report(st(abs(drift) < tol), f"{tag} temperature drift",
           f"{drift:+.3f} K over {int(span)} steps (tol {tol:.2f} K)")
    epeak = max(abs(e) for e in eflow)
    for name, vals in (("W", work), ("Q", heat)):
        mean, se = mean_se(vals[half:])
        tol = 4.0 * se + 1.0e-3 * epeak
        report(st(abs(mean) < tol), f"{tag} steady {name}",
               f"{mean:+.4e} W (tol {tol:.2e} W, transient peak "
               f"|E| {epeak:.3e} W)")

    gridfile = os.path.join(outdir, f"grid.macro_heating.{wv}.txt")
    if os.path.exists(gridfile) and ok:
        ratio = box_rotation(gridfile) / OMEGA
        report(st(abs(ratio - 1.0) < 0.03), f"{tag} rigid rotation",
               f"gas Omega / wall omega = {ratio:.4f}")
    else:
        report("FAIL" if strict else "INFO", f"{tag} rigid rotation",
               "no grid dump (run did not finish)")

    return {"steps": steps, "tgas": tgas, "work": work, "heat": heat,
            "eflow": eflow, "degas": [e - egas[0] for e in egas],
            "pred": cum_pred}


def plot(plotdir, data):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("      matplotlib not available, no plot")
        return
    style = {"yes": {"color": "tab:blue"},
             "no": {"color": "tab:orange", "ls": "--"}}
    fig, ax = plt.subplots(2, 2, figsize=(12, 8))
    for wv, d in data.items():
        if d is None:
            continue
        ax[0, 0].plot(d["steps"][1:], d["tgas"][1:], label=f"wallvel {wv}",
                      lw=0.8, **style[wv])
    ax[0, 0].axhline(TWALL, color="k", lw=0.8, ls=":", label="Twall")
    ax[0, 0].set_xlabel("step")
    ax[0, 0].set_ylabel("gas temperature, thermal (K)")
    ax[0, 0].set_title("temperature: 100 K start, plateau, drift")
    ax[0, 0].legend()

    d = data.get("yes")
    if d:
        ax[0, 1].plot(d["steps"], d["work"], lw=0.8, label="W = omega T")
        ax[0, 1].plot(d["steps"], d["heat"], lw=0.8, label="Q = E - W")
        ax[0, 1].plot(d["steps"], d["eflow"], lw=0.8, label="E (etot)")
        ax[0, 1].set_xlim(0, 2000)
        ax[0, 1].set_xlabel("step")
        ax[0, 1].set_ylabel("gas -> walls (W)")
        ax[0, 1].set_title("wallvel yes: transient (first 2000 steps)")
        ax[0, 1].legend()

        ax[1, 1].plot(d["steps"], d["degas"], lw=1.5,
                      label="Egas(t) - Egas(0)")
        ax[1, 1].plot(d["steps"], d["pred"], lw=1.0, ls="--",
                      label="-integral of E dt")
        ax[1, 1].set_xlim(0, 2000)
        ax[1, 1].set_xlabel("step")
        ax[1, 1].set_ylabel("J")
        ax[1, 1].set_title("wallvel yes: first law for the gas")
        ax[1, 1].legend()

    for wv, d in data.items():
        if d is None:
            continue
        half = len(d["steps"]) // 2
        ax[1, 0].plot(d["steps"][half:], d["work"][half:], lw=0.6,
                      label=f"W, wallvel {wv}", **style[wv])
        ax[1, 0].plot(d["steps"][half:], d["heat"][half:], lw=0.6,
                      alpha=0.6, label=f"Q, wallvel {wv}",
                      color=style[wv]["color"], ls=":")
    ax[1, 0].axhline(0.0, color="k", lw=0.8)
    ax[1, 0].set_xlabel("step")
    ax[1, 0].set_ylabel("gas -> walls (W)")
    ax[1, 0].set_title("second half: W and Q around zero?")
    ax[1, 0].legend(fontsize=8)

    fig.tight_layout()
    fname = os.path.join(plotdir, "macro_heating.png")
    fig.savefig(fname, dpi=110)
    print(f"      wrote {fname}")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)
    outdir = sys.argv[1]
    plotdir = sys.argv[2] if len(sys.argv) > 2 else outdir
    m = species_mass("N")
    print(f"macro energy budget, rotating square container, omega "
          f"{OMEGA:.2f} rad/s, Twall {TWALL:.0f} K, gas starts at 100 K")
    data = {}
    for wv, strict in (("yes", True), ("no", False)):
        data[wv] = analyse(outdir, wv, m, strict)
    plot(plotdir, data)
    nfail = results.count("FAIL")
    npass = results.count("PASS")
    print()
    print(f"{npass} passed, {nfail} failed")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
