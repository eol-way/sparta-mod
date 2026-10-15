#!/usr/bin/env python3
"""
Velocity distribution of particles diffusely reflected by a moving wall
(in.test_translation, in.test_rotation), compared with the kernel theory

usage: python3 verify_distribution.py OUTDIR [PLOTDIR]

OUTDIR holds dump.translation.MODEL_WALLVEL.*.txt and
dump.rotation.MODEL_WALLVEL.*.txt (MODEL = diffuse | vanish).
Statistics are pure Python.  If matplotlib is available, histograms with
the theoretical curves are written to PLOTDIR (default OUTDIR).
Exit code 1 if a check fails.

Theory: a fully diffuse wall moving with velocity u_w re-emits molecules
with a Maxwellian at Twall in the wall frame, flux weighted:
  normal component w_n = (v - u_w).n > 0 : (2 w_n / c^2) exp(-w_n^2 / c^2)
  tangential components                  : Gaussian, mean (u_w).t, sd c/sqrt(2)
  with c = sqrt(2 k Twall / m)
Particles sampled in a thin layer next to the wall are density weighted
(a slow particle stays longer in the layer).  Each sample is weighted by
1 / (time it spends in the layer) to recover the flux-weighted (emitted)
distribution.  For the rotating cylinder, every particle is traced back
along its straight path to the emission point on r = R, so the normal and
tangential components are taken in the frame of that point (the tangential
velocity of a free particle drops as R/r further out, angular momentum).
"""

import glob
import math
import os
import sys

KB = 1.380649e-23
HERE = os.path.dirname(os.path.abspath(__file__))
TWALL = 300.0

# in.test_translation
X0, DX, NLARGE_T, DT = 0.2, 1.0, 1000, 1.0e-5      # piston x = X0 + DX*step/N
U_PISTON = DX / (NLARGE_T * DT)                     # 100 m/s
LAYER_T = 0.05

# in.test_rotation
R_IN = 0.5
OMEGA = 2.0 * math.pi / (1600 * DT)                 # 392.7 rad/s
LAYER_R = 0.05

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


def dumps(outdir, prefix):
    out = []
    for f in glob.glob(os.path.join(outdir, prefix + ".*.txt")):
        out.append((int(os.path.basename(f).split(".")[-2]), f))
    return sorted(out)


def read_dump(fname):
    lines = open(fname).read().splitlines()
    i = [k for k, l in enumerate(lines) if l.startswith("ITEM: ATOMS")][0]
    names = lines[i].split()[2:]
    return [dict(zip(names, map(float, l.split()))) for l in lines[i + 1:]
            if l.strip()]


# ----------------------------------------------------------------------
# samples: (weight, normal, tangential components...) per reflected particle
# ----------------------------------------------------------------------

def samples_translation(outdir, tag):
    """w_n = vx - U in the piston frame, tangential vy, vz
    residence time in the layer moving with the piston = LAYER / w_n"""
    out = []
    for step, f in dumps(outdir, f"dump.translation.{tag}"):
        xp = X0 + DX * step / NLARGE_T
        for p in read_dump(f):
            if not (xp < p["x"] < xp + LAYER_T):
                continue
            wn = p["vx"] - U_PISTON
            if wn <= 0.0:
                continue
            out.append((wn / LAYER_T, wn, p["vy"], p["vz"]))
    return out


def samples_rotation(outdir, tag):
    """trace each outgoing particle back to r = R_IN, components there
    residence time = path length from r = R_IN to R_IN + LAYER / speed"""
    out = []
    nlost = 0
    rout = R_IN + LAYER_R
    for _, f in dumps(outdir, f"dump.rotation.{tag}"):
        for p in read_dump(f):
            x, y, vx, vy = p["x"], p["y"], p["vx"], p["vy"]
            r = math.hypot(x, y)
            if not (R_IN < r < rout) or x * vx + y * vy <= 0.0:
                continue
            v2 = vx * vx + vy * vy
            pv = x * vx + y * vy
            disc = pv * pv - v2 * (r * r - R_IN * R_IN)
            if disc < 0.0:          # straight path never touched r = R_IN
                nlost += 1
                continue
            t = (pv - math.sqrt(disc)) / v2          # back to the wall
            ex, ey = x - t * vx, y - t * vy
            nx, ny = ex / R_IN, ey / R_IN
            vn = vx * nx + vy * ny
            vt = -vx * ny + vy * nx
            # path from r = R_IN to r = rout along v (forward from e)
            ev = ex * vx + ey * vy
            s = (-ev + math.sqrt(ev * ev + v2 * (rout * rout - R_IN * R_IN))) / v2
            out.append((1.0 / s, vn, vt, p["vz"]))
    return out, nlost


# ----------------------------------------------------------------------
# weighted moments
# ----------------------------------------------------------------------

def moments(samples, k):
    w = [s[0] for s in samples]
    x = [s[k] for s in samples]
    sw = sum(w)
    mean = sum(a * b for a, b in zip(w, x)) / sw
    var = sum(a * (b - mean) ** 2 for a, b in zip(w, x)) / sw
    neff = sw * sw / sum(a * a for a in w)
    return mean, math.sqrt(var), math.sqrt(var / neff), neff


def check_case(name, samples, c, ut, strict):
    """normal: Rayleigh mean sqrt(pi) c / 2, rms c
    tangential: mean ut, sd c / sqrt(2); 2nd tangential: mean 0"""
    if len(samples) < 1000:
        report("FAIL", name, f"only {len(samples)} reflected samples")
        return
    tol = 0.01 * c
    checks = []
    mn, sn, en, neff = moments(samples, 1)
    rms = math.sqrt(sum(s[0] * s[1] ** 2 for s in samples) /
                    sum(s[0] for s in samples))
    checks.append(("normal mean", mn, math.sqrt(math.pi) * c / 2.0, en))
    checks.append(("normal rms", rms, c, en))
    mt, st, et, _ = moments(samples, 2)
    checks.append(("tangential mean", mt, ut, et))
    checks.append(("tangential sd", st, c / math.sqrt(2.0), et))
    mz, sz, ez, _ = moments(samples, 3)
    checks.append(("2nd tangential mean", mz, 0.0, ez))
    checks.append(("2nd tangential sd", sz, c / math.sqrt(2.0), ez))
    print(f"      {name}: {len(samples)} samples, effective {neff:.0f}, "
          f"c = {c:.1f} m/s")
    for label, got, want, err in checks:
        ok = abs(got - want) < 4.0 * err + tol
        if strict:
            report("PASS" if ok else "FAIL", f"{name} {label}",
                   f"{got:8.2f} m/s  theory {want:8.2f}  (4 sigma + 1% c = "
                   f"{4.0 * err + tol:.2f})")
        else:
            report("INFO", f"{name} {label}",
                   f"{got:8.2f} m/s  theory {want:8.2f}"
                   f"{'' if ok else '  <- differs'}")


# ----------------------------------------------------------------------
# plots (optional)
# ----------------------------------------------------------------------

def plot(plotdir, fname, title, cases, c, ut, tlabel):
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("      matplotlib not available, no plot")
        return
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2))
    nb = 60
    for ax, k, lo, hi in ((axes[0], 1, 0.0, 3.0 * c),
                          (axes[1], 2, ut - 3.0 * c, ut + 3.0 * c)):
        width = (hi - lo) / nb
        for label, samples, style in cases:
            if not samples:
                continue
            h = [0.0] * nb
            for s in samples:
                j = int((s[k] - lo) / width)
                if 0 <= j < nb:
                    h[j] += s[0]
            norm = sum(s[0] for s in samples) * width
            xs = [lo + (j + 0.5) * width for j in range(nb)]
            ax.step(xs, [v / norm for v in h], where="mid", label=label, **style)
        xs = [lo + (hi - lo) * j / 400.0 for j in range(401)]
        if k == 1:
            ys = [2.0 * x / c ** 2 * math.exp(-(x / c) ** 2) for x in xs]
            ax.set_xlabel("normal velocity in wall frame w_n (m/s)")
        else:
            ys = [math.exp(-((x - ut) / c) ** 2) / (math.sqrt(math.pi) * c)
                  for x in xs]
            ax.axvline(ut, color="k", ls=":", lw=0.8)
            ax.set_xlabel(tlabel)
        ax.plot(xs, ys, "k--", lw=1.2, label="theory (diffuse, Twall)")
        ax.set_ylabel("probability density (flux weighted)")
        ax.legend(fontsize=8)
    fig.suptitle(title, fontsize=11)
    fig.tight_layout()
    os.makedirs(plotdir, exist_ok=True)
    path = os.path.join(plotdir, fname)
    fig.savefig(path, dpi=90)
    plt.close(fig)
    print(f"      wrote {path}")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    outdir = sys.argv[1]
    plotdir = sys.argv[2] if len(sys.argv) > 2 else outdir
    c = math.sqrt(2.0 * KB * TWALL / species_mass("N"))

    # translation: wall velocity (U, 0, 0), normal +x, tangential mean 0

    yes = samples_translation(outdir, "diffuse_yes")
    check_case("translation wallvel yes", yes, c, 0.0, True)
    no = samples_translation(outdir, "diffuse_no")
    check_case("translation wallvel no", no, c, 0.0, False)
    vanish = samples_translation(outdir, "vanish_yes")
    status = "PASS" if len(vanish) <= 0.001 * max(len(yes), 1) else "FAIL"
    report(status, "translation control (vanish)",
           f"{len(vanish)} outgoing particles vs {len(yes)} with diffuse "
           f"(incoming particles mistaken as reflected)")
    plot(plotdir, "velocity_translation.png",
         f"piston moving at {U_PISTON:.0f} m/s, diffuse wall at {TWALL:.0f} K",
         [("wallvel yes", yes, {"color": "tab:blue"}),
          ("wallvel no", no, {"color": "tab:orange", "ls": "--"})],
         c, 0.0, "tangential velocity vy (m/s)")

    # rotation: wall velocity R omega along theta, normal radial

    ut = R_IN * OMEGA
    yes, lost_yes = samples_rotation(outdir, "diffuse_yes")
    check_case("rotation wallvel yes", yes, c, ut, True)
    no, _ = samples_rotation(outdir, "diffuse_no")
    check_case("rotation wallvel no", no, c, ut, False)
    vanish, _ = samples_rotation(outdir, "vanish_yes")
    status = "PASS" if len(vanish) <= 0.001 * max(len(yes), 1) else "FAIL"
    report(status, "rotation control (vanish)",
           f"{len(vanish)} outgoing particles vs {len(yes)} with diffuse, "
           f"{lost_yes} paths not from r = R")
    plot(plotdir, "velocity_rotation.png",
         f"cylinder R = {R_IN} rotating at {OMEGA:.1f} rad/s "
         f"(wall speed {ut:.1f} m/s), diffuse wall at {TWALL:.0f} K",
         [("wallvel yes", yes, {"color": "tab:blue"}),
          ("wallvel no", no, {"color": "tab:orange", "ls": "--"})],
         c, ut, "tangential velocity v_theta at emission point (m/s)")

    nfail = results.count("FAIL")
    print(f"\n{results.count('PASS')} passed, {nfail} failed")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
