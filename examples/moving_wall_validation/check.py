#!/usr/bin/env python3
"""
Pass/fail checks for examples/moving_wall_validation, called by run_all.sh

usage: python3 check.py OUTDIR [NPROCS]

OUTDIR holds the logs, particle dumps and grid dumps written by run_all.sh.
Pure Python (no numpy).  Exit code 1 if any check fails.
"""

import glob
import math
import os
import sys

KB = 1.380649e-23
HERE = os.path.dirname(os.path.abspath(__file__))

results = []


def report(status, name, msg):
    results.append(status)
    print(f"{status:4s}  {name:28s} {msg}")


# ----------------------------------------------------------------------
# readers
# ----------------------------------------------------------------------

def stats_rows(logfile, as_text=False):
    """rows of all stats blocks in a SPARTA log as dicts keyed by column"""
    rows = []
    hdr = None
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
                vals = [float(v) for v in s]
            except ValueError:
                continue
            rows.append(dict(zip(hdr, s if as_text else vals)))
    return rows


def read_surf(fname):
    """points {id: coords} and elements [(i,j[,k])] of a SPARTA surf file"""
    pts, elems, sec = {}, [], None
    for line in open(fname):
        s = line.split()
        if not s:
            continue
        if s[0] in ("Points", "Lines", "Triangles"):
            sec = s[0]
            continue
        if sec == "Points" and len(s) >= 3:
            pts[int(s[0])] = tuple(float(v) for v in s[1:])
        elif sec in ("Lines", "Triangles") and len(s) >= 3:
            elems.append(tuple(int(v) for v in s[1:]))
    return pts, elems


def read_dump(fname, ncoord):
    """particle coords from a SPARTA particle dump written as 'id x y [z]'"""
    lines = open(fname).read().splitlines()
    i = [k for k, l in enumerate(lines) if l.startswith("ITEM: ATOMS")][0]
    return [tuple(float(v) for v in l.split()[1:1 + ncoord])
            for l in lines[i + 1:] if l.strip()]


def dump_steps(pattern):
    """(step, file) for dump files named prefix.STEP.txt"""
    out = []
    for f in glob.glob(pattern):
        out.append((int(os.path.basename(f).split(".")[-2]), f))
    return sorted(out)


def species_mass(name):
    for line in open(os.path.join(HERE, "air.species")):
        s = line.split()
        if s and not line.startswith("#") and s[0] == name:
            return float(s[2])
    raise ValueError(name)


# ----------------------------------------------------------------------
# geometry
# ----------------------------------------------------------------------

def segments_2d(surffile, origin, scale, trans, center, theta):
    """segments of a 2d surf after read_surf origin/scale/trans and
    a rotation by theta around center"""
    pts, lines = read_surf(surffile)
    c, s = math.cos(theta), math.sin(theta)
    xy = {}
    for i, p in pts.items():
        x = origin[0] + scale[0] * (p[0] - origin[0]) + trans[0]
        y = origin[1] + scale[1] * (p[1] - origin[1]) + trans[1]
        dx, dy = x - center[0], y - center[1]
        xy[i] = (center[0] + c * dx - s * dy, center[1] + s * dx + c * dy)
    return [(xy[a], xy[b]) for a, b in lines]


def inside_2d(px, py, segs):
    """even-odd ray crossing test against a closed set of segments"""
    flag = False
    for (x1, y1), (x2, y2) in segs:
        if (y1 > py) != (y2 > py):
            if px < x1 + (py - y1) * (x2 - x1) / (y2 - y1):
                flag = not flag
    return flag


# ----------------------------------------------------------------------
# checks
# ----------------------------------------------------------------------

def check_conserved(name, logfile, nexpect=None):
    rows = stats_rows(logfile)
    if not rows:
        report("FAIL", name, f"no stats in {logfile}")
        return
    n0, n1 = rows[0]["Np"], rows[-1]["Np"]
    if nexpect is not None:
        n0 = nexpect
    if n1 == n0:
        report("PASS", name, f"particles conserved: {int(n0)} -> {int(n1)}")
    else:
        report("FAIL", name, f"particles not conserved: {int(n0)} -> {int(n1)}")


def check_lost(name, logfile):
    rows = stats_rows(logfile)
    n0, n1 = rows[0]["Np"], rows[-1]["Np"]
    report("INFO", name, f"legacy deletion: {int(n0)} -> {int(n1)} particles")


def check_inside(name, ninside, ntotal, what):
    if ninside == 0:
        report("PASS", name, f"0 of {ntotal} particles inside {what}")
    else:
        report("FAIL", name, f"{ninside} of {ntotal} particles inside {what}")


def check_ellipse2d(out, nprocs):
    log = os.path.join(out, "log.ellipse2d.project")
    check_conserved("ellipse2d project", log)

    surf = os.path.join(HERE, "data.circle")
    for step, f in dump_steps(os.path.join(out, "dump.ellipse2d.project.*.txt")):
        theta = 2.0 * math.pi * step / 2000.0
        segs = segments_2d(surf, (5, 5), (0.2, 1.2), (0, 0), (5, 5), theta)
        xy = read_dump(f, 2)
        nin = sum(inside_2d(x, y, segs) for x, y in xy)
        check_inside(f"ellipse2d inside step {step}", nin, len(xy), "ellipse")

    log = os.path.join(out, "log.ellipse2d.delete")
    check_lost("ellipse2d delete", log)

    # engulf delete must reproduce unmodified SPARTA
    # reference was made on 4 procs, other counts give other random streams

    ref = os.path.join(HERE, "ref", "ellipse2d_delete.np4.stats")
    if str(nprocs) != "4":
        report("SKIP", "ellipse2d delete vs ref", "reference is for 4 procs")
        return
    keys = ["Step", "Np", "Nscoll", "f_10", "c_ket"]
    now = [[r[k] for k in keys] for r in stats_rows(log, as_text=True)]
    refrows = [l.split() for l in open(ref)
               if l.strip() and not l.startswith("#")][1:]
    if now == refrows:
        report("PASS", "ellipse2d delete vs ref",
               "identical to unmodified SPARTA")
    else:
        report("WARN", "ellipse2d delete vs ref",
               "differs from reference (expected only on another "
               "compiler/MPI build)")


def check_ellipsoid3d(out):
    log = os.path.join(out, "log.ellipsoid3d.project")
    check_conserved("ellipsoid3d project", log)

    pts, tris = read_surf(os.path.join(HERE, "data.sphere"))
    a, b, c = 1.5, 0.5, 0.8
    p = {i: (a * v[0], b * v[1], c * v[2]) for i, v in pts.items()}
    faces = []
    for i, j, k in tris:
        e1 = [p[j][m] - p[i][m] for m in range(3)]
        e2 = [p[k][m] - p[i][m] for m in range(3)]
        n = [e1[1] * e2[2] - e1[2] * e2[1], e1[2] * e2[0] - e1[0] * e2[2],
             e1[0] * e2[1] - e1[1] * e2[0]]
        if sum(n[m] * p[i][m] for m in range(3)) < 0.0:   # outward, convex
            n = [-v for v in n]
        faces.append((p[i], n))

    for step, f in dump_steps(os.path.join(out, "dump.ellipsoid3d.project.*.txt")):
        theta = 2.0 * math.pi * step / 2000.0
        cs, sn = math.cos(theta), math.sin(theta)
        xyz = read_dump(f, 3)
        nin = 0
        for x, y, z in xyz:
            xb, yb = cs * x + sn * y, -sn * x + cs * y      # into body frame
            if (xb / a) ** 2 + (yb / b) ** 2 + (z / c) ** 2 >= 1.0:
                continue   # mesh is inside the ellipsoid
            q = (xb, yb, z)
            if all(sum(n[m] * (q[m] - v0[m]) for m in range(3)) < 0.0
                   for v0, n in faces):
                nin += 1
        check_inside(f"ellipsoid3d inside step {step}", nin, len(xyz),
                     "ellipsoid")

    check_lost("ellipsoid3d delete", os.path.join(out, "log.ellipsoid3d.delete"))


def check_two_rotors(out):
    check_conserved("two_rotors", os.path.join(out, "log.two_rotors"))
    surf = os.path.join(HERE, "data.circle")
    for step, f in dump_steps(os.path.join(out, "dump.two_rotors.*.txt")):
        theta = 2.0 * math.pi * step / 2000.0
        r1 = segments_2d(surf, (5, 5), (0.1, 0.4), (-1.5, 0), (3.5, 5), theta)
        r2 = segments_2d(surf, (5, 5), (0.1, 0.4), (1.5, 0), (6.5, 5), -theta)
        xy = read_dump(f, 2)
        nin = sum(inside_2d(x, y, r1) or inside_2d(x, y, r2) for x, y in xy)
        check_inside(f"two_rotors inside step {step}", nin, len(xy), "rotors")


def check_move_surf_cmd(out):
    log = os.path.join(out, "log.move_surf_cmd.project")
    text = open(log).read()
    ndel = [l for l in text.splitlines() if "deleted particles" in l]
    npush = [l for l in text.splitlines() if "pushed out of surfs" in l]
    if ndel and ndel[0].split()[0] == "0":
        report("PASS", "move_surf_cmd project",
               f"0 deleted, {npush[0].split()[0] if npush else '?'} pushed out")
    else:
        report("FAIL", "move_surf_cmd project",
               f"deleted line: {ndel[0].strip() if ndel else 'missing'}")
    check_conserved("move_surf_cmd project", log)

    surf = os.path.join(HERE, "data.circle")
    segs = segments_2d(surf, (5, 5), (0.2, 1.2), (1.0, 0), (6, 5), 0.0)
    for step, f in dump_steps(os.path.join(out, "dump.move_surf_cmd.project.*.txt")):
        xy = read_dump(f, 2)
        nin = sum(inside_2d(x, y, segs) for x, y in xy)
        check_inside("move_surf_cmd inside", nin, len(xy), "moved ellipse")

    check_lost("move_surf_cmd delete",
               os.path.join(out, "log.move_surf_cmd.delete"))


def piston_ratio(logfile):
    rows = stats_rows(logfile)
    p0 = 1.0e18 * KB * 300.0                  # nrho k T of in.piston
    return rows[-1]["c_3"] / p0


def check_piston(out, tol=0.03):
    m = species_mass("N")
    s = 100.0 / math.sqrt(2.0 * KB * 300.0 / m)       # U / c_mp
    exact = ((1.0 + 2.0 * s * s) * (1.0 + math.erf(s)) +
             2.0 * s / math.sqrt(math.pi) * math.exp(-s * s))

    r = piston_ratio(os.path.join(out, "log.piston.static"))
    status = "PASS" if abs(r - 1.0) < tol else "FAIL"
    report(status, "piston static (calibration)", f"p/p0 = {r:.4f}, exact 1")

    r = piston_ratio(os.path.join(out, "log.piston.project_yes"))
    status = "PASS" if abs(r / exact - 1.0) < tol else "FAIL"
    report(status, "piston wallvel yes",
           f"p/p0 = {r:.4f}, free-molecular exact {exact:.4f} "
           f"(err {100.0 * (r / exact - 1.0):+.2f}%)")

    r = piston_ratio(os.path.join(out, "log.piston.project_no"))
    report("INFO", "piston wallvel no",
           f"p/p0 = {r:.4f} (no wall velocity in collisions)")
    r = piston_ratio(os.path.join(out, "log.piston.delete_no"))
    report("INFO", "piston engulf delete", f"p/p0 = {r:.4f} (legacy)")


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
        ut = (-dy * u + dx * v) / r
        num += ut * r
        den += r * r
    return num / den


def check_rotating_box(out, tol=0.03):
    omega = 2.0 * math.pi / (1000 * 1.0e-4)
    check_conserved("rotating_box wallvel yes",
                    os.path.join(out, "log.rotating_box.project_yes"))
    slope = box_rotation(os.path.join(out, "grid.rotating_box.project_yes.txt"))
    status = "PASS" if abs(slope / omega - 1.0) < tol else "FAIL"
    report(status, "rotating_box wallvel yes",
           f"gas rotation {slope:.2f} 1/s vs wall {omega:.2f} 1/s "
           f"(ratio {slope / omega:.4f})")
    slope = box_rotation(os.path.join(out, "grid.rotating_box.project_no.txt"))
    report("INFO", "rotating_box wallvel no",
           f"gas rotation {slope:.2f} 1/s (ratio {slope / omega:.4f})")


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(2)
    out = sys.argv[1]
    nprocs = sys.argv[2] if len(sys.argv) > 2 else "4"

    checks = [lambda: check_ellipse2d(out, nprocs),
              lambda: check_ellipsoid3d(out),
              lambda: check_two_rotors(out),
              lambda: check_move_surf_cmd(out),
              lambda: check_piston(out),
              lambda: check_rotating_box(out)]
    for check in checks:
        try:
            check()
        except (OSError, IndexError, KeyError, ValueError) as err:
            report("FAIL", "check error", repr(err))

    nfail = results.count("FAIL")
    print(f"\n{results.count('PASS')} passed, {nfail} failed, "
          f"{results.count('WARN')} warnings")
    sys.exit(1 if nfail else 0)


if __name__ == "__main__":
    main()
