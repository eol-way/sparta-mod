#!/usr/bin/env python3
"""
GIF animations of the moving wall validation runs (generate_all_gifs.sh)

usage: python3 make_gifs.py ANIMDIR GIFDIR [name ...]

ANIMDIR holds the anim.*.particle.STEP.txt and anim.*.surf.STEP.txt dumps
written by the inputs with -var anim N.  One GIF per test is written to
GIFDIR.  Names: ellipse2d two_rotors ellipsoid3d piston rotating_box
(default all).  Needs numpy, matplotlib and Pillow.

Particles inside a surf are drawn in red and counted in the frame title,
so a frame with "inside 0" shows that no particle was left in the wall.
"""

import glob
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                      # noqa: E402
from matplotlib.animation import PillowWriter        # noqa: E402
from matplotlib.collections import LineCollection    # noqa: E402

FPS = 10
DPI = 80
WALL = "0.15"          # wall line color
FILL = "0.80"          # solid fill color
PART = "tab:blue"      # particle color
BAD = "red"            # particles inside a surf


# ----------------------------------------------------------------------
# dump readers
# ----------------------------------------------------------------------

def frames(animdir, prefix):
    """{step: file} for dump files prefix.STEP.txt"""
    out = {}
    for f in glob.glob(os.path.join(animdir, prefix + ".*.txt")):
        out[int(os.path.basename(f).split(".")[-2])] = f
    return out


def read_dump(fname):
    """columns {name: array} of one SPARTA particle or surf dump file"""
    with open(fname) as fp:
        lines = fp.read().splitlines()
    i = [k for k, l in enumerate(lines) if l.startswith("ITEM: ATOMS") or
         l.startswith("ITEM: SURFS")][0]
    names = lines[i].split()[2:]
    rows = [l.split() for l in lines[i + 1:] if l.strip()]
    data = np.array(rows, dtype=float).reshape(-1, len(names))
    return {n: data[:, k] for k, n in enumerate(names)}


def common_steps(*framesets):
    steps = set(framesets[0])
    for f in framesets[1:]:
        steps &= set(f)
    return sorted(steps)


# ----------------------------------------------------------------------
# geometry
# ----------------------------------------------------------------------

def segments_2d(surf):
    return np.column_stack([surf["v1x"], surf["v1y"], surf["v2x"], surf["v2y"]])


def loops_2d(segs, tol=1.0e-6):
    """chain segments v1 -> v2 into closed polygons (traversal order kept)"""
    key = lambda x, y: (round(x / tol), round(y / tol))
    nxt = {key(s[0], s[1]): (s[2], s[3]) for s in segs}
    loops, used = [], set()
    for s in segs:
        k0 = key(s[0], s[1])
        if k0 in used:
            continue
        poly, k = [], k0
        while k in nxt and k not in used:
            used.add(k)
            x, y = nxt[k]
            poly.append((x, y))
            k = key(x, y)
        if len(poly) > 2 and k == k0:
            loops.append(np.array(poly))
    return loops


def signed_area(poly):
    x, y = poly[:, 0], poly[:, 1]
    return 0.5 * np.sum(x * np.roll(y, -1) - np.roll(x, -1) * y)


def solid_2d(px, py, segs):
    """points inside the solid of closed 2d surfs

    SPARTA line normal = z x (v2 - v1) points into the gas, so a loop
    traversed clockwise has solid inside (rotor) and a loop traversed
    counter-clockwise has gas inside (container).  The winding number is
    0 in the gas of a rotor-only geometry and 1 in the gas of a container,
    every other value is solid.
    """
    w = np.zeros(len(px), dtype=int)
    for x1, y1, x2, y2 in segs:
        up = (y1 <= py) & (y2 > py)
        down = (y1 > py) & (y2 <= py)
        cross = up | down
        if not cross.any():
            continue
        xi = x1 + (py[cross] - y1) * (x2 - x1) / (y2 - y1)
        right = px[cross] < xi
        w[cross] += np.where(right, np.where(up[cross], 1, -1), 0)
    gas = 1 if any(signed_area(p) > 0.0 for p in loops_2d(segs)) else 0
    return w != gas


def tri_normals(surf):
    v1 = np.column_stack([surf["v1x"], surf["v1y"], surf["v1z"]])
    v2 = np.column_stack([surf["v2x"], surf["v2y"], surf["v2z"]])
    v3 = np.column_stack([surf["v3x"], surf["v3y"], surf["v3z"]])
    return v1, v2, v3, np.cross(v2 - v1, v3 - v1)       # outward (SPARTA)


def inside_convex_3d(p, surf):
    """points inside a convex closed tri surf (all faces see them inside)"""
    v1, _, _, n = tri_normals(surf)
    if len(p) == 0:
        return np.zeros(0, dtype=bool)
    d = np.einsum("pk,fk->pf", p, n) - np.einsum("fk,fk->f", v1, n)[None, :]
    return np.all(d < 0.0, axis=1)


def section_z0(surf):
    """segments where the tri surf cuts the plane z = 0"""
    v1, v2, v3, _ = tri_normals(surf)
    segs = []
    for a, b, c in zip(v1, v2, v3):
        pts = []
        for p, q in ((a, b), (b, c), (c, a)):
            if (p[2] > 0.0) != (q[2] > 0.0):
                t = p[2] / (p[2] - q[2])
                pts.append(p[:2] + t * (q[:2] - p[:2]))
        if len(pts) == 2:
            segs.append(np.concatenate(pts))
    return np.array(segs).reshape(-1, 4)


# ----------------------------------------------------------------------
# drawing helpers
# ----------------------------------------------------------------------

def draw_wall(ax, segs, fill=True, color=WALL, lw=1.8):
    if fill:
        for poly in loops_2d(segs):
            if signed_area(poly) < 0.0:          # clockwise = solid inside
                ax.fill(poly[:, 0], poly[:, 1], color=FILL, zorder=1)
    ax.add_collection(LineCollection(segs.reshape(-1, 2, 2), colors=color,
                                     linewidths=lw, zorder=3))


def draw_particles(ax, x, y, inside, size=0.6):
    ax.scatter(x[~inside], y[~inside], s=size, c=PART, lw=0, zorder=2)
    if inside.any():
        ax.scatter(x[inside], y[inside], s=12, c=BAD, lw=0, zorder=4)


def write_gif(fig, path, steps, draw):
    writer = PillowWriter(fps=FPS)
    with writer.saving(fig, path, DPI):
        for step in steps:
            draw(step)
            writer.grab_frame()
    plt.close(fig)
    print(f"wrote {path} ({len(steps)} frames, "
          f"{os.path.getsize(path) / 1.0e6:.1f} MB)")


# ----------------------------------------------------------------------
# one GIF per test
# ----------------------------------------------------------------------

def gif_ellipse2d(animdir, gifdir):
    modes = ["delete", "project"]
    pf = {m: frames(animdir, f"anim.ellipse2d.{m}.particle") for m in modes}
    sf = {m: frames(animdir, f"anim.ellipse2d.{m}.surf") for m in modes}
    steps = common_steps(*pf.values(), *sf.values())
    if not steps:
        return
    fig, axes = plt.subplots(1, 2, figsize=(10, 5.4))

    def draw(step):
        for ax, m in zip(axes, modes):
            ax.clear()
            p = read_dump(pf[m][step])
            segs = segments_2d(read_dump(sf[m][step]))
            ins = solid_2d(p["x"], p["y"], segs)
            draw_wall(ax, segs)
            draw_particles(ax, p["x"], p["y"], ins)
            ax.set_xlim(0, 10)
            ax.set_ylim(0, 10)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"engulf {m}\nparticles {len(ins)}  inside wall "
                         f"{int(ins.sum())}", fontsize=11)
        fig.suptitle(f"rotating ellipse, 10% particle sample, step {step}",
                     fontsize=12)

    write_gif(fig, os.path.join(gifdir, "ellipse2d.gif"), steps, draw)


def gif_two_rotors(animdir, gifdir):
    pf = frames(animdir, "anim.two_rotors.particle")
    sf = frames(animdir, "anim.two_rotors.surf")
    steps = common_steps(pf, sf)
    if not steps:
        return
    fig, ax = plt.subplots(figsize=(6, 6.3))

    def draw(step):
        ax.clear()
        p = read_dump(pf[step])
        segs = segments_2d(read_dump(sf[step]))
        ins = solid_2d(p["x"], p["y"], segs)
        draw_wall(ax, segs)
        draw_particles(ax, p["x"], p["y"], ins)
        ax.set_xlim(1, 9)
        ax.set_ylim(1, 9)
        ax.set_aspect("equal")
        ax.set_xticks([])
        ax.set_yticks([])
        ax.set_title(f"2 counter-rotating rotors, engulf project, step {step}\n"
                     f"particles {len(ins)} (10% sample)  inside wall "
                     f"{int(ins.sum())}", fontsize=11)

    write_gif(fig, os.path.join(gifdir, "two_rotors.gif"), steps, draw)


def gif_ellipsoid3d(animdir, gifdir):
    modes = ["delete", "project"]
    pf = {m: frames(animdir, f"anim.ellipsoid3d.{m}.particle") for m in modes}
    sf = {m: frames(animdir, f"anim.ellipsoid3d.{m}.surf") for m in modes}
    steps = common_steps(*pf.values(), *sf.values())
    if not steps:
        return
    fig, axes = plt.subplots(1, 2, figsize=(10, 5.4))

    def draw(step):
        for ax, m in zip(axes, modes):
            ax.clear()
            p = read_dump(pf[m][step])
            surf = read_dump(sf[m][step])
            xyz = np.column_stack([p["x"], p["y"], p["z"]])
            ins = inside_convex_3d(xyz, surf)
            draw_wall(ax, section_z0(surf))
            draw_particles(ax, p["x"], p["y"], ins, size=2.0)
            ax.set_xlim(-2, 2)
            ax.set_ylim(-2, 2)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"engulf {m}\nparticles {len(ins)}  inside wall "
                         f"{int(ins.sum())}", fontsize=11)
        fig.suptitle(f"rotating 3d ellipsoid, slab |z| < 0.2, step {step}",
                     fontsize=12)

    write_gif(fig, os.path.join(gifdir, "ellipsoid3d.gif"), steps, draw)


def gif_piston(animdir, gifdir):
    modes = [("project_yes", "engulf project, wallvel yes", "tab:blue"),
             ("project_no", "engulf project, wallvel no", "tab:orange"),
             ("delete_no", "engulf delete (legacy)", "tab:gray")]
    pf = {m: frames(animdir, f"anim.piston.{m}.particle") for m, _, _ in modes}
    sf = {m: frames(animdir, f"anim.piston.{m}.surf") for m, _, _ in modes}
    steps = common_steps(*pf.values(), *sf.values())
    if not steps:
        return
    bins = np.linspace(0.0, 10.0, 51)
    width = bins[1] - bins[0]
    n0 = {m: len(read_dump(pf[m][steps[0]])["x"]) / 9.0 for m, _, _ in modes}

    fig = plt.figure(figsize=(9, 7))
    grid = fig.add_gridspec(4, 1, height_ratios=[1, 1, 1, 2.2], hspace=0.45)
    axes = [fig.add_subplot(grid[k]) for k in range(4)]

    def draw(step):
        for ax, (m, label, color) in zip(axes[:3], modes):
            ax.clear()
            p = read_dump(pf[m][step])
            segs = segments_2d(read_dump(sf[m][step]))
            ax.scatter(p["x"], p["y"], s=0.3, c=color, lw=0)
            draw_wall(ax, segs, fill=False, color=BAD, lw=3)
            xp = segs[0][0]
            ax.axvspan(0, xp, color=FILL)
            ax.set_xlim(0, 10)
            ax.set_ylim(0, 1)
            ax.set_yticks([])
            ax.set_title(f"{label}: particles {len(p['x'])}", fontsize=9,
                         loc="left")
        ax = axes[3]
        ax.clear()
        for m, label, color in modes:
            x = read_dump(pf[m][step])["x"]
            h, _ = np.histogram(x, bins=bins)
            ax.step(bins[:-1], h / (n0[m] * width), where="post", color=color,
                    label=label)
        ax.axhline(1.0, color="k", lw=0.6, ls=":")
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 1.8)
        ax.set_xlabel("x (m)")
        ax.set_ylabel("n / n0")
        ax.legend(loc="upper right", fontsize=8)
        fig.suptitle(f"specular piston at 100 m/s into collisionless gas, "
                     f"step {step} (10% sample)", fontsize=11)

    write_gif(fig, os.path.join(gifdir, "piston.gif"), steps, draw)


def gif_rotating_box(animdir, gifdir, omega=2.0 * math.pi / (1000 * 1.0e-4)):
    modes = [("project_yes", "wallvel yes", "tab:blue"),
             ("project_no", "wallvel no", "tab:orange")]
    pf = {m: frames(animdir, f"anim.rotating_box.{m}.particle")
          for m, _, _ in modes}
    sf = {m: frames(animdir, f"anim.rotating_box.{m}.surf") for m, _, _ in modes}
    steps = common_steps(*pf.values(), *sf.values())
    if not steps:
        return
    rbins = np.linspace(0.0, 1.75, 8)
    rmid = 0.5 * (rbins[1:] + rbins[:-1])

    fig = plt.figure(figsize=(10, 7.5))
    grid = fig.add_gridspec(2, 2, height_ratios=[1.25, 1], hspace=0.3)
    axes = [fig.add_subplot(grid[0, 0]), fig.add_subplot(grid[0, 1]),
            fig.add_subplot(grid[1, :])]

    def draw(step):
        prof = {}
        for ax, (m, label, color) in zip(axes[:2], modes):
            ax.clear()
            p = read_dump(pf[m][step])
            segs = segments_2d(read_dump(sf[m][step]))
            ins = solid_2d(p["x"], p["y"], segs)
            ax.scatter(p["x"][~ins], p["y"][~ins], s=0.3, c=color, lw=0)
            if ins.any():
                ax.scatter(p["x"][ins], p["y"][ins], s=12, c=BAD, lw=0)
            draw_wall(ax, segs, fill=False)
            ax.set_xlim(2.1, 7.9)
            ax.set_ylim(2.1, 7.9)
            ax.set_aspect("equal")
            ax.set_xticks([])
            ax.set_yticks([])
            ax.set_title(f"{label}: inside wall {int(ins.sum())}", fontsize=10)

            dx, dy = p["x"] - 5.0, p["y"] - 5.0
            r = np.hypot(dx, dy)
            ok = (r > 0.0) & (r < rbins[-1])
            ut = (-dy[ok] * p["vx"][ok] + dx[ok] * p["vy"][ok]) / r[ok]
            k = np.digitize(r[ok], rbins) - 1
            prof[m] = [ut[k == j].mean() if (k == j).any() else np.nan
                       for j in range(len(rmid))]
        ax = axes[2]
        ax.clear()
        ax.plot(rmid, omega * rmid, "k--", label="rigid rotation  u = omega r")
        for m, label, color in modes:
            ax.plot(rmid, prof[m], "o-", color=color, label=label)
        ax.set_xlim(0, 1.75)
        ax.set_ylim(-30, 130)
        ax.set_xlabel("r (m)")
        ax.set_ylabel("mean u_theta (m/s)")
        ax.legend(loc="upper left", fontsize=9)
        fig.suptitle(f"square container rotating at {omega:.1f} rad/s, gas "
                     f"initially at rest, step {step} (25% sample)",
                     fontsize=11)

    write_gif(fig, os.path.join(gifdir, "rotating_box.gif"), steps, draw)


GIFS = {"ellipse2d": gif_ellipse2d, "two_rotors": gif_two_rotors,
        "ellipsoid3d": gif_ellipsoid3d, "piston": gif_piston,
        "rotating_box": gif_rotating_box}


def main():
    if len(sys.argv) < 3:
        print(__doc__)
        sys.exit(2)
    animdir, gifdir = sys.argv[1], sys.argv[2]
    names = sys.argv[3:] or list(GIFS)
    os.makedirs(gifdir, exist_ok=True)
    for name in names:
        GIFS[name](animdir, gifdir)


if __name__ == "__main__":
    main()
