#!/usr/bin/env python3
"""
Write the surf files used by in.test_translation and in.test_rotation

usage: python3 make_surfs.py   (writes into the directory of this script)

data.piston3d   plate x = 0.2 spanning y,z in [0,0.2], normal +x
data.ring_inner 2d circle of radius 0.5, 72 lines, normals outward (gas outside)
data.ring_outer 2d circle of radius 1.0, 72 lines, normals inward (gas inside)

SPARTA tri normal = (p2 - p1) x (p3 - p1) and line normal = z x (p2 - p1)
point into the gas.
The rotation test is 2d: a 3d tube touching the z faces of the box cannot be
rotated, round-off of the rotation matrix moves its end points out of the box.
"""

import math
import os

HERE = os.path.dirname(os.path.abspath(__file__))


def write_surf(fname, title, points, tris):
    with open(os.path.join(HERE, fname), "w") as fp:
        fp.write(f"{title}\n\n{len(points)} points\n{len(tris)} triangles\n\n")
        fp.write("Points\n\n")
        for i, (x, y, z) in enumerate(points, 1):
            fp.write(f"{i} {x:.15g} {y:.15g} {z:.15g}\n")
        fp.write("\nTriangles\n\n")
        for i, (a, b, c) in enumerate(tris, 1):
            fp.write(f"{i} {a} {b} {c}\n")


def write_lines(fname, title, points, lines):
    with open(os.path.join(HERE, fname), "w") as fp:
        fp.write(f"{title}\n\n{len(points)} points\n{len(lines)} lines\n\n")
        fp.write("Points\n\n")
        for i, (x, y) in enumerate(points, 1):
            fp.write(f"{i} {x:.15g} {y:.15g}\n")
        fp.write("\nLines\n\n")
        for i, (a, b) in enumerate(lines, 1):
            fp.write(f"{i} {a} {b}\n")


def ring(radius, nline, outward):
    """counter-clockwise points, z x (p2 - p1) is outward if traversed
    clockwise"""
    pts = [(radius * math.cos(2.0 * math.pi * k / nline),
            radius * math.sin(2.0 * math.pi * k / nline)) for k in range(nline)]
    lines = []
    for k in range(nline):
        a, b = k + 1, (k + 1) % nline + 1
        lines.append((b, a) if outward else (a, b))
    return pts, lines


def main():
    plate = [(0.2, 0.0, 0.0), (0.2, 0.2, 0.0), (0.2, 0.2, 0.2), (0.2, 0.0, 0.2)]
    write_surf("data.piston3d", "piston plate x = 0.2, normal +x, gas on +x side",
               plate, [(1, 2, 3), (1, 3, 4)])

    pts, lines = ring(0.5, 72, outward=True)
    write_lines("data.ring_inner", "inner circle R = 0.5, 72 lines, normals "
                "outward", pts, lines)
    pts, lines = ring(1.0, 72, outward=False)
    write_lines("data.ring_outer", "outer circle R = 1.0, 72 lines, normals "
                "inward", pts, lines)
    print("wrote data.piston3d data.ring_inner data.ring_outer")


if __name__ == "__main__":
    main()
