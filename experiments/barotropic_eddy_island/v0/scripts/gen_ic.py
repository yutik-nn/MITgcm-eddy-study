#!/usr/bin/env python3
"""
Generate IC/bathy for barotropic eddy–island experiments.
Writes big-endian float32 binaries into --outdir (default ../input).
Cases:
  - ideal : 1000x500 km @ 5 km; sloped-rim island (50 km core, 30 km slope);
            Gaussian eddy (R=100 km) to the EAST
  - test  : 600x400 km  @ 5 km; hard-wall cylinder island (70 km);
            Gaussian eddy (R=90 km) to the EAST
"""

import os
import argparse
import numpy as np

def write_bin(path, arr):
    os.makedirs(os.path.dirname(path), exist_ok=True)
    arr.astype('>f4').tofile(path)

def geostrophic_from_eta(eta, dx, dy, g=9.81, f0=1.0e-4):
    # np.gradient returns [d/dy, d/dx] when given (y, x)
    dη_dy, dη_dx = np.gradient(eta, dy, dx)
    u = -(g/f0) * dη_dy
    v =  (g/f0) * dη_dx
    return u, v

def generate_ideal(outdir="../input"):
    # Grid
    nx, ny = 200, 100
    dx = dy = 5.0e3
    x = np.arange(nx) * dx
    y = np.arange(ny) * dy
    X, Y = np.meshgrid(x, y)

    # Bathy: 1000 m with a sloped rim island
    cx, cy = 300e3, 250e3
    radius, slope = 50e3, 30e3
    r = np.hypot(X - cx, Y - cy)
    depth = np.full((ny, nx), 1000.0, dtype=np.float64)
    ring = r < radius + slope
    taper = np.zeros_like(r)
    taper[ring] = np.minimum(1.0, (radius + slope - r[ring]) / slope)
    depth[ring] = 1000.0 * (1.0 - taper[ring])
    depth[r < radius] = 0.0

    # SSH: Gaussian eddy to the EAST of island
    ex, ey, R, eta_max = 800e3, 250e3, 100e3, 0.5
    eta = eta_max * np.exp(-((X - ex)**2 + (Y - ey)**2) / (2.0 * R**2))

    # Balanced velocities; zero over land
    u, v = geostrophic_from_eta(eta, dx, dy)
    u[depth == 0.0] = 0.0
    v[depth == 0.0] = 0.0

    # Write
    write_bin(f"{outdir}/ideal_bathy.bin", depth)
    write_bin(f"{outdir}/ideal_eta.bin",   eta)
    write_bin(f"{outdir}/ideal_u.bin",     u)
    write_bin(f"{outdir}/ideal_v.bin",     v)
    print("Wrote ideal_* to", outdir, "| shape = (ny,nx) = (100,200)")

def generate_test(outdir="../input"):
    # Grid
    nx, ny = 120, 80
    dx = dy = 5.0e3
    x = np.arange(nx) * dx
    y = np.arange(ny) * dy
    X, Y = np.meshgrid(x, y)

    # Bathy: hard-wall cylindrical island
    cx, cy = 200e3, 200e3
    radius = 70e3
    r = np.hypot(X - cx, Y - cy)
    depth = np.full((ny, nx), 1000.0, dtype=np.float64)
    depth[r < radius] = 0.0

    # SSH: Gaussian eddy to the EAST of island
    ex, ey, R, eta_max = 450e3, 200e3, 90e3, 0.5
    eta = eta_max * np.exp(-((X - ex)**2 + (Y - ey)**2) / (2.0 * R**2))

    # Balanced velocities; zero over land
    u, v = geostrophic_from_eta(eta, dx, dy)
    u[depth == 0.0] = 0.0
    v[depth == 0.0] = 0.0

    # Write
    write_bin(f"{outdir}/test_bathy.bin", depth)
    write_bin(f"{outdir}/test_eta.bin",   eta)
    write_bin(f"{outdir}/test_u.bin",     u)
    write_bin(f"{outdir}/test_v.bin",     v)
    print("Wrote test_* to", outdir, "| shape = (ny,nx) = (80,120)")

def main():
    ap = argparse.ArgumentParser(description="Generate barotropic eddy–island inputs (no plotting).")
    ap.add_argument("--case", choices=["ideal","test"], default="ideal")
    ap.add_argument("--outdir", default="../input")
    args = ap.parse_args()

    if args.case == "ideal":
        generate_ideal(args.outdir)
    else:
        generate_test(args.outdir)

if __name__ == "__main__":
    main()