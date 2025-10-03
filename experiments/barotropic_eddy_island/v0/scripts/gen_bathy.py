#!/usr/bin/env python3
import numpy as np, os

# --- domain/grid ---
Lx = Ly = 1024e3      # m
Nx = Ny = 256
H  = 4000.0           # ocean depth (m), positive downward

# --- island geometry (true vertical cylinder) ---
xc, yc   = 600e3, 512e3  # center (m)
R_island = 60e3          # radius (m)

# grid (cell centers)
x = (np.arange(Nx) + 0.5) * (Lx/Nx)
y = (np.arange(Ny) + 0.5) * (Ly/Ny)
X, Y = np.meshgrid(x, y, indexing="xy")

# start ocean everywhere with depth H (positive)
bathy = H * np.ones((Ny, Nx), dtype=">f4")   # big-endian float32

# inside the circle: land (depth = 0)
mask = (X - xc)**2 + (Y - yc)**2 <= R_island**2
bathy[mask] = 0.0

os.makedirs("input", exist_ok=True)
bathy.tofile("input/bathy.bin")
print("wrote input/bathy.bin", bathy.shape, bathy.dtype,
      "ocean cells:", (bathy>0).sum(), "land cells:", (bathy==0).sum())