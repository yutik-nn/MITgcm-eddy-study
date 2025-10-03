#!/usr/bin/env python3
# scripts/gen_vel.py
import numpy as np
from pathlib import Path

# --- paths & grid (adjust if you change your grid) ---
BASE = Path(__file__).resolve().parents[1]      # .../v0
INP  = BASE / "input"
OUTP = BASE / "analysis"
OUTP.mkdir(exist_ok=True)

nx, ny = 200, 100
dx, dy = 5e3, 5e3              # meters

# --- physics ---
g   = 9.81
f0  = 5.0e-5                   # set beta=0 for now (f-plane)
beta= 0.0

# --- load η ---
eta = np.fromfile(INP/"initial_eta.bin", dtype=">f4").reshape(ny, nx)
eta -= eta.mean()              # zero-mean to avoid spurious volume change

# --- C-grid locations & Coriolis on faces ---
yU = (np.arange(ny))   * dy     # U at south faces
yV = (np.arange(ny)+1) * dy     # V at north faces
fU = f0 + beta*(yU - 0.0)[:,None]
fV = f0 + beta*(yV - 0.0)[:,None]

# --- centered gradients mapped to faces (periodic-roll for interior stencil) ---
deta_dy_U = (eta - np.roll(eta, 1, axis=0)) / dy
deta_dx_V = (eta - np.roll(eta, 1, axis=1)) / dx

# --- geostrophic balance on faces ---
uU = -(g/fU) * deta_dy_U
vV =  (g/fV) * deta_dx_V

# optional: tame boundaries
uU[0,:]  = 0.0; uU[-1,:] = 0.0
vV[:,0]  = 0.0; vV[:,-1] = 0.0

# --- write binaries (big-endian float32, shape Ny x Nx) ---
uU.astype(">f4").tofile(INP/"u_init.bin")
vV.astype(">f4").tofile(INP/"v_init.bin")
print(f"Saved u_init.bin, v_init.bin  |  Umax={np.abs(uU).max():.3f}  Vmax={np.abs(vV).max():.3f} m/s")

# --- quick preview (centered for visualization) ---
try:
    depth = np.fromfile(INP/"bathy.bin", dtype=">f4").reshape(ny, nx)
    land  = depth == 0.0
except Exception:
    land = np.zeros((ny,nx), bool)

# average face velocities to centers for plotting
uc = 0.5*(uU + np.roll(uU, -1, axis=0))
vc = 0.5*(vV + np.roll(vV, -1, axis=1))
speed = np.sqrt(uc**2 + vc**2)

import matplotlib
matplotlib.use("Agg")  # non-blocking
import matplotlib.pyplot as plt
xC = (np.arange(nx)+0.5)*dx/1e3
yC = (np.arange(ny)+0.5)*dy/1e3
XC, YC = np.meshgrid(xC, yC, indexing="xy")

plt.figure(figsize=(8,4))
pcm = plt.pcolormesh(XC, YC, speed, shading="nearest")
cb  = plt.colorbar(pcm); cb.set_label("Speed (m s$^{-1}$)")
q = 5
plt.quiver(XC[::q,::q], YC[::q,::q], uc[::q,::q], vc[::q,::q],
           scale=5, width=0.003)
if land.any():
    plt.contour(XC, YC, land.astype(float), levels=[0.5], colors="k", linewidths=1.0)
plt.xlabel("x (km)"); plt.ylabel("y (km)")
plt.title("Initial velocity (geostrophic from η)")
plt.axis("equal"); plt.tight_layout()
plt.savefig(OUTP/"vel_preview.png", dpi=150)
print("Preview saved to analysis/vel_preview.png")