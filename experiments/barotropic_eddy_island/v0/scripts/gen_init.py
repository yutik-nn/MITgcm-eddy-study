#!/usr/bin/env python3
import numpy as np, os

# --- grid/domain ---
Lx = Ly = 1024e3
Nx = Ny = 256
Dx = Lx/Nx; Dy = Ly/Ny
x = (np.arange(Nx) + 0.5) * Dx
y = (np.arange(Ny) + 0.5) * Dy
X, Y = np.meshgrid(x, y, indexing="xy")

# --- physics params ---
g    = 9.81
f0   = 5.0e-5
beta = 2.0e-11
H    = 4000.0  # (only used for context)

# --- eddy parameters (balanced Gaussian SSH) ---
x0, y0 = 250e3, 512e3  # initial eddy center
R      = 60e3          # e-folding radius (m)
eta0   = 0.25          # SSH amplitude (m), >0 => anticyclone in NH

r2   = (X-x0)**2 + (Y-y0)**2
etan = eta0 * np.exp(-0.5 * r2 / R**2)

# centered gradients at cell centers
deta_dx = (np.roll(etan,-1,axis=1) - np.roll(etan,1,axis=1)) / (2*Dx)
deta_dy = (np.roll(etan,-1,axis=0) - np.roll(etan,1,axis=0)) / (2*Dy)

# β-plane Coriolis
f = f0 + beta*(Y - 0.0)

# geostrophic balance: u = -(g/f)*∂η/∂y, v = (g/f)*∂η/∂x
u = -(g/f) * deta_dy
v =  (g/f) * deta_dx

os.makedirs("input", exist_ok=True)
etan.astype(">f4").tofile("input/etan_init.bin")
u.astype(">f4").tofile("input/uvel_init.bin")
v.astype(">f4").tofile("input/vvel_init.bin")
print("wrote etan/uvel/vvel init files", etan.shape, u.dtype,
      "Umax≈", float(np.abs(u).max()))