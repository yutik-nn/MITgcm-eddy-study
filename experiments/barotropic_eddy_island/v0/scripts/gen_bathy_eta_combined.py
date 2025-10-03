import matplotlib
matplotlib.use('TkAgg')  # Enable GUI backend
import numpy as np
import matplotlib.pyplot as plt

# Grid parameters
nx = 200    # 1000 km / 5 km
ny = 100    # 500 km / 5 km
dx = 5.e3   # 5 km
dy = 5.e3

# Create coordinate arrays
X, Y = np.meshgrid(np.arange(nx)*dx, np.arange(ny)*dy)

# Island
cx = 300.e3  # 300 km from west
cy = 250.e3  # centered in y
radius = 50.e3      # 50 km radius
slope_width = 30.e3  # 30 km slope

# Distance from center
r = np.sqrt((X - cx)**2 + (Y - cy)**2)

# Initialize depth: 1000 m everywhere
depth = np.ones((ny, nx)) * 1000.0

# Smooth ramp: depth decreases as r increases
mask = r < radius + slope_width
taper = np.zeros_like(r)
taper[mask] = np.minimum(1.0, (radius + slope_width - r[mask]) / slope_width)
depth[mask] = 1000.0 * (1.0 - taper[mask])

# Inside the island, set depth to 0 (dry land)
depth[r < radius] = 0.0

# Save bathymetry
depth.astype('>f4').tofile('../input/bathy.bin')

# --- Gaussian Eddy ---
ex = 800.e3   # 800 km from west (east of island)
ey = 250.e3   # centered in y
eddy_radius = 100.e3
eta0_max = 0.5  # 50 cm SSH anomaly

# Gaussian SSH
r2 = (X - ex)**2 + (Y - ey)**2
eta0 = eta0_max * np.exp(-r2 / (2 * eddy_radius**2))

# Save initial SSH
eta0.astype('>f4').tofile('../input/initial_eta.bin')

# --- Plot 2D ---
plt.figure(figsize=(12, 6))

# Bathymetry
plt.subplot(1, 2, 1)
cf = plt.contourf(X/1000, Y/1000, depth, levels=50, cmap='Blues_r')
plt.colorbar(cf, label='Depth (m)')
plt.contour(X/1000, Y/1000, depth, levels=[0], colors='k', linewidths=2)
plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.title('Bathymetry: Cylindrical Island')
plt.axis('equal')

# Initial SSH
plt.subplot(1, 2, 2)
cf = plt.contourf(X/1000, Y/1000, eta0, levels=20, cmap='RdYlBu_r')
plt.colorbar(cf, label='SSH (m)')
#plt.contour(X/1000, Y/1000, eta0, levels=np.linspace(0.1, 0.5, 5), colors='k', linewidths=0.5)
plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.title('Initial SSH: Gaussian Eddy')
plt.axis('equal')

plt.tight_layout()
plt.savefig('../analysis/island_and_eddy.png', dpi=150, bbox_inches='tight')
plt.show()