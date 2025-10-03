import matplotlib
matplotlib.use('TkAgg')  # Enable GUI backend

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # Required for 3D plots

# Grid parameters
nx = 200    # 1000 km / 5 km = 200 grid points
ny = 100    # 500 km / 5 km = 100 grid points
dx = 5.e3   # 5 km
dy = 5.e3

# Create coordinate arrays
X, Y = np.meshgrid(np.arange(nx)*dx, np.arange(ny)*dy)

# Center of island (now in the center of the domain)
cx = 500.e3  # 500 km from west → center
cy = 250.e3  # 250 km from south → center

# Island parameters [when I tried to make it stick out]
radius = 50.e3      # 50 km radius
elevation = 50.0     # Island sticks 50 m above sea level
slope_width = 30.e3  # Smooth slope over 30 km (~6 grid cells)

# Distance from center
r = np.sqrt((X - cx)**2 + (Y - cy)**2)

# Initialize depth: 1000 m everywhere
depth = np.ones((ny, nx)) * 1000.0

# Smooth ramp: depth decreases as r increases toward radius
mask = r < radius + slope_width
taper = np.zeros_like(r)
taper[mask] = np.minimum(1.0, (radius + slope_width - r[mask]) / slope_width)
depth[mask] = 1000.0 * (1.0 - taper[mask]) + elevation * taper[mask]

# Inside the island, set depth to negative (land above sea level)
depth[r < radius] = elevation

# Save as MITgcm binary (32-bit float, big-endian)
depth.astype('>f4').tofile('../input/bathy.bin')

# Plot 2D
plt.figure(figsize=(10, 5))
plt.imshow(depth, origin='lower', extent=[0, nx*dx/1000, 0, ny*dy/1000], cmap='terrain')
plt.colorbar(label='Depth (m)')
plt.contour(depth, levels=[0], colors='k', linewidths=1)
plt.xlabel('x (km)')
plt.ylabel('y (km)')
plt.title('Bathymetry: Cylindrical Island (0 m = sea level)')
plt.savefig('../analysis/bathy_2d.png', dpi=150, bbox_inches='tight')
plt.show()

# Plot 3D
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')
surf = ax.plot_surface(X/1000, Y/1000, depth - 1000, cmap='terrain', edgecolor='none', alpha=0.9)
ax.set_xlabel('x (km)')
ax.set_ylabel('y (km)')
ax.set_zlabel('Elevation relative to deep ocean (m)')
ax.set_title('3D Bathymetry: Cylindrical Island')
fig.colorbar(surf, ax=ax, shrink=0.5, pad=0.1, label='Depth (m)')
ax.view_init(elev=20, azim=45)
plt.savefig('../analysis/bathy_3d.png', dpi=150, bbox_inches='tight')
plt.show()  # This forces the 3D window to open

















# ------ old version
# import numpy as np
# import matplotlib.pyplot as plt

# # Grid
# nx, ny = 100, 50
# dx, dy = 10.e3, 10.e3
# X, Y = np.meshgrid(np.arange(nx)*dx, np.arange(ny)*dy)

# # Center of island
# cx, cy = 500.e3, 250.e3
# radius = 100.e3  # 100 km radius

# # Depth: 1000 m everywhere, 0 m (land) inside island
# depth = np.ones((ny, nx)) * 1000.
# depth[(X-cx)**2 + (Y-cy)**2 < radius**2] = 0.

# # Save as MITgcm binary
# depth.astype('>f4').tofile('../input/bathy.bin')

# Plot
# plt.imshow(depth, origin='lower')
# plt.colorbar(label='Depth (m)')
# plt.title('Cylindrical Island Bathymetry')
# plt.savefig('bathy.png')
# plt.show()