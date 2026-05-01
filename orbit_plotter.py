import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from mpl_toolkits.mplot3d.art3d import Line3DCollection
import matplotlib.colors as mcolors
import os


data = np.loadtxt(r"C:\Users\sghob\OneDrive\Desktop\Res_code\orbit_test25.txt", delimiter="\t")
#print(data[:,0])
# Example data: Replace these with your actual (r, phi, z) data
#data = data[30 * len(data) // 32: 8 * len(data) // 8:10]  
#data = data[len(data)//3::10]  
#data = data[0 * len(data) // 32:2 * len(data) // 32:10]  
#data = data[:]
r = data[:,1] / 3.086E19  # Radial distance
phi = data[:,2] # Azimuthal angle (radians)
z = data[:,3] / 3.086E19  # Vertical coordinate
t = data[:,0] / (3.15E7 * 1E9)
# Convert to Cartesian coordinates
x = r * np.cos(phi)
y = r * np.sin(phi)
# Create segments for coloring
points = np.array([x, y, z]).T.reshape(-1, 1, 3)
segments = np.concatenate([points[:-1], points[1:]], axis=1)

# Create a colormap
norm = plt.Normalize(t.min(), t.max())
colors = plt.cm.cividis(norm(t))  # You can change 'plasma' to other colormaps

# Create 3D plot
fig = plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111, projection='3d')

# Add colored line segments
lc = Line3DCollection(segments, cmap='magma', norm=norm, linewidth=1)
lc.set_array(t)
ax.add_collection(lc)

'''
# Color bar
cbar = plt.colorbar(lc, ax=ax)
cbar.set_label("Time (Gyr)")
'''
# Colorbar with padding to prevent overlap
norm = mcolors.Normalize(vmin=0, vmax=2)  # Set the color scale range
lc.set_norm(norm)  # Apply the normalization to the LineCollection

cbar = fig.colorbar(lc, ax=ax, pad=0.1)
cbar.set_label("Time (Gyr)", fontsize=18)
cbar.ax.tick_params(labelsize=18)
# Optional: Adjust layout if needed
plt.tight_layout()

# Labels
ax.set_xlabel("X (kpc)", fontsize=18, labelpad = 12)
ax.set_ylabel("Y (kpc)", fontsize=18, labelpad = 12)
ax.set_zlabel("Z (kpc)", fontsize=18, labelpad = 12)
#ax.set_title(r"3D Orbit at 75$\degree$ Inclination")
ax.set_xlim([-1.5,1.5])
ax.set_ylim([-1.5,1.5])
ax.set_zlim([-1.5,1.5])
ax.tick_params(labelsize = 18)

plt.savefig("orbit45.pdf")
plt.show()