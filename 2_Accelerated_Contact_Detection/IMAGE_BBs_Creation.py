import sys, os, pickle
import numpy as np
from pdb import set_trace
from time import time

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from PyClasses.FEAssembly import *
from PyClasses.Contacts import *
from PyClasses.FEModel import *
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from PyClasses.Utilities import plot_coords
#####################
### Setting model ###
#####################
os.chdir(sys.path[0])

# POTATO
[ptt] = pickle.load(open("PotatoAssembly.dat","rb"))
ptt.isRigid = True     # faster solving when True
ndofs = 3*len(ptt.X)

# MODEL
model = FEModel([ptt], [],[])           # [bodies, contacts, BCs, opts*]
ptt.Translate([-6.0, 0.0, 0.0])
model.X = ptt.X.ravel()
ptt.surf.ComputeGrgPatches(np.zeros(ndofs),range(len(ptt.surf.nodes)))

offset = 0.1
# n_per_side = 200 


model_X = model.X.reshape(-1,3)
xa,xb = min(model_X[:,0]), max(model_X[:,0])
ya,yb = min(model_X[:,1]), max(model_X[:,1])
za,zb = min(model_X[:,2]), max(model_X[:,2])

dx, dy, dz = offset*np.array([xb-xa,yb-ya,zb-za])
xmin, xmax = xa-dx, xb+dx
ymin, ymax = ya-dy, yb+dy
zmin, zmax = za-dz, zb+dz


import matplotlib.pyplot as plt

fig = plt.figure(figsize=(10, 10))

# Define the vertices of the bounding box
vertices = np.array([[xmin, ymin, zmin],
                     [xmin, ymin, zmax],
                     [xmin, ymax, zmin],
                     [xmin, ymax, zmax],
                     [xmax, ymin, zmin],
                     [xmax, ymin, zmax],
                     [xmax, ymax, zmin],
                     [xmax, ymax, zmax]])

# Define the 12 edges of the bounding box
edges = [[vertices[j] for j in [0, 1, 3, 2]],
         [vertices[j] for j in [4, 5, 7, 6]],
         [vertices[j] for j in [0, 1, 5, 4]],
         [vertices[j] for j in [2, 3, 7, 6]],
         [vertices[j] for j in [0, 2, 6, 4]],
         [vertices[j] for j in [1, 3, 7, 5]]]

# Upper left: view from +z
ax1 = fig.add_subplot(223, projection='3d', proj_type='ortho')
model.plot(ax1)
plot_coords(ax1,orig=(-4,-3,0))

ax1.add_collection3d(Poly3DCollection(edges, facecolors='cyan', linewidths=1, edgecolors='r', alpha=.06))
ax1.view_init(elev=90, azim=-90)
ax1.set_box_aspect([1,1,1])  # Aspect ratio is 1:1:1
ax1.set_aspect('equal')
ax1.set_zticks([])  # Remove z-axis ticks and numbers
ax1.set_title('View from +z')

# Upper right: view from +x
ax2 = fig.add_subplot(222, projection='3d', proj_type='ortho')
model.plot(ax2)
plot_coords(ax2,orig=(0,-3,-2))

ax2.add_collection3d(Poly3DCollection(edges, facecolors='cyan', linewidths=1, edgecolors='r', alpha=.06))
ax2.view_init(elev=0, azim=0)
ax2.set_box_aspect([1,1,1])  # Aspect ratio is 1:1:1
ax2.set_aspect('equal')
ax2.set_xticks([])  # Remove x-axis ticks and numbers
ax2.set_title('View from +x')

# Lower left: view from -y
ax3 = fig.add_subplot(221, projection='3d', proj_type='ortho')
model.plot(ax3)
plot_coords(ax3,orig=(-4,0,-2))
ax3.add_collection3d(Poly3DCollection(edges, facecolors='cyan', linewidths=1, edgecolors='r', alpha=.06))
ax3.view_init(elev=0, azim=-90)
ax3.set_box_aspect([1,1,1])  # Aspect ratio is 1:1:1
ax3.set_aspect('equal')
ax3.set_yticks([])  # Remove y-axis ticks and numbers
ax3.set_title('View from -y')

# Lower right: isometric view
ax4 = fig.add_subplot(224, projection='3d', proj_type='ortho')
model.plot(ax4)
plot_coords(ax4,orig=(4,-3,0))
ax4.add_collection3d(Poly3DCollection(edges, facecolors='cyan', linewidths=1, edgecolors='r', alpha=.06))
ax4.view_init(elev=30, azim=45)
ax4.set_box_aspect([1,1,1])  # Aspect ratio is 1:1:1
ax4.set_aspect('equal')
ax4.set_position([0.55, 0.1, 0.35, 0.35])  # Adjust the position and size of the subplot

plt.subplots_adjust(left=0, right=1, top=1, bottom=0, wspace=0, hspace=0)

# Save the figure as .png
fig.savefig('figure.png', format='png')

# Save the figure as .pdf
fig.savefig('figure.pdf', format='pdf')

plt.show()
