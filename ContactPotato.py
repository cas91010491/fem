from PyClasses.FEAssembly import *
from PyClasses.Contacts import *
from PyClasses.FEModel import *

import meshio, sys, os
import numpy as np
import matplotlib.pyplot as plt
from pdb import set_trace
from math import pi

####################
### GETTING DATA ###
####################
os.chdir(sys.path[0])

# BLOCK
# mesh   = meshio.read("Meshes/Block50x50x24.msh")
mesh_blk   = meshio.read("Meshes/Cubes/cube843.msh")
X_blk     = mesh_blk.points
hexas_blk = mesh_blk.cells_dict['hexahedron']

# POTATO
mesh_ptt = meshio.read("Meshes/QuadSpheres/QuadSphere4.msh")
# mesh_ptt = meshio.read("Meshes/Cubes/cube843.msh")
X_ptt = mesh_ptt.points
hexas_ptt = mesh_ptt.cells_dict['hexahedron']


######################
### BUILDING MODEL ###
######################

recOut = False          # if this was the last script run, This can be set to True to skip outers computations
### BODIES ###
blk = FEAssembly(X_blk,hexas_blk, name= "BLOCK",recOuters=recOut)
blk.Youngsmodulus = 0.5
ptt = FEAssembly(X_ptt,hexas_ptt, name = "POTATO",recOuters=recOut)
ptt.isRigid = True     # faster solving when True


### Transformations ###
ptt.Resize(4.0,dir='x')
ptt.Resize(3.0,dir='y')
ptt.Resize(2.0,dir='z')
# ptt.Bulge(2.0,center=(3,0,0))
# ptt.Bulge(2.0,center=(8,0,0))
# ptt.Bulge(-2.0,center=(-2.5,0,0))
# ptt.Rotate(pi/2,"y")                            # 22.5° to avoid double master ambiguity
ptt.RandDistort(0.5)
gap = 2.0                                       # Initial gap between ball and blk. ( \!/ doesn't consider smoothing )
# pttUp = max(blk.X[:,2])-min(ptt.X[:,2])+gap
ptt.Translate([6.0,0.0,0.0])
blk.Resize(0.5)
blk.Translate([0.0,0.0,2.0])

### Selections (nodes) ###
blk_bottom  = blk.SelectFlatSide("-z")
blk_top     = blk.SelectFlatSide("+z")
blk_top_few = blk.SelectNodesBySphere([0.0, 0.0, 2.5],3.0, OnSurface=True)
ptt_all     = ptt.SelectAll()
ptt_top     = ptt.SelectFlatSide("+z")
top_ptt     = max(ptt.X[:,2])
ptt_top_few = ptt.SelectNodesBySphere([-4.0, 0.0, top_ptt],2*0.95, OnSurface=True)


blk.Rotate(pi/8,"y") 

### BOUNDARY CONDITIONS ###  [body, nodes, type, directions, values, times(*)]
cond_bd1 = [ptt, ptt.SelectAll(), "dirichlet", "xyz"  , [0.0, 0.0, 0.0] ]
cond_bd2 = [blk, blk_top   , "dirichlet", "xyz"  , [8.0, 0.0,0.0] ]
# cond_bd2 = [blk, blk_top   , "dirichlet", "xyz"  , [0.0, 0.0,-0.52], [0,0.1] ]
# cond_bd3 = [blk, blk_top   , "dirichlet", "xyz"  , [3.0, 9.0, 0.0], [0.1,1.0] ]

BCs = [cond_bd1, cond_bd2]

### CONTACTS ###            # [body, nodes]
slave   = [blk , blk_bottom]
# masterNodes = list(set(ptt.surf.nodes)-set(ptt_top_few))
# master  = [ptt, masterNodes ]
master = ptt

contact1 = Contact(slave, master, kn=0.01, cubicT=None, C1Edges = False, OutPatchAllowance=1e-6, maxGN = 0.05,f0=0.1)       # (slave, master) inputs can be surfaces as well

### MODEL ###
model = FEModel([blk, ptt], [contact1], BCs)           # [bodies, contacts, BCs, opts*]

ndofs = 3*(len(X_blk)+len(X_ptt))
ptt.surf.ComputeGrgPatches(np.zeros(ndofs),range(len(ptt.surf.nodes)))
model.plotNow()       # Uncomment to see and verify geometry

#############
## Running ##
#############
model.Solve(TimeSteps=100, recover=False, max_iter=20)
