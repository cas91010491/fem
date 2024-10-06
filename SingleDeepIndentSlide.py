from PyClasses.FEAssembly import *
from PyClasses.FEModel import *
from PyClasses.Contacts import *

import meshio
import numpy as np
import matplotlib.pyplot as plt
from pdb import set_trace

####################
### GETTING DATA ###
####################
os.chdir(sys.path[0])

# CUBETOP
# mesh_top  = meshio.read("Meshes/Cubes/cube222.msh")
mesh_top  = meshio.read("Meshes/QuadSpheres/HQuadSphere9.msh")
X_top     = mesh_top.points
hexas_top = mesh_top.cells_dict['hexahedron']
# CUBETBOT
mesh_bot  = meshio.read("Meshes/Cubes/cube20x20x4side5.msh")
# mesh_bot  = meshio.read("Meshes/Cubes/cube111.msh")
X_bot     = mesh_bot.points
hexas_bot = mesh_bot.cells_dict['hexahedron']


######################
### BUILDING MODEL ###
######################

### BODIES ###
recOut = False          # if this was the last script run, This can be set to True to skip outers computations
UpperCube = FEAssembly(X_top, hexas_top, name = "UpperCube",recOuters=recOut) 
LowerCube = FEAssembly(X_bot, hexas_bot, name = "LowerCube",recOuters=recOut)
# LowerCube.isRigid = True
UpperCube.Resize(2)
UpperCube.Resize(0.5,dir="z")

### SELECTIONS ###
# UpperCube_top     = UpperCube.SelectFlatSide("+z")
UpperCube_top     = UpperCube.SelectHigherThan("z",-0.1)
# UpperCube_top     = UpperCube.SelectAll()
# UpperCube_bottom  = [24]
# LowerCube_bottom  = LowerCube.SelectAll()
LowerCube_bottom  = LowerCube.SelectFlatSide("-z")
LowerCube_top     = LowerCube.SelectFlatSide("+z")

### TRANSFORMATIONS ###
# UpperCube.X[24][2] -= 1.0
# UpperCube.X[26][2] -= 0.5
# LowerCube.Resize(5)

gap = 0.05
UpperCube.Translate([ -1.75, 0.0, gap/2-min([z for z in UpperCube.X[:,2]]) ])
LowerCube.Translate([ 0.0, 0.0,-gap/2-max([z for z in LowerCube.X[:,2]]) ])
UpperCube_bottom  = UpperCube.SelectLowerThan("z",0.5,OnSurface=True)


### BOUNDARY CONDITIONS ###  [body, nodes, type, directions, values, times(*)]
cond_bd1 = [LowerCube , LowerCube_bottom , "dirichlet", "xyz"  , [0.0, 0.0, 0.0] ]
# cond_bd2 = [UpperCube , UpperCube_top    , "dirichlet", "xyz"  , [0.0, 0.0,-3.0] ]
cond_bd2 = [UpperCube , UpperCube_top         , "dirichlet", "xyz"  , [0.0, 0.0,-0.5], [0.0, 0.12] ]
cond_bd3 = [UpperCube , UpperCube_top         , "dirichlet", "xyz"  , [4.0, 0.0, 0.0], [0.12, 1.0] ]
BCs = [cond_bd1, cond_bd2, cond_bd3]
# BCs = [cond_bd1, cond_bd2]

### CONTACTS ###
slave  = [UpperCube , UpperCube_bottom]
master = [LowerCube , LowerCube_top   ]

contact1 = Contact(slave, master, kn=1, kt= 100, cubicT=None, C1Edges = False, OutPatchAllowance=1e-3, f0=0, maxGN = 0.001, mu=0.2)       # (slave, master) inputs can be surfaces as well
ctcts = [contact1]

### MODEL ###
model = FEModel([LowerCube, UpperCube], ctcts, BCs, UpdateOnIteration=True,)     # Should include BCs
# model.plotNow()

#############
## Running ##
#############
model.Solve(TimeSteps=1000, recover=True, max_iter=10)
