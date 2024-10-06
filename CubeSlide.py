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
# mesh_top  = meshio.read("Meshes/Cubes/cube5x5x5side4.msh")
mesh_top  = meshio.read("Meshes/Cubes/cube111.msh")
X_top     = mesh_top.points
hexas_top = mesh_top.cells_dict['hexahedron']
# CUBETBOT
mesh_bot  = meshio.read("Meshes/Cubes/cube5x5x1.msh")
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
LowerCube.isRigid = True

### SELECTIONS ###
UpperCube_top     = UpperCube.SelectFlatSide("+z")
UpperCube_bottom  = UpperCube.SelectFlatSide("-z")
LowerCube_bottom  = LowerCube.SelectFlatSide("-z")
LowerCube_top     = LowerCube.SelectFlatSide("+z")

### TRANSFORMATIONS ###
UpperCube.Resize(2.01)

gap = 0.095
UpperCube.Translate([-1.2, 0.0, gap/2-min([z for z in UpperCube.X[:,2]]) ])
LowerCube.Translate([ 0.0, 0.0,-gap/2-max([z for z in LowerCube.X[:,2]]) ])


### BOUNDARY CONDITIONS ###  [body, nodes, type, directions, values, times(*)]
# cond_bd1 = [LowerCube , LowerCube_bottom , "dirichlet", "xyz"  , [0.0, 0.0, 0.0] ]
cond_bd1 = [LowerCube , LowerCube.SelectAll() , "dirichlet", "xyz"  , [0.0, 0.0, 0.0] ]
# cond_bd2 = [UpperCube , UpperCube_top         , "dirichlet", "xyz"  , [0.0, 0.0,-3.0] ]
cond_bd2 = [UpperCube , UpperCube_top         , "dirichlet", "xyz"  , [0.0, 0.0,-0.2], [0.0, 0.12] ]
cond_bd3 = [UpperCube , UpperCube_top         , "dirichlet", "xyz"  , [2.4, 0.0, 0.0], [0.12, 1.0] ]
BCs = [cond_bd1, cond_bd2, cond_bd3]
# BCs = [cond_bd1, cond_bd2]

### CONTACTS ###
slave  = [UpperCube , UpperCube_bottom]
master = [LowerCube , LowerCube_top   ]

# contact1 = Contact(UpperCube_ctct, LowerCube_ctct, kn=100, cubicT=-1e-1, C1Edges = False)       # (slave, master) inputs can be surfaces as well
contact1 = Contact(slave, master, kn=1000, cubicT=None, C1Edges = False, OutPatchAllowance=1e-3, maxGN = 0.001, mu=0.5, f0 = 0.05)       # (slave, master) inputs can be surfaces as well
ctcts = [contact1]

### MODEL ###
model = FEModel([LowerCube, UpperCube], ctcts, BCs, UpdateOnIteration=True,)     # Should include BCs
model.plotNow()

#############
## Running ##
#############
model.Solve(TimeSteps=100, recover=False, max_iter=10, Lars=False)
