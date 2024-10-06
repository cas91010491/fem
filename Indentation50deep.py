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

# CUBE
# mesh   = meshio.read("Meshes/Block50x50x24.msh")
mesh   = meshio.read("Meshes/Block24x24x8.msh")
X     = mesh.points
hexas = mesh.cells_dict['hexahedron']

# INDENTOR
meshIdt = meshio.read("Meshes/Indenter/Indentor3.msh")
XS = meshIdt.points
hexasS = meshIdt.cells_dict['hexahedron']


######################
### BUILDING MODEL ###
######################

recOut = False          # if this was the last script run, This can be set to True to skip outers computations
### BODIES ###
bed = FEAssembly(X,hexas, name= "BED",recOuters=recOut)
bed.Youngsmodulus = 0.5
idt = FEAssembly(XS,hexasS, name = "INDENTOR",recOuters=recOut)
bed.isRigid = True     # faster solving when True


### Transformations ###
idt.Rotate(pi/8,"z")                            # 22.5° to avoid double master ambiguity
idt.Resize(2.0)
gap = 0.02                                       # Initial gap between ball and bed. ( \!/ doesn't consider smoothing )
idtUp = max(bed.X[:,2])-min(idt.X[:,2])+gap
idt.Translate([-4.0,0.0, idtUp])

### Selections (nodes) ###
bed_bottom  = bed.SelectFlatSide("-z")
bed_top     = bed.SelectFlatSide("+z")
bed_top_few = bed.SelectNodesBySphere([0.0, 0.0, 2.5],3.0, OnSurface=True)
idt_all     = idt.SelectAll()
idt_top     = idt.SelectFlatSide("+z")
top_idt     = max(idt.X[:,2])
idt_top_few = idt.SelectNodesBySphere([-4.0, 0.0, top_idt],2*0.95, OnSurface=True)

### BOUNDARY CONDITIONS ###  [body, nodes, type, directions, values, times(*)]
cond_bd1 = [bed, bed.SelectAll(), "dirichlet", "xyz"  , [0.0, 0.0, 0.0] ]
# cond_bd2 = [idt, idt_all   , "dirichlet", "xyz"  , [0.0, 0.0,-4.0] ]
cond_bd2 = [idt, idt_top   , "dirichlet", "xyz"  , [0.0, 0.0,-0.52], [0,0.1] ]
cond_bd3 = [idt, idt_top   , "dirichlet", "xyz"  , [3.0, 9.0, 0.0], [0.1,1.0] ]

BCs = [cond_bd1, cond_bd2,cond_bd3]

### CONTACTS ###            # [body, nodes]
master   = [bed , bed_top]
masterNodes = list(set(idt.surf.nodes)-set(idt_top_few))
slave  = [idt, masterNodes ]

contact1 = Contact(slave, master, kn=0.01, cubicT=None, C1Edges = False, OutPatchAllowance=1e-6, maxGN = 0.05,f0=0.1)       # (slave, master) inputs can be surfaces as well

### MODEL ###
model = FEModel([bed, idt], [contact1], BCs)           # [bodies, contacts, BCs, opts*]
# model.plotNow()       # Uncomment to see and verify geometry

#############
## Running ##
#############
model.Solve(TimeSteps=200, recover=True, max_iter=20)
