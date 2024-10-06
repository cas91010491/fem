from PyClasses.FEAssembly import *
from PyClasses.Contacts import *
from PyClasses.FEModel import *

import meshio, sys, os, pickle
import numpy as np
import matplotlib.pyplot as plt
from pdb import set_trace
from math import pi
from scipy import sparse

from time import time


####################
### GETTING DATA ###
####################
os.chdir(sys.path[0])

# BLOCK
mesh_blk   = meshio.read("../Meshes/Block_pseudo2d.msh")
X_blk     = mesh_blk.points
hexas_blk = mesh_blk.cells_dict['hexahedron']
blk = FEAssembly(X_blk,hexas_blk, name= "BLOCK",recOuters=False)
blk.Youngsmodulus = 0.05
blk.Translate([-3,0.0,2.6])
blk.Resize(1.01)


# # 2d-Base
mesh_base   = meshio.read("../Meshes/Base_pseudo2d.msh")
# mesh_ptt = meshio.read("Meshes/Cubes/cube843.msh")
X_base = mesh_base.points
hexas_base = mesh_base.cells_dict['hexahedron']
base = FEAssembly(X_base,hexas_base, name = "BASE")
# ptt.Resize(4.0,dir='x')
base.Resize(2.0,dir='y')
base.Resize(1.1,dir='x')
# ptt.Resize(2.0,dir='z')
# ptt.RandDistort(0.5)
pickle.dump([base],open("Base2dAssembly.dat","wb"))
base.isRigid = True     # faster solving when True

######################
### BUILDING MODEL ###
######################

### Selections (nodes) ###
blk_bottom  = blk.SelectFlatSide("-z")
blk_top     = blk.SelectFlatSide("+z")
# ptt_highernodes = base.SelectHigherThan("z", val = 0.5, Strict = True,OnSurface = True)
ptt_highernodes = base.SelectFlatSide("+z")
print(ptt_highernodes)
lead_face = blk.SelectFlatSide("-y")

slave_nodes = list(set(blk_bottom).intersection(set(lead_face)))
# slave_nodes = blk_bottom

blk_top = list(set(blk_top).intersection(set(lead_face)))

same_pairs = np.zeros((len(lead_face),2),dtype=int)
for ii,node in enumerate(lead_face):
    x,y,z = blk.X[node]
    node_mirror = np.where((abs(blk.X[:,0]-x)<1e-10) & (abs(blk.X[:,2]-z)<1e-10) & (blk.X[:,1]!=y))[0][0]
    same_pairs[ii,:] = [node,node_mirror]


base.X[21,2] += 0.5 
base.X[39,2] += 0.5 

# r = np.zeros(2*(len(blk.X)))
# c = np.zeros(2*(len(blk.X)))
# v = np.ones(2*(len(blk.X)))
rs = np.zeros(2*(len(blk.X)))
cs = np.zeros(2*(len(blk.X)))
vs = np.ones((2*len(blk.X)))
rt = np.zeros((len(blk.X)))
ct = np.zeros((len(blk.X)))
vt = np.ones((len(blk.X)))

ndofs = 3*(len(X_blk)+len(base.X))

for i_p,pair in enumerate(same_pairs):
    node_1 = pair[0]
    node_2 = pair[1]
    # r[[4*i_p,4*i_p+1,4*i_p+2,4*i_p+3]]=[3*node_1,3*node_1+2,3*node_2,3*node_2+2]
    # c[[4*i_p,4*i_p+1,4*i_p+2,4*i_p+3]]=[2*i_p,2*i_p+1,2*i_p,2*i_p+1]
    rt[[2*i_p,2*i_p+1]]=[3*node_1,3*node_1+2]
    ct[[2*i_p,2*i_p+1]]=[2*i_p,2*i_p+1]
    # rs[[2*i_p,2*i_p+1]]=[3*node_2,3*node_2+2]
    # cs[[2*i_p,2*i_p+1]]=[3*node_1,3*node_1+2]
    rs[[4*i_p,4*i_p+1,4*i_p+2,4*i_p+3]]=[3*node_1,3*node_1+2,3*node_2,3*node_2+2]
    cs[[4*i_p,4*i_p+1,4*i_p+2,4*i_p+3]]=[3*node_1,3*node_1+2,3*node_1,3*node_1+2]


Ns = sparse.csr_matrix((vs,(rs,cs)),shape=(ndofs,ndofs))
Nt = sparse.csr_matrix((vt,(rt,ct)),shape=(ndofs,len(blk.X)))
N = [Ns,Nt]


### BOUNDARY CONDITIONS ###  [body, nodes, type, directions, values, times(*)]
cond_bd1 = [base, base.SelectAll(), "dirichlet", "xyz"  , [0.0, 0.0, 0.0] ]     # Base static
cond_bd2 = [blk, blk_top   , "dirichlet", "xz"  , [0.0, 0.0,-0.25], [0.0, 0.2] ]             # Displacement
cond_bd3 = [blk, blk_top   , "dirichlet", "xz"  , [6.0, 0.0, 0.0], [0.2, 1.0] ]             # Displacement
cond_bd4 = [blk, blk.SelectAll()  , "dirichlet", "y"  , [0.0, 0.0,0.0] ]             # Symmetry


BCs = [cond_bd1, cond_bd2,cond_bd3,cond_bd4]

### CONTACTS ###            # [body, nodes]
slave   = [blk , slave_nodes]
master = [base,ptt_highernodes]


contact1 = Contact(slave, master, kn=5, C1Edges = True, maxGN = 0.001,f0=0.1)       # (slave, master) inputs can be surfaces as well

### MODEL ###
model = FEModel([blk, base], [contact1], BCs,transform_2d=N)           # [bodies, contacts, BCs, opts*]

# base.surf.ComputeGrgPatches(np.zeros(ndofs),range(len(base.surf.nodes)))
base.surf.ComputeGrgPatches(np.zeros(ndofs),ptt_highernodes,exactNodesGiven=True)
# contact1.masterSurf.ComputeGrgPatches(np.zeros(ndofs),[])
# model.plotNow(as2D=True,OnlyMasterSurf=True)       # Uncomment to see and verify geometry

#############
## Running ##
#############

# import cProfile
# import pstats
# import io
# pr = cProfile.Profile()
# pr.enable()

t0 = time()

model.Solve(TimeSteps=100, recover=True, ForcedShift=False,max_iter=100)

print("this took",time()-t0,"seconds to compute")


# pr.disable()
# s = io.StringIO()
# ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
# ps.print_stats()

# print(s.getvalue())


