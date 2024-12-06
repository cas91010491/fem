
import meshio, sys, os, pickle
import numpy as np
import matplotlib.pyplot as plt
from pdb import set_trace
from math import pi

sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from PyClasses.FEAssembly import *
from PyClasses.Contacts import *
from PyClasses.FEModel import *

from time import time
import argparse

from tensorflow.keras.models import load_model
from tensorflow.keras.losses import MeanSquaredError
from tensorflow.keras.models import load_model


# Load the neural network model
model_path = 'model_epoch_850.h5'
MT_model = load_model(model_path, custom_objects={'mse': MeanSquaredError()})

# Argument parsing
parser = argparse.ArgumentParser(description='Process a data for the 2d contact model.')
parser.add_argument('--min_method', type=str, required=True, help='minimization method: BFGS, LBFGSNN')
parser.add_argument('--mesh', type=int, required=True, help='choose mesh 5, 10 or 15')
parser.add_argument('--plastic', type=int, required=True, help='boolean for plastic')
parser.add_argument('--ann', type=int, required=True, help='boolean for ANN')

args = parser.parse_args()

# Calculate subspace bounds
minimization_method = args.min_method
mesh = args.mesh
plastic = args.plastic
useANN = args.ann

# minimization_method = "BFGS"
# mesh = 10
# plastic = 0


####################
### GETTING DATA ###
####################
os.chdir(sys.path[0])

# BLOCK
mesh_blk   = meshio.read("../Meshes/Cubes/cube"+str(mesh)+"x"+str(mesh)+"x"+str(mesh)+".msh")

X_blk     = mesh_blk.points
hexas_blk = mesh_blk.cells_dict['hexahedron']
if plastic:
    blk = FEAssembly(X_blk,hexas_blk, name= "BLOCK",recOuters=False,plastic_param=[0.01,0.05,1.0])
else:
    blk = FEAssembly(X_blk,hexas_blk, name= "BLOCK",recOuters=False)
blk.Youngsmodulus = 0.05
blk.Translate([0.0,0.0,3.5])

# # POTATO
[ptt] = pickle.load(open("PotatoAssembly.dat","rb"))
## OR ##
# mesh_ptt = meshio.read("../Meshes/QuadSpheres/QuadSphere4.msh")
# # mesh_ptt = meshio.read("Meshes/Cubes/cube843.msh")
# X_ptt = mesh_ptt.points
# hexas_ptt = mesh_ptt.cells_dict['hexahedron']
# ptt = FEAssembly(X_ptt,hexas_ptt, name = "POTATO",recOuters=recOut)
# ptt.Resize(4.0,dir='x')
# ptt.Resize(3.0,dir='y')
# ptt.Resize(2.0,dir='z')
# ptt.RandDistort(0.5)
# pickle.dump([ptt],open("PotatoAssembly.dat","wb"))
ptt.isRigid = True     # faster solving when True

######################
### BUILDING MODEL ###
######################

### Selections (nodes) ###
blk_bottom  = blk.SelectFlatSide("-z")
blk_top     = blk.SelectFlatSide("+z")
ptt_highernodes = ptt.SelectHigherThan("z", val = 0.5, Strict = True,OnSurface = True)


### BOUNDARY CONDITIONS ###  [body, nodes, type, directions, values, times(*)]
cond_bd1 = [ptt, ptt.SelectAll(), "dirichlet", "xyz"  , [0.0, 0.0, 0.0] ]
cond_bd2 = [blk, blk_top   , "dirichlet", "xyz"  , [12.0, 0.0,0.0] ]


BCs = [cond_bd1, cond_bd2]

### CONTACTS ###            # [body, nodes]
slave   = [blk , blk_bottom]
master = [ptt,ptt_highernodes]

contact1 = Contact(slave, master, kn=1e2, C1Edges = False, maxGN = 0.001,ANNmodel=MT_model if useANN else None)       # (slave, master) inputs can be surfaces as well

### MODEL ###
subname = "_"+("plastic" if plastic else "elastic")+"_"+minimization_method+"_"+str(mesh)+("_ANN" if useANN else "")
model = FEModel([blk, ptt], [contact1], BCs, subname=subname)           # [bodies, contacts, BCs, opts*]

ndofs = 3*(len(X_blk)+len(ptt.X))
ptt.surf.ComputeGrgPatches(np.zeros(ndofs),range(len(ptt.surf.nodes)))
# model.plotNow()       # Uncomment to see and verify geometry

#############
## Running ##
#############

# import cProfile
# import pstats
# import io
# pr = cProfile.Profile()
# pr.enable()

t0 = time()


recov = "OUTPUT_202412060907ContactPotato_Ex1_elastic_BFGS_10/"+"RecoveryData.dat"
model.Solve(TimeSteps=100,max_iter=20, recover=False ,minimethod=minimization_method,plot=0)
# model.Solve(TimeSteps=100,max_iter=20, recover=False ,minimethod=minimization_method,plot=1)


print("this took",time()-t0,"seconds to compute")


# pr.disable()
# s = io.StringIO()
# ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
# ps.print_stats()

# print(s.getvalue())


