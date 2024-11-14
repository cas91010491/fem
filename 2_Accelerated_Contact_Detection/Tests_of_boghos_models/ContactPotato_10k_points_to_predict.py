from PyClasses.FEAssembly import *
from PyClasses.Contacts import *
from PyClasses.FEModel import *

import meshio, sys, os, pickle
import numpy as np
import matplotlib.pyplot as plt
# from mpl_toolkits.mplot3d import Axes3D
from math import pi

from time import time


# import cProfile
# import pstats
# import io
# pr = cProfile.Profile()
# pr.enable()



####################
### GETTING DATA ###
####################
os.chdir(sys.path[0])

# # POTATO
[ptt] = pickle.load(open("PotatoAssembly.dat","rb"))
# ptt.Translate([-6.0, 0.0, 0.0])
nf = len(ptt.X)
ptt.DoFs = np.array( [ [ 3*n + i for i in range(3) ] for n in range(nf) ] )
ptt.isRigid = True     # faster solving when True

ptt.surf.ComputeGrgPatches(np.zeros(3*len(ptt.X)),range(len(ptt.surf.nodes)))


######################
### BUILDING MODEL ###
######################

# Create 10k points and register their patch, and surface coordinate
n_slaves = int(1e4)
patches_ids = np.random.randint(96,size=(n_slaves))
T1T2 = np.random.rand(n_slaves,2)
GNs = np.random.normal(0.0,0.3,(n_slaves,1))
XS = np.zeros((n_slaves,3))


real_actives = [patches_ids[i] if GNs[i] < 0 else None for i in range(n_slaves)]

# Generate the points in space (x,y,z coords)
for i_x in range(n_slaves):
    i_p = patches_ids[i_x]
    t = T1T2[i_x]
    gn = GNs[i_x]
    patch = ptt.surf.patches[i_p]
    xc = patch.Grg(t)
    nor = patch.D3Grg(t)
    XS[i_x] = xc + gn*nor


# Slave Body made out of only scattered points
slave_body = FEAssembly(XS,[[]],"points_body",outers=[[],list(range(n_slaves))])
slave_body.X = XS

u = np.zeros((3*nf+3*n_slaves,),dtype=np.float64)

# Contact object
contact_old = Contact([slave_body,list(range(n_slaves))],ptt)
contact_new = Contact([slave_body,list(range(n_slaves))],ptt)

# Model
model = FEModel([ptt,slave_body],[contact_old,contact_new],[])



print("######\n# Getting Candidates and Actives...\n######")
# Predict using old method
t0 = time()
contact_old.getCandidates(u,CheckActive=True)
print("old takes:",time()-t0)

# Predict using new method
t0 = time()
contact_new.getCandidatesANN(u,CheckActive=True)
print("new takes:",time()-t0)
print("")



print("#######\n# Getting Projections...\n######")

total_nodes_in = np.sum(GNs < 0)

# Predict using old method
t0 = time()
correct_patch_node_in = 0
T1T2_old = -1*np.ones((n_slaves,2))
for idx in range(n_slaves):
    if GNs[idx]<0:
        if contact_old.actives[idx]==patches_ids[idx]:
            correct_patch_node_in += 1
        patch_obj = ptt.surf.patches[patches_ids[idx]]
        T1T2_old[idx] = patch_obj.findProjection(XS[idx])


print("old takes:",time()-t0)
print(f"Actives correctly predicted: ({correct_patch_node_in}/{total_nodes_in})")

# Predict using new method
t0 = time()
correct_patch_node_in = 0
T1T2_new = -1*np.ones((n_slaves,2))
# correct_in_points = [idx  for idx in range(n_slaves) if (GNs[idx]<0 and (contact_new.actives[idx]==patches_ids[idx]))]
all_patches_obj = contact_new.masterSurf.patches
for i_patch, patch in enumerate(all_patches_obj):
    nodes_in_patch = contact_new.actives_sparse[:,i_patch].nonzero()[0]
    if len(nodes_in_patch)>0:
        predictions_for_t1t2 = patch.MinDistANN( XS[nodes_in_patch] , verbose=0 )
        T1T2_new[nodes_in_patch] = predictions_for_t1t2


print("new takes:",time()-t0)
print(f"Actives correctly predicted: ({correct_patch_node_in}/{total_nodes_in})")


# pr.disable()
# s = io.StringIO()
# ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
# ps.print_stats()

# print(s.getvalue())



import matplotlib.pyplot as plt
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
contact_old.plotContact(u,ax)
ax.set_aspect('equal')
plt.show()




set_trace()
