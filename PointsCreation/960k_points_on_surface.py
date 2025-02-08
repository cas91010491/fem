from PyClasses.FEAssembly import *
from PyClasses.Contacts import *
from PyClasses.FEModel import *

import argparse
import sys, os, pickle
import numpy as np
from pdb import set_trace
from time import time

# Set script's directory as current working directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# Load points from Points_on_surf.ft
def load_points(filename):
    points = []
    with open(filename, 'r') as ftfile:
        for line in ftfile:
            points.append([float(x) for x in line.split()])
    return np.array(points)

points = load_points("Points_on_surf.ft")
import pdb; pdb.set_trace()

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
# model.plotNow()       # Uncomment to see and verify geometry


n_per_side = 400

model_X = model.X.reshape(-1,3)

t0 = time()

points = np.zeros((96*n_per_side*n_per_side,6))

# with open("Points_on_surf.csv", 'w') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
# with open("Points_chacking.csv", 'w') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
    # creating a csv writer object
    # csvwriter = csv.writer(csvfile)


    # Iterate over points generated to get closest patch and surface parameters
patches = model.bodies[0].surf.patches
npatches = len(patches)

for i_p,patch in enumerate(patches):
    # nested loop to iterate over surface coordinates 0<(t1,t2)<1
    print(f"Processing patch {i_p+1} of {npatches}")
    points_in_patch = np.zeros((n_per_side*n_per_side,6))
    for n1,t1 in enumerate(np.linspace(0,1,n_per_side)):
        for n2,t2 in enumerate(np.linspace(0,1,n_per_side)):
            x = patch.Grg0([t1,t2])
            # csvwriter.writerow([x[0],x[1],x[2],patch.iquad,t1,t2])
            points_in_patch[n1*n_per_side+n2]=[x[0],x[1],x[2],patch.iquad,t1,t2]
    points[i_p*n_per_side*n_per_side:(i_p+1)*n_per_side*n_per_side,:]=points_in_patch


# Save points to a .ft file
with open("Points_on_surf.ft", 'w') as ftfile:
    for point in points:
        ftfile.write(f"{point[0]} {point[1]} {point[2]} {int(point[3])} {point[4]} {point[5]}\n")

"""


##################################
# Not doing anything with this yet
##################################
# TODO: capture distribution (type and corresponding parameters to be used in next section)
data = np.loadtxt('ContactData', delimiter=',')
GNs = data[:,3]
GNs_pos = GNs[GNs > 0]
GNs_neg = GNs[GNs < 0]
logGNs = np.log10(np.abs(GNs))
gnmin=min(GNs)
gnmax=max(GNs)

params = johnsonsu.fit(GNs)


#########################
### Points generation ###
#########################
n_points = int(1e5)

# gns = johnsonsu.rvs(params[0], params[1], params[2], params[3], size=n_points)
gns = norm.rvs(0,1e-3,size=n_points)
# set_trace()
gns[::2] = johnsonsu.rvs(params[0], params[1], params[2], params[3], size=int(n_points/2))
# gns = johnsonsu.rvs(params[0], params[1], size=int(n_points/2))

t
with open("Points.csv", 'w') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)

    points = np.zeros((n_points,3))
    for i in range(n_points):
        # choose patch
        # ipatch = int(np.random.rand()*len(ptt.surf.patches))
        ipatch = 50
        patch = ptt.surf.patches[ipatch]
        # choose surface coordinate within patch
        xi = np.array([np.random.uniform(0.,1.),np.random.uniform(0,1.)])
        # choose a normal gap for the point on the surface
        # gn = np.random.normal(loc=0,scale=1e-3)/2.
        gn = gns[i]

        results = patch.get_Dgn_for_rigids(gn,xi)
        # compute the randomly produced point
        if i%1000==0:
            print(i)
        points[i]=results[:3]
        csvwriter.writerow(results)

# fig = plt.figure(tight_layout=True)
# ax = fig.add_subplot(111)
# ax.hist(logGNs,bins=1000,color='b')
# # ax.hist(gns,range=(-0.005,0.005),bins=200,color='r')
# plt.show()


########################################
### Plot of results (Might take forever)
########################################
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
model.plot(ax)
ax.scatter(points[:,0],points[:,1],points[:,2],s=0.1,c='k')
plt.show()"""