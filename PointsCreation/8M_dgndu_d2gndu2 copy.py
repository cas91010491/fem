from PyClasses.FEAssembly import *
from PyClasses.Contacts import *
from PyClasses.FEModel import *

import argparse
import sys, os, pickle
import numpy as np
from pdb import set_trace
from time import time


"""
This code generates a set of points close to the surface of the loaded geometry.
The potato-shaped geometry contains 96 patches on its surface and, for each patch,
a number of points is generated which are distributed on the surface parameters and 
a list of normal gaps.
"""

##############################
### Data Generation Inputs ###
##############################

points_per_side = 20
xi_margin = 1e-8


xis = np.linspace(xi_margin,1-xi_margin,points_per_side)
gns = np.concatenate((np.logspace(-1, -8, num=8, base=10), [0], np.logspace(-8, -1, num=8, base=10)))


filename = "NewData_"+str(points_per_side)+"x"+str(points_per_side)+"x"+str(len(gns))+".csv"  # Replace with your actual file path
headers = ['x', 'y', 'z', 'p_id', 'xi1', 'xi2', 'gn', 'x_surf', 'y_surf','z_surf',
                'dgn1','dgn2','dgn3','ddgn11','ddgn12','ddgn13','ddgn22','ddgn23','ddgn33']  # Assign column names



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



import pandas as pd
import numpy as np
from scipy.spatial import KDTree

# Load the CSV data without headers
ft_file = "Points_on_surf.ft"  # Replace with your actual file path
data = pd.read_csv(ft_file, delim_whitespace=True, header=None)  # No headers in the file
data.columns = ['x', 'y', 'z', 'p_id', 'xi1', 'xi2']  # Assign column names
# set_trace()



# Extract surface points and associated information
surface_points = data[['x', 'y', 'z']].values  # (960000, 3)
metadata = data.values  # Full table with all 7 columns

# Build the k-d tree
tree = KDTree(surface_points)

# Function to query the k closest points and make p_id unique
def find_closest_points_unique(query_points, k=10):
    # Query the k-d tree
    distances, indices = tree.query(query_points, k=k)
    
    results = []
    for i, idx_set in enumerate(indices):  # Loop over each query point
        closest_metadata = metadata[idx_set]  # Get metadata for the closest points
        closest_distances = distances[i]  # Get distances for these points
        
        # Create a DataFrame to manage the points
        # df_closest = pd.DataFrame(closest_metadata, columns=['x', 'y', 'z', 'p_id', 'xi1', 'xi2', 'gn'])
        df_closest = pd.DataFrame(closest_metadata, columns=['x', 'y', 'z', 'p_id', 'xi1', 'xi2'])
        df_closest['distance'] = closest_distances
        
        # Keep the closest point for each unique p_id
        df_unique = df_closest.sort_values('distance').drop_duplicates(subset='p_id', keep='first')
        results.append(df_unique.drop(columns=['distance']).values)  # Drop distance column for the final result
    
    return results








t0 = time()

# for each patch, generate points_per_sidexpoints_per_sidexlen(gns) points
xs = np.zeros((len(ptt.surf.patches)*points_per_side**2*len(gns),3))
idx = 0
for ip,patch in enumerate(ptt.surf.patches):
    for xi1 in xis:
        for xi2 in xis:
            for gn in gns:
                xs[idx] = patch.Grg0([xi1,xi2]) + gn*patch.D3Grg([xi1,xi2])
                
                idx += 1



# xs = np.array([[2.35325983, 1.79146486, 1.22203662]])

CANDIDATES = find_closest_points_unique(xs)

# set_trace()

# with open("Points_"+str(i)+'_'+str(j)+'_'+str(k)+".csv", 'w') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
# with open("NewTrainingData.csv", 'w') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
with open("NewTrainingData.csv", 'a+') as csvfile:        #'a' is for "append". If the file doesn't exists, cretes a new one
    # creating a csv writer object
    csvwriter = csv.writer(csvfile)
    
    # Write the headers
    csvwriter.writerow(headers)

    # Iterate over points generated to get closest patch and surface parameters
    patches = model.bodies[0].surf.patches
    npatches = len(patches)
    raw_distances = np.zeros(npatches)

    xm_matrix = np.zeros((npatches,3))
    for ip,patch in enumerate(patches):
        xm_matrix[ip] = patch.BS.x
    distance_matrix = np.linalg.norm(xm_matrix[:, np.newaxis, :] - xs, axis=2)

    for i,xsi in enumerate(xs):
        # for ip,patch in enumerate(patches):
        #     raw_distances[ip] = distance_matrix[ip,i]
        # raw_distances = distance_matrix[:,i]

        # retrieve ip, t1, t2, gn from index i and store as '_0' for initial values
        ip = int(i/(points_per_side**2*len(gns)))
        t1_idx = (i%(points_per_side**2*len(gns)))//(points_per_side*len(gns))
        t1_0 = xis[t1_idx]
        t2_idx = (i%(points_per_side**2*len(gns)))%points_per_side
        t2_0 = xis[t2_idx]
        gn_idx = i%len(gns)
        gn_0 = gns[gn_idx]


        if i== 123:
            set_trace()


        rec_lvl = 15 # seeding for recursive MinDist search (u_i = np.linspace(x0,x1,seeding+1))
        # # get candidates of closest patch
        # mindist = min(raw_distances)
        # maxdist = max(raw_distances)
        # # patch_size = 2.0        # approximate average patch size (observed from BS radius ~1)
        # # size_factor = 0.03      # offset percentage for search radius where search_radius = mindist + size_factor*(maxdist-mindist)
        # patch_size = 2.5        # approximate average patch size (observed from BS radius ~1)
        # size_factor = 0.15      # offset percentage for search radius where search_radius = mindist + size_factor*(maxdist-mindist)
        candidates = np.array(CANDIDATES[i][:,3],dtype=int)

        if i%100==0:
            print("pt_nr:",i,"\txsi=",xsi, "\ttotal_time:",time()-t0)  # to visualize evolution of computation

        t1t2 = np.zeros((len(candidates),2))
        gns = np.zeros(len(candidates))         # global normal vectors
        for ip,p_idx in enumerate(candidates):
            patch = patches[p_idx]
            # if ip==1:
            #     set_trace()
            t1t2[ip] = patch.findProjection(xsi,recursive=rec_lvl)
            xc = patch.Grg0(t1t2[ip])
            nor = patch.D3Grg(t1t2[ip])
            gns[ip] = (xsi-xc)@nor

        cand_idx =  np.argmin(abs(gns))                   # index of the candidate with smallest angle between x
        pid = candidates[cand_idx]                  # index of closest patch
        t1,t2 = t1t2[cand_idx]

        # Unprojected Case!
        trials = 0
        while not((0<=t1<=1) and (0<=t2<=1)):               # loop until a valid projection is found
            trials += 1
            # size_factor *= 1.2      # offset percentage for search radius where search_radius = mindist + size_factor*(maxdist-mindist)
            # candidates = [patch for patch, dist in zip(patches, raw_distances) if dist < mindist+size_factor*(maxdist-mindist)]
            # rec_lvl += 1
            rec_lvl *= 2
            for ip,p_idx in enumerate(candidates):
                patch = patches[p_idx]
                t1t2[ip] = patch.findProjection(xsi,recursive=rec_lvl)
                xc = patch.Grg0(t1t2[ip])
                nor = patch.D3Grg(t1t2[ip])
                inside = ((0<=t1t2[ip][0]<=1) and (0<=t1t2[ip][1]<=1))
                gns[ip] = (xsi-xc)@nor if inside else 1000

            cand_idx =  np.argmin(abs(gns))                   # index of the candidate with smallest angle between x
            pid = candidates[cand_idx]                  # index of closest patch
            t1,t2 = t1t2[cand_idx]

            if trials>10:
                print("more than 10 trials while looking for projections for point",xsi,".")
                print("last candidates: \n",[cand for cand in candidates])
                print("p_id:",pid,"\tgn:",gns[cand_idx])
                pass

        xc, dxcdt, d2xcd2t = patch.Grg([t1,t2], deriv=2)
        D1p, D2p = dxcdt.T
        D3p = np.cross(D1p,D2p)
        nor = D3p/norm(D3p)
        
        dfdt =  -2*np.tensordot(xsi-xc,d2xcd2t,axes=1) + 2*(dxcdt.T @ dxcdt)
        dfdxs = -2*dxcdt.T
        dtdxs = np.linalg.solve(-dfdt,dfdxs)

        # Vars for K---------------------------------------------------
        dndt = patch.dndt([t1,t2])
        dndxs = dndt@dtdxs

        csvwriter.writerow([xsi[0],xsi[1],xsi[2],pid,t1,t2,gns[cand_idx],
                            xc[0],xc[1],xc[2],
                            nor[0],nor[1],nor[2],
                            dndxs[0,0],dndxs[0,1],dndxs[0,2],
                            dndxs[1,1],dndxs[1,2],
                            dndxs[2,2]])

# pr.disable()
# s = io.StringIO()
# ps = pstats.Stats(pr, stream=s).sort_stats('cumulative')
# ps.print_stats()

# print(s.getvalue())



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