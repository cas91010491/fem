

import os
import glob
import sys
import numpy as np
import pandas as pd 
import scipy as sp 
from scipy.spatial import cKDTree
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.interpolate import griddata

# Set the main directory as the current working directory
main_directory = '/home/diego2/PointsCreation'
os.chdir(main_directory)

pnts_fold = 'csv_files'
csvs = glob.glob(pnts_fold +  "/*.csv")
df_list = [ pd.read_csv(csvf,header=None) for csvf in csvs ]
df = pd.concat(df_list, ignore_index=True)
df.columns = ['x', 'y', 'z', 'p_id', 'xi_1', 'xi_2', 'gn']
df = df.sort_values(by=['x', 'y', 'z']).reset_index(drop=True)
pnts3D = df[['x','y','z']].values

# sampling coordinates of slices 
xmap = np.sort(np.unique(np.round(pnts3D[:,0],6)))
ymap = np.sort(np.unique(np.round(pnts3D[:,1],6)))
zmap = np.sort(np.unique(np.round(pnts3D[:,2],6)))

# computing the sampling distances in each direction : 
Dx = np.mean(np.unique(np.diff(xmap)))
Dy = np.mean(np.unique(np.diff(ymap)))
Dz = np.mean(np.unique(np.diff(zmap)))

grid_x, grid_y, grid_z = np.meshgrid(xmap, ymap, zmap, indexing='ij')
gngrid = df.gn.values.reshape(len(grid_x),len(grid_y),len(grid_z))
idxgrid = np.linspace(0,df.shape[0]-1,df.shape[0], dtype=np.int64).reshape(len(grid_x),len(grid_y),len(grid_z))


# Assuming grid_x, grid_y, grid_z are already defined
# Create 3D indices
i_indices, j_indices, k_indices = np.meshgrid(
    np.arange(len(grid_x)), 
    np.arange(len(grid_y)), 
    np.arange(len(grid_z)), 
    indexing='ij'
)

# Flatten the indices into a list of tuples
listidx = np.array([i_indices.flatten(), j_indices.flatten(), k_indices.flatten()]).T


lists = np.linspace(0, len(df)-1, len(df), dtype=np.int64).reshape(len(grid_x),len(grid_y),len(grid_z))


# def estimate_direction(plus, minus, delta,approach):
#         if approach=='midpoint':
#             derivative_approx = (plus - minus) / (2 * delta)
#         if approach=='ext':
#             derivative_approx = (plus - minus) / (delta)
#         return derivative_approx

# def finite_difference_gradient(idx):
#     def compute_gradient(axis, coord_grid, field_grid, i, j, k):
#         if axis == 'x':
#             dim, plus_idx, minus_idx = i, (i + 1, j, k), (i - 1, j, k)
#         elif axis == 'y':
#             dim, plus_idx, minus_idx = j, (i, j + 1, k), (i, j - 1, k)
#         elif axis == 'z':
#             dim, plus_idx, minus_idx = k, (i, j, k + 1), (i, j, k - 1)

#         axis_idx = {"x":0, "y":0, "z":0}

#         # Handle boundary cases
#         if dim == 0:  # Forward difference
#             value_plus = field_grid[plus_idx]
#             value_minus = field_grid[i, j, k]
#             coord_plus = coord_grid[plus_idx]
#             coord_minus = coord_grid[i, j, k]
#             gradient = estimate_direction(value_plus, value_minus, coord_plus - coord_minus, 'ext')
#         elif dim == field_grid.shape[axis_idx[axis]] - 1:  # Backward difference
#             value_plus = field_grid[i, j, k]
#             value_minus = field_grid[minus_idx]
#             coord_plus = coord_grid[i, j, k]
#             coord_minus = coord_grid[minus_idx]
#             gradient = estimate_direction(value_plus, value_minus, coord_plus - coord_minus, 'ext')
#         else:  # Midpoint (central difference)
#             value_plus = field_grid[plus_idx]
#             value_minus = field_grid[minus_idx]
#             coord_plus = coord_grid[plus_idx]
#             coord_minus = coord_grid[minus_idx]
#             gradient = estimate_direction(value_plus, value_minus, 0.5 * (coord_plus - coord_minus), 'midpoint')
#         return gradient

#     # Unpack indices
#     i, j, k = idx

#     # Compute gradients in each direction
#     grad_x = compute_gradient('x', grid_x, gngrid, i, j, k)
#     grad_y = compute_gradient('y', grid_y, gngrid, i, j, k)
#     grad_z = compute_gradient('z', grid_z, gngrid, i, j, k)

#     # Compute the gradient magnitude
#     grad_mag = np.sqrt(grad_x**2 + grad_y**2 + grad_z**2)

#     # Normalize the gradient to have unit length
#     grad_x_normalized = grad_x / grad_mag
#     grad_y_normalized = grad_y / grad_mag
#     grad_z_normalized = grad_z / grad_mag

#     return grad_x, grad_y, grad_z, grad_mag



# grads = np.array([finite_difference_gradient(idx) for idx in listidx])



# df['dgndx'],df['dgndy'],df['dgndz'], df['ndgndl'] = grads[:,0],grads[:,1],grads[:,2],grads[:,3]


# dgndxgrid = df.dgndx.values.reshape(len(grid_x),len(grid_y),len(grid_z))
# dgndygrid = df.dgndy.values.reshape(len(grid_x),len(grid_y),len(grid_z))
# dgndzgrid = df.dgndz.values.reshape(len(grid_x),len(grid_y),len(grid_z))
# ndgndlgrid = df.ndgndl.values.reshape(len(grid_x),len(grid_y),len(grid_z))



# extent = [df.y.min(), df.y.max(), df.z.min(), df.z.max()]

# plt.imshow(gngrid[10,:,:].T, origin='lower', cmap='viridis', extent=extent)
# plt.colorbar()



# vmin=np.min(gngrid)
# vmax=np.max(gngrid)
# extent = [df.y.min(), df.y.max(), df.z.min(), df.z.max()]
# for i in range(len(xmap)):
#     plt.imshow(gngrid[i,:,:].T, origin='lower', cmap='viridis', extent=extent)#, vmin=vmin, vmax=vmax)
#     plt.colorbar()
#     plt.savefig('gn_cutX/' + 'cut_' + str(i) + '.png')
#     plt.close()


# vmin=np.min(ndgndlgrid)
# vmax=np.max(ndgndlgrid)
# for i in range(len(xmap)):
#     plt.imshow(ndgndlgrid[i,:,:].T, origin='lower', cmap='viridis', extent=extent)#, vmin=vmin, vmax=vmax)
#     plt.colorbar()
#     plt.savefig('ndgndl_cutX/' + 'cut_' + str(i) + '.png')
#     plt.close()

# vmin=np.min(dgndxgrid)
# vmax=np.max(dgndxgrid)
# for i in range(len(xmap)):
#     plt.imshow(dgndxgrid[i,:,:].T, origin='lower', cmap='viridis', extent=extent)#, vmin=vmin, vmax=vmax)
#     plt.colorbar()
#     plt.savefig('dgdx_cutX/' + 'cut_' + str(i) + '.png')
#     plt.close()


# vmin=np.min(dgndygrid)
# vmax=np.max(dgndygrid)
# for i in range(len(xmap)):
#     plt.imshow(dgndygrid[i,:,:].T, origin='lower', cmap='viridis', extent=extent)#, vmin=vmin, vmax=vmax)
#     plt.colorbar()
#     plt.savefig('dgdy_cutX/' + 'cut_' + str(i) + '.png')
#     plt.close()


# min=np.min(dgndzgrid)
# vmax=np.max(dgndzgrid)
# for i in range(len(xmap)):
#     plt.imshow(dgndzgrid[i,:,:].T, origin='lower', cmap='viridis', extent=extent)#, vmin=vmin, vmax=vmax)
#     plt.colorbar()
#     plt.savefig('dgdz_cutX/' + 'cut_' + str(i) + '.png')
#     plt.close()


fstd = 2
pbs = []
for i in range(2,gngrid.shape[0]-2):
    for j in range(2,gngrid.shape[1]-2):
        for k in range(2,gngrid.shape[2]-2):
            gn_family= np.array([
                gngrid[i,j,k],
                gngrid[i-1,j,k], gngrid[i-2,j,k], gngrid[i+1,j,k], gngrid[i+2,j,k],
                gngrid[i,j-1,k], gngrid[i,j-2,k], gngrid[i,j+1,k], gngrid[i,j+2,k],
                gngrid[i,j,k-1], gngrid[i,j,k-2], gngrid[i,j,k+1], gngrid[i,j,k+2]])
            cond = np.logical_and(gn_family[0] < (np.mean(gn_family) + fstd*np.std(gn_family)),
                                  gn_family[0] > (np.mean(gn_family) - fstd*np.std(gn_family)))
            if not cond :
                print("prorblem" + str([i,j,k]))
                prob = (i,j,k)
                gnprob = gn_family[0]
                gnnew = np.mean(gn_family)
                pbs.append([prob,gnprob,gnnew])


# dfc[(dfc.x == grid_x[pbs[10][0]]) & (dfc.y == grid_y[pbs[10][0]]) & (dfc.z == grid_z[pbs[10][0]])][['x','y','z']]


dfc = df 
for i in range(len(pbs)):
    inds = pbs[i][0]
    gnn = pbs[i][2]
    dfc.loc[(dfc.x == grid_x[inds]) & (dfc.y == grid_y[inds]) & (dfc.z == grid_z[inds]),'gn'] = gnn


bound_d =5.0;
# df_bound = dfc[(dfc.gn>-1*bound_d)&(dfc.gn<bound_d)].reset_index(drop=True)
df_bound = dfc[(dfc.gn>-0.5)&(dfc.gn<1.0)].reset_index(drop=True)
pnts3D_bound = df_bound[['x','y','z']].values

if not os.path.exists('Zcuts5'):
    os.makedirs('Zcuts5')

for i in range(len(zmap)):
    cut_pnts = df_bound[df_bound.z==zmap[i]]
    plt.figure(figsize=(10,10))
    plt.scatter(cut_pnts[cut_pnts.gn<0].x, cut_pnts[cut_pnts.gn<0].y, color='red',s=0.1)
    plt.scatter(cut_pnts[cut_pnts.gn>0].x, cut_pnts[cut_pnts.gn>0].y, color='blue', s=0.1)
    plt.xlim(np.min(xmap), np.max(xmap))
    plt.ylim(np.min(ymap), np.max(ymap))
    plt.savefig('Zcuts5/' + 'cut_' + str(i) + '.png')

if not os.path.exists('Ycuts5'):
    os.makedirs('Ycuts5')


for i in range(len(ymap)):
    cut_pnts = df_bound[df_bound.y==ymap[i]]
    plt.figure(figsize=(10,10))
    plt.scatter(cut_pnts[cut_pnts.gn<0].x, cut_pnts[cut_pnts.gn<0].z, color='red',s=0.1)
    plt.scatter(cut_pnts[cut_pnts.gn>0].x, cut_pnts[cut_pnts.gn>0].z, color='blue', s=0.1)
    plt.xlim(np.min(xmap), np.max(xmap))
    plt.ylim(np.min(zmap), np.max(zmap))
    plt.savefig('Ycuts5/' + 'cut_' + str(i) + '.png')

if not os.path.exists('Xcuts5'):
    os.makedirs('Xcuts5')


for i in range(len(zmap)):
    cut_pnts = df_bound[df_bound.x==xmap[i]]
    plt.figure(figsize=(10,10))
    plt.scatter(cut_pnts[cut_pnts.gn<0].y, cut_pnts[cut_pnts.gn<0].z, color='red',s=0.1)
    plt.scatter(cut_pnts[cut_pnts.gn>0].y, cut_pnts[cut_pnts.gn>0].z, color='blue', s=0.1)
    plt.xlim(np.min(ymap), np.max(ymap))
    plt.ylim(np.min(zmap), np.max(zmap))
    plt.savefig('Xcuts5/' + 'cut_' + str(i) + '.png')



diX , idxiX = kd_treeX.query([pnts3D[555555,0],0], k=3)
diX, idxiX = np.array(diX[1::]), idxiX[1::]
diY , idxiY = kd_treeY.query([pnts3D[555555,1],0], k=3)
diY, idxiY = np.array(diY[1::]), idxiY[1::]
diZ , idxiZ = kd_treeZ.query([pnts3D[555555,2],0], k=3)
diZ, idxiZ = np.array(diZ[1::]), idxiZ[1::]



bound_d = 0.5;
df_bound = df[(df.gn>-1*bound_d)&(df.gn<bound_d)].reset_index(drop=True)
pnts3D_bound = df_bound[['x','y','z']].values
kd_tree_bound = cKDTree(pnts3D_bound)



# fig = plt.figure()
# ax = fig.add_subplot(projection='3d')
# ax.scatter(df_bound.x, df_bound.y, df_bound.z,s=0.1)


# plt.figure(figsize=(25,10))
# sns.countplot(x=df_bound.p_id)


# # In[ ]:


# testp = [np.min(pnts3D[:,0])-1.0, np.min(pnts3D[:,1])-1.0, np.min(pnts3D[:,2])-1.0]


# # In[ ]:


# dis, idx = kd_tree.query(testp, k=1)
# dis_bound, idx_bound = kd_tree_bound.query(testp, k=1)


# # In[ ]:


# kd_tree.query(testp, k=20)


# # In[ ]:


# print(df.iloc[idx].p_id)
# print(df_bound.iloc[idx_bound].p_id)


# # In[ ]:


# print("percentage of the data reduction : ", 100*(1.0 - df_bound.shape[0]/df.shape[0]), '%')


# # In[ ]:


# np.sum([dfl.shape[0] for dfl in df_list])


# # In[ ]:


# df.iloc[idx].gn


# # In[ ]:


# pnts3D


# # In[ ]:


# np.median(np.sort(np.unique(pnts3D[:,0])))


# # In[ ]:





# # In[ ]:


# np.unique(np.diff(xmap))


# # In[ ]:


# np.unique(np.diff(ymap))


# # In[ ]:


# np.unique(np.diff(zmap))


# # In[ ]:


# xmap


# # In[ ]:


# len(xmap)


# # In[ ]:
