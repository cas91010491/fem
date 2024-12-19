import argparse
import sys, os, pickle
import numpy as np
from pdb import set_trace
from time import time
import seaborn as sns


sys.path.append(os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from PyClasses.FEAssembly import *
from PyClasses.Contacts import *
from PyClasses.FEModel import *
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
import pandas as pd
from keras.models import load_model
from tensorflow.keras.losses import MeanSquaredError



import pandas as pd
import os
os.chdir(sys.path[0])




# Define headers
points_headers = ['x', 'y', 'z', 'p_id', 'xi1', 'xi2', 'gn']
derivs_headers = ['dgn1', 'dgn2', 'dgn3', 'd2gn11', 'd2gn22', 'd2gn33', 'd2gn12', 'd2gn13', 'd2gn23']
all_headers = points_headers + derivs_headers

# Initialize an empty DataFrame
combined_df = pd.DataFrame(columns=all_headers)

# Loop through the file indices
for i in range(3):
    for j in range(3):
        for k in range(3):
            points_file = f'csv_files/Points_{i}_{j}_{k}.csv'
            derivs_file = f'derivs/derivs_{i}_{j}_{k}.csv'
            
            # Read the points file (no header)
            points_df = pd.read_csv(points_file, header=None, names=points_headers)
            
            # Read the derivs file (with header)
            derivs_df = pd.read_csv(derivs_file)
            
            # Concatenate the dataframes
            merged_df = pd.concat([points_df, derivs_df], axis=1)
            
            # Append to the combined dataframe
            combined_df = pd.concat([combined_df, merged_df], ignore_index=True)

# Now combined_df contains all the data with the desired headers
print(combined_df.head())



#####################
### Setting model ###
#####################

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



# Load the ANN model
ann_model = load_model('model_epoch_1000.h5', custom_objects={'mse': MeanSquaredError()})

# Load test data as DataFrame
df = pd.read_csv('MultiTaskNN/test_data.csv', delimiter=',', skiprows=1, names=['x', 'y', 'z', 'p_id', 'xi1', 'xi2', 'gn'])
df = df[(df['gn'] > -0.25) & (df['gn'] < 0.5)]
# df = df.head(300)

# Reset the index of the DataFrame
df.reset_index(drop=True, inplace=True)


# import pdb; pdb.set_trace()

num_rows = len(df)
print(f"Number of rows in the DataFrame: {num_rows}")



# Make predictions using the ANN model
predictions = ann_model.predict(df[['x', 'y', 'z']].values)
dist_pred = predictions[0]
clas_pred = predictions[1]
proj_pred = predictions[2]


n_candids = 4

possible_actives = np.where(dist_pred<0.1)[0]   # Slave nodes close to the master surface
candids = -1*np.ones((num_rows,n_candids),dtype=int)
candids[possible_actives] = np.argsort(clas_pred[possible_actives], axis=1)[:, ::-1][:, :n_candids]
# self.candids[possible_actives] = np.argsort(predictions[1][possible_actives], axis=1)[:, ::-1][:, :n_candids]
t1t2 = np.array(proj_pred[:,:-1].reshape(-1,96,2,order='F'),dtype=np.float64)


# Round the coordinates to 6 decimal places to account for precision errors
df_rounded = df.copy()
df_rounded[['x', 'y', 'z']] = df_rounded[['x', 'y', 'z']].round(6)

combined_df_rounded = combined_df.copy()
combined_df_rounded[['x', 'y', 'z']] = combined_df_rounded[['x', 'y', 'z']].round(6)

# Merge df_rounded with combined_df_rounded on ['x', 'y', 'z'] to find matching coordinates
merged_df = df_rounded.merge(combined_df_rounded, on=['x', 'y', 'z'], how='inner')

# Extract the corresponding derivative values
derivs = merged_df[['dgn1', 'dgn2', 'dgn3', 'd2gn11', 'd2gn22', 'd2gn33', 'd2gn12', 'd2gn13', 'd2gn23']]

# Now 'derivs' contains the derivative values for the matching coordinates
print(derivs.head())



# import pdb; pdb.set_trace()


# Arrays to store predictions
diff_GN = np.zeros(num_rows)
diff_DGNDU = np.zeros(num_rows)
diff_D2GNDU2 = np.zeros(num_rows)

# Initialize lists to store the rows and differences
no_patch_found_rows = []

# Iterate over points generated to get closest patch and surface parameters
patches = model.bodies[0].surf.patches
npatches = len(patches)


# import pdb; pdb.set_trace()

opa = 0

candids = clas_pred.argsort(axis=1)[:, -n_candids:]
for index, row in df.iterrows():
    xsi = row[['x', 'y', 'z']].values
    t1_true, t2_true = row[['xi1', 'xi2']].values
    pid_true = row['p_id']
    gn_true = row['gn']
    
    is_patch_correct = False        # measures ONLY tangential correspondance
    # changed = False

    # set_trace()
    at_least_one = 0
    more_than_one = 0

    nodedata = 100*np.ones((len(candids[index]),12))
    for ican,candid in enumerate(candids[index]):
        patch = patches[candid]
        t0 = t1t2[index,candid].copy()
        t = patch.findProjection(xsi,ANNapprox=True,t0=t0,recursive=1)
        t, gn ,dgndu,d2gd2u = patch.KC_fless_rigidMaster(xsi,10,cubicT=None,t=t, test_gns=True)
        is_patch_correct = 0-opa<t[0]<1+opa and 0-opa<t[1]<1+opa    #boolean
        if is_patch_correct:
            at_least_one = 1
            more_than_one += 1
            nodedata[ican] = [gn,t[0],t[1],dgndu[0],dgndu[1],dgndu[2],d2gd2u[0,0],d2gd2u[0,1],d2gd2u[0,2],d2gd2u[1,1],d2gd2u[1,2],d2gd2u[2,2]]
    

    if at_least_one>0:

        right_cand = np.argmin(np.abs(nodedata[:,0]))

        # if more_than_one>1 and right_cand>0:
        #     set_trace()

        gn,t[0],t[1],dgndu[0],dgndu[1],dgndu[2] = nodedata[right_cand,0:6]
        d2gd2u[0,0],d2gd2u[0,1],d2gd2u[0,2] = nodedata[right_cand,6:9]
        d2gd2u[1,1],d2gd2u[1,2],d2gd2u[2,2] = nodedata[right_cand,9:12]
        d2gd2u[1,0] = d2gd2u[0,1]
        d2gd2u[2,0] = d2gd2u[0,2]
        d2gd2u[2,1] = d2gd2u[1,2]

        
        patch_id = candids[index][right_cand]

    else:
        print("No correct patch found for point",index)
        no_patch_found_rows.append(row[['x', 'y', 'z', 'p_id', 'xi1', 'xi2', 'gn']])
        continue
        # import pdb; pdb.set_trace()

    # import pdb; pdb.set_trace()

    try:
        diff_GN[index] = abs(gn-gn_true)
        diff_DGNDU[index] = norm(dgndu-derivs.iloc[index][['dgn1', 'dgn2', 'dgn3']].values)
        D2true6 = derivs.iloc[index][['d2gn11', 'd2gn22', 'd2gn33', 'd2gn12', 'd2gn13', 'd2gn23']].values
        D2GN_true = np.array([[D2true6[0], D2true6[3], D2true6[4]],[D2true6[3], D2true6[1], D2true6[5]], [D2true6[4], D2true6[5], D2true6[2]]])
        diff_D2GNDU2[index] = norm(d2gd2u-D2GN_true)
    except:
        import pdb; pdb.set_trace()


import pdb; pdb.set_trace()


error = diff_GN

dgnlog = np.array([np.log10(er) if er != 0 else -15 for er in error])

x_data = dgnlog  # Filtered X data
y_data = df['gn']     # Corresponding Y data








plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 12,  # Adjust axis label size
    "axes.titlesize": 14,  # Adjust title size
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
})


# Create the figure and axis layout
fig = plt.figure(figsize=(10, 8))
grid = plt.GridSpec(4, 100, hspace=0.05, wspace=0.05)  

# Main hexbin plot
main_ax = fig.add_subplot(grid[1:4, 0:75])  # Allocate 75% of width to the main plot
hb = main_ax.hexbin(x_data, y_data, gridsize=50, cmap='jet', mincnt=1)
main_ax.set_xlabel(r'$\log_{10}{e(\hat{g} )}$', fontsize=16)  # LaTeX-style label
main_ax.set_ylabel(r'\textbf{True Signed Distance}', fontsize=14)  # LaTeX-style label

cbar_ax = fig.add_subplot(grid[1:4, 95:98])  
cbar = fig.colorbar(hb, cax=cbar_ax, orientation='vertical')
cbar_ax.set_ylabel(r'\textbf{Frequency}', fontsize=14)

# Top histogram
x_hist = fig.add_subplot(grid[0, 0:75], sharex=main_ax)
sns.histplot(x=x_data, bins=50, kde=True, ax=x_hist, color='blue')
x_hist.set_ylabel(r'\textbf{Frequency}', fontsize=12)
x_hist.tick_params(axis='x', labelbottom=False)  
x_hist.set_xlim(main_ax.get_xlim()) 

# Side histogram
y_hist = fig.add_subplot(grid[1:4, 76:90])  # Side histogram next to main plot
sns.histplot(y=y_data, bins=50, kde=True, ax=y_hist, color='blue', orientation='horizontal')
y_hist.set_xlabel(r'\textbf{Frequency}', fontsize=12)
y_hist.tick_params(axis='y', labelleft=False)  # Hide y-axis labels for side histogram
y_hist.set_ylim(main_ax.get_ylim())  # Align y-axis limits with the main plot

# Adjust layout
plt.tight_layout(rect=[0, 0, 1, 0.95])  # Add space for the title
plt.savefig('MY_GN.pdf')
plt.show()

# At the end of the script, save the no_patch_found_rows and diffs as CSV files
no_patch_found_df = pd.DataFrame(no_patch_found_rows)
no_patch_found_df.to_csv('no_patch_found_rows.csv', index=False)

# import pdb; pdb.set_trace()
diffs = np.array([df['gn'],diff_GN, diff_DGNDU, diff_D2GNDU2]).T
diffs_df = pd.DataFrame(diffs, columns=['gn', 'e_gn', 'e_dgn', 'e_d2gn'])
diffs_df.to_csv('diffs.csv', index=False)

