import os
import sys
import numpy as np
import pandas as pd
import tensorflow as tf
from tensorflow.keras.models import load_model
from tensorflow.keras.utils import register_keras_serializable
from sklearn.metrics import classification_report, confusion_matrix
import seaborn as sns

import matplotlib.pyplot as plt

# Set the script's directory as the current working directory
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

# Add the script's directory to the system path
sys.path.append(script_dir)


# Define the headers
headers = ['x', 'y', 'z', 'p_id', 'xi1', 'xi2', 'gn']

# Load the data from the CSV file with headers
points_on_surf_path = os.path.join(script_dir, '..', 'Points_on_surf.csv')
data = pd.read_csv(points_on_surf_path, header=None, names=headers)

print(data.columns)


# Step 1: Filter the data
edge_data = data[(data['xi1'].isin([0, 1])) | (data['xi2'].isin([0, 1]))]

# Step 2: Create a dictionary to store the points for each patch
patch_points = {}
for _, row in edge_data.iterrows():
    p_id = row['p_id']
    point = (row['x'], row['y'], row['z'])
    if p_id not in patch_points:
        patch_points[p_id] = set()
    patch_points[p_id].add(point)

# Step 3: Initialize the neighborhood matrix
patch_ids = list(patch_points.keys())
num_patches = len(patch_ids)
neighbourhood_matrix = np.zeros((num_patches, num_patches), dtype=int)

# Step 4: Compare the points of each patch pair
for i in range(num_patches):
    for j in range(i + 1, num_patches):
        p_id1 = patch_ids[i]
        p_id2 = patch_ids[j]
        common_points = patch_points[p_id1].intersection(patch_points[p_id2])
        if len(common_points) == 1:
            neighbourhood_matrix[i, j] = 1
            neighbourhood_matrix[j, i] = 1
        elif len(common_points) > 1:
            neighbourhood_matrix[i, j] = 2
            neighbourhood_matrix[j, i] = 2

# Convert the neighborhood matrix to a DataFrame for better readability
neighbourhood_matrix_df = pd.DataFrame(neighbourhood_matrix, index=patch_ids, columns=patch_ids)

# Print the neighborhood matrix
# print(neighbourhood_matrix_df)

# Plot the neighborhood matrix in gray scale
plt.figure(figsize=(10, 8))
sns.heatmap(neighbourhood_matrix_df, cmap='gray', cbar=False, linewidths=.5, linecolor='black')
plt.gca().set_aspect('equal', adjustable='box')
plt.xticks(ticks=np.arange(0, len(patch_ids), 2) + 0.5, labels=np.array(patch_ids[::2], dtype=int), rotation=90)
plt.yticks(ticks=np.arange(0, len(patch_ids), 2) + 0.5, labels=np.array(patch_ids[::2], dtype=int), rotation=0)
plt.title('Neighbourhood Matrix')
plt.xlabel('Patch ID')
plt.ylabel('Patch ID')
plt.show()


