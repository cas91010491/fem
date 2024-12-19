import argparse
import sys, os, pickle
import numpy as np
from pdb import set_trace
from time import time
import seaborn as sns
import matplotlib.pyplot as plt


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


# Load diffs.csv as a dataframe
diffs_df = pd.read_csv('diffs.csv')

# data = pd.read_csv('test_data.csv')
# data = data[(data['gn'] > -0.25) & (data['gn'] < 0.5)]



import pdb; pdb.set_trace()

# Define the columns to plot and their corresponding x-labels
columns_to_plot = ['e_gn', 'e_dgn', 'e_d2gn']
x_labels = [
    r'$\log_{10}{e(\hat{g} )}$', 
    r'$\log_{10}{e(\hat{\mathbf{n}}_\textrm{lead} )}$', 
    r'$\log_{10}{e(\hat{\frac{\textrm{d}\,\mathbf{n}_\textrm{lead}}{\textrm{d}\,\mathbf{u}_i}} )}$'
]

# Update plot parameters
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 12,  # Adjust axis label size
    "axes.titlesize": 14,  # Adjust title size
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
})

# Loop through each column and create a plot
for col, x_label in zip(columns_to_plot, x_labels):
    x_data = np.array([np.log10(er) if er != 0 else -15 for er in diffs_df[col]])
    y_data = diffs_df['gn']

    # Create the figure and axis layout
    fig = plt.figure(figsize=(8, 6))  # Smaller figure size
    grid = plt.GridSpec(4, 100, hspace=0.05, wspace=0.05)  

    # Main hexbin plot
    main_ax = fig.add_subplot(grid[1:4, 0:75])  # Allocate 75% of width to the main plot
    hb = main_ax.hexbin(x_data, y_data, gridsize=50, cmap='jet', mincnt=1)
    main_ax.set_xlabel(x_label, fontsize=16)  # LaTeX-style label
    main_ax.set_ylabel(r'\textbf{True Signed Distance}', fontsize=14)  # LaTeX-style label

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

    # Adjust layout to add more margin at the bottom
    plt.subplots_adjust(bottom=0.3)

    # Save the plot as a PDF
    plt.savefig(f'{col}_vs_gn.pdf')
    plt.close()

# Extra plot for gn distribution where e_gn > 1e-10
filtered_gn = diffs_df[diffs_df['e_gn'] > 1e-10]['gn']

fig = plt.figure(figsize=(8, 6))  # Smaller figure size
sns.histplot(filtered_gn, bins=50, kde=True, color='blue')
plt.xlabel(r'\textbf{True Signed Distance}', fontsize=14)
plt.ylabel(r'\textbf{Frequency}', fontsize=12)
plt.title(r'\textbf{Distribution of gn for $e\_\hat{g} > 1e-10$}', fontsize=16)

# Adjust layout to add more margin at the bottom
plt.subplots_adjust(bottom=0.3)

# Save the plot as a PDF
plt.savefig('gn_distribution_for_e_gn_greater_than_1e-10.pdf')
plt.close()

# Extra plot for gn distribution where e_dgn > 1e-7
filtered_gn = diffs_df[diffs_df['e_dgn'] > 1e-7]['gn']

fig = plt.figure(figsize=(8, 6))  # Smaller figure size
sns.histplot(filtered_gn, bins=50, kde=True, color='blue')
plt.xlabel(r'\textbf{True Signed Distance}', fontsize=14)
plt.ylabel(r'\textbf{Frequency}', fontsize=12)
plt.title(r'\textbf{Distribution of gn for $e\_\hat{\mathbf{n}}_\textrm{lead} > 1e-7$}', fontsize=16)

# Adjust layout to add more margin at the bottom
plt.subplots_adjust(bottom=0.3)

# Save the plot as a PDF
plt.savefig('gn_distribution_for_e_dgn_greater_than_1e-7.pdf')
plt.close()

# Extra plot for gn distribution where e_d2gn > 1e-4
filtered_gn = diffs_df[diffs_df['e_d2gn'] > 1e-4]['gn']

fig = plt.figure(figsize=(8, 6))  # Smaller figure size
sns.histplot(filtered_gn, bins=50, kde=True, color='blue')
plt.xlabel(r'\textbf{True Signed Distance}', fontsize=14)
plt.ylabel(r'\textbf{Frequency}', fontsize=12)
plt.title(r'\textbf{Distribution of gn for $e\_\hat{\frac{\textrm{d}\,\mathbf{n}_\textrm{lead}}{\textrm{d}\,\mathbf{u}_i}} > 1e-4$}', fontsize=16)

# Adjust layout to add more margin at the bottom
plt.subplots_adjust(bottom=0.3)

# Save the plot as a PDF
plt.savefig('gn_distribution_for_e_d2gn_greater_than_1e-4.pdf')
plt.close()

