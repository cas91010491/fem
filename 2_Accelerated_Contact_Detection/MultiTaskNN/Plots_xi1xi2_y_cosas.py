import os
import seaborn as sns
import pandas as pd
from tensorflow.keras.models import load_model
from tensorflow.keras.losses import MeanSquaredError
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 12,  # Adjust axis label size
    "axes.titlesize": 14,  # Adjust title size
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
})

# Set the directory of the main file as the current working directory
os.chdir(os.path.dirname(os.path.abspath(__file__)))

# Load the data from the CSV file
data = pd.read_csv('test_data.csv')

query_patch = 50

data = data[(data['gn'] > -0.25) & (data['gn'] < 0.5)]
data = data[data['p_id'] == query_patch]

# Load the neural network model
model_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'model_epoch_1000.h5')
model = load_model(model_path, custom_objects={'mse': MeanSquaredError()})

# Define 'points' as the first three columns of data (headers 'x', 'y', 'z')
points = data[['x', 'y', 'z']]
xi1_true = data['xi1']
xi2_true = data['xi2']

XI12 = model.predict(points)[2]
XI1 = XI12[:, query_patch]
XI2 = XI12[:, 96 + query_patch]

# Calculate errors
error_xi1 = np.abs(XI1.reshape(-1) - xi1_true)
error_xi2 = np.abs(XI2.reshape(-1) - xi2_true)

dgnlog_xi1 = np.array([np.log10(er) if er != 0 else -15 for er in error_xi1])
dgnlog_xi2 = np.array([np.log10(er) if er != 0 else -15 for er in error_xi2])

x_data_xi1 = dgnlog_xi1  # Filtered X data for xi1
y_data_xi1 = xi1_true    # Corresponding Y data for xi1

x_data_xi2 = dgnlog_xi2  # Filtered X data for xi2
y_data_xi2 = xi2_true    # Corresponding Y data for xi2

# Create the figure and axis layout
fig = plt.figure(figsize=(10, 5))
grid = plt.GridSpec(4, 220, hspace=0.05, wspace=0.1)  # Increase wspace for more separation

# Main hexbin plot for xi1
main_ax_xi1 = fig.add_subplot(grid[1:4, 0:85])  # Allocate 85 columns to the main plot for xi1
hb_xi1 = main_ax_xi1.hexbin(x_data_xi1, y_data_xi1, gridsize=50, cmap='jet', mincnt=1)
main_ax_xi1.set_xlabel(r'$\log_{10}{e(\xi_1)}$', fontsize=16)  # LaTeX-style label
main_ax_xi1.set_ylabel(r'$\xi_1$', fontsize=16)

# Top histogram for xi1
top_hist_xi1 = fig.add_subplot(grid[0, 0:85], sharex=main_ax_xi1)
top_hist_xi1.hist(x_data_xi1, bins=50, color='gray')
top_hist_xi1.axis('off')

# Side histogram for xi1
side_hist_xi1 = fig.add_subplot(grid[1:4, 85:90], sharey=main_ax_xi1)
side_hist_xi1.hist(y_data_xi1, bins=50, orientation='horizontal', color='gray')
side_hist_xi1.axis('off')

# Invisible plot to create space
invisible_ax = fig.add_subplot(grid[1:4, 90:95])
invisible_ax.axis('off')

# Main hexbin plot for xi2
main_ax_xi2 = fig.add_subplot(grid[1:4, 100:185])  # Allocate 85 columns to the main plot for xi2
hb_xi2 = main_ax_xi2.hexbin(x_data_xi2, y_data_xi2, gridsize=50, cmap='jet', mincnt=1)
main_ax_xi2.set_xlabel(r'$\log_{10}{e(\xi_2)}$', fontsize=16)  # LaTeX-style label
main_ax_xi2.set_ylabel(r'$\xi_2$', fontsize=16)

# Top histogram for xi2
top_hist_xi2 = fig.add_subplot(grid[0, 100:185], sharex=main_ax_xi2)
top_hist_xi2.hist(x_data_xi2, bins=50, color='gray')
top_hist_xi2.axis('off')

# Side histogram for xi2
side_hist_xi2 = fig.add_subplot(grid[1:4, 185:190], sharey=main_ax_xi2)
side_hist_xi2.hist(y_data_xi2, bins=50, orientation='horizontal', color='gray')
side_hist_xi2.axis('off')

# Adjust layout
plt.tight_layout()

# Save the plot as a PDF file
fig.savefig('xi1_xi2_plots.pdf')

# Show the plot
plt.show()