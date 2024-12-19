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
data = data[(data['gn'] > -0.25) & (data['gn'] < 0.5)]
# data = data.head(100)

# Load the model from the parent directory
# import pdb; pdb.set_trace()

# Load the neural network model
model_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'model_epoch_1000.h5')
model = load_model(model_path, custom_objects={'mse': MeanSquaredError()})

# Define 'points' as the first three columns of data (headers 'x', 'y', 'z')
points = data[['x', 'y', 'z']]
gn_true = data['gn']


gn_predicted = model.predict(points)[0]


error = np.abs(gn_predicted.reshape(-1)-gn_true)

dgnlog = np.array([np.log10(er) if er != 0 else -15 for er in error])

x_data = dgnlog  # Filtered X data
y_data = gn_true     # Corresponding Y data

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
plt.savefig('SignedDist.pdf')
plt.show()