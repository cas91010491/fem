import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from keras.models import load_model
from keras.losses import MeanSquaredError


# Update plot parameters
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

# Filter the data
data = data[(data['gn'] > -0.25) & (data['gn'] < 0.5)]
data.reset_index(drop=True, inplace=True)

# Load the neural network model
model_path = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'model_epoch_1000.h5')
model = load_model(model_path, custom_objects={'mse': MeanSquaredError()})

# Predict for all points
points = data[['x', 'y', 'z']]
XI12 = model.predict(points)[2]

# Initialize a list to store errors for each patch
all_errors = []

for query_patch in range(96):
    patch_data = data[data['p_id'] == query_patch]
    indices = patch_data.index

    xi1_true = patch_data['xi1']
    xi2_true = patch_data['xi2']

    xi1_pred = XI12[indices, query_patch]
    xi2_pred = XI12[indices, 96 + query_patch]

    error_xi1 = np.abs(xi1_pred - xi1_true)
    error_xi2 = np.abs(xi2_pred - xi2_true)

    error = np.log10(np.sqrt(error_xi1**2 + error_xi2**2))
    all_errors.append(error)

minerr = np.min([np.min(errors) for errors in all_errors])
maxerr = np.max([np.max(errors) for errors in all_errors])

# Create a 2D histogram plot using imshow
hist_data = np.zeros((96, 200))

for i, errors in enumerate(all_errors):
    hist, bins = np.histogram(errors, bins=200)
    hist_data[i, :] = hist

# Apply logarithm to the frequency data
log_hist_data = np.log1p(hist_data)  # Use log1p to avoid log(0)

plt.imshow(log_hist_data.T, aspect='auto', cmap='viridis', origin='lower',extent=[0, 95,minerr,maxerr])
# plt.colorbar(label=r'$\log_{10}(\text{Frequency})$')
plt.colorbar(label=r'$\log_{10}$(Frequency)')
plt.xlabel('Patch Number')
plt.ylabel(r'Euclidean error for ($\hat{\xi_1}$,$\hat{\xi_2}$)')
plt.title('Logarithm of Error Distribution Across Patches')

# Adjust y-axis to represent error range from 0 to 1
# plt.yticks(np.linspace(0, 199, 11), np.round(np.linspace(0, 1, 11), 2))

# Save the plot as a PDF file
output_path = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'error_xix1_all_patches.pdf')
plt.savefig(output_path, format='pdf')

plt.show()

