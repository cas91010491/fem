import os
import pandas as pd
import matplotlib.pyplot as plt

# Set the main directory as the current working directory
main_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(main_dir)

# Load the training history from the CSV file
training_history_path = os.path.join(main_dir, 'training_history.csv')
training_history = pd.read_csv(training_history_path)

# Extract the relevant columns
signed_distance_loss = training_history['signed_distance_output_loss']
val_signed_distance_loss = training_history['val_signed_distance_output_loss']
classification_loss = training_history['classification_output_loss']
val_classification_loss = training_history['val_classification_output_loss']
projection_loss = training_history['projection_output_loss']
val_projection_loss = training_history['val_projection_output_loss']

# Configure matplotlib to use LaTeX for text rendering and set font sizes
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "font.serif": ["Computer Modern Roman"],
    "axes.titlesize": 15,  # Increased by 50%
    "axes.labelsize": 12,  # Increased by 50%
    "xtick.labelsize": 9,  # Increased by 50%
    "ytick.labelsize": 9,  # Increased by 50%
    "legend.fontsize": 12,  # Increased by 50%
    "figure.figsize": (9, 3)
})

# Create a figure with three subplots
fig, axs = plt.subplots(1, 3)

# Plot signed distance loss
axs[0].plot(signed_distance_loss, label='Train')
axs[0].plot(val_signed_distance_loss, label='Validation')
axs[0].set_title('Signed Distance Loss')
axs[0].set_xlabel('Epoch')
axs[0].set_ylabel('Loss')
axs[0].legend()

# Plot classification loss
axs[1].plot(classification_loss, label='Training')
axs[1].plot(val_classification_loss, label='Validation')
axs[1].set_title('Classification Loss')
axs[1].set_xlabel('Epoch')
axs[1].set_ylabel('Loss')
axs[1].legend()

# Plot projection loss
axs[2].plot(projection_loss, label='Training')
axs[2].plot(val_projection_loss, label='Validation')
axs[2].set_title('Projection Loss')
axs[2].set_xlabel('Epoch')
axs[2].set_ylabel('Loss')
axs[2].legend()

# Adjust layout
plt.tight_layout()

# Save the plot as a PDF file
fig.savefig('MTL_TrainingHist.pdf')

# Show the plot
plt.show()

# Create a second figure with log scale on the y-axis
fig_log, axs_log = plt.subplots(1, 3, figsize=(9, 3))

# Plot signed distance loss with log scale
axs_log[0].plot(signed_distance_loss, label='Training')
axs_log[0].plot(val_signed_distance_loss, label='Validation')
axs_log[0].set_title('Signed Distance Loss (Log Scale)')
axs_log[0].set_xlabel('Epoch')
axs_log[0].set_ylabel('Loss')
axs_log[0].set_yscale('log')
axs_log[0].legend()

# Plot classification loss with log scale
axs_log[1].plot(classification_loss, label='Training')
axs_log[1].plot(val_classification_loss, label='Validation')
axs_log[1].set_title('Classification Loss (Log Scale)')
axs_log[1].set_xlabel('Epoch')
axs_log[1].set_ylabel('Loss')
axs_log[1].set_yscale('log')
axs_log[1].legend()

# Plot projection loss with log scale
axs_log[2].plot(projection_loss, label='Training')
axs_log[2].plot(val_projection_loss, label='Validation')
axs_log[2].set_title('Projection Loss (Log Scale)')
axs_log[2].set_xlabel('Epoch')
axs_log[2].set_ylabel('Loss')
axs_log[2].set_yscale('log')
axs_log[2].legend()

# Adjust layout
plt.tight_layout()

# Save the log scale plot as a PDF file
fig_log.savefig('MTL_TrainingHist_LogScale.pdf')

# Show the log scale plot
plt.show()

