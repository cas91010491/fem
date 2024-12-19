import numpy as np
import matplotlib.pyplot as plt

# Update plot parameters
plt.rcParams.update({
    "text.usetex": True,
    "font.family": "serif",
    "axes.labelsize": 12,  # Adjust axis label size
    "axes.titlesize": 14,  # Adjust title size
    "xtick.labelsize": 10,
    "ytick.labelsize": 10,
})

# Create a smaller figure
plt.figure(figsize=(6, 4))  # Adjust the width and height as needed

# Define the activation functions
def sigmoid(x):
    return 1 / (1 + np.exp(-x))

def relu(x):
    return np.maximum(0, x)

def tanh(x):
    return np.tanh(x)

# Generate x values
x = np.linspace(-10, 10, 400)

# Compute y values for each activation function
y_sigmoid = sigmoid(x)
y_relu = relu(x)
y_tanh = tanh(x)

# Create the plot
plt.plot(x, y_relu, label='ReLU', color='orange', linewidth=3)
plt.plot(x, y_tanh, label='Tanh', color='red', linestyle='--')
plt.plot(x, y_sigmoid, label='Sigmoid', color='blue', linestyle=':', linewidth=2)

# Add title and labels
plt.title('Activation Functions')

# Set lims
plt.xlim(-6, 6)
plt.ylim(-2, 4)

# Remove axis numbers
plt.xticks([])
plt.yticks([])

# Add legend
plt.legend()

# Set aspect ratio to square
plt.gca().set_aspect('equal', adjustable='box')

# Add axis labels inside the plot
plt.text(5.5, -0.5, 'z', fontsize=12, verticalalignment='bottom', horizontalalignment='right')
plt.text(-0.5, 3.5, 'y', fontsize=12, verticalalignment='top', horizontalalignment='right')

# Show the plot
plt.axhline(0, color='black', linewidth=0.5)
plt.axvline(0, color='black', linewidth=0.5)
plt.show()