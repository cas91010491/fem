import numpy as np

# Example array for multiple nodes
m = 3  # number of nodes
n = 4  # number of patches per node

# Example array of shape (m, n, 2)
patch_probabilities = np.array([
    [
        [66, 1.0],
        [65, 8.93964461e-20],
        [50, 2.49087664e-32],
        [30, 1.10387883e-32]
    ],
    [
        [45, 0.9],
        [46, 0.05],
        [47, 0.03],
        [48, 0.02]
    ],
    [
        [77, 0.6],
        [76, 0.3],
        [75, 0.1],
        [74, 0.0]
    ]
])

# Select the best patch for each node
best_patches = patch_probabilities[:, 0, :]

# Output the results
for node_index, (best_patch_id, highest_probability) in enumerate(best_patches):
    print(f"Node {node_index}: Best patch ID: {best_patch_id}, with probability: {highest_probability}")
