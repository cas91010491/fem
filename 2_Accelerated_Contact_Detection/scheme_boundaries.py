import numpy as np
import matplotlib.pyplot as plt

# Parameters for the shape
theta = np.linspace(0, 2 * np.pi, 50000)
r = 1 + 0.2 * np.sin(3 * theta) + 0.1 * np.cos(7 * theta) + 0.3 * np.cos(1 + theta)

# Original curve in Cartesian coordinates
x = r * np.cos(theta)
y = r * np.sin(theta)

# Offset distances
offset_distance = -4.5 / 4
in_set_distance = -1.5 / 4

# Compute derivatives for normals
dx = np.gradient(x)
dy = np.gradient(y)
norm = np.sqrt(dx**2 + dy**2)

# Normalized normal vectors (perpendicular)
nx = -dy / norm
ny = dx / norm

# Outer and inner offset curves
x_outer = x + offset_distance * nx
y_outer = y + offset_distance * ny
x_inner = x - in_set_distance * nx
y_inner = y - in_set_distance * ny

# Number of zones
num_zones = 5
zone_indices = np.linspace(0, len(theta) - 1, num_zones + 1).astype(int)

# Colors for each zone
zone_colors = ['lightcoral', 'lightskyblue', 'lightgreen', 'khaki', 'plum']

# Plotting
plt.figure(figsize=(6, 6))
plt.plot(x, y, color='k', lw=2.0)
plt.plot(x_outer, y_outer, '--', color='k')
plt.plot(x_inner, y_inner, '--', color='k')

# Fill zones between inner and outer offset curves
for i in range(num_zones):
    # Indices defining each segment
    start_idx, end_idx = zone_indices[i], zone_indices[i + 1]
    
    # Boundary points within this segment
    x_zone_inner, y_zone_inner = x_inner[start_idx:end_idx], y_inner[start_idx:end_idx]
    x_zone_outer, y_zone_outer = x_outer[start_idx:end_idx], y_outer[start_idx:end_idx]

    # Start and end points on the original curve
    x_start, y_start = x[start_idx], y[start_idx]
    x_end, y_end = x[end_idx], y[end_idx]

    # Normals for perpendicular lines at segment boundaries
    nx_start, ny_start = nx[start_idx], ny[start_idx]
    nx_end, ny_end = nx[end_idx], ny[end_idx]

    # Perpendicular line endpoints
    x_inner_start, y_inner_start = x_zone_inner[0], y_zone_inner[0]
    x_inner_end, y_inner_end = x_zone_inner[-1], y_zone_inner[-1]
    
    x_outer_start, y_outer_start = x_zone_outer[0], y_zone_outer[0]
    x_outer_end, y_outer_end = x_zone_outer[-1], y_zone_outer[-1]

    # Build coordinates for the current segment (inner to outer boundary and two perpendicular radial lines)
    x_zone = np.concatenate([[x_outer_start], x_zone_outer, [x_outer_end], 
                             [x_inner_end], x_zone_inner[::-1], [x_inner_start]])
    y_zone = np.concatenate([[y_outer_start], y_zone_outer, [y_outer_end], 
                             [y_inner_end], y_zone_inner[::-1], [y_inner_start]])

    # Plot the filled zone
    if zone_colors[i % len(zone_colors)]=='lightskyblue': 
        zorder = 1e6
    else:
        zorder = None
    plt.fill(x_zone, y_zone, color=zone_colors[i % len(zone_colors)],zorder=zorder, alpha=0.6)

# Fill the inner offset region in white to cover the center
plt.fill(x_inner, y_inner, color='white', alpha=1.0)

# Formatting
plt.axis('equal')
plt.axis('off')
plt.show()
