import numpy as np
import matplotlib.pyplot as plt

# Define the Rosenbrock function and its gradient
def rosenbrock(x):
    return (1 - x[0])**2 + 100 * (x[1] - x[0]**2)**2

def rosenbrock_gradient(x):
    grad = np.zeros_like(x)
    grad[0] = -2 * (1 - x[0]) - 400 * x[0] * (x[1] - x[0]**2)
    grad[1] = 200 * (x[1] - x[0]**2)
    return grad

def rosenbrock_hessian(x):
    hess = np.zeros((2, 2))
    hess[0, 0] = 2 - 400 * (x[1] - 3 * x[0]**2)
    hess[0, 1] = -400 * x[0]
    hess[1, 0] = -400 * x[0]
    hess[1, 1] = 200
    return hess

# Improved Trust Region Subproblem Solver
def trust_region_subproblem_detailed(grad, hess, delta):
    """
    Solve the trust region subproblem using Newton's method with truncation.
    """
    # Newton's step
    try:
        p_newton = -np.linalg.solve(hess, grad)
    except np.linalg.LinAlgError:
        # Handle singular Hessian by fallback to steepest descent
        p_newton = -grad

    # Check if Newton's step is within the trust region
    if np.linalg.norm(p_newton) <= delta:
        return p_newton

    # If Newton's step is outside, find a point on the trust region boundary
    # Steepest descent direction
    p_steepest = -grad / np.linalg.norm(grad) * delta

    # Compute the intersection point using quadratic equation
    p_boundary = p_steepest  # Default to steepest descent if Newton's direction is invalid
    a = p_newton @ p_newton
    b = 2 * p_newton @ grad
    c = grad @ grad - delta**2
    discriminant = b**2 - 4 * a * c
    if discriminant >= 0:
        tau = (-b + np.sqrt(discriminant)) / (2 * a)
        if 0 <= tau <= 1:
            p_boundary = tau * p_newton

    return p_boundary

# Generate the contour plot for Rosenbrock
x = np.linspace(-2, 2, 400)
y = np.linspace(-1, 3, 400)
X, Y = np.meshgrid(x, y)
Z = rosenbrock([X, Y])

# Reset initial variables for synchronization corrections
xk = np.array([-1.5, 1.5])  # Initial guess
delta = 0.5  # Initial trust region radius
iterations = 30  # Number of iterations

image_counter = 1

for k in range(iterations):
    grad = rosenbrock_gradient(xk)
    hess = rosenbrock_hessian(xk)
    
    # Solve trust region subproblem
    pk = trust_region_subproblem_detailed(grad, hess, delta)

    # Compute new point and reduction ratio
    xk_new = xk + pk
    actual_reduction = rosenbrock(xk) - rosenbrock(xk_new)
    predicted_reduction = -grad @ pk - 0.5 * pk @ hess @ pk
    rho = actual_reduction / predicted_reduction if predicted_reduction != 0 else 0

    # Update trust region radius and decide to accept/reject step
    if rho > 0.75 and np.linalg.norm(pk) >= 0.9 * delta:
        delta = min(2 * delta, 2.0)
    elif rho < 0.25:
        delta = max(0.25 * delta, 0.1)
    if rho > 0.1:
        xk = xk_new  # Accept the step
    else:
        xk_new = xk  # Reject the step, stay at the current point

    # Ensure the new estimate is within the trust region
    if np.linalg.norm(xk - xk_new) > delta:
        xk_new = xk + (delta / np.linalg.norm(xk_new - xk)) * (xk_new - xk)

    # Plot Current Estimate with Newton and Steepest Descent Directions
    newton_step = -np.linalg.solve(hess, grad)
    steepest_descent = -grad / np.linalg.norm(grad) * delta

    # Current estimate
    plt.figure(figsize=(6, 6))
    plt.contour(X, Y, Z, levels=50, cmap='viridis', alpha=0.8)
    plt.plot(xk[0], xk[1], 'ro', label='Current Estimate')  # Current estimate
    plt.arrow(xk[0], xk[1], newton_step[0], newton_step[1],
              color='green', head_width=0.1, length_includes_head=True, label="Newton's Step")
    plt.arrow(xk[0], xk[1], steepest_descent[0], steepest_descent[1],
              color='orange', head_width=0.1, length_includes_head=True, label="Steepest Descent")
    plt.plot(xk_new[0], xk_new[1], 'bo', label='New Estimate')  # New estimate
    circle = plt.Circle((xk[0], xk[1]), delta, color='blue', fill=False, label='Trust Region')
    plt.gca().add_artist(circle)
    plt.title(f"Trust Region Method - Iteration {k + 1}")
    plt.xlabel("x")
    plt.ylabel("y")
    plt.legend(loc='upper right')
    plt.savefig(f"trust_region_step_{image_counter}.png")
    plt.close()
    image_counter += 1

    # Update for next iteration
    xk = xk_new
