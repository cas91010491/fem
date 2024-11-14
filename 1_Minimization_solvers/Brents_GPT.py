import numpy as np

def brents_method(f, a, b, tol=1e-5, max_iter=100):
    """
    Brent's method for finding the root of f within the interval [a, b].
    
    Parameters:
    f       -- The function whose root we want to find (in your case, the derivative of the objective function)
    a, b    -- The initial interval [a, b] where a root is suspected to exist
    tol     -- Tolerance for the stopping criterion
    max_iter -- Maximum number of iterations
    
    Returns:
    The value of x that is the root of f(x) within the interval [a, b].
    """
    
    fa = f(a)
    fb = f(b)
    
    if fa * fb > 0:
        raise ValueError("The function must have different signs at a and b")
    
    if abs(fa) < abs(fb):
        a, b = b, a
        fa, fb = fb, fa
    
    c = a
    fc = fa
    d = e = b - a
    
    for iteration in range(max_iter):
        if fb == 0 or abs(b - a) < tol:
            return b, fb  # Root found
        
        if fa != fc and fb != fc:
            # Inverse quadratic interpolation
            s = (a * fb * fc) / ((fa - fb) * (fa - fc)) + \
                (b * fa * fc) / ((fb - fa) * (fb - fc)) + \
                (c * fa * fb) / ((fc - fa) * (fc - fb))
        else:
            # Secant method
            s = b - fb * (b - a) / (fb - fa)
        
        # Conditions for bisection method
        cond1 = (s < (3 * a + b) / 4 or s > b)
        cond2 = (e < tol or abs(s - b) >= abs(e) / 2)
        cond3 = (abs(s - b) >= abs(d) / 2)
        cond4 = (abs(b - a) < tol)
        cond5 = (abs(fb) >= abs(fa))

        if cond1 or cond2 or cond3 or cond4 or cond5:
            # Perform bisection if conditions are met
            s = (a + b) / 2
            e = d = b - a
        
        else:
            d = e
            e = b - s
        
        fs = f(s)
        a, fa = b, fb
        if fa * fs < 0:
            b, fb = s, fs
        else:
            c, fc = s, fs
        
        if abs(fa) < abs(fb):
            a, b = b, a
            fa, fb = fb, fa
    
    raise RuntimeError(f"Maximum iterations ({max_iter}) reached without convergence")

# Example usage:
# Define the derivative of your objective function
def f_derivative(x):
    return 2 * x - 4  # For example: derivative of f(x) = x^2 - 4x + 4

# Find the root of the derivative (where the function is minimized)
root = brents_method(f_derivative, a=0, b=5, tol=1e-6)
print(f"Root found at x = {root}")