#%%
import numpy as np
from scipy.optimize import minimize

# Define the grid points and their corresponding free energy values
x_values = np.linspace(0, 10, 50)  # Replace with appropriate x range and number of points
y_values = np.linspace(0, 10, 50)  # Replace with appropriate y range and number of points

# Create a grid of x, y values
X, Y = np.meshgrid(x_values, y_values)

# Define a simple quadratic free energy surface for demonstration
A = 1.0
B = 0.5
C = 0.1
free_energy_surface = A * X**2 + B * Y**2 + C * X * Y  # Replace with your actual free energy calculations
print(free_energy_surface)

 # Define your 2D free energy surface as a NumPy array
# Replace this with your actual data

# Define the starting and ending points on the surface
start_point = (1, 1)  # Replace with your actual starting point
end_point = (2, 2)        # Replace with your actual ending point

# Define the number of intermediate points for the string
num_intermediate_points = 9

# Initialize the string as a set of linearly spaced points between start and end points
path = np.linspace(start_point, end_point, num_intermediate_points + 2)
print(path)

# Define the energy function to be minimized (distance along the string)
def energy_function(x):
    # Calculate the sum of distances between consecutive points
    distances = np.linalg.norm(np.diff(x, axis=0), axis=1)
    return np.sum(distances)
res=energy_function(path)
print(res)

# Optimize the string using the minimize function from scipy
result = minimize(energy_function, path, method='L-BFGS-B', options={'disp': True})

# Extract the optimized string points
#optimized_string = result.x.reshape(-1, 2)

# Now `optimized_string` contains the optimized path along the free energy surface



# Now, `free_energy_surface` is a 2D array representing the free energy values at each (x, y) point on the grid


# %%
# Make data
u = np.linspace(0, 2 * np.pi, 100)
v = np.linspace(0, np.pi, 100)
x = 10 * np.outer(np.cos(u), np.sin(v))
y = 10 * np.outer(np.sin(u), np.sin(v))
z = 10 * np.outer(np.ones(np.size(u)), np.cos(v))

print(x)
print(y)
print(x.shape,y.shape,z.shape)
# %%
