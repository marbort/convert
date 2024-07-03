import numpy as np
import matplotlib.pyplot as plt
import sys

with open(sys.argv[1],'r') as ifile:
    lines=ifile.readlines()

# Generate example data (you can replace this with your own data)
x = [float(x.split()[0]) for x in lines if ("&" not in x) and ("@" not in x) and ("#" not in x)]
y = [float(x.split()[1]) for x in lines if ("&" not in x) and  ("@" not in x) and ("#" not in x)]

# Smooth the data (optional)
#smoothed_y = np.convolve(y, np.ones(10)/10, mode='valid')

# Compute the second derivative using central differences
second_derivative = np.gradient(np.gradient(y, x), x)

# Find indices where the second derivative changes sign
inflection_indices = np.where(np.diff(np.sign(second_derivative)))[0]

inflection_x=np.array(x)[inflection_indices]
inflection_y=np.array(y)[inflection_indices]


# Plot the data and inflection points
plt.figure(figsize=(8, 6))
plt.plot(x, y, label='Data')
plt.scatter(inflection_x, inflection_y, color='red', label='Inflection Points')
plt.xlabel('X')
plt.ylabel('Y')
plt.title('Inflection Points Detection')
plt.legend()
plt.grid(True)
plt.show()