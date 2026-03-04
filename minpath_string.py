import numpy as np
import matplotlib.pyplot as plt
from scipy.ndimage import gaussian_filter, sobel
import sys

def interpolate_path(start, end, num_images):
    return np.linspace(start, end, num_images)

def compute_energy_gradient(energy_grid):
    # Use Sobel filters for gradients
    grad_y = sobel(energy_grid, axis=0)
    grad_x = sobel(energy_grid, axis=1)
    print("Gradient shapes:", grad_x.shape, grad_y.shape)
    return -grad_x, -grad_y  # negative gradient = force

def normalize(v):
    return v / (np.linalg.norm(v) + 1e-12)

def move_images(images, grad_x, grad_y, alpha=0.2):
    new_images = []
    for i in range(1, len(images) - 1):  # skip endpoints (fixed)
        print(f"Processing image {i}: {images[i]}")
        y, x = images[i]
        iy, ix = int(round(y)), int(round(x))
        
        # Ensure indices stay inside the grid
        iy = np.clip(iy, 0, grad_x.shape[0] - 1)
        ix = np.clip(ix, 0, grad_x.shape[1] - 1)
        #print(f"{iy}, {ix}")
        # Get local gradient (force)
        fx = grad_x[iy, ix]
        fy = grad_y[iy, ix]
        force = np.array([fy, fx])  # (dy, dx)

        # Tangent to the path (forward - backward)
        tangent = normalize(images[i + 1] - images[i - 1])
        # Project force perpendicular to tangent
        force_perp = force - np.dot(force, tangent) * tangent

        # Move image
        new_pos = images[i] + alpha * force_perp
        new_images.append(new_pos)
    return np.vstack(([images[0]], new_images, [images[-1]]))  # keep endpoints fixed

def redistribute_images(images):
    # Equal arc length redistribution
    cumdist = np.cumsum(np.linalg.norm(np.diff(images, axis=0), axis=1))
    cumdist = np.insert(cumdist, 0, 0)
    total = cumdist[-1]
    interp_dists = np.linspace(0, total, len(images))
    new_images = [images[0]]
    for d in interp_dists[1:-1]:
        idx = np.searchsorted(cumdist, d)
        t = (d - cumdist[idx - 1]) / (cumdist[idx] - cumdist[idx - 1] + 1e-10)
        interp_point = (1 - t) * images[idx - 1] + t * images[idx]
        new_images.append(interp_point)
    new_images.append(images[-1])
    return np.array(new_images)

def string_method(energy_grid, start, end, num_images=20, n_iter=100, alpha=0.2):
    grad_x, grad_y = compute_energy_gradient(energy_grid)
    images = interpolate_path(np.array(start), np.array(end), num_images)
    print(f"Images {images}")

    for _ in range(n_iter):
        images = move_images(images, grad_x, grad_y, alpha=alpha)
        images = redistribute_images(images)

    return images

def extract_data(input):
    data=[]
    with open(input,'r') as ifile:
        file=np.loadtxt(ifile,unpack=True)
    
    cv1=np.unique(file[0])
    cv2=np.unique(file[1])
    val=np.array(file[2])
    free_grid=val.reshape(len(cv1),len(cv2))
    

    return(cv1,len(cv1),cv2,len(cv2),free_grid)


# ==== Example usage ====
# Generate a synthetic PES for demonstration

x,numx,y,numy,Z = extract_data(sys.argv[1])
X, Y = np.meshgrid(x, y)

#x = np.linspace(-3, 3, 101)
#y = np.linspace(-3, 4, 101)
#X, Y = np.meshgrid(x, y)
#Z = X**4 + Y**4 - 4*X**2 - 4*Y**2 + 10*np.exp(-((X+1.5)**2 + (Y+1.5)**2))  # two basins

#assert Z.shape == Z2.shape, "Z must match the grid shape"
#assert X.shape == X2.shape and Y.shape == Y2.shape, "X and Y must match the grid shape"



energy_grid = gaussian_filter(Z, sigma=1)
print("X grid shape:", X.shape)
print("Y grid shape:", Y.shape)
print("Energy grid shape:", energy_grid.shape)

# Map continuous to grid coordinates
def to_grid_coords(xval, yval, xgrid, ygrid):
    x_idx = np.argmin(np.abs(xgrid - xval))
    y_idx = np.argmin(np.abs(ygrid - yval))
    return y_idx, x_idx


start_point = to_grid_coords(0.0, 1.0, x, y)
end_point   = to_grid_coords(2.0,  3.0, x, y)

mep_images = string_method(energy_grid, start_point, end_point,
                           num_images=20, n_iter=200, alpha=0.3)

# ==== Plotting ====
plt.figure(figsize=(8, 6))
plt.contourf(X, Y, energy_grid, levels=40, cmap="viridis")
grid_x = x
grid_y = y
ix = np.clip((mep_images[:,1]).astype(int), 0, len(x) - 1)
iy = np.clip((mep_images[:,0]).astype(int), 0, len(y) - 1)
path_coords_x = x[ix]
path_coords_y = y[iy]
plt.plot(path_coords_x, path_coords_y, 'r.-', label="MEP")
plt.scatter([x[start_point[1]], x[end_point[1]]],
            [y[start_point[0]], y[end_point[0]]],
            c='white', s=100, edgecolor='black', label="Endpoints")
plt.legend()
plt.title("Minimum Energy Path on 2D Grid (String Method)")
plt.xlabel("X")
plt.ylabel("Y")
plt.colorbar(label="Energy")
plt.show()
