import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import os
import glob
import plumed
import argparse
import sys
import numpy as np
from scipy.optimize import minimize
import scipy.ndimage.filters as filters
import scipy.ndimage.morphology as morphology
import scipy.ndimage
from matplotlib.colors import ListedColormap, LinearSegmentedColormap


sys.path.append("/home/marco/SHARED/GitHub/mepfinder/mepfinder")
# from  flooder import Flooder


# Python3 program to illustrate
# Saddle point


# Method to find saddle point
def findSaddlePoint(mat, n):
    # Process all rows one
    # by one
    for i in range(n):
        # Find the minimum element
        # of row i.
        # Also find column index of
        # the minimum element
        min_row = mat[i][0]
        col_ind = 0
        for j in range(1, n):
            if min_row > mat[i][j]:
                min_row = mat[i][j]
                col_ind = j

        # Check if the minimum element
        # of row is also the maximum
        # element of column or not
        k = 0
        for k in range(n):
            # Note that col_ind is fixed
            if min_row < mat[k][col_ind]:
                break
            k += 1

        # If saddle point present in this
        # row then print
        if k == n:
            print("Value of Saddle Point ", min_row)
            return True

    # If Saddle Point found
    return False


def allSaddles(matrix):
    rowmins = []
    rowmaxs = []
    colmins = []
    colmaxs = []

    for i, row in enumerate(matrix):
        m = min(row)
        M = max(row)
        for j, x in enumerate(row):
            if x == m:
                rowmins.append((i, j))
            if x == M:
                rowmaxs.append((i, j))

    t = [list(column) for column in zip(*matrix)]  # transpose of matrix

    for j, col in enumerate(t):
        m = min(col)
        M = max(col)
        for i, x in enumerate(col):
            if x == m:
                colmins.append((i, j))
            if x == M:
                colmaxs.append((i, j))

    return (set(rowmins) & set(colmaxs)) | (set(rowmaxs) & set(colmins))


def detect_local_minima(arr):
    # https://stackoverflow.com/questions/3684484/peak-detection-in-a-2d-array/3689710#3689710
    """
    Takes an array and detects the troughs using the local maximum filter.
    Returns a boolean mask of the troughs (i.e. 1 when
    the pixel's value is the neighborhood maximum, 0 otherwise)
    """
    # define an connected neighborhood
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#generate_binary_structure
    neighborhood = scipy.ndimage.generate_binary_structure(len(arr.shape), 2)
    # apply the local minimum filter; all locations of minimum value
    # in their neighborhood are set to 1
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.filters.html#minimum_filter
    local_min = scipy.ndimage.minimum_filter(arr, footprint=neighborhood) == arr
    # local_min is a mask that contains the peaks we are
    # looking for, but also the background.
    # In order to isolate the peaks we must remove the background from the mask.
    #
    # we create the mask of the background
    background = arr == 0
    #
    # a little technicality: we must erode the background in order to
    # successfully subtract it from local_min, otherwise a line will
    # appear along the background border (artifact of the local minimum filter)
    # http://www.scipy.org/doc/api_docs/SciPy.ndimage.morphology.html#binary_erosion
    eroded_background = scipy.ndimage.binary_erosion(
        background, structure=neighborhood, border_value=1
    )
    #
    # we obtain the final mask, containing only peaks,
    # by removing the background from the local_min mask
    detected_minima = local_min ^ eroded_background
    return np.where(detected_minima)


def extract_data(input):
    data = []

    with open(input, "r") as ifile:
        file = np.loadtxt(ifile)
    cv1_temp = np.array([x[0] for x in file])
    cv2_temp = np.array([x[1] for x in file])
    cv1 = np.unique(cv1_temp)
    cv2 = np.unique(cv2_temp)
    val = np.array([x[2] for x in file])
    free_grid = val.reshape(len(cv1), len(cv2))
    min_pt = detect_local_minima(free_grid)
    return (cv1, cv2, free_grid, min_pt)


def symmetryze(free_grid,cv1,cv2):
    symm_free_grid = np.empty_like(free_grid)
    fes_like=np.array([[0,0,0]])
    for i in range(len(free_grid)):
        for j in range(len(free_grid[i])):
            symm_free_grid[i][j] = (free_grid[i][j] + free_grid[j][i]) / 2
    np.savetxt("data_symm.dat", symm_free_grid, fmt="%6.3f")
    return symm_free_grid


def get_minimum_path(cv1, cv2, free_grid):
    min_free_cv1 = [(min(x), np.where(x == min(x))[0][0]) for x in free_grid]
    min_free_cv2 = [
        (
            min([x[i] for x in free_grid]),
            np.where([x[i] for x in free_grid] == min([x[i] for x in free_grid])),
        )
        for i in range(len(cv2))
    ]
    min_path_cv1 = [(cv1[i], cv2[x[1]], x[0]) for i, x in enumerate(min_free_cv1)]
    min_path_cv2 = [(cv1[x[1]], cv2[i], x[0]) for i, x in enumerate(min_free_cv2)]
    return (min_path_cv1, min_path_cv2)


def minpath(input, x_start, y_start, x_end, y_end):
    # Define your 2D free energy surface as a NumPy array
    # Replace this with your actual data
    free_energy_surface = input

    # Define the starting and ending points on the surface
    start_point = (x_start, y_start)  # Replace with your actual starting point
    end_point = (x_end, y_end)  # Replace with your actual ending point

    # Define the number of intermediate points for the string
    num_intermediate_points = 9

    # Initialize the string as a set of linearly spaced points between start and end points
    path = np.linspace(start_point, end_point, num_intermediate_points + 2)

    # Define the energy function to be minimized (distance along the string)
    def energy_function(x):
        # Calculate the sum of distances between consecutive points
        distances = np.linalg.norm(np.diff(x, axis=0), axis=1)
        return np.sum(distances)

    # Optimize the string using the minimize function from scipy
    # result = minimize(energy_function(free_energy_surface), path, method='L-BFGS-B', options={'disp': True})

    # Extract the optimized string points
    # optimized_string = result.x.reshape(-1, 2)

    # Now `optimized_string` contains the optimized path along the free energy surface
    return optimized_string


def dim_red(cv1, cv2, free_grid, temp, output):
    kbt = 0.0083144621 * temp  # in kj/mol
    prob_grid = np.exp(-free_grid / kbt)
    reduced_prob_2 = [np.trapz(y, cv2) for y in prob_grid]
    reduced_fes_2 = -kbt * np.log(reduced_prob_2)
    offset2 = min(reduced_fes_2)
    reduced_fes_2 -= offset2
    reduced_prob_1 = [np.trapz([y[i] for y in prob_grid], cv1) for i in range(len(cv1))]
    reduced_fes_1 = -kbt * np.log(reduced_prob_1)
    offset1 = min(reduced_fes_1)
    reduced_fes_1 -= offset1
    with open(f"{output}_cv1.dat", "w") as ofile:
        for i, item in enumerate(reduced_fes_1):
            ofile.write(f"{cv1[i]}  {item}\n")
    with open(f"{output}_cv2.dat", "w") as ofile2:
        for i, item in enumerate(reduced_fes_2):
            ofile2.write(f"{cv2[i]}  {item}\n")
    return (reduced_fes_1, reduced_fes_2)


def dim_red_state(cv1, cv2, free_grid, temp):
    kbt = 0.0083144621 * temp  # in kj/mol
    prob_grid = np.exp(-free_grid / kbt)
    reduced_prob_2 = [np.trapz(y, cv2) for y in prob_grid]
    reduced_fes_2 = -kbt * np.log(reduced_prob_2)
    offset2 = min(reduced_fes_2)
    reduced_fes_2 -= offset2
    reduced_prob_1 = [np.trapz([y[i] for y in prob_grid], cv1) for i in range(len(cv1))]
    reduced_fes_1 = -kbt * np.log(reduced_prob_1)
    offset1 = min(reduced_fes_1)
    reduced_fes_1 -= offset1
    return (reduced_fes_1, reduced_fes_2)


def plot2d(x, y, maxz, value, file, labx, laby, cmap, minima, min_pt, symm):
    fig = plt.figure(figsize=(16, 12), dpi=150)
    font = {"family": "Formular", "weight": "normal", "size": 46}
    mpl.rc("font", **font)
    mpl.rcParams["axes.linewidth"] = 3
    mpl.rcParams["lines.linewidth"] = 3
    # lev=int(round(np.max(np.ma.masked_invalid(value))/10,0))
    MAX = int(maxz)
    print(f"grid[0][0]={value[0][0]}")

    # plt.imshow(np.rot90(value),extent=(min(x),max(x),min(y),max(y)))
    # kjmol/plot
    lev = range(0, MAX + 5, 5)
    # print(len(x),len(y))
    value[value > MAX] = MAX + 10
    CLines = plt.contour(
        x,
        y,
        value,
        levels=range(0, MAX, 20),
        vmin=0,
        vmax=MAX,
        linewidths=1.5,
        colors="black",
    )
    """
    if max(x) > max(y):
        max_ax=max(x)
    else:
        max_ax=max(y)
    if min(x) > min(y):
        min_ax=min(y)
    else:
        min_ax=min(x)
    """
    print(max(x), max(y))
    tix = np.linspace(0, round(max(x)), round(round(max(x)) / 0.5) + 1)
    tiy = np.linspace(0, round(max(y)), round(round(max(y)) / 0.5) + 1)
    print(tix,tiy,max(x), max(y))
    # plt.clabel(CLines,levels=range(0,MAX,20), inline=True, fontsize=16,colors='black')
    print(f"Plotting {file} with xlim {plt.xlim()} and ylim {plt.ylim()}")
    plt.contourf(x, y, value, lev, vmin=0, vmax=MAX, cmap=cmap)
    plt.xticks(tix, labels=[str(round(x, 1)) for x in tix])
    plt.yticks(tiy)

    ###

    # kcal/mol plot
    # lev=[x/2 for x in range(0,21,1)]
    # val_kcal=[x/4.184 for x in value]
    # plt.contourf(x, y,val_kcal,lev,vmin=0,vmax=10,cmap=cmap)

    plt.xlabel(labx)
    plt.ylabel(laby)
    # plt.xticks(np.arange(0,max(x),0.5))
    # bounds=[1,2,3,4]
    # cbarkcal
    # cbar=plt.colorbar(label="$\Delta A\ (kcal\ mol^{-1})$",ticks=range(0,11,1))
    # cbar.ax.set_ylim(0,10)
    # cbar kj/mol

    cbar = plt.colorbar(
        label="$\Delta A\ (kJ\ mol^{-1})$", ticks=range(0, MAX + 20, 20)
    )
    cbar.ax.set_ylim(0, MAX)

    # for i in minpath:
    #    plt.scatter(x[i[0]],y[i[1]],color='black')
    # plt.scatter(x[minima[0]],y[minima[1]])
    # plt.scatter(x[maxima[0]],y[maxima[1]],color='red')
    # plt.xlim([0.75,3.25])
    # plt.ylim([0.0,1.5])
    if minima:
        min_crd = []
        with open(minima, "r") as ifile:
            lines = ifile.readlines()
        for line in lines:
            min_crd.append(
                (int(line.split()[1]), float(line.split()[5]), float(line.split()[-1]))
            )
        for pt in min_crd:
            plt.scatter(pt[1], pt[2], color="white")
    good_minima = []
    try:
        for i, item in enumerate(min_pt[0]):
            if value[item, min_pt[1][i]] > 60:
                pass
            else:
                good_minima.append(
                    (x[min_pt[1][i]], y[item], value[item][min_pt[1][i]])
                )
                # plt.scatter(x[min_pt[1][i]],y[item],color='red')
        with open("minima.dat", "w") as ofile:
            for i in good_minima:
                ofile.write(" ".join([f"{j:10.4f}" for j in i]) + "\n")
    except:
        print("No minima found. Skipping")
        pass
    ax = fig.gca()

    # circle1=plt.Circle(xy=(.95,.95), radius=0.15, color='r',fill=False )
    # circle2=plt.Circle(xy=(1.85,1.85), radius=0.15, color='r',fill=False )
    # ax.add_patch(circle1)
    # ax.add_patch(circle2)

    plt.tight_layout()
    if symm:
        plt.savefig("{}_symm.png".format(file), format="png")
    else:
        plt.savefig("{}.png".format(file), format="png")


def plot3d(x, y, value, input, labx, laby):
    with open(input, "r") as ifile:
        file = np.loadtxt(ifile)
    fig = plt.figure(figsize=(16, 10), dpi=150)
    font = {"family": "sans", "weight": "normal", "size": 32}
    mpl.rc("font", **font)
    ax = fig.add_subplot(projection="3d")
    x = np.array([x[0] for x in file]).reshape(len(x), len(y), order="F")
    y = np.array([x[1] for x in file]).reshape(len(x), len(y), order="F")
    # x=np.unique(np.array([x[0] for x in file]))
    # y=np.unique(np.array([x[1] for x in file]))
    ax.plot_surface(x, y, value)
    ax.axes.set_zlim3d(bottom=0, top=50)
    plt.savefig("{}_3D.png".format(input), format="png")


def plotminpath(min_path_cv1, min_path_cv2, file):
    fig = plt.figure(figsize=(16, 10), dpi=150)
    font = {"family": "sans", "weight": "normal", "size": 32}
    mpl.rc("font", **font)
    mpl.rcParams["axes.linewidth"] = 3
    try:
        cv1 = np.loadtxt("fes-rew_square_sparse_cv1.dat")
        # print(cv1)
    except:
        print("No cv1 reweight found")
    try:
        cv2 = np.loadtxt("fes-rew_square_sparse_cv2.dat")
    except:
        print("No cv2 reweight found")
    # lev=int(round(np.max(np.ma.masked_invalid(value))/10,0))
    # lev=range(0,1000,5)
    # plt.imshow(np.rot90(value),extent=(min(x),max(x),min(y),max(y)))
    # plt.plot([x[0] for x in min_path_cv1],[x[2] for x in min_path_cv1],label="CV1")
    plt.plot([x[0] for x in cv1], [x[1] for x in cv1], label="CV1 minpath")
    # plt.plot([x[1] for x in min_path_cv2],[x[2] for x in min_path_cv2],label="CV2")
    plt.plot([x[0] for x in cv2], [x[1] for x in cv2], label="CV2 minpath")
    plt.xlabel("CV")
    plt.ylabel("Free Energy ($kJ\ mol^{-1})$")
    # plt.xlim([0.5,2.0])
    # plt.ylim([-1.,50.0])
    plt.legend()
    # plt.colorbar(label="Free Energy ($kJ\ mol^{-1})$")
    plt.savefig("{}_minpath.png".format(file), format="png")


def plot_reduced(cv1, cv2, min_path_cv1, min_path_cv2, file):
    fig = plt.figure(figsize=(16, 10), dpi=150)
    font = {"family": "sans", "weight": "normal", "size": 32}
    mpl.rc("font", **font)
    plt.plot(cv1, min_path_cv1, label="CV1")
    plt.plot(cv2, min_path_cv2, label="CV2")
    plt.xlabel("CV")
    plt.ylabel("Free Energy ($kJ\ mol^{-1})$")
    # plt.xlim([0.5,4.0])
    plt.ylim([-1.0, 220.0])
    plt.legend()
    # plt.colorbar(label="Free Energy ($kJ\ mol^{-1})$")
    plt.savefig("{}_reduced.png".format(file), format="png")


def plot_convergence(inputs):
    font = {"family": "sans", "weight": "normal", "size": 32}
    mpl.rc("font", **font)
    for input in inputs:
        fig = plt.figure(figsize=(16, 10), dpi=300)
        files = glob.glob(f"{input}_*.dat")
        colors = plt.cm.coolwarm(np.linspace(0, 1, len(files)))
        if len(files) == 0:
            print(f"No files found for input {input}")
            continue
        else:
            for i, file in enumerate(files):
                with open(file, "r") as ifile:
                    data = np.loadtxt(ifile, unpack=True)
                plt.plot(data[0], data[1], label=i, linewidth=2, color=colors[i])
            plt.ylim([-5, 125])
            plt.legend(ncols=5)
            plt.savefig(f"{input}_convergence.png", format="png")


def plot_horizontal_cbar(x, y, maxz, value, cmap):
    MAX = int(maxz)
    lev = range(0, MAX + 5, 5)
    fig = plt.figure(figsize=(16, 12), dpi=150)
    plt.contourf(x, y, value, lev, vmin=0, vmax=MAX, cmap=cmap)

    cbar = plt.colorbar(
        label="$\Delta A\ (kJ\ mol^{-1})$",
        ticks=range(0, MAX + 20, 20),
        orientation="horizontal",
    )
    cbar.ax.set_xlim(0, MAX)
    plt.savefig("cbar.png", format="png")


def main():
    parser = argparse.ArgumentParser(description="Plot data")
    parser.add_argument("--input", dest="input", type=str, help="input data")
    parser.add_argument(
        "--output",
        dest="output",
        type=str,
        help="output file",
        default="2D_reduced_compare.png",
    )
    parser.add_argument(
        "--title", dest="title", type=str, help="plot title", default=""
    )
    parser.add_argument(
        "--labels",
        dest="labels",
        default=None,
        type=str,
        nargs="+",
        help="Plot labels label",
    )
    parser.add_argument(
        "--xlab", dest="xlab", type=str, help="X axis label", default="CV1"
    )
    parser.add_argument(
        "--ylab", dest="ylab", type=str, help="Y axis label", default="CV2"
    )
    parser.add_argument(
        "--minima",
        dest="minima",
        type=str,
        help="file containining minima values",
        default=None,
    )
    parser.add_argument("--limy", dest="limy", type=float, help="y MAX value")
    parser.add_argument("--max", dest="max", type=int, help="z MAX value", default=200)
    parser.add_argument(
        "--symm",
        dest="symm",
        help="symmetrize plot",
        action="store_true",
        default=False,
    )
    parser.add_argument("--p1", dest="p1", help="p1 coords", nargs="+", type=str)
    parser.add_argument("--p2", dest="p2", help="p2 coords", nargs="+", type=str)
    parser.add_argument(
        "--tol",
        dest="tol",
        help="tolearnce to find min around point",
        type=float,
        default=0.1,
    )

    args = parser.parse_args()

    cv1, cv2, free_grid, min_pt = extract_data(args.input)

    # cv1_state,cv2_state,free_grid_state,min_pt_state=extract_data(sys.argv[2])
    # print(cv1[min_pt[1][0]],cv2[min_pt[0][0]],free_grid[min_pt[0][0],min_pt[1][0]])
    file = os.path.splitext(args.input)[0]
    labx = args.xlab
    laby = args.ylab
    MAX = args.max
    minima = args.minima

    cmap_active = "rainbow"
    saddle = allSaddles(free_grid)
    # print(saddle)
    # cmap_active=ListedColormap(np.linspace([0.16862745098, 0.219607843137,1,1],[1, 1,1,1],12)) #blue
    # cmap_active=ListedColormap(np.linspace([1.0, 0.40784313725490196,0,1],[1, 1,1,1],12)) #orange
    # cmap_active=ListedColormap(np.linspace([0.1960784313725490, 0.7686274509803922,0.4980392156862745,1],[1, 1,1,1],12)) #green
    0.1960784313725490

    # 0.16862745098, 0.219607843137,1
    # cmap_active=LinearSegmentedColormap.from_list("mycmap",["#2b38ff","#FFFFFF"])#"#24BC99","#D6BB61","#E18F26","#FF6800"])
    # oldcmap"#000000","#2b38ff","#17d9ff","#f7059b"])

    # cmap_active='Blues_r'
    # path=minpath(free_grid,1,1,2,2)
    min_path_cv1, min_path_cv2 = get_minimum_path(cv1, cv2, free_grid)
    reduced_fes_1, reduced_fes_2 = dim_red(cv1, cv2, free_grid, 300, args.input)
    # reduced_fes_1_state,reduced_fes_2_state=dim_red_state(cv1_state,cv2_state,free_grid_state,300)
    # minima=detect_local_minima(free_grid)
    # maxima=detect_local_minima(-free_grid)
    if args.symm:
        symm_grid = symmetryze(free_grid,cv1,cv2)
        plot2d(
            cv1,
            cv2,
            MAX,
            symm_grid,
            file,
            labx,
            laby,
            cmap_active,
            minima,
            min_pt,
            args.symm,
        )
        args.symm = False
    plot2d(
        cv1,
        cv2,
        MAX,
        free_grid,
        file,
        labx,
        laby,
        cmap_active,
        minima,
        min_pt,
        args.symm,
    )
    plot_horizontal_cbar(cv1, cv2, MAX, free_grid, cmap_active)
    plot3d(cv1, cv2, free_grid, file + ".dat", labx, laby)
    # plotminpath(min_path_cv1,min_path_cv2,file)
    plot_reduced(cv1, cv2, reduced_fes_1, reduced_fes_2, file)
    # plot_reduced_state(cv1_state,cv2_state,reduced_fes_1_state,reduced_fes_2_state,file)
    plot_convergence(
        [
            "fes-rew_square_sparse_cv1_walls",
            "fes-rew_square_sparse_cv2_walls",
            "fes-rew_square_sparse_cv1",
            "fes-rew_square_sparse_cv2",
        ]
    )


if __name__ == "__main__":
    main()
