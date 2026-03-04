import dpdata as dp
import glob
import os
import subprocess
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import regex as re
from natsort import natsorted, ns
import json
import argparse
import plumed


def count_2D_dataset(path,grid_1_size,grid_2_size,cv1,cv2):
    "Count the number of points in a 2D dataset clustering them in tiles of size grid_1_size x grid_2_size"
    with open(path, 'r') as f:
        data=plumed.read_as_pandas(f)
    x_tiles = np.arange(data[cv1].min(),data[cv1].max()+grid_1_size, grid_1_size)
    y_tiles = np.arange(data[cv2].min(),data[cv2].max()+grid_2_size, grid_2_size)

    counts, _, _ = np.histogram2d(data[cv1], data[cv2], bins=[x_tiles, y_tiles])
    return counts,x_tiles, y_tiles

def plot_2D_dataset(path, grid_1_size, grid_2_size, output_path,cv1, cv2, MAX=5000):
    fig = plt.figure(figsize=(16, 12), dpi=150)
    font = {"family": "Formular", "weight": "normal", "size": 46}
    mpl.rc("font", **font)
    mpl.rcParams["axes.linewidth"] = 3
    mpl.rcParams["lines.linewidth"] = 3
    "Plot a 2D dataset clustering them in tiles of size grid_1_size x grid_2_size with 3D bars"
    counts,x,y = count_2D_dataset(path, grid_1_size, grid_2_size,cv1,cv2)

    plt.imshow(counts.T, origin='lower', cmap='rainbow', aspect='auto',extent=[x[0], x[-1], y[0], y[-1]],
               vmax= MAX  , vmin=000, interpolation='nearest')
    
    #plt.contourf(x[:-1], y[:-1], counts.T, range(0,5000,10), vmin=0, vmax=5000, cmap='rainbow', alpha=0.8)
    plt.colorbar(label='Counts')
    plt.xlabel(f"{cv1} (binned by {grid_1_size})")
    plt.ylabel(f"{cv2} (binned by {grid_2_size})")
    plt.title('2D Dataset Counts')
    plt.grid(True)
    plt.xticks(np.arange(-0.5,2.5,0.5))
    plt.savefig(output_path)
    plt.close()

def main():
    parser = argparse.ArgumentParser(description="Process and plot 2D datasets.")
    parser.add_argument('input_path', type=str, help='Path to the input dataset file')
    parser.add_argument('--grid_1_size', type=float, default=0.1, help='Size of the first grid dimension')
    parser.add_argument('--grid_2_size', type=float, default=0.1, help='Size of the second grid dimension')
    parser.add_argument('--output_path', type=str, default='output_plot.png', help='Path to save the output plot')
    parser.add_argument('--cv1', type=str, default='CV1', help='Name of the first CV')
    parser.add_argument('--cv2', type=str, default='CV2', help='Name of the second CV')
    parser.add_argument('--max', type=int, default=5000, help='Maximum value for color scaling in the plot')

    args = parser.parse_args()

    plot_2D_dataset(args.input_path, args.grid_1_size, args.grid_2_size, args.output_path,args.cv1, args.cv2,args.max)
    
if __name__ == "__main__":
    main()
    
    