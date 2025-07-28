import matplotlib.pyplot as plt
import numpy as np
import plumed
import argparse

argparser = argparse.ArgumentParser(description='Plot CV1 vs CV2 from COLVAR file')
argparser.add_argument('--colvar', type=str, help='Path to the COLVAR file', default='COLVAR')
argparser.add_argument('--stride', type=int, help='Stride for reading frames', default=1)
argparser.add_argument('--col1', type=str, help='Column name for CV1', default='CV1')
argparser.add_argument('--col2', type=str, help='Column name for CV2', default='CV2')
argparser.add_argument('--output', type=str, help='Output file for the plot', default='scatter_CV1_CV2.png')
args = argparser.parse_args()

colvar=plumed.read_as_pandas(args.colvar)

fig= plt.figure(figsize=(8, 6))
ax = fig.add_subplot(111)
ax.scatter(colvar[args.col1][::args.stride], colvar[args.col2][::args.stride], s=1, alpha=0.5)
ax.set_xlabel(args.col1)
ax.set_ylabel(args.col2)
ax.set_title(f'Scatter plot of {args.col1} vs {args.col2}')
ax.grid(True)
plt.savefig(args.output)
plt.close(fig)
print(f"Scatter plot saved to {args.output}")
print(f"Total frames plotted: {len(colvar[args.col1]) // args.stride}")
