#!/usr/bin/env python3
"""
Bootstrap error analysis for PLUMED COLVAR files from metadynamics runs.

Reads a PLUMED COLVAR file, performs Boltzmann reweighting using the bias
potential, and estimates the free energy surface (FES) with bootstrap error bars.

Usage examples:
  # Basic 1D FES with default parameters
  python bootstrap_error.py --input COLVAR --cv-col cv1 --bias-col metad.bias

  # Custom bootstrap parameters
  python bootstrap_error.py --input COLVAR --cv-col cv1 --bias-col metad.bias \
      --nsamples 500 --npoints 10000 --nbins 80 --temp 350

  # Block bootstrap with equilibration skip and stride
  python bootstrap_error.py --input COLVAR --cv-col cv1 --bias-col metad.bias \
      --block-size 50 --skip 1000 --stride 5

  # 2D FES
  python bootstrap_error.py --input COLVAR --cv-col cv1 cv2 --bias-col metad.bias \
      --nbins 40 --output fes2d_bootstrap.dat

  # No reweighting (e.g. unbiased simulation or COLVAR without bias column)
  python bootstrap_error.py --input COLVAR --cv-col cv1 --no-reweight
"""

import numpy as np
import argparse
import sys


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

def parse_colvar(filename):
    """
    Parse a PLUMED COLVAR file.

    Returns
    -------
    fields : list of str
        Column names extracted from the ``#! FIELDS`` header line.
    data : np.ndarray, shape (N, ncols)
        Numerical data as a 2-D array.
    """
    fields = None
    data_lines = []

    with open(filename, "r") as fh:
        for line in fh:
            stripped = line.strip()
            if not stripped:
                continue
            # Header lines
            if stripped.startswith("#"):
                if "FIELDS" in stripped:
                    # e.g.  #! FIELDS time cv1 metad.bias metad.rbias ...
                    # Also handles  ### FIELDS ...  or  ###! FIELDS ...
                    idx = stripped.index("FIELDS")
                    fields = stripped[idx + len("FIELDS"):].split()
                continue
            # Data line
            data_lines.append([float(x) for x in stripped.split()])

    if fields is None:
        # Fallback: generate generic column names
        ncols = len(data_lines[0]) if data_lines else 0
        fields = [f"col{i}" for i in range(ncols)]
        print(f"WARNING: No FIELDS header found. Using generic names: {fields}",
              file=sys.stderr)

    data = np.array(data_lines)
    return fields, data


# ---------------------------------------------------------------------------
# Bootstrap helpers
# ---------------------------------------------------------------------------

def _weighted_histogram_1d(cv, weights, bins):
    """Weighted histogram for 1D CV."""
    hist, bin_edges = np.histogram(cv, bins=bins, weights=weights, density=False)
    centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    return centres, hist


def _weighted_histogram_2d(cv1, cv2, weights, bins):
    """Weighted histogram for 2D CV."""
    hist, xedges, yedges = np.histogram2d(
        cv1, cv2, bins=bins, weights=weights, density=False
    )
    xcentres = 0.5 * (xedges[:-1] + xedges[1:])
    ycentres = 0.5 * (yedges[:-1] + yedges[1:])
    return xcentres, ycentres, hist


def _hist_to_fes(hist, kbt):
    """Convert a histogram to a free energy surface (set min to 0)."""
    fes = np.full_like(hist, np.nan, dtype=float)
    mask = hist > 0
    fes[mask] = -kbt * np.log(hist[mask])
    fes[mask] -= np.nanmin(fes[mask])  # shift minimum to zero
    return fes


def _block_resample_indices(n_total, n_points, block_size, rng):
    """
    Draw *n_points* data-point indices using block bootstrap.

    Whole blocks of *block_size* consecutive indices are drawn with
    replacement until at least *n_points* indices are collected; the
    result is then trimmed to exactly *n_points*.
    """
    n_blocks_needed = int(np.ceil(n_points / block_size))
    # Possible block starting positions
    max_start = n_total - block_size
    if max_start < 0:
        raise ValueError(
            f"block_size ({block_size}) is larger than the dataset ({n_total})"
        )
    starts = rng.integers(0, max_start + 1, size=n_blocks_needed)
    indices = np.concatenate(
        [np.arange(s, s + block_size) for s in starts]
    )
    return indices[:n_points]


# ---------------------------------------------------------------------------
# Main bootstrap engine
# ---------------------------------------------------------------------------

def bootstrap_fes_1d(cv, weights, n_samples, n_points, n_bins, kbt,
                     block_size, rng, cv_min=None, cv_max=None):
    """
    Bootstrap estimate of the 1D FES with error bars.

    Returns
    -------
    centres : np.ndarray
        Bin centres.
    fes_mean : np.ndarray
        Mean FES across bootstrap samples.
    fes_std : np.ndarray
        Standard deviation (error) of the FES.
    """
    n_total = len(cv)
    all_fes = []

    # Fixed bin edges so all samples share the same grid
    lo = cv_min if cv_min is not None else cv.min()
    hi = cv_max if cv_max is not None else cv.max()
    bin_edges = np.linspace(lo, hi, n_bins + 1)

    for _ in range(n_samples):
        if block_size > 1:
            idx = _block_resample_indices(n_total, n_points, block_size, rng)
        else:
            idx = rng.integers(0, n_total, size=n_points)

        hist, _ = np.histogram(cv[idx], bins=bin_edges,
                               weights=weights[idx], density=False)
        fes = _hist_to_fes(hist, kbt)
        all_fes.append(fes)

    all_fes = np.array(all_fes)
    centres = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    fes_mean = np.nanmean(all_fes, axis=0)
    fes_std = np.nanstd(all_fes, axis=0)

    # Shift mean to zero at the minimum
    valid = np.isfinite(fes_mean)
    if valid.any():
        fes_mean[valid] -= np.nanmin(fes_mean[valid])

    return centres, fes_mean, fes_std


def bootstrap_fes_2d(cv1, cv2, weights, n_samples, n_points, n_bins, kbt,
                     block_size, rng, cv_min=None, cv_max=None):
    """
    Bootstrap estimate of the 2D FES with error bars.

    Parameters
    ----------
    cv_min, cv_max : list of float or None
        [xmin, ymin] and [xmax, ymax] for the bin edges.

    Returns
    -------
    xcentres, ycentres : np.ndarray
        Bin centres along each axis.
    fes_mean : np.ndarray, shape (nbins, nbins)
    fes_std  : np.ndarray, shape (nbins, nbins)
    """
    n_total = len(cv1)
    all_fes = []

    xlo = cv_min[0] if cv_min is not None else cv1.min()
    xhi = cv_max[0] if cv_max is not None else cv1.max()
    ylo = cv_min[1] if cv_min is not None else cv2.min()
    yhi = cv_max[1] if cv_max is not None else cv2.max()
    xedges = np.linspace(xlo, xhi, n_bins + 1)
    yedges = np.linspace(ylo, yhi, n_bins + 1)

    for _ in range(n_samples):
        if block_size > 1:
            idx = _block_resample_indices(n_total, n_points, block_size, rng)
        else:
            idx = rng.integers(0, n_total, size=n_points)

        hist, _, _ = np.histogram2d(cv1[idx], cv2[idx],
                                    bins=[xedges, yedges],
                                    weights=weights[idx], density=False)
        fes = _hist_to_fes(hist, kbt)
        all_fes.append(fes)

    all_fes = np.array(all_fes)
    xcentres = 0.5 * (xedges[:-1] + xedges[1:])
    ycentres = 0.5 * (yedges[:-1] + yedges[1:])
    fes_mean = np.nanmean(all_fes, axis=0)
    fes_std = np.nanstd(all_fes, axis=0)

    valid = np.isfinite(fes_mean)
    if valid.any():
        fes_mean[valid] -= np.nanmin(fes_mean[valid])

    return xcentres, ycentres, fes_mean, fes_std


# ---------------------------------------------------------------------------
# I/O
# ---------------------------------------------------------------------------

def write_fes_1d(filename, centres, fes_mean, fes_std, cv_name):
    """Write 1D FES to a PLUMED-style output file."""
    with open(filename, "w") as fh:
        fh.write(f"#! FIELDS {cv_name} fes fes_error\n")
        fh.write(f"#! SET min_{cv_name} {centres[0]:.6f}\n")
        fh.write(f"#! SET max_{cv_name} {centres[-1]:.6f}\n")
        fh.write(f"#! SET nbins_{cv_name} {len(centres)}\n")
        for c, m, s in zip(centres, fes_mean, fes_std):
            if np.isfinite(m):
                fh.write(f" {c:15.8f} {m:15.8f} {s:15.8f}\n")
            else:
                fh.write(f" {c:15.8f} {'inf':>15s} {'inf':>15s}\n")
    print(f"1D FES written to {filename}")


def write_fes_2d(filename, xc, yc, fes_mean, fes_std, cv_names):
    """Write 2D FES to a PLUMED-style grid file."""
    with open(filename, "w") as fh:
        fh.write(f"#! FIELDS {cv_names[0]} {cv_names[1]} fes fes_error\n")
        fh.write(f"#! SET min_{cv_names[0]} {xc[0]:.6f}\n")
        fh.write(f"#! SET max_{cv_names[0]} {xc[-1]:.6f}\n")
        fh.write(f"#! SET nbins_{cv_names[0]} {len(xc)}\n")
        fh.write(f"#! SET min_{cv_names[1]} {yc[0]:.6f}\n")
        fh.write(f"#! SET max_{cv_names[1]} {yc[-1]:.6f}\n")
        fh.write(f"#! SET nbins_{cv_names[1]} {len(yc)}\n")
        for i, x in enumerate(xc):
            for j, y in enumerate(yc):
                m = fes_mean[i, j]
                s = fes_std[i, j]
                if np.isfinite(m):
                    fh.write(f" {x:15.8f} {y:15.8f} {m:15.8f} {s:15.8f}\n")
                else:
                    fh.write(f" {x:15.8f} {y:15.8f} {'inf':>15s} {'inf':>15s}\n")
            fh.write("\n")  # blank line between x-blocks (gnuplot convention)
    print(f"2D FES written to {filename}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def build_parser():
    parser = argparse.ArgumentParser(
        description="Bootstrap error analysis for PLUMED COLVAR files "
                    "(metadynamics / OPES).",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=__doc__,
    )
    # Required
    parser.add_argument("--input", "-i", required=True,
                        help="Path to the PLUMED COLVAR file.")
    # Column selection
    parser.add_argument("--cv-col", nargs="+", default=None,
                        help="Name(s) of the CV column(s) in the FIELDS "
                             "header.  1 name → 1D FES, 2 names → 2D FES.  "
                             "Default: auto-detect first non-time, non-bias "
                             "column.")
    parser.add_argument("--bias-col", nargs="+", default=["metad.bias"],
                        help="Name(s) of the bias-potential column(s). "
                             "Multiple columns are summed together "
                             "(e.g. --bias-col metad.bias wall.bias). "
                             "Default: metad.bias.")
    parser.add_argument("--no-reweight", action="store_true",
                        help="Skip Boltzmann reweighting (uniform weights). "
                             "Useful for unbiased simulations or if the "
                             "COLVAR has no bias column(s).")
    # Bootstrap parameters
    parser.add_argument("--nsamples", "-n", type=int, default=200,
                        help="Number of bootstrap resamples (default: 200).")
    parser.add_argument("--npoints", "-p", type=int, default=None,
                        help="Number of data points drawn per bootstrap "
                             "sample.  Default: total length of the "
                             "(filtered) dataset.")
    parser.add_argument("--block-size", "-b", type=int, default=1,
                        help="Block size for block bootstrap to account for "
                             "autocorrelation (default: 1 = standard "
                             "bootstrap).")
    # Histogram / FES
    parser.add_argument("--nbins", type=int, default=50,
                        help="Number of bins for the FES histogram "
                             "(default: 50).")
    parser.add_argument("--cv-min", nargs="+", type=float, default=None,
                        help="Minimum CV value(s) for the FES grid. "
                             "One value per CV dimension. "
                             "Default: data minimum.")
    parser.add_argument("--cv-max", nargs="+", type=float, default=None,
                        help="Maximum CV value(s) for the FES grid. "
                             "One value per CV dimension. "
                             "Default: data maximum.")
    parser.add_argument("--temp", "-T", type=float, default=300.0,
                        help="Temperature in Kelvin (default: 300).")
    # Data filtering
    parser.add_argument("--stride", type=int, default=1,
                        help="Read every Nth frame (default: 1, read all).")
    parser.add_argument("--skip", type=int, default=0,
                        help="Number of initial frames to discard as "
                             "equilibration (default: 0).")
    # Misc
    parser.add_argument("--output", "-o", default="fes_bootstrap.dat",
                        help="Output filename (default: fes_bootstrap.dat).")
    parser.add_argument("--seed", type=int, default=None,
                        help="Random seed for reproducibility.")
    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    # ---- Constants ----
    kB = 8.314462618e-3  # kJ/(mol·K)
    kbt = kB * args.temp

    # ---- Parse COLVAR ----
    fields, data = parse_colvar(args.input)
    print(f"Read {data.shape[0]} frames, {data.shape[1]} columns")
    print(f"  FIELDS: {' '.join(fields)}")

    # ---- Apply skip + stride ----
    data = data[args.skip::args.stride]
    n_frames = data.shape[0]
    print(f"After skip={args.skip}, stride={args.stride}: {n_frames} frames")

    if n_frames == 0:
        print("ERROR: No data remaining after skip/stride filtering.",
              file=sys.stderr)
        sys.exit(1)

    # ---- Identify columns ----
    def col_index(name):
        if name in fields:
            return fields.index(name)
        raise ValueError(
            f"Column '{name}' not found in FIELDS: {fields}"
        )

    # Bias column names for skip-set
    bias_col_names = set(args.bias_col)

    # CV columns
    if args.cv_col is not None:
        cv_names = args.cv_col
    else:
        # Auto-detect: first column that is not 'time' and not the bias col
        skip_names = {"time"} | bias_col_names
        cv_names = [f for f in fields if f not in skip_names][:1]
        if not cv_names:
            print("ERROR: Cannot auto-detect CV column. "
                  "Use --cv-col explicitly.", file=sys.stderr)
            sys.exit(1)
        print(f"  Auto-detected CV column: {cv_names[0]}")

    cv_indices = [col_index(n) for n in cv_names]
    ndim = len(cv_names)
    if ndim > 2:
        print("ERROR: Only 1D or 2D FES is supported.", file=sys.stderr)
        sys.exit(1)

    # ---- Weights ----
    if args.no_reweight:
        weights = np.ones(n_frames)
        print("  Reweighting: DISABLED (uniform weights)")
    else:
        try:
            # Sum all bias columns into a total bias
            total_bias = np.zeros(n_frames)
            for bcol in args.bias_col:
                bias_idx = col_index(bcol)
                total_bias += data[:, bias_idx]
            # Numerical stabilisation: subtract max before exp
            weights = np.exp((total_bias - total_bias.max()) / kbt)
            bias_str = " + ".join(args.bias_col)
            print(f"  Reweighting with column(s): {bias_str} "
                  f"(kBT = {kbt:.4f} kJ/mol)")
        except ValueError as e:
            print(f"WARNING: {e}. "
                  f"Proceeding with uniform weights.", file=sys.stderr)
            weights = np.ones(n_frames)

    # ---- Bootstrap parameters ----
    n_points = args.npoints if args.npoints is not None else n_frames
    rng = np.random.default_rng(args.seed)

    print(f"\nBootstrap parameters:")
    print(f"  nsamples   = {args.nsamples}")
    print(f"  npoints    = {n_points}")
    print(f"  block_size = {args.block_size}")
    print(f"  nbins      = {args.nbins}")
    print(f"  seed       = {args.seed}")

    # ---- Validate cv-min / cv-max dimensions ----
    for label, vals in [("--cv-min", args.cv_min), ("--cv-max", args.cv_max)]:
        if vals is not None and len(vals) != ndim:
            print(f"ERROR: {label} expects {ndim} value(s) (one per CV), "
                  f"got {len(vals)}.", file=sys.stderr)
            sys.exit(1)

    # ---- Run bootstrap ----
    if ndim == 1:
        cv = data[:, cv_indices[0]]
        cmin = args.cv_min[0] if args.cv_min is not None else None
        cmax = args.cv_max[0] if args.cv_max is not None else None
        centres, fes_mean, fes_std = bootstrap_fes_1d(
            cv, weights, args.nsamples, n_points, args.nbins, kbt,
            args.block_size, rng, cv_min=cmin, cv_max=cmax
        )
        write_fes_1d(args.output, centres, fes_mean, fes_std, cv_names[0])
    else:
        cv1 = data[:, cv_indices[0]]
        cv2 = data[:, cv_indices[1]]
        xc, yc, fes_mean, fes_std = bootstrap_fes_2d(
            cv1, cv2, weights, args.nsamples, n_points, args.nbins, kbt,
            args.block_size, rng, cv_min=args.cv_min, cv_max=args.cv_max
        )
        write_fes_2d(args.output, xc, yc, fes_mean, fes_std, cv_names)

    print("Done.")


if __name__ == "__main__":
    main()
