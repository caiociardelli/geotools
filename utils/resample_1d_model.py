#!/usr/bin/python
# -*- coding: utf-8 -*-

"""
 GeoTools

 Author: Caio Ciardelli, Northwestern University, October 2025

 This program is free software; you can redistribute it and/or modify
 it under the terms of the GNU General Public License as published by
 the Free Software Foundation; either version 3 of the License, or
 (at your option) any later version.

 This program is distributed in the hope that it will be useful,
 but WITHOUT ANY WARRANTY; without even the implied warranty of
 MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 GNU General Public License for more details.

 You should have received a copy of the GNU General Public License along
 with this program; if not, write to the Free Software Foundation, Inc.,
 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.

-----------------------------------------------------------------------------------------------

This script resamples a 1-D Earth model in two stages:
  1. Linear collapse: Remove points in constant-gradient intervals (no error).
  2. Importance-based pruning: Remove least important points until target ratio.

It preserves the first four header lines verbatim and maintains original data line formatting.
Outputs reduction percentage, max interpolation errors, and a comparison plot.

Usage: python resample_1d_model.py --input INPUT --output OUTPUT [options]

Requirements:
- argparse library
- math library
- heapq library
- numpy library
- matplotlib library

-----------------------------------------------------------------------------------------------
"""

import argparse
import math
import heapq
import numpy as np
import matplotlib.pyplot as plt


def parse_arguments():
  # Parse command-line arguments for script configuration
  parser = argparse.ArgumentParser(
      description="Resample a 1-D Earth model with linear collapse and importance-based pruning."
  )
  parser.add_argument("-i", "--input", required=True, help="Input model file path")
  parser.add_argument("-o", "--output", required=True, help="Output resampled model file path")
  parser.add_argument("--depth-col", type=int, default=0, help="Index of depth column (default: 0)")
  parser.add_argument("--check-cols", type=int, nargs="*", default=None,
                      help="Columns to check for resampling (default: all except depth)")
  parser.add_argument("--rtol", type=float, default=1e-6, help="Relative slope tolerance for linear collapse (default: 1e-6)")
  parser.add_argument("--atol", type=float, default=1e-8, help="Absolute slope tolerance for linear collapse (default: 1e-8)")
  parser.add_argument("--no-linear-collapse", action="store_false", dest="enable_linear",
                      help="Disable linear collapse stage")
  parser.add_argument("--target-keep-ratio", type=float, default=0.50,
                      help="Fraction of data rows to keep (default: 0.50)")
  parser.add_argument("--normalize", action="store_true", default=True,
                      help="Normalize errors by column range (default: True)")
  parser.add_argument("--error-agg", choices=["max", "l2", "l1"], default="max",
                      help="Error aggregation method: max, l2, l1 (default: max)")
  parser.add_argument("--preserve-zero", action="store_true", default=True,
                      help="Preserve duplicate depths as segment breaks (default: True)")
  parser.add_argument("--save-fig", action="store_true", default=False,
                      help="Save comparison plot as PNG (default: False)")
  parser.add_argument("--max-cols", type=int, default=3, help="Max columns in plot grid (default: 3)")
  parser.add_argument("--header-lines", type=int, default=4, help="Number of header lines to preserve (default: 4)")
  parser.add_argument("--fortran-d", action="store_true", default=True,
                      help="Accept Fortran-style 'D' in floating-point numbers (default: True)")
  return parser.parse_args()


def to_float(token, fortran_d):
  # Convert token to float, handling Fortran 'D' notation if enabled
  if fortran_d:
    token = token.replace('D', 'E').replace('d', 'e')
  return float(token)


def parse_numeric_line(line, fortran_d):
  # Parse a line into numeric values if possible
  tokens = line.strip().split()
  if not tokens:
    return None
  try:
    return [to_float(t, fortran_d) for t in tokens]
  except ValueError:
    return None


def parse_file(path, header_lines, fortran_d):
  # Read and parse the input file, separating headers and data
  with open(path, 'r', encoding='utf-8') as f:
    lines = f.readlines()
  data_vals, data_raw, is_data_line = [], [], []
  for line in lines[header_lines:]:
    vals = parse_numeric_line(line, fortran_d)
    is_data_line.append(vals is not None)
    if vals:
      data_vals.append(vals)
      data_raw.append(line.rstrip('\n'))
  return lines, is_data_line, data_vals, data_raw


def write_resampled(in_lines, is_data_line, data_raw, keep_mask, out_path, header_lines):
  # Write the resampled model to the output file, preserving headers
  out_lines = [in_lines[i] for i in range(header_lines) if i != 2]
  n_data_lines = sum(1 for k in keep_mask if k)
  out_lines[1] = out_lines[1].split()[0] + f" {n_data_lines}\n"
  for i, raw_line in enumerate(data_raw):
    if keep_mask[i]:
      out_lines.append(raw_line + '\n')
  with open(out_path, 'w', encoding='utf-8') as f:
    f.writelines(out_lines)


def get_segments(z, preserve_zero):
  # Identify segments in depth values, handling zero differences if preserved
  if not z:
    return []
  segments = []
  start = 0
  for i in range(len(z) - 1):
    if preserve_zero and math.isclose(z[i], z[i + 1], rel_tol=0.0, abs_tol=0.0):
      segments.append((start, i))
      start = i + 1
  segments.append((start, len(z) - 1))
  return [(a, b) for a, b in segments if b >= a]


def slopes_equal(a, b, rtol, atol):
  # Check if two slopes are approximately equal within tolerances
  m = max(1.0, abs(a), abs(b))
  return abs(a - b) <= atol + rtol * m


def vectors_equal(a, b, rtol, atol):
  # Check if two vectors of slopes are equal within tolerances
  return len(a) == len(b) and all(slopes_equal(x, y, rtol, atol) for x, y in zip(a, b))


def compute_slopes(x0, x1, y0, y1):
  # Compute slopes between two points for multiple y values
  dx = x1 - x0 if x1 != x0 else 1e-10
  return [(b - a) / dx for a, b in zip(y0, y1)]


def select_columns(ncols, depth_col, check_cols):
  # Select columns to check for resampling, excluding depth if not specified
  if check_cols:
    return sorted(set(check_cols))
  return [j for j in range(ncols) if j != depth_col]


def linear_collapse(data_vals, depth_col, check_cols, rtol, atol, preserve_zero):
  # Perform linear collapse to remove points in constant-gradient intervals
  n = len(data_vals)
  if n == 0:
    return []
  keep = [False] * n
  keep[0] = True
  i = 0
  while i < n - 1:
    xi, xip1 = data_vals[i][depth_col], data_vals[i + 1][depth_col]
    if preserve_zero and math.isclose(xi, xip1, rel_tol=0.0, abs_tol=0.0):
      keep[i + 1] = True
      i += 1
      continue
    y0 = [data_vals[i][j] for j in check_cols]
    y1 = [data_vals[i + 1][j] for j in check_cols]
    m_ref = compute_slopes(xi, xip1, y0, y1)
    j = i + 1
    while j < n - 1:
      xj, xjp1 = data_vals[j][depth_col], data_vals[j + 1][depth_col]
      if preserve_zero and math.isclose(xj, xjp1, rel_tol=0.0, abs_tol=0.0):
        break
      yj = [data_vals[j][k] for k in check_cols]
      yjp1 = [data_vals[j + 1][k] for k in check_cols]
      if not vectors_equal(m_ref, compute_slopes(xj, xjp1, yj, yjp1), rtol, atol):
        break
      j += 1
    keep[j] = True
    i = j
  return keep


def compute_scales(data_vals, cols):
  # Compute scaling factors (ranges) for each column to normalize errors
  mins = [float('inf')] * len(cols)
  maxs = [float('-inf')] * len(cols)
  for row in data_vals:
    for idx, c in enumerate(cols):
      v = row[c]
      mins[idx] = min(mins[idx], v)
      maxs[idx] = max(maxs[idx], v)
  return [hi - lo if hi - lo > 0 and math.isfinite(hi - lo) else 1.0 for lo, hi in zip(mins, maxs)]


def interpolate_vector(t, yL, yR):
  # Linearly interpolate between two vectors
  return [a + t * (b - a) for a, b in zip(yL, yR)]


def point_error(i, iL, iR, z, data_vals, cols, scales, normalize, agg):
  # Compute interpolation error at a point using specified aggregation
  zi, zL, zR = z[i], z[iL], z[iR]
  t = (zi - zL) / (zR - zL) if zR != zL else 0.0
  y = [data_vals[i][c] for c in cols]
  yL = [data_vals[iL][c] for c in cols]
  yR = [data_vals[iR][c] for c in cols]
  yhat = interpolate_vector(t, yL, yR)
  errs = [abs(a - b) / (s if normalize else 1.0) for a, b, s in zip(y, yhat, scales)]
  if agg == 'l2':
    return math.sqrt(sum(e * e for e in errs))
  elif agg == 'l1':
    return sum(errs)
  return max(errs) if errs else 0.0


def importance_prune(data_vals, depth_col, cols, init_keep, target_keep_count, normalize, agg, preserve_zero):
  # Prune points based on importance (error contribution) using a priority queue
  n = len(data_vals)
  if n <= 2:
    return [True] * n
  z = [row[depth_col] for row in data_vals]
  segs = get_segments(z, preserve_zero)
  keep = init_keep[:]
  scales = compute_scales(data_vals, cols) if normalize else [1.0] * len(cols)

  prev_k, next_k = [-1] * n, [-1] * n
  for a, b in segs:
    last = -1
    for i in range(a, b + 1):
      if keep[i]:
        if last != -1:
          next_k[last], prev_k[i] = i, last
        last = i

  heap, version = [], [0] * n

  def push_if_interior(i):
    if i < 0 or prev_k[i] == -1 or next_k[i] == -1:
      return
    imp = point_error(i, prev_k[i], next_k[i], z, data_vals, cols, scales, normalize, agg)
    version[i] += 1
    heapq.heappush(heap, (imp, version[i], i))

  for a, b in segs:
    keep[a] = keep[b] = True
    i = next_k[a]
    while i != -1 and i != b:
      push_if_interior(i)
      i = next_k[i]

  current_kept = sum(keep)
  goal = max(2 * len(segs), min(target_keep_count, n))

  while current_kept > goal and heap:
    imp, ver, i = heapq.heappop(heap)
    if ver != version[i]:
      continue
    L, R = prev_k[i], next_k[i]
    if L == -1 or R == -1:
      continue
    keep[i] = False
    next_k[L], prev_k[R] = R, L
    current_kept -= 1
    push_if_interior(L)
    push_if_interior(R)

  return keep


def evaluate_max_error(data_vals, depth_col, cols, keep, normalize):
  # Evaluate maximum interpolation errors (absolute and normalized) after pruning
  z = [row[depth_col] for row in data_vals]
  segs = get_segments(z, True)
  kept_in_seg = [[i for i in range(a, b + 1) if keep[i]] or [a, b] for a, b in segs]
  for k in kept_in_seg:
    if k[0] != segs[kept_in_seg.index(k)][0]:
      k.insert(0, segs[kept_in_seg.index(k)][0])
    if k[-1] != segs[kept_in_seg.index(k)][1]:
      k.append(segs[kept_in_seg.index(k)][1])
  scales = compute_scales(data_vals, cols) if normalize else [1.0] * len(cols)
  max_abs = [0.0] * len(cols)
  max_norm = [0.0] * len(cols)

  for (a, b), seg_kept in zip(segs, kept_in_seg):
    for j in range(len(seg_kept) - 1):
      L, R = seg_kept[j], seg_kept[j + 1]
      if z[R] == z[L]:
        continue
      for i in range(L + 1, R):
        if keep[i]:
          continue
        t = (z[i] - z[L]) / (z[R] - z[L]) if z[R] != z[L] else 0.0
        for idx, c in enumerate(cols):
          yi = data_vals[i][c]
          yL = data_vals[L][c]
          yR = data_vals[R][c]
          yhat = yL + t * (yR - yL)
          eabs = abs(yi - yhat)
          enorm = eabs / scales[idx]
          max_abs[idx] = max(max_abs[idx], eabs)
          max_norm[idx] = max(max_norm[idx], enorm)

  return max_abs, max_norm


def plot_comparison(data_vals, keep_mask, depth_col, check_cols, out_png, max_cols):
  # Generate and optionally save a comparison plot of original vs. resampled data
  arr = np.asarray(data_vals, dtype=float)
  z, kept = arr[:, depth_col], np.asarray(keep_mask, dtype=bool)
  z_k = z[kept]
  nplots, ncols = len(check_cols), min(max_cols, max(1, len(check_cols)))
  nrows = (nplots + ncols - 1) // ncols
  fig, axes = plt.subplots(nrows, ncols, figsize=(5 * ncols, 8 * nrows), squeeze=False)
  axes = axes.ravel()

  for ip, col in enumerate(check_cols):
    ax = axes[ip]
    y, y_k = arr[:, col], arr[kept, col]
    ax.plot(y, z, lw=1.3, label="original")
    ax.plot(y_k, z_k, lw=2.0, label="resampled")
    ax.set_xlabel(f'Column {col}')
    ax.set_ylabel(f'Depth (col {depth_col})')
    ax.invert_yaxis()
    ax.grid(True, ls=':', alpha=0.5)
    ax.set_title(f'Column {col}')
    ax.legend()

  for j in range(ip + 1, len(axes)):
    fig.delaxes(axes[j])

  fig.tight_layout()
  if out_png:
    fig.savefig(out_png, dpi=200, bbox_inches='tight')
  plt.show()


def main():
  # Main function to orchestrate resampling process
  args = parse_arguments()
  lines, is_data_line, data_vals, data_raw = parse_file(args.input, args.header_lines, args.fortran_d)
  if not data_vals:
    write_resampled(lines, is_data_line, data_raw, [], args.output, args.header_lines)
    return

  check_cols = select_columns(len(data_vals[0]), args.depth_col, args.check_cols)
  keep = linear_collapse(data_vals, args.depth_col, check_cols, args.rtol, args.atol, args.preserve_zero) \
      if args.enable_linear else [True] * len(data_vals)
  target_keep = max(2, int(math.ceil(args.target_keep_ratio * len(data_vals))))
  keep_final = importance_prune(data_vals, args.depth_col, check_cols, keep, target_keep,
                                args.normalize, args.error_agg, args.preserve_zero)

  print(f"Debug: Number of kept points: {sum(keep_final)} out of {len(data_vals)}")
  write_resampled(lines, is_data_line, data_raw, keep_final, args.output, args.header_lines)

  n_in, n_out = len(data_vals), sum(keep_final)
  reduction = 100.0 * (1.0 - n_out / max(1, n_in))
  print(f"Resampled: kept {n_out}/{n_in} data rows ({reduction:.1f}% reduction)")

  max_abs, max_norm = evaluate_max_error(data_vals, args.depth_col, check_cols, keep_final, args.normalize)
  print("Max absolute error per column:", ", ".join(f"c{c}:{e:.4g}" for c, e in zip(check_cols, max_abs)))
  if args.normalize:
    print("Max normalized error per column:", ", ".join(f"c{c}:{e:.4g}" for c, e in zip(check_cols, max_norm)))

  out_png = f"{args.output}.png" if args.save_fig else None
  plot_comparison(data_vals, keep_final, args.depth_col, check_cols, out_png, args.max_cols)


if __name__ == "__main__":
  main()