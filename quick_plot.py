#!/usr/bin/env python3

import argparse
from plotting_functions import save_quick_plot


def main():
    parser = argparse.ArgumentParser(description='Quick observed-vs-predicted trench motion plot (no GMT).')
    parser.add_argument('--predicted', required=True, help='Predicted file path (lat lon vt_cm_yr azim).')
    parser.add_argument('--observed', required=True, help='Observed file path (tnew.*.dat).')
    parser.add_argument('--output', required=True, help='Output PNG path.')
    parser.add_argument('--title', default='Quick Trench Motion Check', help='Figure title.')
    parser.add_argument('--tol', type=float, default=0.1, help='Coordinate rounding tolerance in degrees.')
    parser.add_argument('--neutral', type=float, default=0.3, help='Neutral sign threshold in cm/yr.')
    args = parser.parse_args()

    save_quick_plot(
        predicted_path=args.predicted,
        observed_path=args.observed,
        output_path=args.output,
        title=args.title,
        tol=args.tol,
        neutral=args.neutral,
    )


if __name__ == '__main__':
    main()
