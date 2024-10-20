import argparse
import sys

import numpy as np
import mrcfile

from slabify.utils import measure_thickness
from slabify import slabify


def parse_args() -> argparse.Namespace:
    """
    Parses command-line arguments for the lamella boundary mask generator.

    Returns:
        argparse.Namespace: Parsed command-line arguments.
    """
    parser = argparse.ArgumentParser(
        description="Lamella boundary mask generator.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    general_opts = parser.add_argument_group(title="General options")
    general_opts.add_argument(
        "-i",
        "--input",
        type=str,
        help="Path to input tomogram in MRC format.",
        required=True,
    )
    general_opts.add_argument(
        "-o",
        "--output",
        type=str,
        help="Path to output mask in MRC format.",
        required=True,
    )
    general_opts.add_argument(
        "--output-masked", "--output_masked", type=str, help="Path to output masked volume in MRC format."
    )
    general_opts.add_argument(
        "--border", type=int, default=0, help="Voxels to exclude from the border in XY."
    )
    general_opts.add_argument(
        "--offset",
        type=int,
        default=0,
        help="Voxels to offset from the border along Z, i.e. to make a thicker (positive value) or thinner (negative value) slab mask.",
    )
    general_opts.add_argument(
        "--angpix", type=float, help="Pixel size (in Angstroms) for output mask."
    )
    general_opts.add_argument(
        "--measure",
        default=False,
        action="store_true",
        help="Measure thickness of the lamella slab mask at the center of the tomogram.",
    )

    manual_opts = parser.add_argument_group(
        title="Manual masking options (using model with points to define the boundary)"
    )
    manual_opts.add_argument(
        "--points",
        type=str,
        help="Path to IMOD model containing boundary control points. Both .mod and .txt formats are accepted, but if a .mod file is supplied, IMOD is expected to be present in your $PATH.",
    )

    auto_opts = parser.add_argument_group(
        title="Automatic masking options (using local analysis of variance)"
    )
    auto_opts.add_argument(
        "--n-samples", "--n_samples",
        type=int,
        default=50000,
        help="Number of points to sample at random within the tomogram.",
    )
    auto_opts.add_argument(
        "--boxsize",
        type=int,
        default=32,
        help="Box size (in pixels) to analyze variance around each sampled point.",
    )
    auto_opts.add_argument(
        "--z-min", "--z_min",
        type=int,
        default=1,
        help="Minimum Z slice to consider in variance analysis, starting from 1. If not specified, will use the first slice.",
    )
    auto_opts.add_argument(
        "--z-max", "--z_max",
        type=int,
        default=None,
        help="Maximum Z slice to consider in variance analysis. If not specified, will use the last slice.",
    )
    auto_opts.add_argument(
        "--iterations",
        type=int,
        default=3,
        help="Number of iterations for automatic fitting of top and bottom planes defining the lamella boundaries.",
    )
    auto_opts.add_argument(
        "--simple",
        default=False,
        action="store_true",
        help="Fit a single plane through the center of the lamella and assume fixed lamella thickness.",
    )
    auto_opts.add_argument(
        "--thickness",
        type=int,
        help="Total lamella thickness (in pixels) when doing simple automatic boundary masking. By default the program will use half the tomogram size in Z.",
    )
    auto_opts.add_argument(
        "--percentile",
        type=float,
        default=95,
        help="Percentile of highest variance locations to select for fitting.",
    )
    auto_opts.add_argument("--seed", type=int, default=4056, help="Random seed.")

    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)

    return parser.parse_args()


def main():
    """
    Main function to generate the lamella boundary mask based on the parsed arguments.
    """
    args = parse_args()
    input = args.input
    output = args.output
    output_masked = args.output_masked
    points = args.points
    border = args.border
    offset = args.offset
    n_samples = args.n_samples
    boxsize = args.boxsize
    z_min = args.z_min
    z_max = args.z_max
    iterations = args.iterations
    simple = args.simple
    measure = args.measure
    thickness = args.thickness
    percentile = args.percentile
    seed = args.seed
    angpix = args.angpix

    inmrc = mrcfile.open(input, permissive=True)
    tomo = np.array(inmrc.data)

    bmask = slabify(
        tomo=tomo,
        points=points,
        border=border,
        offset=offset,
        n_samples=n_samples,
        boxsize=boxsize,
        z_min=z_min,
        z_max=z_max,
        iterations=iterations,
        simple=simple,
        thickness=thickness,
        percentile=percentile,
        seed=seed,
    )

    if not angpix:
        print(
            f"Warning: --angpix was not provided. Using the pixel size from the input tomogram's header: {inmrc.voxel_size.x:.3f}"
        )
        angpix = inmrc.voxel_size.x

    if measure:
        measure_thickness(mask=bmask, z_offset=offset, angpix=angpix, verbose=True)

    # Save output mask:
    with mrcfile.new(output, overwrite=True) as outmrc:
        outmrc.set_data(bmask)
        outmrc.voxel_size = angpix

    if output_masked:
        # Save output masked volume:
        with mrcfile.new(output_masked, overwrite=True) as outmrc:
            outmrc.set_data(tomo * bmask)
            outmrc.voxel_size = angpix
