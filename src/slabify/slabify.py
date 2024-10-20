import os
import subprocess

import numpy as np

from slabify.boundary_mask_auto import create_boundary_mask_auto
from slabify.stopgap_tm import sg_tm_create_boundary_mask
from slabify.utils import apply_mask_border


def slabify(
    tomo: np.ndarray,
    points: str = None,
    border: int = 0,
    offset: float = 0.0,
    n_samples: int = 50000,
    boxsize: int = 32,
    z_min: int = 1,
    z_max: int = None,
    iterations: int = 3,
    simple: bool = False,
    thickness: int = None,
    percentile: float = 95,
    seed: int = 4056,
):
    # Check if a points file was provided:
    if isinstance(points, str):
        points_basename, points_ext = os.path.splitext(points)

        # If the points file is provided in binary format (.mod) we call IMOD's model2point to convert it. IMOD must be already loaded for this to work:
        if points_ext == ".mod":
            subprocess.run(["model2point", points, points_basename + ".txt"])
            points = points_basename + ".txt"

        boundary = np.loadtxt(points)

        # Create boundary mask using given points:
        bmask = sg_tm_create_boundary_mask(
            mask_size=tomo.shape, boundary=boundary, z_offset=offset
        )

    else:
        # Create a boundary mask automatically:
        bmask = create_boundary_mask_auto(
            tomo=tomo,
            N=n_samples,
            boxsize=boxsize,
            z_min=z_min,
            z_max=z_max,
            z_offset=offset,
            simple=simple,
            thickness=thickness,
            iterations=iterations,
            percentile=percentile,
            seed=seed,
        )

    # Mask out some voxels away from the border in XY:
    bmask = apply_mask_border(mask=bmask, xy_border=border)

    return bmask
