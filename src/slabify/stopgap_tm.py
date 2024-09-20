import numpy as np
from slabify.utils import affine_fit, make_boundary_mask, plane_equation_Z_from_XY


def sg_tm_create_boundary_mask(
    mask_size: np.ndarray, boundary: np.ndarray, z_offset: float = 0.0
) -> np.ndarray:
    """
    Create a boundary mask for a tomogram assuming a slab geometry.
    Reimplementation in Python of the MATLAB function:
    https://github.com/wan-lab-vanderbilt/STOPGAP/blob/master/sg_toolbox/sg_tm_create_boundary_mask.m

    NOTE: the Python version of this function assumes volume dimensions in the mrcfile convention: [Z, Y, X].

    Args:
        mask_size (np.ndarray): The size of the tomogram to be masked in format (Z, Y, X).
        boundary (np.ndarray): A plain-text file containing points that define the slab geometry.
        xy_border (int, optional): Number of voxels to mask on each edge. Defaults to 0.

    Returns:
        mask (np.ndarray): The boundary mask.

    Original author: William Wan
    Converted to Python by: Ricardo D. Righetto (with some Copilot help)

    """

    dims = mask_size

    # Parse top and bottom points
    n_points = boundary.shape[0]
    n_sets = n_points // 4
    if n_points % 4 != 0:
        raise ValueError("Input boundary coordinates must be supplied in sets of 4")

    top = np.zeros((n_sets * 2, 3))
    bottom = np.zeros((n_sets * 2, 3))

    # Parse point list:
    for i in range(n_sets):
        top[2 * i] = boundary[4 * i]
        top[2 * i + 1] = boundary[4 * i + 1]
        bottom[2 * i] = boundary[4 * i + 2]
        bottom[2 * i + 1] = boundary[4 * i + 3]

    # Generate boundaries
    X, Y = np.meshgrid(np.arange(dims[2]), np.arange(dims[1]))

    n_1, _, p_1 = affine_fit(top)
    n_2, _, p_2 = affine_fit(bottom)

    Z_top = plane_equation_Z_from_XY(X=X, Y=Y, n=n_1, p=p_1, z_offset=z_offset)
    Z_bottom = plane_equation_Z_from_XY(X=X, Y=Y, n=n_2, p=p_2, z_offset=-z_offset)

    Z_top[Z_top > dims[0]] = dims[0]
    Z_bottom[Z_bottom < 0] = 0

    mask = make_boundary_mask(mask_size=dims, Z_top=Z_top, Z_bottom=Z_bottom)

    return mask
