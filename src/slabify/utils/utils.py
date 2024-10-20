import numpy as np


def measure_thickness(
    mask: np.ndarray,
    yx_coords: tuple[int, int] = None,
    z_offset: float = 0.0,
    angpix: float = None,
    verbose=False,
) -> tuple[float, float]:
    """
    Measure thickness in Z of binary slab mask at given point.

    Args:
        mask (np.ndarray):  The lamella slab mask.
        yx_coords (tuple, optional): the y,x coordinates where to measure the Z thickness.
        z_offset (float, optional): Offset in the Z direction for the mask. Defaults to 0.0.
        angpix (float, optional): The pixel size in Angstroms for displaying the measured thickness, only relevant if verbose is True.
        verbose (bool, optional): Whether to print the measured thickness with and without Z offset to the screen.

    Returns:
        tuple: Slab mask thickness with and without Z offset.
    """
    # If no coordinates are passed, we use the center of the volume:
    if not yx_coords:
        yx_coords = [
            mask.shape[1] // 2,
            mask.shape[2] // 2,
        ]  # Remember mrcfile convention is [Z, Y, X]

    t_with_offset = np.sum(mask[:, yx_coords[0], yx_coords[1]])
    t_no_offset = t_with_offset - 2 * z_offset

    if verbose:

        if not angpix:
            units = "px"
        else:
            units = "Å"
            t_with_offset *= angpix
            t_no_offset *= angpix

        print(f"Slab mask thickness with Z offset: {t_with_offset:.2f} {units}")
        print(f"Slab mask thickness without Z offset: {t_no_offset:.2f} {units}")

    return t_with_offset, t_no_offset


def sample_points(
    mask_size: np.ndarray, N: int = 50000, boxsize: int = 32, z_min: int = 1, z_max: int = None, seed: int = 42
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Sample points at random within the tomogram, staying away from the borders.

    Args:
        mask_size (np.ndarray): The dimensions of the tomogram.
        N (int, optional): Number of points to sample. Defaults to 50000.
        boxsize (int, optional): Box size in pixels around each sampled point. Defaults to 32.
        z_min (int, optional): Minimum Z slice to sample, starting from 1. Defaults to 1.
        z_max (int, optional): Maximum Z slice to sample. Defaults to None, which corresponds to the highest slice.
        seed (int, optional): Random seed for reproducibility. Defaults to 42.

    Returns:
        tuple: Arrays of sampled Z, Y, X coordinates.
    """
    dims = mask_size
    half_box = boxsize // 2
    if not z_max:
        z_max = dims[0]

    np.random.seed(seed=seed)
    Z_rand = np.random.randint(z_min + half_box - 1, z_max - half_box, size=(N,))
    Y_rand = np.random.randint(half_box, dims[1] - half_box, size=(N,))
    X_rand = np.random.randint(half_box, dims[2] - half_box, size=(N,))

    return Z_rand, Y_rand, X_rand


def variance_at_points(
    tomo: np.ndarray,
    Z: np.ndarray,
    Y: np.ndarray,
    X: np.ndarray,
    N: int = 50000,
    boxsize: int = 32,
) -> np.ndarray:
    """
    Calculate the variance around each sampled point within the specified box size.

    Args:
        tomo (np.ndarray): The tomogram array.
        Z (np.ndarray): Array of sampled Z coordinates.
        Y (np.ndarray): Array of sampled Y coordinates.
        X (np.ndarray): Array of sampled X coordinates.
        N (int, optional): Number of sampled points. Defaults to 50000.
        boxsize (int, optional): Box size in pixels around each sampled point. Defaults to 32.

    Returns:
        variances (np.ndarray): Array containing the variance at each sampled point.
    """

    half_box = boxsize // 2

    variances = np.zeros((N,))

    for i in np.arange(N):

        x, y, z = X[i], Y[i], Z[i]
        x_min, x_max = x - half_box, x + half_box
        y_min, y_max = y - half_box, y + half_box
        z_min, z_max = z - half_box, z + half_box

        box = tomo[z_min:z_max, y_min:y_max, x_min:x_max]

        variances[i] = np.var(box)

    return variances


def make_boundary_mask(mask_size: np.ndarray, Z_top: np.ndarray, Z_bottom: np.ndarray):
    """
    Creates the boundary mask volume given the top and bottom planes.

    Args:
        mask_size (np.ndarray): The dimensions of the tomogram.
        Z_top (np.ndarray): Volume containing the Z value of each voxel, according to top plane equation.
        Z_bottom (np.ndarray): Volume containing the Z value of each voxel, according to bottom plane equation.

    Returns:
        mask (np.ndarray): The binary boundary mask, with value of 1 within the boundary and 0 outside.
    """
    dims = mask_size

    # Generate mask
    mask = np.zeros(dims, dtype=np.int8)

    # Loop over X,Y and fill each column in Z with the mask where appropriate:
    for i in range(dims[2]):  # Loop over X
        for j in range(dims[1]):  # Loop over Y
            z1 = int(round(Z_top[j, i]))
            z2 = int(round(Z_bottom[j, i]))
            mask[z2:z1, j, i] = 1

    return mask


def apply_mask_border(mask: np.ndarray, xy_border: int = 0) -> np.ndarray:
    """
    Applies a border exclusion to the mask.

    Args:
        mask (np.ndarray): The binary mask to be adjusted.
        xy_border (int, optional): Number of pixels to exclude from the XY border. Defaults to 0.

    Returns:
        mask (np.ndarray): The adjusted binary mask.
    """
    if xy_border > 0:
        border = np.zeros(mask.shape, dtype=np.int8)
        border[:, xy_border:-xy_border, xy_border:-xy_border] = 1
        mask *= border

    return mask


def affine_fit(X: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Computes the plane that best fits a set of sample points.

    Args:
        X (np.ndarray): N by 3 matrix where each row is a sample point.

    Returns:
        tuple: (n, V, p), where
            n (np.ndarray): Unit vector normal to the plane.
            V (np.ndarray): 3 by 2 matrix with orthonormal basis vectors of the plane.
            p (np.ndarray): A point belonging to the plane.
    """
    # Compute mean of samples
    p = np.mean(X, axis=0)

    # Reduce samples
    R = X - p

    # Compute principal directions
    V, _, _ = np.linalg.svd(R.T)
    n = V[:, -1]
    V = V[:, 0:-1]

    return n, V, p


def plane_equation_Z_from_XY(
    X: np.ndarray, Y: np.ndarray, n: np.ndarray, p: np.ndarray, z_offset: float = 0
) -> np.ndarray:
    """
    Calculates the Z-coordinate of a plane given X and Y coordinates and plane coefficients.

    Args:
        X (np.ndarray): Array of X-coordinates.
        Y (np.ndarray): Array of Y-coordinates.
        coeffs (np.ndarray): Coefficients of the plane equation.
        intercept (float): Intercept of the plane equation.

    Returns:
        Z (np.ndarray): The Z-coordinate values corresponding to the input X and Y coordinates.
    """
    return -(n[0] * X + n[1] * Y - np.dot(n, p)) / n[2] + z_offset


def distance_of_points_to_plane(
    x: np.ndarray, y: np.ndarray, z: np.ndarray, n: np.ndarray, p: np.ndarray
) -> np.ndarray:
    """
    Calculates the distance of points to a plane.

    Args:
        points (np.ndarray): Array of points.
        coeffs (np.ndarray): Coefficients of the plane equation.
        intercept (float): Intercept of the plane equation.

    Returns:
        D (np.ndarray): The distances of the points from the plane.
    """
    d = -np.dot(n, p)

    return (n[2] * z + n[1] * y + n[0] * x + d) / np.sqrt(
        n[2] ** 2 + n[1] ** 2 + n[0] ** 2
    )
