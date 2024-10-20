import numpy as np
from sklearn.linear_model import RANSACRegressor

from slabify.utils import (
    distance_of_points_to_plane,
    make_boundary_mask,
    plane_equation_Z_from_XY,
    sample_points,
    variance_at_points,
)


def create_boundary_mask_auto(
    tomo: np.ndarray,
    N: int = 50000,
    boxsize: int = 32,
    z_min: int = 1,
    z_max: int = None,
    z_offset: float = 0.0,
    simple: bool = False,
    thickness: int = None,
    iterations: int = 3,
    percentile: float = 95,
    seed: int = 42,
) -> np.ndarray:
    """
    Automatically create a slab (boundary) mask by fitting one or two planes to enclose the points with high variance.

    Args:
        tomo (np.ndarray): The tomogram array.
        N (int, optional): Number of points to sample. Defaults to 50000.
        boxsize (int, optional): Box size in pixels to analyze variance around each sampled point. Defaults to 32.
        z_min (int, optional): Minimum Z slice to sample, starting from 1. Defaults to 1.
        z_max (int, optional): Maximum Z slice to sample. Defaults to None, which corresponds to the highest slice.
        z_offset (float, optional): Offset in the Z direction for the mask. Defaults to 0.0.
        simple (bool, optional): Whether to use the simple masking method (single plane). Defaults to False.
        thickness (int, optional): Total thickness of the lamella in pixels (used in simple mode). Defaults to None.
        iterations (int, optional): Number of iterations for plane fitting. Defaults to 3.
        percentile (float, optional): Percentile of highest variance locations to use for fitting. Defaults to 95.
        seed (int, optional): Random seed for reproducibility. Defaults to 42.

    Returns
    -------
        mask (np.ndarray): Binary mask representing the lamella slab.
    """
    dims = tomo.shape
    # Sample N points at random:
    Z_rand, Y_rand, X_rand = sample_points(
        mask_size=dims, N=N, boxsize=boxsize, z_min=z_min, z_max=z_max, seed=seed
    )
    # Calculate the variance around each point:
    variances = variance_at_points(
        tomo=tomo, Z=Z_rand, Y=Y_rand, X=X_rand, N=N, boxsize=boxsize
    )

    variance_thr = np.percentile(variances, percentile)
    idx = variances[:] > variance_thr
    idx = idx.squeeze()
    # We now threshold to only work with variances and coordinates for points above the threshold:
    # Hopefully this represents points with "interesting" density, i.e. within the lamella:
    variances = variances[idx]
    Z_rand, Y_rand, X_rand = Z_rand[idx], Y_rand[idx], X_rand[idx]

    # Robust linear fit using RANSAC:
    fit = RANSACRegressor()
    fit.fit(np.transpose([X_rand, Y_rand]), Z_rand)
    n = np.array([fit.estimator_.coef_[0], fit.estimator_.coef_[1], -1])
    p = np.array([0, 0, fit.estimator_.intercept_])
    n_top = n_bottom = np.zeros(n.shape)  # Initialize
    p_top = p_bottom = np.zeros(p.shape)  # Initialize

    X, Y = np.meshgrid(np.arange(dims[2]), np.arange(dims[1]))

    if simple:
        # Fit a single plane through the ~center of the lamella slab:
        if not thickness:
            thickness = float(dims[0]) / 2
        half_thickness = thickness / 2
        # The top and bottom planes are now defined
        Z_top = plane_equation_Z_from_XY(
            X=X, Y=Y, n=n, p=p, z_offset=half_thickness + z_offset
        )
        Z_bottom = plane_equation_Z_from_XY(
            X=X, Y=Y, n=n, p=p, z_offset=-half_thickness - z_offset
        )

    else:
        # Fit two planes to the high-variance points: one for the top and one for the bottom.

        # We start with a single plane at the center:
        D = distance_of_points_to_plane(x=X_rand, y=Y_rand, z=Z_rand, n=n, p=p)
        print(f"D.mean = {D.mean():.3f}; D.min = {D.min():.3f}; D.max = {D.max():.3f}")

        # 'above' and 'below' definitions must be inverted because of internal coordinate conventions:
        # above = D < 0
        # below = D > 0
        percentile_dist = 98
        thr_above = -np.percentile(-D, percentile_dist)
        print(f"thr_above = {thr_above:.3f}")
        above = D < thr_above
        print(f"D[above].mean = {D[above].mean():.3f}; D[above].min = {D[above].min():.3f}; D[above].max = {D[above].max():.3f}")
        # above = above.squeeze()
        thr_below = np.percentile(D, percentile_dist)
        print(f"thr_below = {thr_below:.3f}")
        below = D > thr_below
        print(f"D[below].mean = {D[below].mean():.3f}; D[below].min = {D[below].min():.3f}; D[below].max = {D[below].max():.3f}")
        # below = below.squeeze()

        X_above, Y_above, Z_above = X_rand[above], Y_rand[above], Z_rand[above]
        X_below, Y_below, Z_below = X_rand[below], Y_rand[below], Z_rand[below]

        fit_top = RANSACRegressor()
        fit_top.fit(np.transpose([X_above, Y_above]), Z_above)
        fit_bottom = RANSACRegressor()
        fit_bottom.fit(np.transpose([X_below, Y_below]), Z_below)
        n_top = np.array(
            [fit_top.estimator_.coef_[0], fit_top.estimator_.coef_[1], -1]
        )
        p_top = np.array([0, 0, fit_top.estimator_.intercept_])
        n_bottom = np.array(
            [fit_bottom.estimator_.coef_[0], fit_bottom.estimator_.coef_[1], -1]
        )
        p_bottom = np.array([0, 0, fit_bottom.estimator_.intercept_])

    Z_top = plane_equation_Z_from_XY(
        X=X, Y=Y, n=n_top, p=p_top, z_offset=+z_offset
    )
    Z_bottom = plane_equation_Z_from_XY(
        X=X, Y=Y, n=n_bottom, p=p_bottom, z_offset=-z_offset
    )

    Z_top[Z_top > dims[0]] = dims[0]
    Z_bottom[Z_bottom < 0] = 0

    mask = make_boundary_mask(mask_size=dims, Z_top=Z_top, Z_bottom=Z_bottom)

    return mask
