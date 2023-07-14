"""
winds.util.common
=================

A module containing functions common to single and multi-Doppler radar wind
retrievals.

"""

import numpy as np


def parse_grid_resolution(
        grid, check_uniform=True, debug=False, verbose=False):
    """
    Parse the resolution of each dimension found in a Grid.

    Parameters
    ----------
    grid : Grid
        Py-ART Grid containing axes data.
    check_uniform : bool, optional
        True to check if all grid dimensions have uniform resolution, and if so
        return a scalar value for each dimension. If False, the resolutions
        between each grid point for each dimension are returned.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print progress and identification information, False to
        suppress.

    Returns
    -------
    dx : float or ndarray
        The grid resolution in the x dimension.
    dy : float or ndarray
        The grid resolution in the y dimension.
    dz : float or ndarray
        The grid resolution in the z dimension.

    """

    # Compute coordinate differences
    dx = np.diff(grid.axes['x_disp']['data'], n=1)
    dy = np.diff(grid.axes['y_disp']['data'], n=1)
    dz = np.diff(grid.axes['z_disp']['data'], n=1)
    #- added by oue for RHI grid (X-Z plane)
    if (len(grid.axes['y_disp']['data'])==1):
	dy = np.diff(grid.axes['y_disp']['data'], n=0)
    #---
    if check_uniform:
        # Check for variable grid resolutions
        dx_uniform = np.isclose(dx.std(), 0.0, atol=2.0e-3)
        dy_uniform = np.isclose(dy.std(), 0.0, atol=2.0e-3)
        dz_uniform = np.isclose(dz.std(), 0.0, atol=2.0e-3)
        if not dx_uniform or not dy_uniform or not dz_uniform:
            raise ValueError('One or more grid dimensions are non-uniform')

        dx = dx[0]
        dy = dy[0]
        dz = dz[0]

        if verbose:
            print 'x-dimension resolution: {:.2f} m'.format(dx)
            print 'y-dimension resolution: {:.2f} m'.format(dy)
            print 'z-dimension resolution: {:.2f} m'.format(dz)

    return dx, dy, dz
