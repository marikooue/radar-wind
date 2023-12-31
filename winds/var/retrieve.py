"""
winds.var.retrieve
==================

Retrieve the 3-D wind field from a single or multiple scanning Doppler radars
using a variational algorithm with multiple analysis constraints. Typical
constraints include Doppler velocity observations, mass continiuty, and spatial
smoothness.

"""

import os
import getpass
import platform
import numpy as np
import pyproj
from netCDF4 import num2date
from datetime import datetime
from scipy.optimize import minimize

from pyart.core import Grid
from pyart.config import get_field_name, get_metadata

from ..core import Profile, Solver
from ..util import standard
from . import common, _cost


def solve_wind_field(
        grids, solver=None, profile=None, background=None, first_guess='zero',
        fall_speed='Caya', atol=0.1, legacy=False, debug=False, verbose=False, **kwargs):
    """
    Retrieve the 3-D wind field from a single or multiple scanning Doppler
    radars.

    Parameters
    ----------
    grids : list-like
        Input grids used to derive the wind field.
    solver : Solver, optional
        Solver defining the necessary parameters for the wind retrieval
        and functional minimization.
    profile : Profile, optional
        Profile containing valid state variables such as altitude above mean
        sea level, temperature, and pressure. If None, a standard atmosphere
        defined by the U.S. Standard Model is used.
    background : 'zero' or 'profile', list or tuple, optional
        The background wind field. A background wind field of zero everywhere
        can be specified, as well as using the input Profile. If a list or
        tuple, then it is assumed the (u, v, w) wind components of the
        background field are supplied, in that order. If None, no background
        wind field constraint will be used.
    first_guess : 'zero' or 'profile', list or tuple, optional
        Similar to the background parameter except this specifies the initial
        (first) guess wind field. Typically this is set to zero everywhere.
    fall_speed : 'Caya' or 'supplied', optional
        The hydrometeor fall speed relation to use.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    Returns
    -------
    conv : pyart.core.Grid
        Grid containing retrieved (u, v, w) wind components, coordinate
        information, and metadata.

    """

    # Parse input grids
    if isinstance(grids, Grid):
        grids = [grids]
    if verbose:
        print 'Number of input grids: {}'.format(len(grids))

    # Parse solver
    if solver is None:
        solver = Solver(grids)

    # Parse analysis domain parameters
    nz, ny, nx = solver.shape
    N = solver.N

    # Parse atmospheric profile
    if profile is None:
        profile = standard.atmosphere(
            grid=grids[0], debug=debug, verbose=verbose)

    # Parse first (initial) guess field for control variables,
    #
    # u0 = (u1, u2, ... , uN)
    # v0 = (v1, v2, ... , vN)
    # w0 = (w1, w2, ... , wN).
    #
    # Finally, pack the first guess such that it is in vector form,
    #
    # xo = (u1, ... , uN, v1, ... , vN, w1, ... , wN).
    #
    # This is also how the analysis vector defined by the control variables is
    # to be structured
    if isinstance(first_guess, str):
        if first_guess.upper() == 'ZERO':
            u0 = np.zeros(N, dtype=np.float64)
            v0 = np.zeros(N, dtype=np.float64)
            w0 = np.zeros(N, dtype=np.float64)
        elif first_guess.upper() == 'PROFILE':
            u0 = profile.u.repeat(ny * nx).reshape(nz, ny, nx).flatten()
            v0 = profile.v.repeat(ny * nx).reshape(nz, ny, nx).flatten()
            w0 = profile.w.repeat(ny * nx).reshape(nz, ny, nx).flatten()
        else:
            raise ValueError('Unsupported first guess field')

    elif isinstance(first_guess, (list, tuple)):
        u0, v0, w0 = first_guess

    else:
        raise ValueError('Unsupported first guess field')

    x0 = np.concatenate((u0, v0, w0), axis=0)

    # Parse background field
    solver.add_background(background=background, profile=profile)

    # Compute radar line-of-sight components
    common._add_grid_line_sight_components(
        grids, solver, debug=debug, verbose=verbose)

    # Compute hydrometeor fall velocities
    if fall_speed.upper() == 'CAYA':
        common.add_fall_speed_caya(
            grids, solver, profile=profile, debug=debug, verbose=verbose)
    #else:
    #    raise ValueError('Unsupported hydrometeor fall speed relation')
    common._check_fall_speeds(grids, solver)

    # Compute radar observation weights
    if not np.all([solver.vdop_weight_field in grid.fields for grid in grids]):
	print 'Computing grid_radial_velocity_weight'
        common.add_grid_radial_velocity_weight(
            grids, solver, use_weights=solver.use_weights, debug=debug, verbose=verbose, **kwargs)

    # Determine available radar coverage in the analysis grid
    if not np.all([solver.cover_field in grid.fields for grid in grids]):
        common.add_grid_coverage(
            grids, solver, atol=atol,debug=debug, verbose=verbose, **kwargs)

    # Extrapolate profile throughout domain
    profile.extrapolate_to_grid(grids[0])

    # Parse the cost function and its gradient (Jacobian)
    fun = _cost._compute_wind_cost
    jac = _cost._compute_wind_grad

    # Peform heavily-smoothed first pass. This is designed to improve the
    # initial conditions by retrieving the mean horizontal wind flow
    if solver.first_pass:

        if verbose:
            print 'Performing heavily-smoothed first pass'

        # The first pass is meant to retrieve the large-scale horizontal wind
        # flow. This means that we want to ignore mass continuity and set the
        # spatial smoothness constraint weights to relatively "large" values
        solver_first = solver.copy()
        solver_first.continuity = None
        if solver_first.smoothness.upper() == 'POTVIN':
            solver_first.lambda_s1 = 100.0 * solver.lambda_s1
            solver_first.lambda_s2 = 100.0 * solver.lambda_s2
            solver_first.lambda_s3 = 100.0 * solver.lambda_s3
            solver_first.lambda_s4 = 100.0 * solver.lambda_s4

        # Call the SciPy solver
        res = minimize(
            fun, x0, args=(grids, solver_first, profile, debug, verbose),
            jac=jac, method=solver_first.method, hess=None, hessp=None,
            bounds=None, constraints=None,
            options=solver_first.options_first_pass)

        # Set the vertical velocity to 0 m/s everywhere after the first pass
        # since this pass does not attempt to properly retrieve vertical air
        # motion
        x0 = res.x.astype(np.float64)
        x0[2*N:3*N] = 0.0

    # Retrieve the full (u, v, w) wind field
    res = minimize(
        fun, x0, args=(grids, solver, profile, debug, verbose), jac=jac,
        method=solver.method, hess=None, hessp=None, bounds=None,
        constraints=None, options=solver.options)

    # Parse control variables from analysis vector. The analysis vector is
    # assumed to be packed,
    #
    # x = x(u1, u2, ... , uN, v1, v2, ... , vN, w1, w2, ... , wN)
    #
    # So u is packed first, then v, and finally w
    xopt = res.x.astype(np.float64)
    u = xopt[0:N].reshape(nz, ny, nx)
    v = xopt[N:2*N].reshape(nz, ny, nx)
    w = xopt[2*N:3*N].reshape(nz, ny, nx)

    # Parse wind component metadata from Py-ART configuration file and wind
    # field data
    u_wind = get_metadata(solver.u_field)
    v_wind = get_metadata(solver.v_field)
    w_wind = get_metadata(solver.w_field)
    u_wind['data'] = u.astype(np.float32)
    v_wind['data'] = v.astype(np.float32)
    w_wind['data'] = w.astype(np.float32)

    # Populate fields dictionary
    fields = {
        solver.u_field: u_wind,
        solver.v_field: v_wind,
        solver.w_field: w_wind,
        solver.cover_field: grids[0].fields[solver.cover_field],
        solver.refl_network_field: grids[0].fields[solver.refl_network_field],
        solver.fall_field: grids[0].fields[solver.fall_field],
        }

    # Populate grid metadata
    metadata = common._populate_metadata()

    # Create wind retrieval grid
    if legacy:
        axes = common._populate_legacy_axes(grids)
        conv = Grid.from_legacy_parameters(fields, axes, metadata)

    # Compute the horizontal wind divergence field
    common.add_wind_divergence(conv, solver, debug=debug, verbose=verbose)

    return conv
