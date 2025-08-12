"""
winds.var.common
================

Routines common to the variational wind retrieval module.

"""

import pyproj
import getpass
import platform
import numpy as np
from datetime import datetime
from scipy import ndimage
from scipy.interpolate import UnivariateSpline

from pyart.config import get_fillvalue, get_field_name, get_metadata
from pyart.util.datetime_utils import datetimes_from_dataset

from ..core import Profile
from ..util import standard, grid_time
from . import _operators

# Necessary and/or potential future improvements to the common submodule:
#
# * Revisit radar line-of-sight function to properly handle analysis grid
#   altitude and radar altitude.



def add_grid_radial_velocity_weight(
        grids, solver, use_weights=False, debug=False, verbose=False):
    """
    Compute radial velocity observation weights. Radial velocity observations
    can be treated as equal throughout the grid or can be weighted by their
    distance-dependent weights.

    Parameters
    ----------
    grids : list-like
        Input grids containing radial velocity observations.
    solver : Solver
        Solver used to retrieve the wind field.
    use_weights : bool, optional
        True to use distance-dependent weights, False to treat all reflectivity
        observations equal.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    """

    for grid in grids:

        if verbose:
            grid_name = grid.metadata['radar_0_instrument_name']
            print 'Radial velocities from grid: {}'.format(grid_name)

        if use_weights and solver.weight_field in grid.fields:
            print  'Weighting radial velocity: ',solver.weight_field
            weight = grid.fields[solver.weight_field]['data']
            #- modified by oue; for crsim (crsim might not have gqi_field)
	    weight[weight<solver.min_weight]=0.0
            #weight *= solver.lambda_o * grid.fields[solver.gqi_field]['data']
            if solver.gqi_field in grid.fields:
                weight *= solver.lambda_o * grid.fields[solver.gqi_field]['data']
            else:
                weight *= solver.lambda_o 
            #---
        else:
            weight = np.zeros(solver.shape, dtype=np.float32)
            #weight = np.squeeze(weight) # added by oue for RHI grid X-Z plane
            weight.fill(solver.lambda_o)
            print use_weights

        # Set radial velocity weight to zero where no radial velocity data is
        # available
        vdop = grid.fields[solver.vdop_field]['data']
	if (len(vdop.shape)==2): # added by oue for RHI grid X-Z plane
	    vdop=np.expand_dims(vdop,axis=1)
        weight[np.ma.getmaskarray(vdop)] = 0.0
	weight[vdop<-200] = 0.0

        if debug:
            print 'Min radial velocity weight: {:.3f}'.format(weight.min())
            print 'Max radial velocity weight: {:.3f}'.format(weight.max())

        # Create new radial velocity weight field
        weight_dict = get_metadata(solver.vdop_weight_field)
        weight_dict['data'] = weight.astype(np.float32)
        weight_dict['valid_max'] = solver.lambda_o
        grid.add_field(
            solver.vdop_weight_field, weight_dict, replace_existing=True)

    return


def add_grid_network_reflectivity(
        grids, solver, ignore_radars=None, use_weights=True, debug=False,
        verbose=False):
    """
    Compute the reflectivity field from a network of scanning radars.
    Reflectivity observations from multiple radars can be treated as equal or
    as a weighted average using their distance-dependent weights.

    Parameters
    ----------
    grids : list-like
        Input grids which contain reflectivity fields to compute the
        reflectivity of the network of radars.
    solver : Solver
        Solver used to retrieve the wind field.
    ignore_radars : list or tuple, optional
        Names of radars as listed in the grid metadata which should be ignored
        when computing the radar network reflectivity field. This is to allow
        for the exclusion of certain radars that may have degraded and/or
        compromised reflectivity observations (e.g., X-band and other shorter
        wavelength radars).
    use_weights : bool, optional
        True to use distance-dependent weights, False to treat all reflectivity
        observations equal.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print progress and identification information, False to
        suppress.

    """

    refl_data = []
    weight_data = []

    for grid in grids:

        # Ignore grid if name matches any of those listed
        grid_name = grid.metadata['radar_0_instrument_name']
        if ignore_radars is not None and grid_name in ignore_radars:
            continue

        if verbose:
            print('Reflectivity observations from grid: {}'.format(grid_name))

        # Convert reflectivity data to linear units
        refl = np.ma.power(
            10.0, (grid.fields[solver.refl_field]['data'] / 10.0))
	if(len(refl.shape)==2): # added by oue for RHI grid X-Z plane
	    refl=np.expand_dims(refl,axis=1)
        refl_data.append(refl)

        if use_weights and solver.weight_field in grid.fields:
            weight = grid.fields[solver.weight_field]['data']
        else:
            weight = np.ones(solver.shape, dtype=np.float32)
	    #weight = np.squeeze(weight) # added by oue for RHI grid X-Z plane
        weight_data.append(weight)

        if debug:
            print 'Min reflectivity weight: {:.3f}'.format(weight.min())
            print 'Max reflectivity weight: {:.3f}'.format(weight.max())

    # Compute weighted average of reflectivity data and convert reflectivity
    # back to decibel
    refl = np.ma.average(refl_data, weights=weight_data, axis=0)
    refl = 10.0 * np.ma.log10(refl)
    refl = np.ma.masked_invalid(refl, copy=False)
    refl.set_fill_value(solver.fill_value)

    if debug:
        print 'Min network reflectivity: {:.3f} dBZ'.format(refl.min())
        print 'Max network reflectivity: {:.3f} dBZ'.format(refl.max())

    # Create radar network reflectivity field dictionary
    refl_dict = get_metadata(solver.refl_field)
    refl_dict['data'] = refl.astype(np.float32)
    refl_dict['_FillValue'] = refl.fill_value
    for grid in grids:
        grid.add_field(
            solver.refl_network_field, refl_dict, replace_existing=True)

    return


def add_fall_speed_caya(
        grids, solver, profile=None, debug=False, verbose=False):
    """
    Compute hydrometeor fall speeds using the empirical relationship from Caya
    (2001).

    Parameters
    ----------
    grids : list-like
        Input grids which contain radar network reflectivity field.
    solver : Solver
        Solver used to retrieve the wind field.
    profile : Profile
        Profile containing temperature data in Celsius. If no air temperature
        data is available or no profile is provided, the U.S. Standard Model is
        used.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    References
    ----------

    """

    # Parse vertical grid coordinate
    nz, ny, nx = solver.shape
    Z = np.repeat(solver.z, ny * nx, axis=0).reshape(nz, ny, nx)

    # Define liquid and ice relations
    liquid = lambda M, Z: 5.94 * M**(1.0 / 8.0) * np.exp(Z / 20000.0)
    ice = lambda M, Z: 1.15 * M**(1.0 / 12.0) * np.exp(Z / 20000.0)

    # Compute precipitation concentration from reflectivity data
    refl = grids[0].fields[solver.refl_network_field]['data']
    M = np.ma.exp((refl - 43.1) / 7.6)

    # Parse air temperature data
    if solver.tdry_field in grids[0].fields:
        tdry = grids[0].fields[solver.tdry_field]['data']
    else:
        if profile is None:
            profile = standard.atmosphere(grid=grids)
        tdry = profile.tdry.repeat(ny * nx, axis=0).reshape(nz, ny, nx)

    # Compute hydrometeor fall speed
    # Positive values correspond to downwards motion
    # Set fall speeds to zero where reflectivity observations are unavailable
    vfall = np.ma.where(tdry >= 0.0, liquid(M, Z), ice(M, Z))
    vfall = np.ma.filled(vfall, 0.0)

    if debug:
        print 'Min hydrometeor fall speed: {:.3f} m/s'.format(vfall.min())
        print 'Max hydrometeor fall speed: {:.3f} m/s'.format(vfall.max())

    # Create hydrometeor fall speed field dictionary
    vfall_dict = get_metadata(solver.fall_field)
    vfall_dict['data'] = vfall.astype(np.float32)
    for grid in grids:
        grid.add_field(solver.fall_field, vfall_dict, replace_existing=True)

    return


def add_wind_divergence(grid, solver, debug=False, verbose=False):
    """
    Compute wind divergence field.

    Parameters
    ----------
    grid : pyart.core.Grid
        Grid containing retrieved wind fields.
    solver : Solver
        Solver used to retrieve the wind field.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    """

    # Parse wind field components
    u = grid.fields[solver.u_field]['data']
    v = grid.fields[solver.v_field]['data']
    w = grid.fields[solver.w_field]['data']

    # Prepare wind field data for ingest into Fortran wrapper
    u = np.asfortranarray(u, dtype=np.float64)
    v = np.asfortranarray(v, dtype=np.float64)
    w = np.asfortranarray(w, dtype=np.float64)

    # Compute the wind field divergence and its horizontal component
    div, dudx, dvdy, dwdz = _operators.div_wind(
        u, v, w, dx=solver.dx, dy=solver.dy, dz=solver.dz, proc=solver.proc,
        finite_order=solver.finite_order, fill_value=solver.fill_value)
    hdiv = dudx + dvdy

    # Create full wind divergence field dictionary
    div_dict = get_metadata(solver.div_field)
    div_dict['data'] = div.astype(np.float32)
    grid.add_field(solver.div_field, div_dict, replace_existing=True)

    # Create horizontal wind divergence field dictionary
    hdiv_dict = get_metadata(solver.hdiv_field)
    hdiv_dict['data'] = hdiv.astype(np.float32)
    grid.add_field(solver.hdiv_field, hdiv_dict, replace_existing=True)

    return


def add_grid_coverage(grids, solver, atol=0.1, debug=False, verbose=False):
    """
    Determine coverage of the input grids, i.e., the total number of available
    radar observations at each grid point.

    Parameters
    ----------
    grids : list-like
        Input grids containing Doppler velocity observation weights used to
        determine the coverage within the analysis domain.
    solver : Solver
        Solver used to derive the wind field.
    atol : float, optional
        The absolute tolerance used to determine whether a weight is considered
        zero.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    """

    nobs = np.zeros(solver.shape, dtype=np.int32)
    #nobs=np.squeeze(nobs) #added by oue for RHI grid X-Z plane

    for grid in grids:

        if verbose:
            grid_name = grid.metadata['radar_0_instrument_name']
            print 'Coverage by grid: {}'.format(grid_name)

        # Parse radial velocity weight data
        weight = grid.fields[solver.vdop_weight_field]['data']

        # Determine where weight is significant enough to be considered a
        # valid observation
        significant = ~np.isclose(weight, 0.0, rtol=1.0e-15, atol=atol)
        nobs += significant.astype(np.int32)

    if debug:
        print 'Min radar coverage: {}'.format(nobs.min())
        print 'Max radar coverage: {}'.format(nobs.max())

    # Create new radar coverage field dictionary
    cover_dict = get_metadata(solver.cover_field)
    cover_dict['data'] = nobs.astype(np.int32)
    cover_dict['valid_min'] = 0
    cover_dict['valid_max'] = len(grids)
    for grid in grids:
        grid.add_field(solver.cover_field, cover_dict, replace_existing=True)

    return


def arm_mergesonde_profile(
        solver, sonde, rate_method='cubic', debug=False, verbose=False):
    """
    Create Profile from ARM merged sounding product.

    Parameters
    ----------
    solver : Solver
        Solver used to retrieve the wind field.
    sonde : netCDF4.Dataset
        ARM merged sounding dataset.
    rate_method : 'cubic' or 'finite', optional
        Method used to compute rates of change with altitude.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    Returns
    -------
    profile : Profile
        Profile derived from ARM merged sounding dataset.

    """

    # Parse vertical coordinate
    # The heights defined in the ARM merged sounding product are recorded in
    # kilometers above mean sea level
    alt_grid = solver.z + solver.alt_0
    alt_sonde = sonde.variables['height'][:] * 1000.0

    # Parse time data
    time_sonde = datetimes_from_dataset(sonde)
    idx = np.abs(time_sonde - solver.time_target).argmin()

    if debug:
        print 'Solver target time: {}'.format(solver.time_target)
        print 'Sounding target time: {}'.format(time_sonde[idx])

    # Parse merged sounding data
    # Dry bulb temperature -> degrees Celsius
    # Dew point temperature -> degrees Celsius
    # Pressure -> kPa
    # Wind speed -> m/s
    # Relative humidity -> percent
    # Specific humidity ->
    tdry_sonde = sonde.variables['temp'][idx,:]
    dp_sonde = sonde.variables['dp'][idx,:]
    pres_sonde = sonde.variables['bar_pres'][idx,:]
    vap_pres_sonde = sonde.variables['vap_pres'][idx,:]
    rh_sonde = sonde.variables['rh'][idx,:]
    sh_sonde = sonde.variables['sh'][idx,:]
    wspd_sonde = sonde.variables['wspd'][idx,:]
    wdir_sonde = sonde.variables['wdir'][idx,:]
    u_sonde = sonde.variables['u_wind'][idx,:]
    v_sonde = sonde.variables['v_wind'][idx,:]
    w_sonde = np.zeros_like(tdry_sonde, subok=False)

    # Compute the density of dry air from the equation of state
    rho_sonde = (pres_sonde * 1000.0) / (287.058 * (tdry_sonde + 273.15))

    # Interpolate merged sounding profile using cubic spline
    tdry = UnivariateSpline(alt_sonde, tdry_sonde).__call__(alt_grid)
    dp = UnivariateSpline(alt_sonde, dp_sonde).__call__(alt_grid)
    pres = UnivariateSpline(alt_sonde, pres_sonde).__call__(alt_grid)
    e = UnivariateSpline(alt_sonde, vap_pres_sonde).__call__(alt_grid)
    rho = UnivariateSpline(alt_sonde, rho_sonde).__call__(alt_grid)
    rh = UnivariateSpline(alt_sonde, rh_sonde).__call__(alt_grid)
    sh = UnivariateSpline(alt_sonde, sh_sonde).__call__(alt_grid)
    wspd = UnivariateSpline(alt_sonde, wspd_sonde).__call__(alt_grid)
    wdir = UnivariateSpline(alt_sonde, wdir_sonde).__call__(alt_grid)
    u = UnivariateSpline(alt_sonde, u_sonde).__call__(alt_grid)
    v = UnivariateSpline(alt_sonde, v_sonde).__call__(alt_grid)
    w = UnivariateSpline(alt_sonde, w_sonde).__call__(alt_grid)

    return Profile(alt=alt_grid, tdry=tdry, dp=dp, pres=pres, rho=rho, e=e,
                   wspd=wspd, wdir=wdir, rh=rh, sh=sh, u=u, v=v, w=w,
                   rate_method=rate_method)


def _add_grid_line_sight_components(
        grids, solver, proj='lcca', datum='WGS84', ellps='WGS84', debug=False,
        verbose=False):
    """
    Compute line-of-sight components using specified map projection.

    Parameters
    ----------
    grids : list-like
        Input grids which geolocation information.
    solver : Solver
        Solver used to retrieve the wind field.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    """

    for grid in grids:

        if verbose:
            grid_name = grid.metadata['radar_0_instrument_name']
            print 'Line-of-sight components for grid: {}'.format(grid_name)

        # Parse geolocation data
        lat_0 = grid.origin_latitude['data'][0]
        lon_0 = grid.origin_longitude['data'][0]
        alt_0 = grid.origin_altitude['data'][0]
        lat_radar = grid.metadata['radar_0_latitude']
        lon_radar = grid.metadata['radar_0_longitude']
        alt_radar = grid.metadata['radar_0_altitude']

        # Create map projection
        pj = pyproj.Proj(
            proj=proj, lat_0=lat_0, lon_0=lon_0, x_0=0.0, y_0=0.0, datum=datum,
            ellps=ellps)

        # Compute (x, y, z) offset of radar in the grid from map projection
        x_radar, y_radar = pj(lon_radar, lat_radar)
        z_radar = alt_radar - alt_0

        if debug:
            print 'Radar x offset in grid: {:.3f} km'.format(x_radar / 1000.0)
            print 'Radar y offset in grid: {:.3f} km'.format(y_radar / 1000.0)
            print 'Radar z offset in grid: {:.3f} km'.format(z_radar / 1000.0)

        # Parse grid coordinates and change reference frame
        x = grid.x['data'] - x_radar
        y = grid.y['data'] - y_radar
        z = grid.z['data'] - z_radar
        Z, Y, X = np.meshgrid(z, y, x, indexing='ij')

        # Compute line-of-sight components
        R = np.sqrt(X**2 + Y**2 + Z**2)
        xhat = X / R
        yhat = Y / R
        zhat = Z / R

        # Create new x-component field dictionary
        xhat_dict = get_metadata(solver.xhat_field)
        xhat_dict['data'] = xhat.astype(np.float32)
        grid.add_field(solver.xhat_field, xhat_dict, replace_existing=True)

        # Create new y-component field dictionary
        yhat_dict = get_metadata(solver.yhat_field)
        yhat_dict['data'] = yhat.astype(np.float32)
        grid.add_field(solver.yhat_field, yhat_dict, replace_existing=True)

        # Create new z-component field dictionary
        zhat_dict = get_metadata(solver.zhat_field)
        zhat_dict['data'] = zhat.astype(np.float32)
        grid.add_field(solver.zhat_field, zhat_dict, replace_existing=True)

    return


def _check_fall_speeds(grids, solver):
    """ Check no missing data in hydrometeor fall speed field. """

    for grid in grids:
        vfall = grid.fields[solver.fall_field]['data']
        grid.fields[solver.fall_field]['data'] = np.ma.filled(vfall, 0.0)

    return


def _populate_legacy_axes(grids):
    """ Populate legacy axes data and metadata. """

    return {
        'time': grids[0].axes['time'],
        'time_start': grids[0].axes['time_start'],
        'time_end': grids[0].axes['time_end'],
        'x_disp': grids[0].axes['x_disp'],
        'y_disp': grids[0].axes['y_disp'],
        'z_disp': grids[0].axes['z_disp'],
        'alt': grids[0].axes['alt'],
        'lat': grids[0].axes['lat'],
        'lon': grids[0].axes['lon'],
        }


def _populate_metadata():
    """ Populate grid metadata. """

    # Datastreams attribute
    datastream_description = (
        'A string consisting of the datastream(s), datastream version(s), '
        'and datastream date (range).')

    # History attribute
    history = 'created by user {} on {} at {}'.format(
        getpass.getuser(), platform.node(),
        datetime.now().strftime('%Y-%m-%dT%H:%M:%S'))

    return {
        'process_version': '',
        'references': '',
        'Conventions': '',
        'site_id': '',
        'site': '',
        'facility_id': '',
        'project': '',
        'state': '',
        'comment': '',
        'institution': '',
        'country': '',
        'description': '',
        'title': '',
        'project': '',
        'input_datastreams_description': datastream_description,
        'input_datastreams_num': '',
        'input_datastreams': '',
        'history': history,
        }
