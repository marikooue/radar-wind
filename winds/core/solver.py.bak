"""
winds.core.solver
=================

A class for defining a solver to be used for variational scanning radar wind
retrievals. Parameters such as the analysis constraint weights, number of
minimization iterations, and what type of mass continuity constraint to use,
among many others, are addressed with the solver class.

"""


import numpy as np
from copy import deepcopy
from scipy.optimize import minimize
from netCDF4 import num2date

from pyart.core import Grid
from pyart.config import get_field_name, get_fillvalue

from ..util import common
from ..util.grid_time import datetimes_from_grid


_BUILTIN_METHODS = [
        'Nelder-Mead',
        'Powell',
        'CG',
        'BFGS',
        'Newton-CG',
        'L-BFGS-B',
        'TNC',
        'COBYLA',
        'SLSQP',
        'dogleg',
        'trust-ncg',
        ]


class Solver(object):
    """
    A class for defining the parameters necessary for retrieving the 3-D wind
    field from scanning Doppler radars.

    Attributes
    ----------
    constant_resolution : bool, optional
        True if analysis domain has constant resolution in all dimensions,
        False otherwise. Note that variable resolution grids are currently not
        supported by 3D-VAR technique.
    method : str or callable, optional
        If string is given then must be one compatible with SciPy's
        minimization routines. See SciPy documentation for list of available
        methods. The default specifies SciPy's conjugate-gradient method.
    controls : int, optional
        The number of control variables. For typical multi-Doppler wind
        retrievals, this is the three wind components.
    first_pass : bool, optional
        True to perform a heavily-smoothed first pass. This is designed to
        retrieve the large-scale horizontal flow which is then used as the
        initial conditions for the full 3-D wind retrieval.
    maxiter : int, optional
        Maximum number of minimization iterations to perform. Results from
        sensitivity testing indicate that this should typically be greater or
        equal to 200.
    maxiter_first_pass : int, optional
        Maximum number of minimization iterations to perform during a heavily-
        smoothed first pass. Typically this should be set between 50 to 100
        iterations.
    lambda_o : float, optional
        The radar Doppler velocity observation weight. If the observations are
        treated as variable, then this weight represents the maximum radar
        observation weight possible. In the literature this is generally set
        to 1.
    lambda_s1, lambda_s2, lambda_s3, lambda_s4 : float, optional
        Smoothness constraint weights as in equation (6) of Potvin et al.
        (2012).
    length_scale : float, optional
        Length scale in meters for relevant cost functions. The length scale is
        designed to make the dimensionality of each cost uniform, as well as
        bring each cost within 1-3 orders of magnitude of each other.
    fill_value : float, optional
        Value indicating missing or bad data. Default value uses Py-ART
        configuration file.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print relevant information, False to suppress.

    Examples
    --------

    """
    def __init__(self, grids, constant_resolution=True, method='CG',
                 controls=3, maxiter=200, continuity='Potvin',
                 smoothness='Potvin', impermeability='weak',
                 length_scale=500.0, first_pass=True, maxiter_first_pass=50,
                 finite_order='low', proc=1, use_weights=False,min_weight=0.001,adv_flag=None, fill_value=None, debug=False,
                 verbose=False, **kwargs):
        """ Initialize. """

        # Default solver attributes
        self.constant_resolution = constant_resolution
        self.method = method
        self.controls = controls
        self.maxiter = maxiter
        self.continuity = continuity
        self.smoothness = smoothness
        self.length_scale = length_scale
        self.first_pass = first_pass
        self.maxiter_first_pass = maxiter_first_pass
        self.finite_order = finite_order
        self.proc = proc
        self.adv_flag = adv_flag
	self.use_weights = use_weights
	self.min_weight = min_weight

        # Fill value attribute
        if fill_value is None:
            self.fill_value = get_fillvalue()

        # Parse grid field names
        self._parse_field_names(**kwargs)

        # Group all relevant solver options together. This is primarily for the
        # use of the SciPy optimize routines
        self.options = {
            'maxiter': self.maxiter,
            'gtol': kwargs.get('gtol', 1.0),
            'disp': verbose,
            }
        self.options_first_pass = {
            'maxiter': self.maxiter_first_pass,
            'gtol': kwargs.get('gtol', 1.0),
            'disp': verbose,
            }

        # Parse input grids
        if isinstance(grids, Grid):
            grids = [grids]

        # Parse grid coordinates and shape
        self.x = grids[0].x['data']
        self.y = grids[0].y['data']
        self.z = grids[0].z['data']
        self.nz, self.ny, self.nx = self.z.size, self.y.size, self.x.size
        self.shape = (self.nz, self.ny, self.nx)

        # Parse grid resolutions
        self.dx, self.dy, self.dz = common.parse_grid_resolution(
            grids[0], check_uniform=self.constant_resolution, debug=debug,
            verbose=verbose)

        # Parse grid origin
        self.lat_0 = grids[0].origin_latitude['data'][0]
        self.lon_0 = grids[0].origin_longitude['data'][0]
        self.alt_0 = grids[0].origin_altitude['data'][0]

        self._compute_functional_size(verbose=verbose)

        # Parse solver target time, defined as the half-way time between the
        # earliest grid start and latest grid end
        self.time_target = datetimes_from_grid(grids[0]).min()
        self.time_start = datetimes_from_grid(grids[0]).min()
        self.time_end = datetimes_from_grid(grids[0]).max()
        if verbose:
            print 'Solver target time: {}'.format(self.time_target)

        # Parse minimizaton method
        if isinstance(self.method, str):
            if self.method in _BUILTIN_METHODS:
                self.solver = None
        else:
            raise ValueError('Unsupported minimization method')

        # Parse radar radial velocity observational constraint parameters
        self.lambda_o = kwargs.get('lambda_o', 1.0)

        # Parse air mass continuity constraint parameters
        if self.continuity.upper() == 'POTVIN':
            self._parse_potvin_continuity(**kwargs)

        elif self.continuity.upper() in ['PROTAT', 'BOTTOM-UP', 'TOP-DOWN']:
            self.lambda_c = kwargs.get('lambda_c', 1.0)

        else:
            raise ValueError('Unsupported mass continuity constraint')

        # Parse smoothness constraint parameters
        if self.smoothness.upper() == 'POTVIN':
            self._parse_potvin_smoothness(**kwargs)
        else:
            raise ValueError('Unsupported smoothness constraint')

        # Parse background constraint parameters
        # Typically an atmospheric sounding is used as the background field,
        # and since a measurement of vertical air motion is not available from
        # these soundings no weight should be given to the vertical component
        self.lambda_ub = kwargs.get('lambda_ub', 0.01)
        self.lambda_vb = kwargs.get('lambda_vb', 0.01)
        self.lambda_wb = kwargs.get('lambda_wb', 0.0)

        # Parse impermeability constraint parameters
        # The default constraint weight value implies this constraint behaves
        # like a strong constraint rather than a weak constraint
        self.impermeability = impermeability
        if self.impermeability == 'WEAK' or 'weak': #modified by oue
            self.lambda_p = kwargs.get('lambda_p', 1000.0)
            #self.lambda_p = kwargs.get('lambda_p', 0.0001)
        else:
            self.lambda_p = kwargs.get('lambda_p', impermeability)
            self.impermeability.upper=impermeability

        # Parse minimization parameters, e.g., iteration number, value of
        # the cost function, and cost function gradient magnitude
        self.cost_iter = 0
        self.cost_grad_iter = 0
        self.cost_value = []
        self.cost_grad_mag = []


    def copy(self):
        """ Return a deep copy of self. """
        return deepcopy(self)


    def add_background(self, background=None, profile=None):
        """
        Add a background wind field to the Solver.

        Parameters
        ----------
        background : 'zero' or 'profile' or list or tuple, optional
            The background field. If a list or tuple is specified, it is
            assumed to contain the (u, v, w) background wind components, in
            that order.
        profile : Profile, optional
            An atmospheric profile containing estimates of (u, v, w) at each
            grid height.

        """

        if background is None:
            ub, vb, wb = None, None, None

        elif isinstance(background, str):
            nz, ny, nx = self.nz, self.ny, self.nx
            if background.upper() == 'ZERO':
                ub = np.zeros((nz, ny, nx), dtype=np.float64)
                vb = np.zeros_like(ub, subok=True)
                wb = np.zeros_like(ub, subok=True)
            elif background.upper() == 'PROFILE' and profile is not None:
                ub = profile.u.repeat(ny * nx).reshape(nz, ny, nx)
                vb = profile.v.repeat(ny * nx).reshape(nz, ny, nx)
                wb = profile.w.repeat(ny * nx).reshape(nz, ny, nx)
            else:
                raise ValueError('Unsupported background field')

        elif isinstance(background, (list, tuple)):
            ub, vb, wb = background
            background = 'supplied'

        else:
            raise ValueError('Unsupported background field')

        self.background = background
        self.ub, self.vb, self.wb = ub, vb, wb

        return


    def _compute_functional_size(self, verbose=True):
        """ Compute functional size. """

        self.N = self.nz * self.ny * self.nx
        self.functional_size = self.controls * self.N
        if verbose:
            print 'Functional size: {}'.format(self.functional_size)
        return


    def _parse_potvin_continuity(self, **kwargs):
        """ Parse mass continuity weights similar to Potvin et al. (2012). """
        self.lambda_c = kwargs.get('lambda_c', 500.0)
        return


    def _parse_potvin_smoothness(self, **kwargs):
        """ Parse smoothness weights similar to Potvin et al. (2012). """
        self.lambda_s1 = kwargs.get('lambda_s1', 1.0)
        self.lambda_s2 = kwargs.get('lambda_s2', 1.0)
        self.lambda_s3 = kwargs.get('lambda_s3', 1.0)
        self.lambda_s4 = kwargs.get('lambda_s4', 0.1)
        return


    def _parse_field_names(self, **kwargs):
        """ Parse grid field names. """

        # Radar observations and retrieved variable field names
        self.refl_field = kwargs.get(
            'refl_field', get_field_name('corrected_reflectivity'))
        self.vdop_field = kwargs.get(
            'vdop_field', get_field_name('corrected_velocity'))
        self.dist_field = kwargs.get(
            'dist_field', get_field_name('nearest_neighbor_distance'))
        self.weight_field = kwargs.get(
            'weight_field', get_field_name('nearest_neighbor_weight'))
        self.gqi_field = kwargs.get(
            'gqi_field', get_field_name('grid_quality_index'))
        self.fall_field = kwargs.get(
            'fall_field', get_field_name('hydrometeor_fall_speed'))
        self.refl_network_field = kwargs.get(
            'refl_network_field', 'radar_network_reflectivity')
        self.vdop_weight_field = kwargs.get(
            'vdop_weight_field', 'corrected_velocity_weight')

        # Wind retrieval field names
        self.u_field = kwargs.get(
            'u_field', get_field_name('eastward_wind_component'))
        self.v_field = kwargs.get(
            'v_field', get_field_name('northward_wind_component'))
        self.w_field = kwargs.get(
            'w_field', get_field_name('vertical_wind_component'))
        self.div_field = kwargs.get(
            'div_field', get_field_name('wind_divergence'))
        self.hdiv_field = kwargs.get(
            'hdiv_field', get_field_name('horizontal_wind_divergence'))

        # Atmospheric state variable field names
        self.tdry_field = kwargs.get('tdry_field', 'temperature')
        self.rho_field = kwargs.get('rho_field', 'density')

        # Radar pointing direction field names
        self.range_field = kwargs.get('range_field', 'range')
        self.azimuth_field = kwargs.get('azimuth_field', 'azimuth')
        self.elevation_field = kwargs.get('elevation_field', 'elevation')
        self.xhat_field = kwargs.get('xhat_field', 'x_component')
        self.yhat_field = kwargs.get('yhat_field', 'y_component')
        self.zhat_field = kwargs.get('zhat_field', 'z_component')

        # Miscellaneous radar field names
        self.cover_field = kwargs.get('cover_field', 'radar_coverage')

        return
