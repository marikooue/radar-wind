"""
winds.core.profile
==================

A class for handling profile data, e.g., radiosonde measurements. This class
can be used to store basic atmospheric profile measurements such as
temperature, pressure, and horizontal wind speed and direction, but there are
attributes for numerous other measurements such as water vapour pressure,
relative humidity, radio refractivity, and the speed of sound, just to name a
few.

"""

import numpy as np

from scipy.interpolate import UnivariateSpline

from pyart.config import get_fillvalue


class Profile(object):
    """
    An object for storing atmospheric profiles, e.g., radiosonde measurements.
    Rates of change with respect to height or altitude are automatically
    computed for all available measurements.

    Attributes
    ----------
    alt : ndarray
        Altitude above mean sea level (AMSL) in meters.
    tdry : ndarray
        Dry bulb temperature (air temperature) in degrees Celsius.
    dp : ndarray
        Dew point temperature in degrees Celsius.
    bar_pres : ndarray
        Barometric (atmospheric) pressure in kilopascals.
    rho : ndarray
        Air density in kilograms per cubic meter.
    wspd : ndarray
        Wind speed in meters per second.
    wdir : ndarray
        Wind direction in degrees assuming meteorological convention.
    u : ndarray
        Eastward wind component in meters per second.
    v : ndarray
        Northward wind component in meters per second.
    w : ndarray
        Vertical wind component in meters per second. This measurement is
        typically not available from radiosonde and other profiling
        instruments or techniques.
    rh : ndarray
        Relative humidity in percent. Relative humidity is the ratio of the
        partial pressure of water vapour to the equilibrium vapour pressure of
        water at the same temperature.
    sh : ndarray
        Specific humidity in grams of water vapour per grams of total air mass
        including dry mass. Specific humidity is also known as the moisture
        content.
    e : ndarray
        Water vapor pressure in hectopascals.
    cs : ndarray
        Speed of sound in meters per second. In dry air at 20 deg C, the speed
        of sound is approximately 343.2 meters per second.
    rate_method : 'cubic_spline' or 'finite', str
        The method used to compute rates of change with respect to altitude.
        If 'cubic_spline', a spline of degree 3 is first used to fit the
        sampled data, whereby the derivative of the fit is used to compute
        rates of change with height. If 'finite', first and second order finite
        differences are used to compute rates of change with height.
    finite_order : 'low' or 'high', str

    Examples
    --------

    """
    def __init__(self, alt=None, tdry=None, dp=None, pres=None, rho=None,
                 wspd=None, wdir=None, u=None, v=None, w=None, rh=None,
                 sh=None, e=None, cs=None, rate_method='cubic',
                 finite_order='low'):
        """ Initialize. """

        # Parse default profile measurements
        self.alt = alt
        self.tdry = tdry
        self.dp = dp
        self.pres = pres
        self.rho = rho
        self.wspd = wspd
        self.wdir = wdir
        self.u = u
        self.v = v
        self.w = w
        self.rh = rh
        self.sh = sh
        self.e = e
        self.cs = cs
        self.rate_method = rate_method
        self.finite_order = finite_order

        # Compute rates of change with respect to altitude
        self.compute_rates()


    def compute_rates(self):
        """ Compute the rate of change with respect to altitude. """

        for attr, data in self.__dict__.items():

            # Check if attribute is a measurement
            if attr != 'alt' and isinstance(data, np.ndarray):

                # Compute rate of change of variable with respect to altitude
                rate = _rate_of_change(self.alt, data, self.rate_method)

                # Define new attribute to reflect the rate of change of
                # the variable with respect to altitude
                setattr(self, 'd{}dz'.format(attr), rate)

        return


    def add_grid_measurement(
            self, grid, field, attr, inverse=False, fill_value=None,
            debug=False, verbose=False):
        """
        Add a measurement field from a grid.

        Parameters
        ----------
        grid : pyart.core.Grid
            Grid containing input measurement field and axes data.
        field : str
            The measurement field name.
        attr : str
            The corresponding profile attribute for the grid measurement.
        inverse : bool, optional
            True to invert the grid measurement field.
        fill_value : float, optional
            The value indicating bad or missing data in the grid field.
        debug : bool, optional
            True to print debugging information, False to suppress.
        verbose : bool, optional
            True to print relevant information, False to suppress.

        """

        from ..var import _operators
        from ..util import common

        # Parse fill value
        if fill_value is None:
            fill_value = get_fillvalue()

        # Parse grid data and prepare for ingest into Fortran wrapper
        data = grid.fields[field]['data']
        if inverse:
            data = 1.0 / data

        # Parse grid resolutions
        dx, dy, dz = common.parse_grid_resolution(
            grid, check_uniform=True, debug=debug, verbose=verbose)

        # Prepare data for ingest into Fortran wrapper
        data = np.asfortranarray(data, dtype=np.float64)

        # Compute the gradient of the grid field
        gradient = _operators.gradient(
            data, dx=dx, dy=dy, dz=dz, finite_order=self.finite_order,
            fill_value=fill_value, proc=1)

        # Save the field and its gradient components to appropriate attributes
        setattr(self, attr, data)
        setattr(self, 'd{}dx'.format(attr), gradient[0])
        setattr(self, 'd{}dy'.format(attr), gradient[1])
        setattr(self, 'd{}dz'.format(attr), gradient[2])

        return


    def extrapolate_to_grid(
            self, grid, fill_value=None, debug=False, verbose=False):
        """ Extrapolate profile measurements to grid. """

        from ..var import _operators
        from ..util import common

        # Parse fill value
        if fill_value is None:
            fill_value = get_fillvalue()

        # Parse grid dimensions
	if (len(grid.fields[grid.fields.keys()[0]]['data'].shape)==2): #added by oue RHI grid X-Z plane
		nz,  nx = grid.fields[grid.fields.keys()[0]]['data'].shape
		ny=1
	else:
        	nz, ny, nx = grid.fields[grid.fields.keys()[0]]['data'].shape

        # Parse grid resolutions
        dx, dy, dz = common.parse_grid_resolution(
            grid, check_uniform=True, debug=debug, verbose=verbose)

        for attr, data in self.__dict__.items():

            # Ignore altitude and rate of change attributes
            if attr == 'alt' or 'dz' in attr:
                continue

            # Ignore non-array attributes
            if not isinstance(data, np.ndarray):
                continue

            # Ignore arrays with more than one dimension
            if isinstance(data, np.ndarray) and data.ndim > 1:
                continue

            # Extend profile data throughout grid and prepare data for
            # ingest into Fortran wrapper
            data = data.repeat(ny * nx).reshape(nz, ny, nx)
            data = np.asfortranarray(data, dtype=np.float64)

            # Compute gradient of extended profile
            gradient = _operators.gradient(
                data, dx=dx, dy=dy, dz=dz, finite_order=self.finite_order,
                fill_value=fill_value, proc=1)

            # Save the field and its gradient components to appropriate
            # attributes
            setattr(self, attr, data)
            setattr(self, 'd{}dx'.format(attr), gradient[0])
            setattr(self, 'd{}dy'.format(attr), gradient[1])
            setattr(self, 'd{}dz'.format(attr), gradient[2])

        return


def _rate_of_change(alt, data, method='cubic'):
    """
    Compute the rate of change with respect to altitude of the sampled data.
    This routine supports variable sample spacings.

    Parameters
    ----------
    alt : ndarray
        Altitude in meters of the sampled data.
    data : ndarray
        Sampled data at the specified altitudes.
    method : 'cubic' or 'finite', optional
        Method used to compute rates of change with altitude. Default method
        uses a cubic spline to estimate the function describing the sampled
        data as a function of altitude, which then provides an analytical
        solution to the rate of change of the sampled data with respect to
        altitude.

    Returns
    -------
    rate : ndarray
        The rate of change with respect to altitude of the sampled data.

    """

    # Case where no data is available
    if alt is None or data is None:
        return None

    if method.upper() == 'CUBIC':
        rate = UnivariateSpline(alt, data, k=3, s=None).__call__(alt, nu=1)

    elif method.upper() == 'FINITE':
        rate = np.zeros_like(data, subok=False)
        for k in range(data.size):
            if k > 0 and k < data.size - 2:
                rate[k] = (data[k+1] - data[k-1]) / (alt[k+1] - alt[k-1])
            elif k == 0:
                rate[k] = (data[k+1] - data[k]) / (alt[k+1] - alt[k])
            else:
                rate[k] = (data[k] - data[k-1]) / (alt[k] - alt[k-1])

    else:
        raise ValueError('Unsupported rate of change method')

    return rate
