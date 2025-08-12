"""
winds.util.standard
===================

A module for creating standard atmospheric and meteorological profiles. Most of
this module relies on the U.S. Standard Atmosphere model.

"""

import numpy as np
import pandas as pd

from StringIO import StringIO
from scipy.interpolate import UnivariateSpline

from ..core import Profile


US_STANDARD = """
altitude  temperature   pressure     density    speed_sound    viscosity
0.00      15.000        101.325      1.225000   340.294        1.81206e-5
0.50      11.750        95.4608      1.167270   338.370        1.79579e-5
1.00      8.5000        89.8746      1.111640   336.434        1.77943e-5
1.50      5.2500        84.5560      1.058070   334.487        1.76298e-5
2.00      2.0000        79.4952      1.006490   332.529        1.74645e-5
2.50      -1.250        74.6825      0.956859   330.560        1.72983e-5
3.00      -4.500        70.1085      0.909122   328.578        1.71311e-5
3.50      -7.750        65.7641      0.863229   326.584        1.69630e-5
4.00      -11.00        61.6402      0.819129   324.579        1.67940e-5
4.50      -14.25        57.7283      0.776775   322.560        1.66241e-5
5.00      -17.50        54.0199      0.736116   320.529        1.64531e-5
5.50      -20.75        50.5068      0.697106   318.486        1.62813e-5
6.00      -24.00        47.1810      0.659679   316.428        1.61084e-5
6.50      -27.25        44.0348      0.623844   314.358        1.59345e-5
7.00      -30.50        41.0607      0.589501   312.274        1.57596e-5
7.50      -33.75        38.2514      0.556624   310.175        1.55837e-5
8.00      -37.00        35.5998      0.525168   308.063        1.54068e-5
8.50      -40.25        33.0990      0.495090   305.935        1.52288e-5
9.00      -43.50        30.7425      0.466348   303.793        1.50498e-5
9.50      -46.75        28.5236      0.438901   301.636        1.48696e-5
10.0      -50.00        26.4363      0.412707   299.463        1.46884e-5
10.5      -53.25        24.4744      0.387725   297.275        1.45061e-5
11.0      -56.50        22.6321      0.363918   295.070        1.43226e-5
11.5      -56.50        20.9162      0.336327   295.070        1.43226e-5
12.0      -56.50        19.3304      0.310828   295.070        1.43226e-5
12.5      -56.50        17.8648      0.287262   295.070        1.43226e-5
13.0      -56.50        16.5104      0.265483   295.070        1.43226e-5
13.5      -56.50        15.2587      0.245355   295.070        1.43226e-5
14.0      -56.50        14.1018      0.226753   295.070        1.43226e-5
14.5      -56.50        13.0327      0.209562   295.070        1.43226e-5
15.0      -56.50        12.0446      0.193674   295.070        1.43226e-5
"""


def atmosphere(grid=None, debug=False, verbose=False):
    """
    Define standard atmospheric temperature, pressure, and density profiles.
    These profiles are derived from the U.S. Standard Atmosphere model. Output
    wind fields are simple zero wind fields.

    Note that the discontinuity in the temperature lapse rate, i.e., from the
    dry adiabatic lapse rate to zero lapse rate, produces poor results when
    using more sophisticated methods (e.g., cubic splines) to compute rates of
    change with altitude.

    Parameters
    ----------
    grid : Grid, optional
        Py-ART Grid defining the altitudes to output atmospheric profiles. If
        no grid is available, altitude data up to 15 km AMSL is used.
    debug : bool, optional
        True to print debugging information, False to suppress.
    verbose : bool, optional
        True to print progress and identification information, False to
        suppress.

    Returns
    -------
    profile : Profile
        Profile containing altitude in meters AMSL, temperature in Celsius,
        pressure in kilopascals, air density in kilograms per cubic meter, and
        a few other standard measurements.

    """

    # Parse U.S. Standard Atmosphere data
    data = pd.read_csv(
        StringIO(US_STANDARD), header='infer', delim_whitespace=True)
    alt = data.altitude.values * 1000.0
    tdry = data.temperature.values
    pres = data.pressure.values
    rho = data.density.values
    cs = data.speed_sound.values

    # Interpolate standard atmosphere data to grid heights using a cubic spline
    if grid is not None:
        alt_grid = grid.axes['z_disp']['data'] + grid.axes['alt']['data'][0]
        tdry = UnivariateSpline(alt, tdry, k=3, s=None).__call__(alt_grid)
        pres = UnivariateSpline(alt, pres, k=3, s=None).__call__(alt_grid)
        rho = UnivariateSpline(alt, rho, k=3, s=None).__call__(alt_grid)
        cs = UnivariateSpline(alt, cs, k=3, s=None).__call__(alt_grid)
        alt = alt_grid

    if debug:
        print 'First level temperature: {:2f} deg C'.format(tdry[0])
        print 'First level pressure: {:.2f} kPa'.format(pres[0])
        print 'First level air density: {:.2f} kg/m^3'.format(rho[0])

    # Define zero wind fields as default
    wspd = np.zeros_like(tdry, subok=False)
    u = np.zeros_like(tdry, subok=False)
    v = np.zeros_like(tdry, subok=False)
    w = np.zeros_like(tdry, subok=False)

    return Profile(alt=alt, tdry=tdry, pres=pres, rho=rho, wspd=wspd, u=u, v=v,
                   w=w, cs=cs, rate_method='finite')
