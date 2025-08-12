"""
winds.core.metric
=================

A class for storing basic checks and metrics on wind retrievals derived from
one or more scanning Doppler radars. This class is designed to give the user
immediate feedback on wind retrieval quality.

"""

import numpy as np

from pyart.core import Grid
from pyart.config import get_fillvalue, get_field_name, get_metadata

from ..util import standard
from ..var import _operators, _continuity


class Metric(object):
    """
    A class for storing and reporting basic checks and metrics on wind
    retrievals derived from one or more scanning Doppler radars.

    Attributes
    ----------
    solver : Solver
        Solver and its associated attributes, including analysis domain
        parameters and functional minimization parameters.
    u : ndarray
        Eastward wind component in meters per second.
    v : ndarray
        Northward wind component in meters per second.
    w : ndarray
        Vertical wind component in meters per second.
    anelastic : bool
        True if the anelastic assumption for mass continuity has been used,
        False otherwise.
    alpha : float
        Relative mass continuity residual computed over the entire analysis
        domain. This provides a metric for the degree to which the wind field
        globally satisfies mass continuity.
    alpha_profile : ndarray
        Relative mass continuity residual as a function of height. This
        provides a metric for the degree to which the wind field satisfies mass
        continuity at each height level.
    impermeability_mbe : float
        Vertical velocity mean bias error (MBE) to surface impermeability.
    impermeability_mae : float
        Vertical velocity mean absolute error (MAE) to surface impermeability.
    impermeability_rmse : float
        Vertical velocity root mean squared error (RMSE) to surface
        impermeability.

    """
    def __init__(self, conv, grids, solver, profile, debug=False,
                 verbose=False):
        """ Initialize.

        Parameters
        ----------
        conv : Grid
            Py-ART grid containing (u, v, w) wind components.
        grids : list or tuple
            Input grids used to retrieve the wind field.
        solver : Solver
            Solver used to retrieve the wind field.
        profile : Profile, optional
            Profile used to retrieve the wind field. If None, the U.S. Standard
            Atmosphere model is used.
        debug : bool, optional
            True to print debugging information, False to suppress.
        verbose : bool, optional
            True to print progress and indentification information, False to
            suppress.

        """

        # Default metric attributes
        self.solver = solver

        # Parse input grids
        if isinstance(grids, Grid):
            grids = [grids]

        # Parse atmospheric profile
        if profile is None:
            profile = standard.atmosphere(
                grid=conv, debug=debug, verbose=verbose)

        # Wind component attributes
        self.u = conv.fields[self.solver.u_field]['data']
        self.v = conv.fields[self.solver.v_field]['data']
        self.w = conv.fields[self.solver.w_field]['data']

        # Parse mass continuity equation, e.g., determine if the anelastic
        # approximation was used
        methods = ['POTVIN', 'PROTAT', 'BOTTOM-UP', 'TOP-DOWN']
        if self.solver.continuity.upper() in methods:
            self.anelastic = True
        else:
            self.anelastic = False

        # Parse functional minimization results
        self.functional_minimization(debug=debug, verbose=verbose)

        # Compute relative continuity residual
        self.continuity_residual(profile, debug=debug, verbose=verbose)

        # Compute radar observation residuals
        self.radar_observation_residual(grids, debug=debug, verbose=verbose)

        # Compute surface impermeability residuals
        self.impermeability_residual(debug=debug, verbose=verbose)


    def functional_minimization(self, debug=False, verbose=False):
        """
        Parse and store functional minimization information, such as the values
        of the cost function and the magnitude of its gradient as a function of
        iteration number.

        Parameters
        ----------
        debug : bool, optional
            True to print debugging information, False to suppress.
        verbose : bool, optional
            True to print progress information, False to suppress.

        """

        # Parse value of cost function as a function of iteration
        self.cost_iter = np.arange(self.solver.cost_iter)
        self.cost_value = np.asarray(self.solver.cost_value)

        if verbose:
            min_idx = self.cost_value.argmin()
            max_idx = self.cost_value.argmax()
            min_cost = self.cost_value.min()
            max_cost = self.cost_value.max()
            print 'Min cost function index: {}'.format(min_idx)
            print 'Max cost function index: {}'.format(max_idx)
            print 'Min cost function value: {:.3e}'.format(min_cost)
            print 'Max cost function value: {:.3e}'.format(max_cost)

        # Value of cost function gradient magnitude as a function of iteration
        self.cost_grad_iter = np.arange(0, self.solver.cost_grad_iter)
        self.cost_grad_mag = np.asarray(self.solver.cost_grad_mag)

        if verbose:
            min_idx = self.cost_grad_mag.argmin()
            max_idx = self.cost_grad_mag.argmax()
            min_grad = self.cost_grad_mag.min()
            max_grad = self.cost_grad_mag.max()
            print 'Min cost gradient magnitude index: {}'.format(min_idx)
            print 'Max cost gradient magnitude index: {}'.format(max_idx)
            print 'Min cost gradient magnitude: {:.3e}'.format(min_grad)
            print 'Max cost gradient magnitude: {:.3e}'.format(max_grad)

        return


    def continuity_residual(self, profile, debug=False, verbose=False):
        """
        Compute mass continuity residuals, e.g., the relative mass continuity
        residual defined as the ratio of squared mass continuity residual with
        the sum of squares of each term in mass continuity.

        Parameters
        ----------
        profile : Profile
            Profile used to retrieve the wind field.
        debug : bool, optional
            True to print debugging information, False to suppress.
        verbose : bool, optional
            True to print progress and indentification information, False to
            suppress.

        """

        if verbose:
            print 'Computing mass continuity residual'

        # Parse atmospheric state data
        rho = profile.rho
        drhodx = profile.drhodx
        drhody = profile.drhody
        drhodz = profile.drhodz

        # Parse wind field and prepare for ingest into Fortran wrapper
        u = np.asfortranarray(self.u, dtype=np.float64)
        v = np.asfortranarray(self.v, dtype=np.float64)
        w = np.asfortranarray(self.w, dtype=np.float64)

        # Compute wind divergence field
        div, dudx, dvdy, dwdz = _operators.div_wind(
            u, v, w, dx=self.solver.dx, dy=self.solver.dy, dz=self.solver.dz,
            finite_order=self.solver.finite_order,
            fill_value=self.solver.fill_value, proc=self.solver.proc)

        if self.anelastic:
            # The anelastic mass continuity equation is given by,
            #
            # D = du/dx + dv/dy + dw/dz + ...
            #   + (u * drho/dx + v * drho/dy + w * drho/dz) / rho = 0
            #
            # The first 3 terms represent the mass flux due to wind divergence,
            # and the final 3 terms represent the mass flux due to advection.
            #
            # We define the relative mass continuity residual as the ratio
            # of D squared with the square of each term in the mass continuity
            # equation
            D = dudx + dvdy + dwdz + \
                u * drhodx / rho + \
                v * drhody / rho + \
                w * drhodz / rho
            norm = dudx**2 + dvdy**2 + dwdz**2 + \
                   (u * drhodx / rho)**2 + \
                   (v * drhody / rho)**2 + \
                   (w * drhodz / rho)**2

            alpha = np.sqrt(np.mean(D**2 / norm))
            alpha_profile = np.sqrt(np.mean(D**2 / norm, axis=(1, 2)))

        if verbose:
            print 'Relative mass continuity residual: {:.3f}'.format(alpha)

        self.alpha = alpha
        self.alpha_profile = alpha_profile

        return


    def radar_observation_residual(self, grids, debug=False, verbose=False):
        """
        Compute radar Doppler (radial) velocity observation residuals, e.g.,
        bias, absolute, and squared errors. Wind retrievals that diverge too
        far from input radar observations should be considered suspect.

        Parameters
        ----------
        grids : list or tuple
            Input grids used to derive the wind field.
        debug : bool, optional
            True to print debugging information, False to suppress.
        verbose : bool, optional
            True to print progress and indentification information, False to
            suppress.

        """

        if verbose:
            print 'Computing radar observation residuals'

        # Parse analysis domain dimensions
        nz, ny, nx = self.solver.shape

        for i, grid in enumerate(grids):

            # Parse grid instrument name
            grid_name = grid.metadata['radar_0_instrument_name']
            if not grid_name:
                grid_name = 'Grid_{}'.format(i)

            # Parse radar data
            vdop_obs = grid.fields[self.solver.vdop_field]['data']
            vfall = grid.fields[self.solver.fall_field]['data']
	    if(len(vdop_obs.shape)==2): #added by oue RHI grid X-Z plane
		vdop_obs = np.expand_dims(vdop_obs,axis=1)
		vfall = np.expand_dims(vfall,axis=1)
            ic = grid.fields[self.solver.xhat_field]['data']
            jc = grid.fields[self.solver.yhat_field]['data']
            kc = grid.fields[self.solver.zhat_field]['data']

            # Compute the projected Doppler (radial) velocity of the wind field
            vdop = self.u * ic + self.v * jc + (self.w - vfall) * kc

            # Compute the bias between wind field Doppler velocity and radar
            # Doppler velocity observations
            bias = vdop - vdop_obs

            # Compute residuals over entire analysis domain
            mbe = np.ma.mean(bias)
            mae = np.ma.mean(np.abs(bias))
            rmse = np.ma.sqrt(np.ma.mean(bias**2))

            if verbose:
                print 'Grid instrument name: {}'.format(grid_name)
                print 'Observation MBE: {:.3f} m/s'.format(mbe)
                print 'Observation MAE: {:.3f} m/s'.format(mae)
                print 'Observation RMSE: {:.3f} m/s'.format(rmse)

            # Compute residuals as a function of height
            bias = bias.reshape(nz, ny * nx)
            mbe_profile = np.ma.mean(bias, axis=1)
            mae_profile = np.ma.mean(np.abs(bias), axis=1)
            rmse_profile = np.ma.sqrt(np.ma.mean(bias**2, axis=1))

            results = {
                'vdop_mbe': mbe,
                'vdop_mae': mae,
                'vdop_rmse': rmse,
                'vdop_mbe_profile': mbe_profile,
                'vdop_mae_profile': mae_profile,
                'vdop_rmse_profile': rmse_profile,
                }

            # Store results as a new attribute
            setattr(self, grid_name.replace('-','').replace(' ', '_'), results)

        return


    def impermeability_residual(self, debug=False, verbose=False):
        """
        Compute residuals, e.g., bias, absolute, and squared errors of the
        surface impermeability condition, i.e., vertical velocity vanishing at
        the surface.

        Parameters
        ----------
        debug : bool, optional
            True to print debugging information, False to suppress.
        verbose : bool, optional
            True to print progress information, False to suppress.

        """

        if verbose:
            print 'Computing surface impermeability residual'

        # Compute surface impermeability residuals
        bias = self.w[0]
        mbe = np.mean(bias)
        mae = np.mean(np.abs(bias))
        rmse = np.sqrt(np.mean(bias**2))

        if verbose:
            print 'Impermeability MBE: {:.3f} m/s'.format(mbe)
            print 'Impermeability MAE: {:.3f} m/s'.format(mae)
            print 'Impermeability RMSE: {:.3f} m/s'.format(rmse)

        self.impermeability_mbe = mbe
        self.impermeability_mae = mae
        self.impermeability_rmse = rmse

        return
