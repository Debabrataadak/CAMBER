# claude: helpers for time-varying cosmic birefringence (whole file added)
r"""
Time-varying cosmic birefringence from an axion-like field.

A spectator pseudoscalar :math:`\phi` with potential :math:`V=m^2\phi^2/2`, coupled to photons by
:math:`-1/4 g\phi F\tilde F`, rotates the polarization of a photon emitted at conformal time :math:`\eta` by

.. math:: \beta(\eta) = \frac{g}{2}\left[\phi(\eta_0)-\phi(\eta)\right] = \beta_{\rm ini}\,[f(\eta)-f(\eta_0)],
          \qquad f=\phi/\phi_{\rm ini},\quad \beta_{\rm ini}=-g\phi_{\rm ini}/2

(the convention of Namikawa, Murai & Naokawa, arXiv:2506.20824). The field obeys
:math:`\ddot\phi + 3H\dot\phi + m^2\phi = 0` in the background of the given cosmology, starting frozen at
:math:`\phi_{\rm ini}` deep in the radiation era; its energy density is neglected.

Use :func:`axion_birefringence_table` to get :math:`\beta(z)` and pass it to
:meth:`.model.CAMBparams.set_birefringence_table`.
"""

import numpy as np
from scipy.integrate import solve_ivp
from scipy.interpolate import CubicSpline

from .camb import get_background

# hbar * c in eV Mpc, and c in km/s: H[eV] = H[km/s/Mpc] / c * hbar_c
_HBAR_C_EV_MPC = 6.39488e-30
_C_KM_S = 299792.458


def hubble_ev(results, z):
    """Hubble rate H(z) in eV from a CAMB results (or background) object."""
    return results.hubble_parameter(np.asarray(z)) / _C_KM_S * _HBAR_C_EV_MPC


def solve_axion_field(results, mass_ev, a_ini=1e-8, wkb_ratio=1000.0, rtol=1e-10):
    r"""
    Solve :math:`f'' + (3 + d\ln H/d\ln a) f' + (m/H)^2 f = 0` in :math:`x = \ln a`, with :math:`f = 1, f' = 0` at
    ``a_ini``. Once :math:`m/H` exceeds ``wkb_ratio`` the field oscillates with an amplitude decaying as
    :math:`a^{-3/2}`; integration stops there and :math:`f` is taken as zero afterwards. The neglected oscillation
    has amplitude :math:`\sim (H/m)^{3/4}` in the radiation era and :math:`\sim (H/m)` in the matter era (8e-3 and
    1.5e-3 for the default); it is fast and averages out in the line-of-sight integral (raising ``wkb_ratio`` to
    1e4 changes EB by < 1e-7 for m <= 1e-27 eV).

    :return: (dense solution of (f, f') in x, x where integration stopped, f today)
    """
    x_grid = np.linspace(np.log(a_ini), 0.0, 6000)
    lnH = CubicSpline(x_grid, np.log(hubble_ev(results, np.exp(-x_grid) - 1)))
    dlnH = lnH.derivative()

    def rhs(x, y):
        return [y[1], -(3 + dlnH(x)) * y[1] - (mass_ev / np.exp(lnH(x))) ** 2 * y[0]]

    def heavy(x, y):
        return mass_ev / np.exp(lnH(x)) - wkb_ratio

    heavy.terminal = True
    heavy.direction = 1
    sol = solve_ivp(
        rhs,
        (x_grid[0], 0.0),
        [1.0, 0.0],
        method="DOP853",
        rtol=rtol,
        atol=1e-13,
        dense_output=True,
        events=heavy,
        max_step=0.05,
    )
    if not sol.success:
        raise RuntimeError(sol.message)
    x_end = sol.t[-1]
    f0 = sol.y[0, -1] if x_end >= -1e-12 else 0.0
    return sol.sol, x_end, f0


def axion_birefringence_table(pars, mass_ev, beta_ini, z=None, results=None):
    r"""
    Rotation angle :math:`\beta(z)` in degrees for an axion-like field of mass ``mass_ev`` (eV) and amplitude
    ``beta_ini`` (degrees, :math:`=-g\phi_{\rm ini}/2`), on a redshift grid for
    :meth:`.model.CAMBparams.set_birefringence_table`.

    :param pars: :class:`.model.CAMBparams` giving the background cosmology
    :param mass_ev: field mass in eV
    :param beta_ini: amplitude in degrees
    :param z: increasing redshift grid; default is 0 to 1e5 with :math:`\Delta\ln(1+z) = 5\times10^{-4}`, enough to
              resolve oscillations during recombination for masses up to about 1e-26 eV
    :param results: optional precomputed :func:`.camb.get_background` results for ``pars`` (saves time in scans)
    :return: (z, beta) arrays
    """
    if results is None:
        results = get_background(pars)
    if z is None:
        z = np.expm1(np.arange(0.0, np.log(1e5 + 1), 5e-4))
    z = np.asarray(z, dtype=np.float64)
    sol, x_end, f0 = solve_axion_field(results, mass_ev)
    x = -np.log1p(z)
    f = np.zeros_like(z)
    ok = x <= x_end
    if np.any(ok):
        f[ok] = sol(x[ok])[0]
    # f(z) = 0 after the field starts fast oscillations (z below the switch); before a_ini it is frozen at 1
    f[x < np.log(1e-8)] = 1.0
    return z, beta_ini * (f - f0)
