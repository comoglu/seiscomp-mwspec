Mw(spec) is a per-station moment magnitude derived from a Brune
:math:`\omega^2` fit of the ground-displacement spectrum. It is a port of the
Seisan SPEC/AUTOMAG spectral-magnitude algorithm (L. Ottemoeller) and fills the
one magnitude type that Seisan provides but SeisComP did not: spectral
:math:`M_w` computed independently at each station, as opposed to scaling from
mB/Mwp or full moment-tensor inversion.

The plugin contributes both an amplitude processor and a magnitude processor of
type ``Mw(spec)``.

Amplitude
=========

For each contributing pick the phase window is deconvolved to ground
displacement and Fourier transformed. The displacement amplitude spectrum of the
signal and of a preceding noise window are corrected for path attenuation,

.. math::

   A_\mathrm{corr}(f) = A(f) \, \big/ \,
   \left[ e^{-\pi f\, t\, / Q(f)} \cdot e^{-\pi \kappa f} \right]

with :math:`Q(f) = Q_0 f^{\alpha}` (or :math:`Q_0(1+(f/f_c)^{\alpha})` when a Q
corner frequency is configured), travel time :math:`t` and near-surface
attenuation :math:`\kappa`. A usable frequency band is selected from the
signal-to-noise spectrum, and the Brune model

.. math::

   \log_{10} A(f) = \log_{10}\Omega_0 - \log_{10}\!\left(1 + (f/f_0)^2\right)

is fitted by grid search for the low-frequency plateau :math:`\Omega_0` and the
corner frequency :math:`f_0`. The amplitude carries :math:`\Omega_0` (unit
``nm*s``); the corner frequency is carried as the amplitude period.

Magnitude
=========

The seismic moment and magnitude follow Seisan automag:

.. math::

   M_0 &= \frac{4\pi\,\rho\,c^3\,R\,\Omega_0}{F \cdot R_{\theta\phi}} \\
   M_w &= \tfrac{2}{3}\log_{10} M_0 - 6.06

where :math:`\rho` and :math:`c` are the density and phase velocity at the
source depth (interpolated from the configured layered model), :math:`R` is the
equivalent geometric-spreading distance (hypocentral for P; Herrmann-Kijko for
regional S), :math:`F` the free-surface factor and :math:`R_{\theta\phi}` the
average radiation pattern. As :math:`M_w` is already a moment magnitude, the
network Mw estimation is the identity.

Configuration
=============

Add ``mwspec`` to the ``plugins`` parameter and enable the ``Mw(spec)``
amplitude and magnitude in :ref:`scamp` / :ref:`scmag` (or :ref:`scolv`). The
shared velocity/attenuation/density model and the analysed phase are set once in
the magnitude binding (``magnitudes.Mw(spec).model`` and ``.phase``); spectral
measurement parameters are set in the amplitude binding. See the binding
parameter descriptions for details.
