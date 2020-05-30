r""" """

from triqs.gf.descriptor_base import Base, Function
from triqs.gf.meshes import MeshImFreq, MeshReFreq

from mesh_refreq_pts import MeshReFreqPts

from math import copysign, pi

import warnings
import numpy

##################################################

class SemiCircular (Base):
    r"""Hilbert transform of a semicircular density of states, i.e.

     .. math::
        g(z) = \int \frac{A(\omega)}{z-\omega} d\omega

    where :math:`A(\omega) = \theta( D - |\omega|) 2 \sqrt{ D^2 - \omega^2}/(\pi D^2)`.

    (Only works in combination with frequency Green's functions.)
    """
    def __init__ (self, half_bandwidth, chem_potential=0.):
        """:param half_bandwidth: :math:`D`, the half bandwidth of the
semicircular density of states
        :param chem_potential: :math:`\mu`, the chemical potential of the    |
semicircular density of states, corresponds to minus the center of the
semicircle
"""
        Base.__init__(self, half_bandwidth=half_bandwidth, chem_potential=chem_potential)

    def __str__(self): return "SemiCircular(%s, %s)"%self.half_bandwidth, chem_potential

    def __call__(self,G):
        D = self.half_bandwidth
        mu = self.chem_potential
        Id = complex(1,0) if len(G.target_shape) == 0 else numpy.identity(G.target_shape[0],numpy.complex_)
        from cmath import sqrt
        if type(G.mesh) == MeshImFreq:
            def f(om_):
                om = om_ + mu
                return (om - 1j*copysign(1,om.imag)*sqrt(D*D - om**2))/D/D*2*Id
        elif type(G.mesh) in [MeshReFreq, MeshReFreqPts]:
            def f(om_):
              om = om_.real + mu
              if (om > -D) and (om < D):
                return (2.0/D**2) * (om - 1j* sqrt(D**2 - om**2))
              else:
                return (2.0/D**2) * (om - copysign(1,om) * sqrt(om**2 - D**2))
        else:
            raise TypeError, "This initializer is only correct in frequency"

        Id = 1. if len(G.target_shape) == 0 else numpy.identity(G.target_shape[0])
        Function(f)(G)

        return G

##################################################

class Flat (Base):
    r"""The Hilbert transform of a flat density of states, with cut-off

    .. math::
        g(z) = \int \frac{A(\omega)}{z-\omega} d\omega

    where :math:`A(\omega) = \theta( D^2 - \omega^2)/(2D)`.

    (Only works in combination with frequency Green's functions.)
    """
    def __init__ (self, half_bandwidth):
        """:param half_bandwidth: :math:`D`, the half bandwidth """
        Base.__init__(self, half_bandwidth=half_bandwidth)

    def __str__(self): return "Flat(%s)"%self.half_bandwidth

    def __call__(self,G):

        D = self.half_bandwidth
        Id = 1. if len(G.target_shape) == 0 else numpy.identity(G.target_shape[0], numpy.complex_)

        if type(G.mesh) == MeshImFreq:
            f = lambda om: (-1/(2.0*D)) * numpy.log(numpy.divide(om-D,om+D)) * Id
        elif type(G.mesh) in [MeshReFreq, MeshReFreqPts]:
            def f(om):
              if (om.real > -D) and (om.real < D):
                return -numpy.log(numpy.divide(abs(om-D),abs(om+D)))*Id/(2*D) - 1j*pi*Id/(2*D)
              else:
                return -numpy.log(numpy.divide(abs(om-D),abs(om+D)))*Id/(2*D)
        else:
            raise TypeError, "This initializer is only correct in frequency"

        # Silence "RuntimeWarning: divide by zero encountered in divide"
        old_err = numpy.seterr(divide='ignore')

        Function(f)(G)
        numpy.seterr(**old_err)
        return G

#########################################################################

class Omega_(Base):
    r"""The function:math:`\omega \rightarrow \omega` """
    def __str__(self): return "Omega"
    def __call__(self,G):
        if G.mesh.__class__.__name__ not in ['MeshImFreq', 'MeshReFreq', 'MeshReFreqPts']:
            raise TypeError, "This initializer is only correct in frequency"

        Id = 1. if G.target_rank == 0 else numpy.identity(G.target_shape[0])

        for n,om in enumerate(G.mesh): G.data[n,...] = om*Id
        return G

Omega = Omega_()
