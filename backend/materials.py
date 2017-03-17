"""
    materials.py is a subroutine of NumBAT that defines Material objects,
    these represent dispersive lossy refractive indices and possess
    methods to interpolate n from tabulated data.

    Copyright (C) 2016  Bjorn Sturmberg, Kokou Dossou
"""

import numpy as np
from scipy.interpolate import interp1d
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

data_location = '../backend/data/'


class Material(object):
    """ Represents a material with a refractive index n.

        TODO

        Args:
            todo

        Currently included materials are;

        .. tabularcolumns:: |c|

        +--------------------+
        | **Semiconductors** |
        +--------------------+
        |    Si              |
        +--------------------+
        |    SiO2            |
        +--------------------+
        |    As2S3           |
        +--------------------+
        |    GaAs            |
        +--------------------+

    """
    def __init__(self, mat_props):
        self.n = mat_props[0]
        self.s = mat_props[1]
        self.c_11 = mat_props[2]
        self.c_12 = mat_props[3]
        self.c_44 = mat_props[4]
        self.p_11 = mat_props[5]
        self.p_12 = mat_props[6]
        self.p_44 = mat_props[7]
        self.eta_11 = mat_props[8]
        self.eta_12 = mat_props[9]
        self.eta_44 = mat_props[10]


def isotropic_stiffness(E, v):
    """
    Calculate the stiffness matrix components of isotropic 
    materials, given the two free parameters:
    E: Youngs_modulus
    v: Poisson_ratio

    Ref: http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_isotropic.cfm
    """
    c_11 = E*(1-v)/((1+v)*(1-2*v))
    c_12 = E*(v)/((1+v)*(1-2*v))
    c_44 = (E*(1-2*v)/((1+v)*(1-2*v)))/2

    return c_11, c_12, c_44


Air = Material([1,None,None,None,None,None,None,None,None,None,None])

# Silicon - http://dx.doi.org/10.1364/OL.41.002338
# Refractive index
n = 3.48
# Density
s = 2329  # kg/m3
# Stiffness tensor components.
c_11 = 165.6e9; c_12 = 63.9e9; c_44 = 79.5e9  # Pa
# Photoelastic tensor components
p_11 = -0.094; p_12 = 0.017; p_44 = -0.051
# Acoustic loss tensor components.
eta_11 = 5.9e-3 ; eta_12 = 5.16e-3 ; eta_44 = 0.62e-3  # Pa s
# Put acoustic parameters together for convenience.
Si = Material([n, s, c_11, c_12, c_44, p_11, p_12, p_44,
               eta_11, eta_12, eta_44])

# Silica - http://dx.doi.org/10.1364/OL.41.002338
n = 1.45
s = 2200  # kg/m3
c_11 = 78.6e9; c_12 = 16.1e9; c_44 = 31.2e9
p_11 = 0.12; p_12 = 0.27; p_44 = -0.075
eta_11 = 1.6e-3 ; eta_12 = 1.29e-3 ; eta_44 = 0.16e-3  # Pa s
SiO2 = Material([n, s, c_11, c_12, c_44, p_11, p_12, p_44,
               eta_11, eta_12, eta_44])

# As2S3 - http://dx.doi.org/10.1364/OL.41.002338
n = 2.37
s = 3200  # kg/m3
c_11 = 18.7e9; c_12 = 6.1e9; c_44 =6.4e9 # Pa
p_11 = 0.25; p_12 = 0.24; p_44 = 0.005
eta_11 = 1.8e-3 ; eta_12 = 1.45e-3 ; eta_44 = 0.18e-3  # Pa s
As2S3 = Material([n, s, c_11, c_12, c_44, p_11, p_12, p_44,
               eta_11, eta_12, eta_44])

# GaAs - http://dx.doi.org/10.1364/OL.41.002338
n = 3.37
s = 5320  # kg/m3
c_11 = 119e9; c_12 = 53.4e9; c_44 = 59.6e9
p_11 = -0.165; p_12 = -0.14; p_44 = -0.072
eta_11 = 7.49e-3 ; eta_12 = 0.72e-3 ; eta_44 = 0.72e-3  # Pa s
GaAs = Material([n, s, c_11, c_12, c_44, p_11, p_12, p_44,
               eta_11, eta_12, eta_44])
