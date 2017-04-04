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
        if self.structure.symmetry_flag == 1:
            self.c_11 = mat_props[2]
            self.c_12 = mat_props[3]
            self.c_44 = mat_props[4]
            self.p_11 = mat_props[5]
            self.p_12 = mat_props[6]
            self.p_44 = mat_props[7]
            self.eta_11 = mat_props[8]
            self.eta_12 = mat_props[9]
            self.eta_44 = mat_props[10]
        # Fully anisotropic case
        else:
            self.c_11 = mat_props[2]
            self.c_12 = mat_props[3]
            self.c_13 = mat_props[4]
            self.c_14 = mat_props[5]
            self.c_15 = mat_props[6]
            self.c_16 = mat_props[7]
            self.c_21 = mat_props[8]
            self.c_22 = mat_props[9]
            self.c_23 = mat_props[10]
            self.c_24 = mat_props[11]
            self.c_25 = mat_props[12]
            self.c_26 = mat_props[13]
            self.c_31 = mat_props[14]
            self.c_32 = mat_props[15]
            self.c_33 = mat_props[16]
            self.c_34 = mat_props[17]
            self.c_35 = mat_props[18]
            self.c_36 = mat_props[19]
            self.c_41 = mat_props[20]
            self.c_42 = mat_props[21]
            self.c_43 = mat_props[22]
            self.c_44 = mat_props[23]
            self.c_45 = mat_props[24]
            self.c_46 = mat_props[25]
            self.c_51 = mat_props[26]
            self.c_52 = mat_props[27]
            self.c_53 = mat_props[28]
            self.c_54 = mat_props[29]
            self.c_55 = mat_props[30]
            self.c_56 = mat_props[31]
            self.c_61 = mat_props[32]
            self.c_62 = mat_props[33]
            self.c_63 = mat_props[34]
            self.c_64 = mat_props[35]
            self.c_65 = mat_props[36]
            self.c_66 = mat_props[37]

            self.p_11 = mat_props[36+2]
            self.p_12 = mat_props[36+3]
            self.p_13 = mat_props[36+4]
            self.p_14 = mat_props[36+5]
            self.p_15 = mat_props[36+6]
            self.p_16 = mat_props[36+7]
            self.p_21 = mat_props[36+8]
            self.p_22 = mat_props[36+9]
            self.p_23 = mat_props[36+10]
            self.p_24 = mat_props[36+11]
            self.p_25 = mat_props[36+12]
            self.p_26 = mat_props[36+13]
            self.p_31 = mat_props[36+14]
            self.p_32 = mat_props[36+15]
            self.p_33 = mat_props[36+16]
            self.p_34 = mat_props[36+17]
            self.p_35 = mat_props[36+18]
            self.p_36 = mat_props[36+19]
            self.p_41 = mat_props[36+20]
            self.p_42 = mat_props[36+21]
            self.p_43 = mat_props[36+22]
            self.p_44 = mat_props[36+23]
            self.p_45 = mat_props[36+24]
            self.p_46 = mat_props[36+25]
            self.p_51 = mat_props[36+26]
            self.p_52 = mat_props[36+27]
            self.p_53 = mat_props[36+28]
            self.p_54 = mat_props[36+29]
            self.p_55 = mat_props[36+30]
            self.p_56 = mat_props[36+31]
            self.p_61 = mat_props[36+32]
            self.p_62 = mat_props[36+33]
            self.p_63 = mat_props[36+34]
            self.p_64 = mat_props[36+35]
            self.p_65 = mat_props[36+36]
            self.p_66 = mat_props[36+37]

            self.eta_11 = mat_props[36+36+2]
            self.eta_12 = mat_props[36+36+3]
            self.eta_13 = mat_props[36+36+4]
            self.eta_14 = mat_props[36+36+5]
            self.eta_15 = mat_props[36+36+6]
            self.eta_16 = mat_props[36+36+7]
            self.eta_21 = mat_props[36+36+8]
            self.eta_22 = mat_props[36+36+9]
            self.eta_23 = mat_props[36+36+10]
            self.eta_24 = mat_props[36+36+11]
            self.eta_25 = mat_props[36+36+12]
            self.eta_26 = mat_props[36+36+13]
            self.eta_31 = mat_props[36+36+14]
            self.eta_32 = mat_props[36+36+15]
            self.eta_33 = mat_props[36+36+16]
            self.eta_34 = mat_props[36+36+17]
            self.eta_35 = mat_props[36+36+18]
            self.eta_36 = mat_props[36+36+19]
            self.eta_41 = mat_props[36+36+20]
            self.eta_42 = mat_props[36+36+21]
            self.eta_43 = mat_props[36+36+22]
            self.eta_44 = mat_props[36+36+23]
            self.eta_45 = mat_props[36+36+24]
            self.eta_46 = mat_props[36+36+25]
            self.eta_51 = mat_props[36+36+26]
            self.eta_52 = mat_props[36+36+27]
            self.eta_53 = mat_props[36+36+28]
            self.eta_54 = mat_props[36+36+29]
            self.eta_55 = mat_props[36+36+30]
            self.eta_56 = mat_props[36+36+31]
            self.eta_61 = mat_props[36+36+32]
            self.eta_62 = mat_props[36+36+33]
            self.eta_63 = mat_props[36+36+34]
            self.eta_64 = mat_props[36+36+35]
            self.eta_65 = mat_props[36+36+36]
            self.eta_66 = mat_props[36+36+37]

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
