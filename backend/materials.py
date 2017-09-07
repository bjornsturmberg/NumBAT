"""
    materials.py is a subroutine of NumBAT that defines Material objects,
    these represent dispersive lossy refractive indices and possess
    methods to interpolate n from tabulated data.

    Copyright (C) 2017  Bjorn Sturmberg, Kokou Dossou
"""


data_location = '../backend/material_data/'


class Material(object):
    """ Represents a material with:

            Refractive index []
            Density [kg/m3]
            Stiffness tensor component [Pa]
            Photoelastic tensor component []
            Acoustic loss tensor component [Pa s]

    """
 

def isotropic_stiffness(E, v):
    """
    Calculate the stiffness matrix components of isotropic 
    materials, given the two free parameters.

    Ref: www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_isotropic.cfm

    Args:
        E  (float): Youngs_modulus

        v  (float): Poisson_ratio
    """
    c_11 = E*(1-v)/((1+v)*(1-2*v))
    c_12 = E*(v)/((1+v)*(1-2*v))
    c_44 = (E*(1-2*v)/((1+v)*(1-2*v)))/2

    return c_11, c_12, c_44
