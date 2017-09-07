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
 


Vacuum = Material("Vacuum")

Si_2016_Smith = Material("Si_2016_Smith")
Si_2015_Van_Laer = Material("Si_2015_Van_Laer")
Si_2013_Laude = Material("Si_2013_Laude")
Si_test_anisotropic = Material("Si_test_anisotropic")

SiO2_2016_Smith = Material("SiO2_2016_Smith")
SiO2_2015_Van_Laer = Material("SiO2_2015_Van_Laer")
SiO2_2013_Laude = Material("SiO2_2013_Laude")

As2S3_2017_Morrison = Material("As2S3_2017_Morrison")
As2S3_2016_Smith = Material("As2S3_2016_Smith")

GaAs_2016_Smith = Material("GaAs_2016_Smith")
