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

import json
import re

data_location = '../backend/material_data/'


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
    # def __init__(self, mat_props,
    #              doi=None, date=None, author=None):
    #     self.doi = doi
    #     self.date = date
    #     self.author = author
    #     self.n = mat_props[0]
    #     self.s = mat_props[1]
    #     if len(mat_props) == 11:
    #         self.c_11 = mat_props[2]
    #         self.c_12 = mat_props[3]
    #         self.c_44 = mat_props[4]
    #         self.p_11 = mat_props[5]
    #         self.p_12 = mat_props[6]
    #         self.p_44 = mat_props[7]
    #         self.eta_11 = mat_props[8]
    #         self.eta_12 = mat_props[9]
    #         self.eta_44 = mat_props[10]
    #     # Fully anisotropic case
    #     else:
    #         self.c_11 = mat_props[2]
    #         self.c_12 = mat_props[3]
    #         self.c_13 = mat_props[4]
    #         self.c_14 = mat_props[5]
    #         self.c_15 = mat_props[6]
    #         self.c_16 = mat_props[7]
    #         self.c_21 = mat_props[8]
    #         self.c_22 = mat_props[9]
    #         self.c_23 = mat_props[10]
    #         self.c_24 = mat_props[11]
    #         self.c_25 = mat_props[12]
    #         self.c_26 = mat_props[13]
    #         self.c_31 = mat_props[14]
    #         self.c_32 = mat_props[15]
    #         self.c_33 = mat_props[16]
    #         self.c_34 = mat_props[17]
    #         self.c_35 = mat_props[18]
    #         self.c_36 = mat_props[19]
    #         self.c_41 = mat_props[20]
    #         self.c_42 = mat_props[21]
    #         self.c_43 = mat_props[22]
    #         self.c_44 = mat_props[23]
    #         self.c_45 = mat_props[24]
    #         self.c_46 = mat_props[25]
    #         self.c_51 = mat_props[26]
    #         self.c_52 = mat_props[27]
    #         self.c_53 = mat_props[28]
    #         self.c_54 = mat_props[29]
    #         self.c_55 = mat_props[30]
    #         self.c_56 = mat_props[31]
    #         self.c_61 = mat_props[32]
    #         self.c_62 = mat_props[33]
    #         self.c_63 = mat_props[34]
    #         self.c_64 = mat_props[35]
    #         self.c_65 = mat_props[36]
    #         self.c_66 = mat_props[37]

    #         self.p_11 = mat_props[36+2]
    #         self.p_12 = mat_props[36+3]
    #         self.p_13 = mat_props[36+4]
    #         self.p_14 = mat_props[36+5]
    #         self.p_15 = mat_props[36+6]
    #         self.p_16 = mat_props[36+7]
    #         self.p_21 = mat_props[36+8]
    #         self.p_22 = mat_props[36+9]
    #         self.p_23 = mat_props[36+10]
    #         self.p_24 = mat_props[36+11]
    #         self.p_25 = mat_props[36+12]
    #         self.p_26 = mat_props[36+13]
    #         self.p_31 = mat_props[36+14]
    #         self.p_32 = mat_props[36+15]
    #         self.p_33 = mat_props[36+16]
    #         self.p_34 = mat_props[36+17]
    #         self.p_35 = mat_props[36+18]
    #         self.p_36 = mat_props[36+19]
    #         self.p_41 = mat_props[36+20]
    #         self.p_42 = mat_props[36+21]
    #         self.p_43 = mat_props[36+22]
    #         self.p_44 = mat_props[36+23]
    #         self.p_45 = mat_props[36+24]
    #         self.p_46 = mat_props[36+25]
    #         self.p_51 = mat_props[36+26]
    #         self.p_52 = mat_props[36+27]
    #         self.p_53 = mat_props[36+28]
    #         self.p_54 = mat_props[36+29]
    #         self.p_55 = mat_props[36+30]
    #         self.p_56 = mat_props[36+31]
    #         self.p_61 = mat_props[36+32]
    #         self.p_62 = mat_props[36+33]
    #         self.p_63 = mat_props[36+34]
    #         self.p_64 = mat_props[36+35]
    #         self.p_65 = mat_props[36+36]
    #         self.p_66 = mat_props[36+37]

    #         self.eta_11 = mat_props[36+36+2]
    #         self.eta_12 = mat_props[36+36+3]
    #         self.eta_13 = mat_props[36+36+4]
    #         self.eta_14 = mat_props[36+36+5]
    #         self.eta_15 = mat_props[36+36+6]
    #         self.eta_16 = mat_props[36+36+7]
    #         self.eta_21 = mat_props[36+36+8]
    #         self.eta_22 = mat_props[36+36+9]
    #         self.eta_23 = mat_props[36+36+10]
    #         self.eta_24 = mat_props[36+36+11]
    #         self.eta_25 = mat_props[36+36+12]
    #         self.eta_26 = mat_props[36+36+13]
    #         self.eta_31 = mat_props[36+36+14]
    #         self.eta_32 = mat_props[36+36+15]
    #         self.eta_33 = mat_props[36+36+16]
    #         self.eta_34 = mat_props[36+36+17]
    #         self.eta_35 = mat_props[36+36+18]
    #         self.eta_36 = mat_props[36+36+19]
    #         self.eta_41 = mat_props[36+36+20]
    #         self.eta_42 = mat_props[36+36+21]
    #         self.eta_43 = mat_props[36+36+22]
    #         self.eta_44 = mat_props[36+36+23]
    #         self.eta_45 = mat_props[36+36+24]
    #         self.eta_46 = mat_props[36+36+25]
    #         self.eta_51 = mat_props[36+36+26]
    #         self.eta_52 = mat_props[36+36+27]
    #         self.eta_53 = mat_props[36+36+28]
    #         self.eta_54 = mat_props[36+36+29]
    #         self.eta_55 = mat_props[36+36+30]
    #         self.eta_56 = mat_props[36+36+31]
    #         self.eta_61 = mat_props[36+36+32]
    #         self.eta_62 = mat_props[36+36+33]
    #         self.eta_63 = mat_props[36+36+34]
    #         self.eta_64 = mat_props[36+36+35]
    #         self.eta_65 = mat_props[36+36+36]
    #         self.eta_66 = mat_props[36+36+37]


    def __init__(self,data_file):

        self.load_data_file(data_file)

    def load_data_file(self, data_file, alt_path=''):  
        """

        Args:
            data_file  (str): name of data file located in NumBAT/backend/material_data
            alt_path  (str): non standard path to data_file
        
        """
        with open(data_location+data_file+'.json','r') as fin:
            s_in = ''.join(fin.readlines())
            s_in = re.sub(r'//.*\n','\n', s_in)

            self._params = json.loads(s_in)

            self.file_name = self._params['file_name']  # Name of this file, will be used as identifier
            self.chemical = self._params['chemical']  # Chemical composition
            self.author = self._params['author']  # Author of data
            self.date = self._params['date']  # Year of data publication/measurement
            self.institution = self._params['institution']  # Source institution
            self.doi = self._params['doi']  # doi or, failing that, the http address

            self.n = self._params['n']  # Refractive index []
            self.s = self._params['s']  # Density [kg/m3]
            self.c_11 = self._params['c_11']  # Stiffness tensor component [Pa]
            self.c_12 = self._params['c_12']  # Stiffness tensor component [Pa]
            self.c_44 = self._params['c_44']  # Stiffness tensor component [Pa]
            self.p_11 = self._params['p_11']  # Photoelastic tensor component []
            self.p_12 = self._params['p_12']  # Photoelastic tensor component []
            self.p_44 = self._params['p_44']  # Photoelastic tensor component []
            self.eta_11 = self._params['eta_11']  # Acoustic loss tensor component [Pa s]
            self.eta_12 = self._params['eta_12']  # Acoustic loss tensor component [Pa s]
            self.eta_44 = self._params['eta_44']  # Acoustic loss tensor component [Pa s]

            try:  # full anisotropic tensor components
                self.c_11 = self._params['c_11']
                self.c_12 = self._params['c_12']
                self.c_13 = self._params['c_13']
                self.c_14 = self._params['c_14']
                self.c_15 = self._params['c_15']
                self.c_16 = self._params['c_16']
                self.c_21 = self._params['c_21']
                self.c_22 = self._params['c_22']
                self.c_23 = self._params['c_23']
                self.c_24 = self._params['c_24']
                self.c_25 = self._params['c_25']
                self.c_26 = self._params['c_26']
                self.c_31 = self._params['c_31']
                self.c_32 = self._params['c_32']
                self.c_33 = self._params['c_33']
                self.c_34 = self._params['c_34']
                self.c_35 = self._params['c_35']
                self.c_36 = self._params['c_36']
                self.c_41 = self._params['c_41']
                self.c_42 = self._params['c_42']
                self.c_43 = self._params['c_43']
                self.c_44 = self._params['c_44']
                self.c_45 = self._params['c_45']
                self.c_46 = self._params['c_46']
                self.c_51 = self._params['c_51']
                self.c_52 = self._params['c_52']
                self.c_53 = self._params['c_53']
                self.c_54 = self._params['c_54']
                self.c_55 = self._params['c_55']
                self.c_56 = self._params['c_56']
                self.c_61 = self._params['c_61']
                self.c_62 = self._params['c_62']
                self.c_63 = self._params['c_63']
                self.c_64 = self._params['c_64']
                self.c_65 = self._params['c_65']
                self.c_66 = self._params['c_66']
                self.p_11 = self._params['p_11']
                self.p_12 = self._params['p_12']
                self.p_13 = self._params['p_13']
                self.p_14 = self._params['p_14']
                self.p_15 = self._params['p_15']
                self.p_16 = self._params['p_16']
                self.p_21 = self._params['p_21']
                self.p_22 = self._params['p_22']
                self.p_23 = self._params['p_23']
                self.p_24 = self._params['p_24']
                self.p_25 = self._params['p_25']
                self.p_26 = self._params['p_26']
                self.p_31 = self._params['p_31']
                self.p_32 = self._params['p_32']
                self.p_33 = self._params['p_33']
                self.p_34 = self._params['p_34']
                self.p_35 = self._params['p_35']
                self.p_36 = self._params['p_36']
                self.p_41 = self._params['p_41']
                self.p_42 = self._params['p_42']
                self.p_43 = self._params['p_43']
                self.p_44 = self._params['p_44']
                self.p_45 = self._params['p_45']
                self.p_46 = self._params['p_46']
                self.p_51 = self._params['p_51']
                self.p_52 = self._params['p_52']
                self.p_53 = self._params['p_53']
                self.p_54 = self._params['p_54']
                self.p_55 = self._params['p_55']
                self.p_56 = self._params['p_56']
                self.p_61 = self._params['p_61']
                self.p_62 = self._params['p_62']
                self.p_63 = self._params['p_63']
                self.p_64 = self._params['p_64']
                self.p_65 = self._params['p_65']
                self.p_66 = self._params['p_66']
                self.eta_11 = self._params['eta_11']
                self.eta_12 = self._params['eta_12']
                self.eta_13 = self._params['eta_13']
                self.eta_14 = self._params['eta_14']
                self.eta_15 = self._params['eta_15']
                self.eta_16 = self._params['eta_16']
                self.eta_21 = self._params['eta_21']
                self.eta_22 = self._params['eta_22']
                self.eta_23 = self._params['eta_23']
                self.eta_24 = self._params['eta_24']
                self.eta_25 = self._params['eta_25']
                self.eta_26 = self._params['eta_26']
                self.eta_31 = self._params['eta_31']
                self.eta_32 = self._params['eta_32']
                self.eta_33 = self._params['eta_33']
                self.eta_34 = self._params['eta_34']
                self.eta_35 = self._params['eta_35']
                self.eta_36 = self._params['eta_36']
                self.eta_41 = self._params['eta_41']
                self.eta_42 = self._params['eta_42']
                self.eta_43 = self._params['eta_43']
                self.eta_44 = self._params['eta_44']
                self.eta_45 = self._params['eta_45']
                self.eta_46 = self._params['eta_46']
                self.eta_51 = self._params['eta_51']
                self.eta_52 = self._params['eta_52']
                self.eta_53 = self._params['eta_53']
                self.eta_54 = self._params['eta_54']
                self.eta_55 = self._params['eta_55']
                self.eta_56 = self._params['eta_56']
                self.eta_61 = self._params['eta_61']
                self.eta_62 = self._params['eta_62']
                self.eta_63 = self._params['eta_63']
                self.eta_64 = self._params['eta_64']
                self.eta_65 = self._params['eta_65']
                self.eta_66 = self._params['eta_66']
                self.anisotropic = True
            except KeyError:
                self.anisotropic = False


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


Vacuum = Material("Vacuum")
Si_2016_Smith = Material("Si_2016_Smith")
SiO2_2016_Smith = Material("SiO2_2016_Smith")
SiO2_2013_Laude = Material("SiO2_2013_Laude")
As2S3_2016_Smith = Material("As2S3_2016_Smith")
As2S3_2017_Morrison = Material("As2S3_2017_Morrison")
GaAs_2016_Smith = Material("GaAs_2016_Smith")

# wl = 1550 nm
