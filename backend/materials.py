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
    def __init__(self,data_file):

        self.load_data_file(data_file)


    def load_data_file(self, data_file, alt_path=''):  
        """
        Load data from json file.
        
        Args:
            data_file  (str): name of data file located in NumBAT/backend/material_data
            
            alt_path  (str): non standard path to data_file
        
        """
