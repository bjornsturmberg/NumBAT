"""
    objects.py is a subroutine of NumBAT. It contains the Struct
    objects that represent the structure being simulated.

    Copyright (C) 2016  Bjorn Sturmberg, Kokou Dossou, Christian Wolff.

"""

import os
import numpy as np
from fortran import NumBAT

msh_location = '../backend/fortran/msh/'

# # Acknowledgements
# print '\n##################################################################\n'\
#     + 'NumBAT is brought to you by Bjorn Sturmberg, Kokou Dossou, \n' \
#     + 'and Christian Wolff, with support from the ARC \n' \
#     + 'Starting NumBAT calculation ...\n' + \
#       '##################################################################\n'


class Struct(object):
    """ Represents a structured layer.

        Args:
            unitcell_x  (float): The horizontal period of the unit cell \
                in nanometers.

            inc_x  (float): The horizontal diameter of the inclusion in nm.

        Keyword Args:
            unitcell_y  (float): The vertical period of the unit cell \
                in nanometers. If None, unitcell_y = unitcell_x.

            inc_y  (float): The vertical diameter of the inclusion in nm.

            inc_shape  (str): Shape of inclusions that have template mesh, \
                currently: 'circular', 'rectangular'. Rectangular is default.

            diameter2-16  (float): The diameters of further inclusions in nm. \
                Implemented up to diameter6 for 1D_arrays.

            inclusion_a  : A :Material: instance for first inclusion, \
                specified as dispersive refractive index (eg. materials.Si_c) \
                or nondispersive complex number (eg. Material(1.0 + 0.0j)).

            inclusion_b  : A :Material: instance for the second \
                inclusion medium.

            inclusion_c  : A :Material: instance for the third \
                inclusion medium.

            inclusion_d  : A :Material: instance for the fourth \
                inclusion medium.

            inclusion_e  : A :Material: instance for the fifth \
                inclusion medium.

            background  : A :Material: instance for the background medium.

            loss  (bool): If False, Im(n) = 0, if True n as in \
                :Material: instance.

            make_mesh_now  (bool): If True, program creates a FEM mesh with \
                provided :NanoStruct: parameters. If False, must provide \
                mesh_file name of existing .mail that will be run despite \
                :NanoStruct: parameters.

            force_mesh  (bool): If True, a new mesh is created despite \
                existence of mesh with same parameter. This is used to make \
                mesh with equal period etc. but different lc refinement.

            mesh_file  (str): If using a set premade mesh give its name \
                including .mail if 2D_array (eg. 600_60.mail), or .txt if \
                1D_array. It must be located in backend/fortran/msh/

            lc_bkg  (float): Length constant of meshing of background medium \
                (smaller = finer mesh)

            lc2  (float): factor by which lc_bkg will be reduced on inclusion \
                surfaces; lc_surface = cl_bkg / lc2.

            lc3-6'  (float): factor by which lc_bkg will be reduced at center \
                of inclusions.

            plotting_fields  (bool): Unless set to true field data deleted.\
                Also plots modes (ie. FEM solutions) in gmsh format. \
                Plots epsilon*|E|^2 & choice of real/imag/abs of \
                x,y,z components & field vectors. Fields are saved as gmsh \
                files, but can be converted by running the .geo file found in \
                Bloch_fields/PNG/

            plot_real  (bool): Choose to plot real part of modal fields.

            plot_imag  (bool): Choose to plot imaginary part of modal fields.

            plot_abs  (bool): Choose to plot absolute value of modal fields.
    """

    def __init__(self, unitcell_x, inc_x,
                 unitcell_y=None, inc_y=None, inc_shape='rectangular',
                 background=materials.Material(1.0 + 0.0j),
                 inclusion_a=materials.Material(1.0 + 0.0j),
                 inclusion_b=materials.Material(1.0 + 0.0j),
                 inclusion_c=materials.Material(1.0 + 0.0j),
                 inclusion_d=materials.Material(1.0 + 0.0j),
                 inclusion_e=materials.Material(1.0 + 0.0j),
                 inclusion_f=materials.Material(1.0 + 0.0j),
                 loss=True, height_nm=100.0,
                 diameter2=0,  diameter3=0, diameter4=0, diameter5=0,
                 diameter6=0, diameter7=0, diameter8=0, diameter9=0,
                 diameter10=0, diameter11=0, diameter12=0, diameter13=0,
                 diameter14=0, diameter15=0, diameter16=0,
                 make_mesh_now=True, force_mesh=True,
                 mesh_file='NEED_FILE.mail',
                 lc_bkg=0.09, lc2=1.0, lc3=1.0, lc4=1.0, lc5=1.0, lc6=1.0,
                 plotting_fields=False, plot_real=1, plot_imag=0, plot_abs=0,
                 plot_field_conc=False):
        self.unitcell_x = float(unitcell_x)
        self.diameter1 = diameter1
        if unitcell_y is None:
            self.unitcell_y = float(unitcell_x)
        else:
            self.unitcell_y = float(unitcell_y)
        self.inc_shape = inc_shape
        self.background = background
        self.inclusion_a = inclusion_a
        self.inclusion_b = inclusion_b
        self.inclusion_c = inclusion_c
        self.inclusion_d = inclusion_d
        self.inclusion_e = inclusion_e
        self.inclusion_f = inclusion_f
        self.loss = loss
        self.diameter2 = diameter2
        self.diameter3 = diameter3
        self.diameter4 = diameter4
        self.diameter5 = diameter5
        self.diameter6 = diameter6
        self.diameter7 = diameter7
        self.diameter8 = diameter8
        self.diameter9 = diameter9
        self.diameter10 = diameter10
        self.diameter11 = diameter11
        self.diameter12 = diameter12
        self.diameter13 = diameter13
        self.diameter14 = diameter14
        self.diameter15 = diameter15
        self.diameter16 = diameter16
        if diameter3 != 0:
            self.nb_typ_el = 4
        elif diameter2 != 0:
            self.nb_typ_el = 3
        else:
            self.nb_typ_el = 2
        self.lc = lc_bkg
        self.lc2 = lc2
        self.lc3 = lc3
        self.lc4 = lc4
        self.lc5 = lc5
        self.lc6 = lc6
        self.force_mesh = force_mesh
        if make_mesh_now is True:
            self.make_mesh()
        else:
            self.mesh_file = mesh_file
        if plotting_fields is True:
            self.plotting_fields = 1
            if not os.path.exists("Bloch_fields"):
                os.mkdir("Bloch_fields")
            if not os.path.exists("Bloch_fields/PDF"):
                os.mkdir("Bloch_fields/PDF")
        else: self.plotting_fields = 0
        self.plot_real = plot_real
        self.plot_imag = plot_imag
        self.plot_abs = plot_abs
        self.plot_field_conc = plot_field_conc


    def make_mesh(self):
        if self.inc_shape in ['circle', 'ellipse', 'square']:
            if self.diameter10 > 0:
                supercell = 16
                msh_name = '%(d)s_%(dy)s_%(dia)s_%(dias)s_%(diass)s_%(diasss)s_%(diassss)s' % {
               'd': dec_float_str(self.period),
               'dy': dec_float_str(self.period_y),
               'dia': dec_float_str(self.diameter1),
               'dias': dec_float_str(self.diameter2),
               'dias': dec_float_str(self.diameter2),
               'diass': dec_float_str(self.diameter3),
               'diasss': dec_float_str(self.diameter4),
               'diassss': dec_float_str(self.diameter5)}
            elif self.diameter5 > 0:
                supercell = 9
                msh_name = '%(d)s_%(dy)s_%(dia)s_%(dias)s_%(diass)s_%(diasss)s_%(diassss)s' % {
               'd': dec_float_str(self.period),
               'dy': dec_float_str(self.period_y),
               'dia': dec_float_str(self.diameter1),
               'dias': dec_float_str(self.diameter2), 'diass': dec_float_str(self.diameter3),
               'diasss': dec_float_str(self.diameter4), 'diassss': dec_float_str(self.diameter5)}
            elif self.diameter4 > 0:
                supercell = 4
                msh_name = '%(d)s_%(dy)s_%(dia)s_%(dias)s_%(diass)s_%(diasss)s' % {
               'd': dec_float_str(self.period),
               'dy': dec_float_str(self.period_y),
               'dia': dec_float_str(self.diameter1),
               'dias': dec_float_str(self.diameter2), 'diass': dec_float_str(self.diameter3),
               'diasss': dec_float_str(self.diameter4)}
            elif self.diameter3 > 0:
                supercell = 3
                msh_name = '%(d)s_%(dy)s_%(dia)s_%(dias)s_%(diass)s' % {
               'd': dec_float_str(self.period),
               'dy': dec_float_str(self.period_y),
               'dia': dec_float_str(self.diameter1),
               'dias': dec_float_str(self.diameter2), 'diass': dec_float_str(self.diameter3)}
            elif self.diameter2 > 0:
                supercell = 2
                msh_name = '%(d)s_%(dy)s_%(dia)s_%(dias)s' % {'d': dec_float_str(self.period),
               'dy': dec_float_str(self.period_y),
               'dia': dec_float_str(self.diameter1),
               'diameters': dec_float_str(self.diameter2)}
            elif self.diameter1 > 0:
                supercell = 1
                msh_name = '%(d)s_%(dy)s_%(dia)s' % {
                           'd': dec_float_str(self.period),
                           'dy': dec_float_str(self.period_y),
                           'dia': dec_float_str(self.diameter1)}
            else:
                raise ValueError, "must have at least one cylinder of nonzero diameter."

            if self.ellipticity != 0:
                msh_name = msh_name + '_e_%(e)s' % {'e': dec_float_str(self.ellipticity),}
            if self.inc_shape == 'square':
                msh_name = msh_name + '_sq'

            if not os.path.exists(msh_location + msh_name + '.mail') or self.force_mesh is True:
                geo_tmp = open(msh_location + '%s_msh_template.geo' % supercell, "r").read()
                geo = geo_tmp.replace('ff = 0;', "ff = %f;" % self.ff)
                geo = geo.replace('d_in_nm = 0;', "d_in_nm = %f;" % self.period)
                geo = geo.replace('dy_in_nm = 0;', "dy_in_nm = %f;" % self.period_y)
                geo = geo.replace('a1 = 0;', "a1 = %f;" % self.diameter1)
                geo = geo.replace('ellipticity = 0;', "ellipticity = %f;" % self.ellipticity)
                if self.inc_shape == 'square': geo = geo.replace('square = 0;', "square = 1;")
                geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)
                if self.posx != 0:
                    # appropriate for old definition of fraction of distance to touching
                    geo = geo.replace('posx = 0;', "posx = %f;" % (self.posx/self.period*(self.period/(2*np.sqrt(supercell)) - self.diameter1/2.0)))
                    # appropriate for % shift of distance of centre point to (ind) unitcell boundary (ie d/2)
                    # geo = geo.replace('posx = 0;', "posx = %f;" % float(self.posx/supercell))
                if self.posy != 0:
                    geo = geo.replace('posy = 0;', "posy = %f;" % (self.posy/self.period*(self.period/(2*np.sqrt(supercell)) - self.diameter1/2.0)))
                    # geo = geo.replace('posy = 0;', "posy = %f;" % float(self.posy/supercell))
                if supercell > 1:
                    geo = geo.replace('a2 = 0;', "a2 = %f;" % self.diameter2)
                    geo = geo.replace('lc4 = lc/1;', "lc4 = lc/%f;" % self.lc4)
                if supercell > 2:
                    geo = geo.replace('a3 = 0;', "a3 = %f;" % self.diameter3)
                    geo = geo.replace('lc5 = lc/1;', "lc5 = lc/%f;" % self.lc5)
                if supercell > 3:
                    geo = geo.replace('a4 = 0;', "a4 = %f;" % self.diameter4)
                    geo = geo.replace('lc6 = lc/1;', "lc6 = lc/%f;" % self.lc6)
                if supercell > 4:
                    geo = geo.replace('a5 = 0;', "a5 = %f;" % self.diameter5)
                    geo = geo.replace('a6 = 0;', "a6 = %f;" % self.diameter6)
                    geo = geo.replace('a7 = 0;', "a7 = %f;" % self.diameter7)
                    geo = geo.replace('a8 = 0;', "a8 = %f;" % self.diameter8)
                    geo = geo.replace('a9 = 0;', "a9 = %f;" % self.diameter9)
                if supercell > 9:
                    geo = geo.replace('a10 = 0;', "a10 = %f;" % self.diameter10)
                    geo = geo.replace('a11 = 0;', "a11 = %f;" % self.diameter11)
                    geo = geo.replace('a12 = 0;', "a12 = %f;" % self.diameter12)
                    geo = geo.replace('a13 = 0;', "a13 = %f;" % self.diameter13)
                    geo = geo.replace('a14 = 0;', "a14 = %f;" % self.diameter14)
                    geo = geo.replace('a15 = 0;', "a15 = %f;" % self.diameter15)
                    geo = geo.replace('a16 = 0;', "a16 = %f;" % self.diameter16)

        else:
            raise NotImplementedError, "\n Selected inc_shape = '%s' \n \
            is not currently implemented. Please make a mesh with gmsh, & \n \
            consider contributing this to NumBAT via github." % self.inc_shape

        self.mesh_file = msh_name + '.mail'
        open(msh_location + msh_name + '.geo', "w").write(geo)
        NumBAT.conv_gmsh(msh_location+msh_name)

        # # Automatically show created mesh in gmsh.
        # gmsh_cmd = 'gmsh '+ msh_location + msh_name + '.msh'
        # os.system(gmsh_cmd)
        # gmsh_cmd = 'gmsh '+ msh_location + msh_name + '.geo'
        # os.system(gmsh_cmd)


