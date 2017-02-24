"""
    objects.py is a subroutine of NumBAT. It contains the Struct
    objects that represent the structure being simulated.

    Copyright (C) 2016  Bjorn Sturmberg, Kokou Dossou, Christian Wolff.

"""

import os
import numpy as np
import materials
from mode_calcs import Simmo
from fortran import NumBAT

msh_location = '../backend/fortran/msh/'

# # Acknowledgments
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

            inc_a_x  (float): The horizontal diameter of the inclusion in nm.

        Keyword Args:
            unitcell_y  (float): The vertical period of the unit cell \
                in nanometers. If None, unitcell_y = unitcell_x.

            inc_a_y  (float): The vertical diameter of the inclusion in nm.

            inc_shape  (str): Shape of inclusions that have template mesh, \
                currently: 'circular', 'rectangular', 'slot'. \
                Rectangular is default.

            slab_a_x  (float): The horizontal diameter in nm of the slab \
                directly below the inclusion.

            slab_a_y  (float): The vertical diameter in nm of the slab \
                directly below the inclusion.

            slab_b_x  (float): The horizontal diameter in nm of the slab \
                separated from the inclusion by slab_a.

            slab_b_y  (float): The vertical diameter in nm of the slab \
                separated from the inclusion by slab_a.

            two_inc_sep  (float): Separation between edges of inclusions in nm.

            incs_y_offset  (float): Vertical offset between centers of \
                inclusions in nm.

            coat_y  (float): The thickness of any coat layer around \
                the inclusion.

            SPECIFY WHAT IS RADIUS WHAT IS DIAMETER

            background  : A :Material: instance for the background medium.

            inc_a_material  : A :Material: instance for the

            inc_b_material  : A :Material: instance for the

            slab_a_material  : A :Material: instance for the

            slab_a_bkg_material  : A :Material: instance for the

            slab_b_material  : A :Material: instance for the

            slab_b_bkg_material  : A :Material: instance for the

            coat_material  : A :Material: instance for the

            bkg_AC  : A list of acoustic parameters

            inc_a_AC  : A list of acoustic parameters

            slab_a_AC  : A list of acoustic parameters

            slab_a_bkg_AC  : A list of acoustic parameters

            slab_b_AC  : A list of acoustic parameters

            slab_b_bkg_AC  : A list of acoustic parameters

            loss  (bool): If False, Im(n) = 0, if True n as in \
                :Material: instance.

            make_mesh_now  (bool): If True, program creates a FEM mesh with \
                provided :NanoStruct: parameters. If False, must provide \
                mesh_file name of existing .mail that will be run despite \
                :NanoStruct: parameters.

            force_mesh  (bool): If True, a new mesh is created despite \
                existence of mesh with same parameter. This is used to make \
                mesh with equal period etc. but different lc refinement.

            mesh_file  (str): If using a set pre-made mesh give its name \
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

    def __init__(self, unitcell_x, inc_a_x,
                 unitcell_y=None, inc_a_y=None, inc_shape='rectangular',
                 slab_a_x=None, slab_a_y=None, slab_b_x=None, slab_b_y=None,
                 coat_y=None, inc_b_x=None, inc_b_y=None,
                 two_inc_sep=None, incs_y_offset=None,
                 bkg_material=materials.Material(1.0 + 0.0j),
                 inc_a_material=materials.Material(1.0 + 0.0j),
                 slab_a_material=materials.Material(0.0 + 0.0j),
                 slab_a_bkg_material=materials.Material(0.0 + 0.0j),
                 slab_b_material=materials.Material(0.0 + 0.0j),
                 slab_b_bkg_material=materials.Material(0.0 + 0.0j),
                 coat_material=materials.Material(0.0 + 0.0j),
                 inc_c_x=None, inc_d_x=None, inc_e_x=None, inc_f_x=None,
                 inc_g_x=None, inc_h_x=None, inc_i_x=None, inc_j_x=None,
                 inc_k_x=None, inc_l_x=None, inc_m_x=None, inc_n_x=None,
                 inc_o_x=None,
                 inc_b_material=materials.Material(0.0 + 0.0j),
                 inc_c_material=materials.Material(0.0 + 0.0j),
                 inc_d_material=materials.Material(0.0 + 0.0j),
                 inc_e_material=materials.Material(0.0 + 0.0j),
                 inc_f_material=materials.Material(0.0 + 0.0j),
                 inc_g_material=materials.Material(0.0 + 0.0j),
                 inc_h_material=materials.Material(0.0 + 0.0j),
                 inc_i_material=materials.Material(0.0 + 0.0j),
                 inc_j_material=materials.Material(0.0 + 0.0j),
                 inc_k_material=materials.Material(0.0 + 0.0j),
                 inc_l_material=materials.Material(0.0 + 0.0j),
                 inc_m_material=materials.Material(0.0 + 0.0j),
                 inc_n_material=materials.Material(0.0 + 0.0j),
                 inc_o_material=materials.Material(0.0 + 0.0j),
                 inc_b_AC=None,
                 inc_c_AC=None, inc_d_AC=None,inc_e_AC=None,inc_f_AC=None,inc_g_AC=None,
                 inc_h_AC=None,inc_i_AC=None,inc_j_AC=None,inc_k_AC=None,inc_l_AC=None,
                 inc_m_AC=None,inc_n_AC=None,inc_o_AC=None,
                 loss=True, bkg_AC=None, inc_a_AC=None, slab_a_AC=None,
                 slab_a_bkg_AC=None, slab_b_AC=None, slab_b_bkg_AC=None, coat_AC=None,
                 make_mesh_now=True, force_mesh=True,
                 mesh_file='NEED_FILE.mail', check_msh=False,
                 lc_bkg=0.09, lc2=1.0, lc3=1.0, lc4=1.0, lc5=1.0, lc6=1.0,
                 plotting_fields=False, plot_real=1, plot_imag=0, plot_abs=0,
                 plot_field_conc=False):
        self.unitcell_x = float(unitcell_x)
        self.inc_a_x = inc_a_x
        if unitcell_y is None:
            self.unitcell_y = float(unitcell_x)
        else:
            self.unitcell_y = float(unitcell_y)
        if inc_a_y is None:
            self.inc_a_y = float(inc_a_x)
        else:
            self.inc_a_y = float(inc_a_y)
        self.inc_b_x = inc_b_x
        if inc_b_x is not None:
            if inc_b_y is None:
                self.inc_b_y = float(inc_b_x)
            else:
                self.inc_b_y = float(inc_b_y)
        self.inc_shape = inc_shape
        self.slab_a_x = slab_a_x
        self.slab_a_y = slab_a_y
        self.slab_b_x = slab_b_x
        self.slab_b_y = slab_b_y
        self.inc_c_x = inc_c_x
        self.inc_d_x = inc_d_x
        self.inc_e_x = inc_e_x
        self.inc_f_x = inc_f_x
        self.inc_g_x = inc_g_x
        self.inc_h_x = inc_h_x
        self.inc_i_x = inc_i_x
        self.inc_j_x = inc_j_x
        self.inc_k_x = inc_k_x
        self.inc_l_x = inc_l_x
        self.inc_m_x = inc_m_x
        self.inc_n_x = inc_n_x
        self.inc_o_x = inc_o_x
        self.coat_y = coat_y
        self.two_inc_sep = two_inc_sep
        self.bkg_material = bkg_material
        self.inc_a_material = inc_a_material
        self.inc_b_material = inc_b_material
        self.slab_a_material = slab_a_material
        self.slab_a_bkg_material = slab_a_bkg_material
        self.slab_b_material = slab_b_material
        self.slab_b_bkg_material = slab_b_bkg_material
        self.inc_c_material = inc_c_material
        self.inc_d_material = inc_d_material
        self.inc_e_material = inc_e_material
        self.inc_f_material = inc_f_material
        self.inc_g_material = inc_g_material
        self.inc_h_material = inc_h_material
        self.inc_i_material = inc_i_material
        self.inc_j_material = inc_j_material
        self.inc_k_material = inc_k_material
        self.inc_l_material = inc_l_material
        self.inc_m_material = inc_m_material
        self.inc_n_material = inc_n_material
        self.inc_o_material = inc_o_material
        self.inc_b_AC = inc_b_AC
        self.inc_c_AC = inc_c_AC
        self.inc_d_AC = inc_d_AC
        self.inc_e_AC = inc_e_AC
        self.inc_f_AC = inc_f_AC
        self.inc_g_AC = inc_g_AC
        self.inc_h_AC = inc_h_AC
        self.inc_i_AC = inc_i_AC
        self.inc_j_AC = inc_j_AC
        self.inc_k_AC = inc_k_AC
        self.inc_l_AC = inc_l_AC
        self.inc_m_AC = inc_m_AC
        self.inc_n_AC = inc_n_AC
        self.inc_o_AC = inc_o_AC
        self.coat_AC = coat_AC
        self.coat_material = coat_material
        self.loss = loss
        if slab_b_x is not None:
            if coat_y is None and inc_b_x is None:
                self.nb_typ_el = 6
            elif coat_y != None and inc_b_x != None:
                self.nb_typ_el = 8
                raise NotImplementedError, "Have not implemented 2 coated " \
                    "inclusions."
            else:
                self.nb_typ_el = 7
        elif slab_a_x is not None:
            if coat_y is None and inc_b_x is None:
                self.nb_typ_el = 4
            elif coat_y != None and inc_b_x != None:
                self.nb_typ_el = 7
                raise NotImplementedError, "Have not implemented 2 coated " \
                    "inclusions."
            else:
                self.nb_typ_el = 5
        else:
            if coat_y is None and inc_b_x is None:
                self.nb_typ_el = 2
            else:
                self.nb_typ_el = 3
        self.check_msh = check_msh
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
            if not os.path.exists("AC_fields"):
                os.mkdir("AC_fields")
        else: self.plotting_fields = 0
        self.plot_real = plot_real
        self.plot_imag = plot_imag
        self.plot_abs = plot_abs
        self.plot_field_conc = plot_field_conc

        # Order must match msh templates!
        self.acoustic_props_tmp = [bkg_AC, inc_a_AC, slab_a_AC, slab_a_bkg_AC, 
                                   slab_b_AC, slab_b_bkg_AC, inc_b_AC, coat_AC]
        # el_conv_table = {}
        # i = 1; j = 1
        # for matter in acoustic_props:
        #     if matter != None:
        #         el_conv_table[i] = j
        #         j += 1
        #     i += 1
        # self.typ_el_AC = el_conv_table
        # print el_conv_table
        acoustic_props = [x for x in self.acoustic_props_tmp if x is not None]
        self.nb_typ_el_AC = len(acoustic_props)
        # Any material not given acoustic_props assumed to be vacuum.
        rho = np.zeros(self.nb_typ_el_AC)
        c_tensor = np.zeros((6,6,self.nb_typ_el_AC))
        c_tensor_z = np.zeros((3,3,3,self.nb_typ_el_AC))
        p_tensor = np.zeros((3,3,3,3,self.nb_typ_el_AC))
        eta_tensor = np.zeros((3,3,3,3,self.nb_typ_el_AC))
        for k_typ in range(self.nb_typ_el_AC):
            if acoustic_props[k_typ]:
                rho[k_typ] = acoustic_props[k_typ][0]
                c_tensor[0,0,k_typ] = acoustic_props[k_typ][1]
                c_tensor[1,1,k_typ] = acoustic_props[k_typ][1]
                c_tensor[2,2,k_typ] = acoustic_props[k_typ][1]
                c_tensor[0,1,k_typ] = acoustic_props[k_typ][2]
                c_tensor[0,2,k_typ] = acoustic_props[k_typ][2]
                c_tensor[1,0,k_typ] = acoustic_props[k_typ][2]
                c_tensor[1,2,k_typ] = acoustic_props[k_typ][2]
                c_tensor[2,0,k_typ] = acoustic_props[k_typ][2]
                c_tensor[2,1,k_typ] = acoustic_props[k_typ][2]
                c_tensor[3,3,k_typ] = acoustic_props[k_typ][3]
                c_tensor[4,4,k_typ] = acoustic_props[k_typ][3]
                c_tensor[5,5,k_typ] = acoustic_props[k_typ][3]

                c_tensor_z[2,2,2,k_typ] = acoustic_props[k_typ][1]
                c_tensor_z[2,0,0,k_typ] = acoustic_props[k_typ][2]
                c_tensor_z[2,1,1,k_typ] = acoustic_props[k_typ][2]
                c_tensor_z[1,1,2,k_typ] = acoustic_props[k_typ][3]
                c_tensor_z[1,2,1,k_typ] = acoustic_props[k_typ][3]
                c_tensor_z[0,0,2,k_typ] = acoustic_props[k_typ][3]
                c_tensor_z[0,2,0,k_typ] = acoustic_props[k_typ][3]

                p_tensor[0,0,0,0,k_typ] = acoustic_props[k_typ][4]
                p_tensor[1,1,1,1,k_typ] = acoustic_props[k_typ][4]
                p_tensor[2,2,2,2,k_typ] = acoustic_props[k_typ][4]
                p_tensor[0,0,1,1,k_typ] = acoustic_props[k_typ][5]
                p_tensor[0,0,2,2,k_typ] = acoustic_props[k_typ][5]
                p_tensor[1,1,0,0,k_typ] = acoustic_props[k_typ][5]
                p_tensor[1,1,2,2,k_typ] = acoustic_props[k_typ][5]
                p_tensor[2,2,0,0,k_typ] = acoustic_props[k_typ][5]
                p_tensor[2,2,1,1,k_typ] = acoustic_props[k_typ][5]
                p_tensor[1,2,1,2,k_typ] = acoustic_props[k_typ][6]
                p_tensor[1,2,2,1,k_typ] = acoustic_props[k_typ][6]
                p_tensor[2,1,1,2,k_typ] = acoustic_props[k_typ][6]
                p_tensor[2,1,2,1,k_typ] = acoustic_props[k_typ][6]
                p_tensor[0,2,0,2,k_typ] = acoustic_props[k_typ][6]
                p_tensor[0,2,2,0,k_typ] = acoustic_props[k_typ][6]
                p_tensor[2,0,0,2,k_typ] = acoustic_props[k_typ][6]
                p_tensor[2,0,2,0,k_typ] = acoustic_props[k_typ][6]
                p_tensor[0,1,0,1,k_typ] = acoustic_props[k_typ][6]
                p_tensor[0,1,1,0,k_typ] = acoustic_props[k_typ][6]
                p_tensor[1,0,0,1,k_typ] = acoustic_props[k_typ][6]
                p_tensor[1,0,1,0,k_typ] = acoustic_props[k_typ][6]

                eta_tensor[0,0,0,0,k_typ] = acoustic_props[k_typ][7]
                eta_tensor[1,1,1,1,k_typ] = acoustic_props[k_typ][7]
                eta_tensor[2,2,2,2,k_typ] = acoustic_props[k_typ][7]
                eta_tensor[0,0,1,1,k_typ] = acoustic_props[k_typ][8]
                eta_tensor[0,0,2,2,k_typ] = acoustic_props[k_typ][8]
                eta_tensor[1,1,0,0,k_typ] = acoustic_props[k_typ][8]
                eta_tensor[1,1,2,2,k_typ] = acoustic_props[k_typ][8]
                eta_tensor[2,2,0,0,k_typ] = acoustic_props[k_typ][8]
                eta_tensor[2,2,1,1,k_typ] = acoustic_props[k_typ][8]
                eta_tensor[1,2,1,2,k_typ] = acoustic_props[k_typ][9]
                eta_tensor[1,2,2,1,k_typ] = acoustic_props[k_typ][9]
                eta_tensor[2,1,1,2,k_typ] = acoustic_props[k_typ][9]
                eta_tensor[2,1,2,1,k_typ] = acoustic_props[k_typ][9]
                eta_tensor[0,2,0,2,k_typ] = acoustic_props[k_typ][9]
                eta_tensor[0,2,2,0,k_typ] = acoustic_props[k_typ][9]
                eta_tensor[2,0,0,2,k_typ] = acoustic_props[k_typ][9]
                eta_tensor[2,0,2,0,k_typ] = acoustic_props[k_typ][9]
                eta_tensor[0,1,0,1,k_typ] = acoustic_props[k_typ][9]
                eta_tensor[0,1,1,0,k_typ] = acoustic_props[k_typ][9]
                eta_tensor[1,0,0,1,k_typ] = acoustic_props[k_typ][9]
                eta_tensor[1,0,1,0,k_typ] = acoustic_props[k_typ][9]
        self.rho = rho
        self.c_tensor = c_tensor
        self.c_tensor_z = c_tensor_z
        self.p_tensor = p_tensor
        self.eta_tensor = eta_tensor


    def make_mesh(self):
        """ Take the parameters specified in python and make a Gmsh FEM mesh.
            Creates a .geo and .msh file, then uses Fortran conv_gmsh routine
            to convert .msh into .mail, which is used in NumBAT FEM routine.
        """
        if self.inc_shape in ['circular', 'rectangular']:
            if self.slab_b_x is not None:
                if self.coat_y is None and self.inc_b_x is None:
                    msh_template = '1_on_2slabs'
                    msh_name = '1_on_2s%(d)s_%(dy)s_%(dia)s_%(dias)s_%(diass)s_%(diasss)s_%(diassss)s' % {
                   'd': dec_float_str(self.unitcell_x),
                   'dy': dec_float_str(self.unitcell_y),
                   'dia': dec_float_str(self.inc_a_x),
                   'dias': dec_float_str(self.inc_a_y),
                   'dias': dec_float_str(self.slab_a_x),
                   'diass': dec_float_str(self.slab_a_y),
                   'diasss': dec_float_str(self.slab_b_x),
                   'diassss': dec_float_str(self.slab_b_y)}
                elif self.coat_y is None and self.inc_b_x is not None:
                    msh_template = '2_on_2slabs'
                    msh_name = '2_on_2s%(d)s_%(dy)s_%(dia)s_%(dias)s_%(diab)s_%(diasb)s_%(diass)s_%(diasss)s_%(diassss)s' % {
                   'd': dec_float_str(self.unitcell_x),
                   'dy': dec_float_str(self.unitcell_y),
                   'dia': dec_float_str(self.inc_a_x),
                   'dias': dec_float_str(self.inc_a_y),
                   'diab': dec_float_str(self.inc_a_x),
                   'diasb': dec_float_str(self.inc_a_y),
                   'dias': dec_float_str(self.slab_a_x),
                   'diass': dec_float_str(self.slab_a_y),
                   'diasss': dec_float_str(self.slab_b_x),
                   'diassss': dec_float_str(self.slab_b_y)}
                elif self.coat_y is not None and self.inc_b_x is not None:
                    raise NotImplementedError, "Have not implemented 2 coated inclusions."
                elif self.coat_y is not None and self.inc_b_x is None:
                        raise NotImplementedError, "Have not implemented 1 coated inclusions."
                else:
                    raise ValueError, "NumBAT doesn't understand you geometry."
            elif self.slab_a_x is not None:
                if self.coat_y is None and self.inc_b_x is None:
                    msh_template = '1_on_slab'
                    msh_name = '1_on_s%(d)s_%(dy)s_%(dia)s_%(dias)s_%(diass)s' % {
                   'd': dec_float_str(self.unitcell_x),
                   'dy': dec_float_str(self.unitcell_y),
                   'dia': dec_float_str(self.inc_a_x),
                   'dias': dec_float_str(self.inc_a_y),
                   'dias': dec_float_str(self.slab_a_x),
                   'diass': dec_float_str(self.slab_a_y)}
                elif self.coat_y is None and self.inc_b_x is not None:
                    msh_template = '2_on_slab'
                    msh_name = '2_on_s%(d)s_%(dy)s_%(dia)s_%(dias)s_%(diab)s_%(diasb)s_%(diass)s' % {
                   'd': dec_float_str(self.unitcell_x),
                   'dy': dec_float_str(self.unitcell_y),
                   'dia': dec_float_str(self.inc_a_x),
                   'dias': dec_float_str(self.inc_a_y),
                   'diab': dec_float_str(self.inc_a_x),
                   'diasb': dec_float_str(self.inc_a_y),
                   'dias': dec_float_str(self.slab_a_x),
                   'diass': dec_float_str(self.slab_a_y)}
                elif self.coat_y is not None and self.inc_b_x is not None:
                    raise NotImplementedError, "Have not implemented 2 coated inclusions."
                elif self.coat_y is not None and self.inc_b_x is None:
                        raise NotImplementedError, "Have not implemented 1 coated inclusions."
                else:
                    raise ValueError, "NumBAT doesn't understand you geometry."
            elif self.inc_a_x is not None:
                if self.coat_y is None and self.inc_b_x is None:
                    msh_template = '1'
                    msh_name = '1_%(d)s_%(dy)s_%(dia)s' % {
                   'd': dec_float_str(self.unitcell_x),
                   'dy': dec_float_str(self.unitcell_y),
                   'dia': dec_float_str(self.inc_a_x),
                   'dias': dec_float_str(self.inc_a_y)}
                elif self.coat_y is None and self.inc_b_x is not None:
                    msh_template = '2'
                    msh_name = '2_%(d)s_%(dy)s_%(dia)s_%(dias)s_%(diab)s' % {
                   'd': dec_float_str(self.unitcell_x),
                   'dy': dec_float_str(self.unitcell_y),
                   'dia': dec_float_str(self.inc_a_x),
                   'dias': dec_float_str(self.inc_a_y),
                   'diab': dec_float_str(self.inc_a_x),
                   'diasb': dec_float_str(self.inc_a_y)}
                elif self.coat_y is not None and self.inc_b_x is not None:
                    raise NotImplementedError, "Have not implemented 2 coated inclusions."
                elif self.coat_y is not None and self.inc_b_x is None:
                        raise NotImplementedError, "Have not implemented 1 coated inclusions."
                else:
                    raise ValueError, "NumBAT doesn't understand you geometry."
            else:
                raise ValueError, "must have at least one nonzero inclusion."
            if self.inc_shape == 'circular': msh_name+='-c'

            if not os.path.exists(msh_location + msh_name + '.mail') or self.force_mesh is True:
                geo_tmp = open(msh_location + '%s_msh_template.geo' % msh_template, "r").read()
                geo = geo_tmp.replace('d_in_nm = 100;', "d_in_nm = %f;" % self.unitcell_x)
                geo = geo.replace('dy_in_nm = 50;', "dy_in_nm = %f;" % self.unitcell_y)
                geo = geo.replace('a1 = 20;', "a1 = %f;" % self.inc_a_x)
                geo = geo.replace('a1y = 10;', "a1y = %f;" % self.inc_a_y)
                if self.inc_shape == 'circular': geo = geo.replace('rect = 1;', "rect = 0;")
                geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)
                if msh_template == '2' or msh_template == '2_on_s' or msh_template == '2_on_2s':
                    geo = geo.replace('a2 = 10;', "a2 = %f;" % self.inc_b_x)
                    geo = geo.replace('a2y = 20;', "a2y = %f;" % self.inc_b_y)
                    geo = geo.replace('sep = 10;', "sep = %f;" % self.two_inc_sep)
                    geo = geo.replace('lc4 = lc/1;', "lc4 = lc/%f;" % self.lc4)
                if msh_template == '2':
                    geo = geo.replace('yoff = -5;;', "yoff = %f;" % self.incs_y_offset)
                if msh_template == '1_on_slab' or msh_template == '1_on_2slabs' or msh_template == '1_on_slab' or msh_template == '2_on_2slabs':
                    geo = geo.replace('slab_width = d_in_nm;', "slab_width = %f;" % self.slab_a_x)
                    geo = geo.replace('slab_height = 10;', "slab_height = %f;" % self.slab_a_y)
                    geo = geo.replace('lc4 = lc/1;', "lc4 = lc/%f;" % self.lc4)
                    geo = geo.replace('lc5 = lc/1;', "lc5 = lc/%f;" % self.lc5)
                if msh_template == '1_on_2slabs' or msh_template == '2_on_2slabs':
                    geo = geo.replace('slab2_width = d_in_nm;', "slab2_width = %f;" % self.slab_b_x)
                    geo = geo.replace('slab2_height = 5;', "slab2_height = %f;" % self.slab_b_y)
                    geo = geo.replace('lc4 = lc/1;', "lc4 = lc/%f;" % self.lc4)
                    geo = geo.replace('lc5 = lc/1;', "lc5 = lc/%f;" % self.lc5)

        elif self.inc_shape in ['slot']:
            if self.coat_AC:
                msh_name = 'slot_c_%(d)s_%(dy)s_%(a)s_%(b)s_%(c)s_%(cc)s_%(ccc)s' % {
                'd': dec_float_str(self.unitcell_x),
                'dy': dec_float_str(self.unitcell_y),
                'a': dec_float_str(self.inc_a_x),
                'b': dec_float_str(self.inc_a_y),
                'c': dec_float_str(self.inc_b_x),
                'cc': dec_float_str(self.slab_a_y),
                'ccc': dec_float_str(self.coat_y)}
                msh_template = 'slot_coated'
                if not os.path.exists(msh_location + msh_name + '.mail') or self.force_mesh is True:
                    geo_tmp = open(msh_location + '%s_msh_template.geo' % msh_template, "r").read()
                    geo = geo_tmp.replace('d_in_nm = 100;', "d_in_nm = %f;" % self.unitcell_x)
                    geo = geo.replace('dy_in_nm = 50;', "dy_in_nm = %f;" % self.unitcell_y)
                    geo = geo.replace('a1 = 20;', "a1 = %f;" % self.inc_a_x)
                    geo = geo.replace('a1y = 10;', "a1y = %f;" % self.inc_a_y)
                    geo = geo.replace('a2 = 20;', "a2 = %f;" % self.inc_b_x)
                    geo = geo.replace('s1y = 10;', "s1y = %f;" % self.slab_a_y)
                    geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                    geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                    geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)
                    geo = geo.replace('c1y = 10;', "c1y = %f;" % self.coat_y)
                self.nb_typ_el = 5
            else:
                msh_name = 'slot_%(d)s_%(dy)s_%(a)s_%(b)s_%(c)s_%(cc)s_%(ccc)s' % {
                'd': dec_float_str(self.unitcell_x),
                'dy': dec_float_str(self.unitcell_y),
                'a': dec_float_str(self.inc_a_x),
                'b': dec_float_str(self.inc_a_y),
                'c': dec_float_str(self.inc_b_x),
                'cc': dec_float_str(self.slab_a_y)}
                msh_template = 'slot'
                if not os.path.exists(msh_location + msh_name + '.mail') or self.force_mesh is True:
                    geo_tmp = open(msh_location + '%s_msh_template.geo' % msh_template, "r").read()
                    geo = geo_tmp.replace('d_in_nm = 100;', "d_in_nm = %f;" % self.unitcell_x)
                    geo = geo.replace('dy_in_nm = 50;', "dy_in_nm = %f;" % self.unitcell_y)
                    geo = geo.replace('a1 = 20;', "a1 = %f;" % self.inc_a_x)
                    geo = geo.replace('a1y = 10;', "a1y = %f;" % self.inc_a_y)
                    geo = geo.replace('a2 = 20;', "a2 = %f;" % self.inc_b_x)
                    geo = geo.replace('s1y = 10;', "s1y = %f;" % self.slab_a_y)
                    geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                    geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                    geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)
                self.nb_typ_el = 4




        elif self.inc_shape in ['SMF']:
            msh_name = 'SMF_%(d)s_%(dy)s_%(a)s_%(b)s_%(c)s_%(d)s_%(e)s_%(f)s_%(g)s' % {
            'd': dec_float_str(self.unitcell_x),
            'dy': dec_float_str(self.unitcell_y),
            'a': dec_float_str(self.inc_a_x),
            'b': dec_float_str(self.inc_b_x),
            'c': dec_float_str(self.inc_c_x),
            'd': dec_float_str(self.inc_d_x),
            'e': dec_float_str(self.inc_e_x),
            'f': dec_float_str(self.inc_f_x),
            'g': dec_float_str(self.inc_g_x)}
            msh_name = msh_name + '_%(h)s_%(i)s_%(j)s_%(k)s_%(l)s_%(m)s_%(n)s_%(o)s' % {
            'h': dec_float_str(self.inc_h_x),
            'i': dec_float_str(self.inc_i_x),
            'j': dec_float_str(self.inc_j_x),
            'k': dec_float_str(self.inc_k_x),
            'l': dec_float_str(self.inc_l_x),
            'm': dec_float_str(self.inc_m_x),
            'n': dec_float_str(self.inc_n_x),
            'o': dec_float_str(self.inc_o_x)}
            msh_template = 'SMF'
            if not os.path.exists(msh_location + msh_name + '.mail') or self.force_mesh is True:
                geo_tmp = open(msh_location + '%s_msh_template.geo' % msh_template, "r").read()
                geo = geo_tmp.replace('d_in_nm = 100;', "d_in_nm = %f;" % self.unitcell_x)
                geo = geo.replace('dy_in_nm = 50;', "dy_in_nm = %f;" % self.unitcell_y)
                geo = geo.replace('a1 = 2;', "a1 = %f;" % self.inc_a_x)
                geo = geo.replace('a2 = 2;', "a2 = %f;" % self.inc_b_x)
                geo = geo.replace('a3 = 2;', "a3 = %f;" % self.inc_c_x)
                geo = geo.replace('a4 = 2;', "a4 = %f;" % self.inc_d_x)
                geo = geo.replace('a5 = 2;', "a5 = %f;" % self.inc_e_x)
                geo = geo.replace('a6 = 2;', "a6 = %f;" % self.inc_f_x)
                geo = geo.replace('a7 = 2;', "a7 = %f;" % self.inc_g_x)
                geo = geo.replace('a8 = 2;', "a8 = %f;" % self.inc_h_x)
                geo = geo.replace('a9 = 2;', "a9 = %f;" % self.inc_i_x)
                geo = geo.replace('a10 = 2;', "a10 = %f;" % self.inc_j_x)
                geo = geo.replace('a11 = 2;', "a11 = %f;" % self.inc_k_x)
                geo = geo.replace('a12 = 2;', "a12 = %f;" % self.inc_l_x)
                geo = geo.replace('a13 = 2;', "a13 = %f;" % self.inc_m_x)
                geo = geo.replace('a14 = 2;', "a14 = %f;" % self.inc_n_x)
                geo = geo.replace('a15 = 2;', "a15 = %f;" % self.inc_o_x)
                geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)

        else:
            raise NotImplementedError, "\n Selected inc_shape = '%s' \n " \
            "is not currently implemented. Please make a mesh with gmsh, & \n " \
            "consider contributing this to NumBAT via gitlab." % self.inc_shape

        self.mesh_file = msh_name + '.mail'
        open(msh_location + msh_name + '.geo', "w").write(geo)
        NumBAT.conv_gmsh(msh_location+msh_name)

        if self.check_msh is True:
            # Automatically show created mesh in gmsh.
            gmsh_cmd = 'gmsh '+ msh_location + msh_name + '.msh'
            os.system(gmsh_cmd)
            gmsh_cmd = 'gmsh '+ msh_location + msh_name + '.geo'
            os.system(gmsh_cmd)


    def calc_EM_modes(self, wl_nm, num_modes, n_eff, **args):
        """ Run a simulation to find the Struct's EM modes.

            Args:
                wl_nm  (float): Wavelength of EM wave in vacuum.

                num_modes  (int): Number of EM modes to solve for.

                n_eff  (float): Guesstimated effective index of
                    fundamental mode, will be origin of FEM search.

            Returns:
                :Simmo: object
        """
        simmo = Simmo(self, wl_nm, num_modes=num_modes, n_eff=n_eff)

        simmo.calc_EM_modes(**args)
        return simmo


    def calc_AC_modes(self, wl_nm, num_modes, k_AC,
                      shift_Hz=None, EM_sim=None, **args):
        """ Run a simulation to find the Struct's acoustic modes.

            Args:
                wl_nm  (float): Wavelength of AC wave in vacuum.

                num_modes  (int): Number of AC modes to solve for.

            Keyword Args:
                shift_Hz  (float): Guesstimated frequency of modes,
                    will be origin of FEM search. NumBAT will make
                    an educated guess if shift_Hz=None.
                    (Technically the shift and invert parameter).

                EM_sim  (:Simmo: object): Typically an acoustic
                    simulation follows on from an optical one.
                    Supply the EM :Simmo: object so the AC FEM mesh
                    can be constructed from this.
                    This is done by removing vacuum regions.

            Returns:
                :Simmo: object
        """
        simmo_AC = Simmo(self, wl_nm, num_modes=num_modes,
                         k_AC=k_AC, shift_Hz=shift_Hz,
                         EM_sim=EM_sim)

        simmo_AC.calc_AC_modes(**args)
        return simmo_AC


def dec_float_str(dec_float):
    """ Convert float with decimal point into string with '_' in place of '.' """
    # string = str(dec_float)
    string = '%8.4f' %dec_float
    string = string.replace('.', '_')
    string = string.replace(" ", "")
    return string
