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

            SPECIFY WHAT IS RADIUS WHAT IS DIAMETER

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

            symmetry_flag  (bool): True if materials all have sufficient \
                symmetry that their tensors contain only 3 unique values.
                If False must specify full [3,3,3,3] tensors.

            material_a  : A :Material: instance - check backend/msh_type_lib

            material_b  : A :Material: instance - check backend/msh_type_lib

            material_c-r  : A :Material: instance - check backend/msh_type_lib

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
                 coat_x=None, coat_y=None, inc_b_x=None, inc_b_y=None,
                 two_inc_sep=None, incs_y_offset=None,
                 inc_c_x=None, inc_d_x=None, inc_e_x=None, inc_f_x=None,
                 inc_g_x=None, inc_h_x=None, inc_i_x=None, inc_j_x=None,
                 inc_k_x=None, inc_l_x=None, inc_m_x=None, inc_n_x=None,
                 inc_o_x=None,
                 material_a=materials.Air, material_b=materials.Air, 
                 material_c=materials.Air, material_d=materials.Air, 
                 material_e=materials.Air, material_f=materials.Air, 
                 material_g=materials.Air, material_h=materials.Air, 
                 material_i=materials.Air, material_j=materials.Air, 
                 material_k=materials.Air, material_l=materials.Air, 
                 material_m=materials.Air, material_n=materials.Air, 
                 material_o=materials.Air, material_p=materials.Air, 
                 material_q=materials.Air, material_r=materials.Air, 
                 loss=True, symmetry_flag=True,
                 make_mesh_now=True, force_mesh=True,
                 mesh_file='NEED_FILE.mail', check_msh=False,
                 lc_bkg=0.09, lc2=1.0, lc3=1.0, lc4=1.0, lc5=1.0, lc6=1.0,
                 plotting_fields=False, plot_real=1, plot_imag=0, plot_abs=0,
                 plot_field_conc=False):
        # Structures geometric shapes
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
        self.coat_x = coat_x
        self.coat_y = coat_y
        self.two_inc_sep = two_inc_sep
        # Structures material properties - need to check geometry definition 
        # to ensure connecting material type with correct surface of geometry
        self.material_a = material_a
        self.material_b = material_b
        self.material_c = material_c
        self.material_d = material_d
        self.material_e = material_e
        self.material_f = material_f
        self.material_g = material_g
        self.material_h = material_h
        self.material_i = material_i
        self.material_j = material_j
        self.material_k = material_k
        self.material_l = material_l
        self.material_m = material_m
        self.material_n = material_n
        self.material_o = material_o
        self.material_p = material_p
        self.material_q = material_q
        self.material_r = material_r

        self.loss = loss
        if slab_b_x is not None:
            if coat_y is None and inc_b_x is None:
                self.nb_typ_el = 6
            elif coat_y != None and inc_b_x != None:
                self.nb_typ_el = 8
                raise NotImplementedError("Have not implemented 2 coated " \
                    "inclusions.")
            else:
                self.nb_typ_el = 7
        elif slab_a_x is not None:
            if coat_y is None and inc_b_x is None:
                self.nb_typ_el = 4
            elif coat_y != None and inc_b_x != None:
                self.nb_typ_el = 7
                raise NotImplementedError("Have not implemented 2 coated " \
                    "inclusions.")
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

        if symmetry_flag is True:
            self.symmetry_flag = 1
        else:
            self.symmetry_flag = 0

        # Order must match msh templates!
        # self.acoustic_props_tmp = [bkg_material, inc_a_material, 
        #                            slab_a_material, slab_a_bkg_material, 
        #                            slab_b_material, slab_b_bkg_material, 
        #                            inc_b_material, coat_material]

        self.acoustic_props_tmp = [material_a, material_b, material_c,
                                   material_d, material_e, material_f, 
                                   material_g, material_h, material_i, 
                                   material_j, material_k, material_l, 
                                   material_m, material_n, material_o, 
                                   material_p, material_q, material_r]
        # el_conv_table = {}
        # i = 1; j = 1
        # for matter in acoustic_props:
        #     if matter != None:
        #         el_conv_table[i] = j
        #         j += 1
        #     i += 1
        # self.typ_el_AC = el_conv_table
        # print el_conv_table
        acoustic_props = [x for x in self.acoustic_props_tmp if x.s is not None]
        self.nb_typ_el_AC = len(acoustic_props)
        # Any material not given acoustic_props assumed to be vacuum.
        rho = np.zeros(self.nb_typ_el_AC)
        c_tensor = np.zeros((6,6,self.nb_typ_el_AC))
        c_tensor_z = np.zeros((3,3,3,self.nb_typ_el_AC))
        p_tensor = np.zeros((3,3,3,3,self.nb_typ_el_AC))
        eta_tensor = np.zeros((3,3,3,3,self.nb_typ_el_AC))
        for k_typ in range(self.nb_typ_el_AC):
            if acoustic_props[k_typ]:
                rho[k_typ] = acoustic_props[k_typ].s

                if symmetry_flag is True:
                    c_tensor[0,0,k_typ] = acoustic_props[k_typ].c_11
                    c_tensor[1,1,k_typ] = acoustic_props[k_typ].c_11
                    c_tensor[2,2,k_typ] = acoustic_props[k_typ].c_11
                    c_tensor[0,1,k_typ] = acoustic_props[k_typ].c_12
                    c_tensor[0,2,k_typ] = acoustic_props[k_typ].c_12
                    c_tensor[1,0,k_typ] = acoustic_props[k_typ].c_12
                    c_tensor[1,2,k_typ] = acoustic_props[k_typ].c_12
                    c_tensor[2,0,k_typ] = acoustic_props[k_typ].c_12
                    c_tensor[2,1,k_typ] = acoustic_props[k_typ].c_12
                    c_tensor[3,3,k_typ] = acoustic_props[k_typ].c_44
                    c_tensor[4,4,k_typ] = acoustic_props[k_typ].c_44
                    c_tensor[5,5,k_typ] = acoustic_props[k_typ].c_44

                    c_tensor_z[2,2,2,k_typ] = acoustic_props[k_typ].c_11
                    c_tensor_z[2,0,0,k_typ] = acoustic_props[k_typ].c_12
                    c_tensor_z[2,1,1,k_typ] = acoustic_props[k_typ].c_12
                    c_tensor_z[1,1,2,k_typ] = acoustic_props[k_typ].c_44
                    c_tensor_z[1,2,1,k_typ] = acoustic_props[k_typ].c_44
                    c_tensor_z[0,0,2,k_typ] = acoustic_props[k_typ].c_44
                    c_tensor_z[0,2,0,k_typ] = acoustic_props[k_typ].c_44

                    p_tensor[0,0,0,0,k_typ] = acoustic_props[k_typ].p_11
                    p_tensor[1,1,1,1,k_typ] = acoustic_props[k_typ].p_11
                    p_tensor[2,2,2,2,k_typ] = acoustic_props[k_typ].p_11
                    p_tensor[0,0,1,1,k_typ] = acoustic_props[k_typ].p_12
                    p_tensor[0,0,2,2,k_typ] = acoustic_props[k_typ].p_12
                    p_tensor[1,1,0,0,k_typ] = acoustic_props[k_typ].p_12
                    p_tensor[1,1,2,2,k_typ] = acoustic_props[k_typ].p_12
                    p_tensor[2,2,0,0,k_typ] = acoustic_props[k_typ].p_12
                    p_tensor[2,2,1,1,k_typ] = acoustic_props[k_typ].p_12
                    p_tensor[1,2,1,2,k_typ] = acoustic_props[k_typ].p_44
                    p_tensor[1,2,2,1,k_typ] = acoustic_props[k_typ].p_44
                    p_tensor[2,1,1,2,k_typ] = acoustic_props[k_typ].p_44
                    p_tensor[2,1,2,1,k_typ] = acoustic_props[k_typ].p_44
                    p_tensor[0,2,0,2,k_typ] = acoustic_props[k_typ].p_44
                    p_tensor[0,2,2,0,k_typ] = acoustic_props[k_typ].p_44
                    p_tensor[2,0,0,2,k_typ] = acoustic_props[k_typ].p_44
                    p_tensor[2,0,2,0,k_typ] = acoustic_props[k_typ].p_44
                    p_tensor[0,1,0,1,k_typ] = acoustic_props[k_typ].p_44
                    p_tensor[0,1,1,0,k_typ] = acoustic_props[k_typ].p_44
                    p_tensor[1,0,0,1,k_typ] = acoustic_props[k_typ].p_44
                    p_tensor[1,0,1,0,k_typ] = acoustic_props[k_typ].p_44

                    eta_tensor[0,0,0,0,k_typ] = acoustic_props[k_typ].eta_11
                    eta_tensor[1,1,1,1,k_typ] = acoustic_props[k_typ].eta_11
                    eta_tensor[2,2,2,2,k_typ] = acoustic_props[k_typ].eta_11
                    eta_tensor[0,0,1,1,k_typ] = acoustic_props[k_typ].eta_12
                    eta_tensor[0,0,2,2,k_typ] = acoustic_props[k_typ].eta_12
                    eta_tensor[1,1,0,0,k_typ] = acoustic_props[k_typ].eta_12
                    eta_tensor[1,1,2,2,k_typ] = acoustic_props[k_typ].eta_12
                    eta_tensor[2,2,0,0,k_typ] = acoustic_props[k_typ].eta_12
                    eta_tensor[2,2,1,1,k_typ] = acoustic_props[k_typ].eta_12
                    eta_tensor[1,2,1,2,k_typ] = acoustic_props[k_typ].eta_44
                    eta_tensor[1,2,2,1,k_typ] = acoustic_props[k_typ].eta_44
                    eta_tensor[2,1,1,2,k_typ] = acoustic_props[k_typ].eta_44
                    eta_tensor[2,1,2,1,k_typ] = acoustic_props[k_typ].eta_44
                    eta_tensor[0,2,0,2,k_typ] = acoustic_props[k_typ].eta_44
                    eta_tensor[0,2,2,0,k_typ] = acoustic_props[k_typ].eta_44
                    eta_tensor[2,0,0,2,k_typ] = acoustic_props[k_typ].eta_44
                    eta_tensor[2,0,2,0,k_typ] = acoustic_props[k_typ].eta_44
                    eta_tensor[0,1,0,1,k_typ] = acoustic_props[k_typ].eta_44
                    eta_tensor[0,1,1,0,k_typ] = acoustic_props[k_typ].eta_44
                    eta_tensor[1,0,0,1,k_typ] = acoustic_props[k_typ].eta_44
                    eta_tensor[1,0,1,0,k_typ] = acoustic_props[k_typ].eta_44

                elif symmetry_flag is False:
                    c_tensor[0,0,k_typ] = acoustic_props[k_typ].c_11
                    c_tensor[0,1,k_typ] = acoustic_props[k_typ].c_12
                    c_tensor[0,2,k_typ] = acoustic_props[k_typ].c_13
                    c_tensor[0,3,k_typ] = acoustic_props[k_typ].c_14
                    c_tensor[0,4,k_typ] = acoustic_props[k_typ].c_15
                    c_tensor[0,5,k_typ] = acoustic_props[k_typ].c_16
                    c_tensor[1,0,k_typ] = acoustic_props[k_typ].c_21
                    c_tensor[1,1,k_typ] = acoustic_props[k_typ].c_22
                    c_tensor[1,2,k_typ] = acoustic_props[k_typ].c_23
                    c_tensor[1,3,k_typ] = acoustic_props[k_typ].c_24
                    c_tensor[1,4,k_typ] = acoustic_props[k_typ].c_25
                    c_tensor[1,5,k_typ] = acoustic_props[k_typ].c_26
                    c_tensor[2,0,k_typ] = acoustic_props[k_typ].c_31
                    c_tensor[2,1,k_typ] = acoustic_props[k_typ].c_32
                    c_tensor[2,2,k_typ] = acoustic_props[k_typ].c_33
                    c_tensor[2,3,k_typ] = acoustic_props[k_typ].c_34
                    c_tensor[2,4,k_typ] = acoustic_props[k_typ].c_35
                    c_tensor[2,5,k_typ] = acoustic_props[k_typ].c_36
                    c_tensor[3,0,k_typ] = acoustic_props[k_typ].c_41
                    c_tensor[3,1,k_typ] = acoustic_props[k_typ].c_42
                    c_tensor[3,2,k_typ] = acoustic_props[k_typ].c_43
                    c_tensor[3,3,k_typ] = acoustic_props[k_typ].c_44
                    c_tensor[3,4,k_typ] = acoustic_props[k_typ].c_45
                    c_tensor[3,5,k_typ] = acoustic_props[k_typ].c_46
                    c_tensor[4,0,k_typ] = acoustic_props[k_typ].c_51
                    c_tensor[4,1,k_typ] = acoustic_props[k_typ].c_52
                    c_tensor[4,2,k_typ] = acoustic_props[k_typ].c_53
                    c_tensor[4,3,k_typ] = acoustic_props[k_typ].c_54
                    c_tensor[4,4,k_typ] = acoustic_props[k_typ].c_55
                    c_tensor[4,5,k_typ] = acoustic_props[k_typ].c_56
                    c_tensor[5,0,k_typ] = acoustic_props[k_typ].c_61
                    c_tensor[5,1,k_typ] = acoustic_props[k_typ].c_62
                    c_tensor[5,2,k_typ] = acoustic_props[k_typ].c_63
                    c_tensor[5,3,k_typ] = acoustic_props[k_typ].c_64
                    c_tensor[5,4,k_typ] = acoustic_props[k_typ].c_65
                    c_tensor[5,5,k_typ] = acoustic_props[k_typ].c_66

                    c_tensor_z[0,0,0,k_typ] = acoustic_props[k_typ].c_15
                    c_tensor_z[0,0,1,k_typ] = acoustic_props[k_typ].c_14
                    c_tensor_z[0,0,2,k_typ] = acoustic_props[k_typ].c_13
                    c_tensor_z[0,1,0,k_typ] = acoustic_props[k_typ].c_65
                    c_tensor_z[0,1,1,k_typ] = acoustic_props[k_typ].c_64
                    c_tensor_z[0,1,2,k_typ] = acoustic_props[k_typ].c_63
                    c_tensor_z[0,2,0,k_typ] = acoustic_props[k_typ].c_55
                    c_tensor_z[0,2,1,k_typ] = acoustic_props[k_typ].c_54
                    c_tensor_z[0,2,2,k_typ] = acoustic_props[k_typ].c_53
                    c_tensor_z[1,0,0,k_typ] = acoustic_props[k_typ].c_65
                    c_tensor_z[1,0,1,k_typ] = acoustic_props[k_typ].c_64
                    c_tensor_z[1,0,2,k_typ] = acoustic_props[k_typ].c_63
                    c_tensor_z[1,1,0,k_typ] = acoustic_props[k_typ].c_25
                    c_tensor_z[1,1,1,k_typ] = acoustic_props[k_typ].c_24
                    c_tensor_z[1,1,2,k_typ] = acoustic_props[k_typ].c_23
                    c_tensor_z[1,2,0,k_typ] = acoustic_props[k_typ].c_45
                    c_tensor_z[1,2,1,k_typ] = acoustic_props[k_typ].c_44
                    c_tensor_z[1,2,2,k_typ] = acoustic_props[k_typ].c_43
                    c_tensor_z[2,0,0,k_typ] = acoustic_props[k_typ].c_55
                    c_tensor_z[2,0,1,k_typ] = acoustic_props[k_typ].c_54
                    c_tensor_z[2,0,2,k_typ] = acoustic_props[k_typ].c_53
                    c_tensor_z[2,1,0,k_typ] = acoustic_props[k_typ].c_45
                    c_tensor_z[2,1,1,k_typ] = acoustic_props[k_typ].c_44
                    c_tensor_z[2,1,2,k_typ] = acoustic_props[k_typ].c_43
                    c_tensor_z[2,2,0,k_typ] = acoustic_props[k_typ].c_35
                    c_tensor_z[2,2,1,k_typ] = acoustic_props[k_typ].c_34
                    c_tensor_z[2,2,2,k_typ] = acoustic_props[k_typ].c_33

                    p_tensor[0,0,0,0,k_typ] = acoustic_props[k_typ].p_11
                    p_tensor[0,0,0,1,k_typ] = acoustic_props[k_typ].p_16
                    p_tensor[0,0,0,2,k_typ] = acoustic_props[k_typ].p_15
                    p_tensor[0,0,1,0,k_typ] = acoustic_props[k_typ].p_16
                    p_tensor[0,0,1,1,k_typ] = acoustic_props[k_typ].p_12
                    p_tensor[0,0,1,2,k_typ] = acoustic_props[k_typ].p_14
                    p_tensor[0,0,2,0,k_typ] = acoustic_props[k_typ].p_15
                    p_tensor[0,0,2,1,k_typ] = acoustic_props[k_typ].p_14
                    p_tensor[0,0,2,2,k_typ] = acoustic_props[k_typ].p_13
                    p_tensor[0,1,0,0,k_typ] = acoustic_props[k_typ].p_61
                    p_tensor[0,1,0,1,k_typ] = acoustic_props[k_typ].p_66
                    p_tensor[0,1,0,2,k_typ] = acoustic_props[k_typ].p_65
                    p_tensor[0,1,1,0,k_typ] = acoustic_props[k_typ].p_66
                    p_tensor[0,1,1,1,k_typ] = acoustic_props[k_typ].p_62
                    p_tensor[0,1,1,2,k_typ] = acoustic_props[k_typ].p_64
                    p_tensor[0,1,2,0,k_typ] = acoustic_props[k_typ].p_65
                    p_tensor[0,1,2,1,k_typ] = acoustic_props[k_typ].p_64
                    p_tensor[0,1,2,2,k_typ] = acoustic_props[k_typ].p_63
                    p_tensor[0,2,0,0,k_typ] = acoustic_props[k_typ].p_51
                    p_tensor[0,2,0,1,k_typ] = acoustic_props[k_typ].p_56
                    p_tensor[0,2,0,2,k_typ] = acoustic_props[k_typ].p_55
                    p_tensor[0,2,1,0,k_typ] = acoustic_props[k_typ].p_56
                    p_tensor[0,2,1,1,k_typ] = acoustic_props[k_typ].p_52
                    p_tensor[0,2,1,2,k_typ] = acoustic_props[k_typ].p_54
                    p_tensor[0,2,2,0,k_typ] = acoustic_props[k_typ].p_55
                    p_tensor[0,2,2,1,k_typ] = acoustic_props[k_typ].p_54
                    p_tensor[0,2,2,2,k_typ] = acoustic_props[k_typ].p_53
                    p_tensor[1,0,0,0,k_typ] = acoustic_props[k_typ].p_61
                    p_tensor[1,0,0,1,k_typ] = acoustic_props[k_typ].p_66
                    p_tensor[1,0,0,2,k_typ] = acoustic_props[k_typ].p_65
                    p_tensor[1,0,1,0,k_typ] = acoustic_props[k_typ].p_66
                    p_tensor[1,0,1,1,k_typ] = acoustic_props[k_typ].p_62
                    p_tensor[1,0,1,2,k_typ] = acoustic_props[k_typ].p_64
                    p_tensor[1,0,2,0,k_typ] = acoustic_props[k_typ].p_65
                    p_tensor[1,0,2,1,k_typ] = acoustic_props[k_typ].p_64
                    p_tensor[1,0,2,2,k_typ] = acoustic_props[k_typ].p_63
                    p_tensor[1,1,0,0,k_typ] = acoustic_props[k_typ].p_21
                    p_tensor[1,1,0,1,k_typ] = acoustic_props[k_typ].p_26
                    p_tensor[1,1,0,2,k_typ] = acoustic_props[k_typ].p_25
                    p_tensor[1,1,1,0,k_typ] = acoustic_props[k_typ].p_26
                    p_tensor[1,1,1,1,k_typ] = acoustic_props[k_typ].p_22
                    p_tensor[1,1,1,2,k_typ] = acoustic_props[k_typ].p_24
                    p_tensor[1,1,2,0,k_typ] = acoustic_props[k_typ].p_25
                    p_tensor[1,1,2,1,k_typ] = acoustic_props[k_typ].p_24
                    p_tensor[1,1,2,2,k_typ] = acoustic_props[k_typ].p_23
                    p_tensor[1,2,0,0,k_typ] = acoustic_props[k_typ].p_41
                    p_tensor[1,2,0,1,k_typ] = acoustic_props[k_typ].p_46
                    p_tensor[1,2,0,2,k_typ] = acoustic_props[k_typ].p_45
                    p_tensor[1,2,1,0,k_typ] = acoustic_props[k_typ].p_46
                    p_tensor[1,2,1,1,k_typ] = acoustic_props[k_typ].p_42
                    p_tensor[1,2,1,2,k_typ] = acoustic_props[k_typ].p_44
                    p_tensor[1,2,2,0,k_typ] = acoustic_props[k_typ].p_45
                    p_tensor[1,2,2,1,k_typ] = acoustic_props[k_typ].p_44
                    p_tensor[1,2,2,2,k_typ] = acoustic_props[k_typ].p_43
                    p_tensor[2,0,0,0,k_typ] = acoustic_props[k_typ].p_51
                    p_tensor[2,0,0,1,k_typ] = acoustic_props[k_typ].p_56
                    p_tensor[2,0,0,2,k_typ] = acoustic_props[k_typ].p_55
                    p_tensor[2,0,1,0,k_typ] = acoustic_props[k_typ].p_56
                    p_tensor[2,0,1,1,k_typ] = acoustic_props[k_typ].p_52
                    p_tensor[2,0,1,2,k_typ] = acoustic_props[k_typ].p_54
                    p_tensor[2,0,2,0,k_typ] = acoustic_props[k_typ].p_55
                    p_tensor[2,0,2,1,k_typ] = acoustic_props[k_typ].p_54
                    p_tensor[2,0,2,2,k_typ] = acoustic_props[k_typ].p_53
                    p_tensor[2,1,0,0,k_typ] = acoustic_props[k_typ].p_41
                    p_tensor[2,1,0,1,k_typ] = acoustic_props[k_typ].p_46
                    p_tensor[2,1,0,2,k_typ] = acoustic_props[k_typ].p_45
                    p_tensor[2,1,1,0,k_typ] = acoustic_props[k_typ].p_46
                    p_tensor[2,1,1,1,k_typ] = acoustic_props[k_typ].p_42
                    p_tensor[2,1,1,2,k_typ] = acoustic_props[k_typ].p_44
                    p_tensor[2,1,2,0,k_typ] = acoustic_props[k_typ].p_45
                    p_tensor[2,1,2,1,k_typ] = acoustic_props[k_typ].p_44
                    p_tensor[2,1,2,2,k_typ] = acoustic_props[k_typ].p_43
                    p_tensor[2,2,0,0,k_typ] = acoustic_props[k_typ].p_31
                    p_tensor[2,2,0,1,k_typ] = acoustic_props[k_typ].p_36
                    p_tensor[2,2,0,2,k_typ] = acoustic_props[k_typ].p_35
                    p_tensor[2,2,1,0,k_typ] = acoustic_props[k_typ].p_36
                    p_tensor[2,2,1,1,k_typ] = acoustic_props[k_typ].p_32
                    p_tensor[2,2,1,2,k_typ] = acoustic_props[k_typ].p_34
                    p_tensor[2,2,2,0,k_typ] = acoustic_props[k_typ].p_35
                    p_tensor[2,2,2,1,k_typ] = acoustic_props[k_typ].p_34
                    p_tensor[2,2,2,2,k_typ] = acoustic_props[k_typ].p_33

                    eta_tensor[0,0,0,0,k_typ] = acoustic_props[k_typ].eta_11
                    eta_tensor[0,0,0,1,k_typ] = acoustic_props[k_typ].eta_16
                    eta_tensor[0,0,0,2,k_typ] = acoustic_props[k_typ].eta_15
                    eta_tensor[0,0,1,0,k_typ] = acoustic_props[k_typ].eta_16
                    eta_tensor[0,0,1,1,k_typ] = acoustic_props[k_typ].eta_12
                    eta_tensor[0,0,1,2,k_typ] = acoustic_props[k_typ].eta_14
                    eta_tensor[0,0,2,0,k_typ] = acoustic_props[k_typ].eta_15
                    eta_tensor[0,0,2,1,k_typ] = acoustic_props[k_typ].eta_14
                    eta_tensor[0,0,2,2,k_typ] = acoustic_props[k_typ].eta_13
                    eta_tensor[0,1,0,0,k_typ] = acoustic_props[k_typ].eta_61
                    eta_tensor[0,1,0,1,k_typ] = acoustic_props[k_typ].eta_66
                    eta_tensor[0,1,0,2,k_typ] = acoustic_props[k_typ].eta_65
                    eta_tensor[0,1,1,0,k_typ] = acoustic_props[k_typ].eta_66
                    eta_tensor[0,1,1,1,k_typ] = acoustic_props[k_typ].eta_62
                    eta_tensor[0,1,1,2,k_typ] = acoustic_props[k_typ].eta_64
                    eta_tensor[0,1,2,0,k_typ] = acoustic_props[k_typ].eta_65
                    eta_tensor[0,1,2,1,k_typ] = acoustic_props[k_typ].eta_64
                    eta_tensor[0,1,2,2,k_typ] = acoustic_props[k_typ].eta_63
                    eta_tensor[0,2,0,0,k_typ] = acoustic_props[k_typ].eta_51
                    eta_tensor[0,2,0,1,k_typ] = acoustic_props[k_typ].eta_56
                    eta_tensor[0,2,0,2,k_typ] = acoustic_props[k_typ].eta_55
                    eta_tensor[0,2,1,0,k_typ] = acoustic_props[k_typ].eta_56
                    eta_tensor[0,2,1,1,k_typ] = acoustic_props[k_typ].eta_52
                    eta_tensor[0,2,1,2,k_typ] = acoustic_props[k_typ].eta_54
                    eta_tensor[0,2,2,0,k_typ] = acoustic_props[k_typ].eta_55
                    eta_tensor[0,2,2,1,k_typ] = acoustic_props[k_typ].eta_54
                    eta_tensor[0,2,2,2,k_typ] = acoustic_props[k_typ].eta_53
                    eta_tensor[1,0,0,0,k_typ] = acoustic_props[k_typ].eta_61
                    eta_tensor[1,0,0,1,k_typ] = acoustic_props[k_typ].eta_66
                    eta_tensor[1,0,0,2,k_typ] = acoustic_props[k_typ].eta_65
                    eta_tensor[1,0,1,0,k_typ] = acoustic_props[k_typ].eta_66
                    eta_tensor[1,0,1,1,k_typ] = acoustic_props[k_typ].eta_62
                    eta_tensor[1,0,1,2,k_typ] = acoustic_props[k_typ].eta_64
                    eta_tensor[1,0,2,0,k_typ] = acoustic_props[k_typ].eta_65
                    eta_tensor[1,0,2,1,k_typ] = acoustic_props[k_typ].eta_64
                    eta_tensor[1,0,2,2,k_typ] = acoustic_props[k_typ].eta_63
                    eta_tensor[1,1,0,0,k_typ] = acoustic_props[k_typ].eta_21
                    eta_tensor[1,1,0,1,k_typ] = acoustic_props[k_typ].eta_26
                    eta_tensor[1,1,0,2,k_typ] = acoustic_props[k_typ].eta_25
                    eta_tensor[1,1,1,0,k_typ] = acoustic_props[k_typ].eta_26
                    eta_tensor[1,1,1,1,k_typ] = acoustic_props[k_typ].eta_22
                    eta_tensor[1,1,1,2,k_typ] = acoustic_props[k_typ].eta_24
                    eta_tensor[1,1,2,0,k_typ] = acoustic_props[k_typ].eta_25
                    eta_tensor[1,1,2,1,k_typ] = acoustic_props[k_typ].eta_24
                    eta_tensor[1,1,2,2,k_typ] = acoustic_props[k_typ].eta_23
                    eta_tensor[1,2,0,0,k_typ] = acoustic_props[k_typ].eta_41
                    eta_tensor[1,2,0,1,k_typ] = acoustic_props[k_typ].eta_46
                    eta_tensor[1,2,0,2,k_typ] = acoustic_props[k_typ].eta_45
                    eta_tensor[1,2,1,0,k_typ] = acoustic_props[k_typ].eta_46
                    eta_tensor[1,2,1,1,k_typ] = acoustic_props[k_typ].eta_42
                    eta_tensor[1,2,1,2,k_typ] = acoustic_props[k_typ].eta_44
                    eta_tensor[1,2,2,0,k_typ] = acoustic_props[k_typ].eta_45
                    eta_tensor[1,2,2,1,k_typ] = acoustic_props[k_typ].eta_44
                    eta_tensor[1,2,2,2,k_typ] = acoustic_props[k_typ].eta_43
                    eta_tensor[2,0,0,0,k_typ] = acoustic_props[k_typ].eta_51
                    eta_tensor[2,0,0,1,k_typ] = acoustic_props[k_typ].eta_56
                    eta_tensor[2,0,0,2,k_typ] = acoustic_props[k_typ].eta_55
                    eta_tensor[2,0,1,0,k_typ] = acoustic_props[k_typ].eta_56
                    eta_tensor[2,0,1,1,k_typ] = acoustic_props[k_typ].eta_52
                    eta_tensor[2,0,1,2,k_typ] = acoustic_props[k_typ].eta_54
                    eta_tensor[2,0,2,0,k_typ] = acoustic_props[k_typ].eta_55
                    eta_tensor[2,0,2,1,k_typ] = acoustic_props[k_typ].eta_54
                    eta_tensor[2,0,2,2,k_typ] = acoustic_props[k_typ].eta_53
                    eta_tensor[2,1,0,0,k_typ] = acoustic_props[k_typ].eta_41
                    eta_tensor[2,1,0,1,k_typ] = acoustic_props[k_typ].eta_46
                    eta_tensor[2,1,0,2,k_typ] = acoustic_props[k_typ].eta_45
                    eta_tensor[2,1,1,0,k_typ] = acoustic_props[k_typ].eta_46
                    eta_tensor[2,1,1,1,k_typ] = acoustic_props[k_typ].eta_42
                    eta_tensor[2,1,1,2,k_typ] = acoustic_props[k_typ].eta_44
                    eta_tensor[2,1,2,0,k_typ] = acoustic_props[k_typ].eta_45
                    eta_tensor[2,1,2,1,k_typ] = acoustic_props[k_typ].eta_44
                    eta_tensor[2,1,2,2,k_typ] = acoustic_props[k_typ].eta_43
                    eta_tensor[2,2,0,0,k_typ] = acoustic_props[k_typ].eta_31
                    eta_tensor[2,2,0,1,k_typ] = acoustic_props[k_typ].eta_36
                    eta_tensor[2,2,0,2,k_typ] = acoustic_props[k_typ].eta_35
                    eta_tensor[2,2,1,0,k_typ] = acoustic_props[k_typ].eta_36
                    eta_tensor[2,2,1,1,k_typ] = acoustic_props[k_typ].eta_32
                    eta_tensor[2,2,1,2,k_typ] = acoustic_props[k_typ].eta_34
                    eta_tensor[2,2,2,0,k_typ] = acoustic_props[k_typ].eta_35
                    eta_tensor[2,2,2,1,k_typ] = acoustic_props[k_typ].eta_34
                    eta_tensor[2,2,2,2,k_typ] = acoustic_props[k_typ].eta_33
                else:
                    raise ValueError("symmetry_flag must be True or False.")

        self.rho = rho
        self.c_tensor = c_tensor
        self.c_tensor_z = c_tensor_z
        self.p_tensor = p_tensor
        self.eta_tensor = eta_tensor

        self.linear_element_shapes = ['rectangular', 'slot', 'rib_coated']
        self.curvilinear_element_shapes = ['circular', 'onion']


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
                    raise NotImplementedError("Have not implemented 2 coated inclusions.")
                elif self.coat_y is not None and self.inc_b_x is None:
                        raise NotImplementedError("Have not implemented 1 coated inclusions.")
                else:
                    raise ValueError("NumBAT doesn't understand you geometry.")
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
                    raise NotImplementedError("Have not implemented 2 coated inclusions.")
                elif self.coat_y is not None and self.inc_b_x is None:
                        raise NotImplementedError("Have not implemented 1 coated inclusions.")
                else:
                    raise ValueError("NumBAT doesn't understand you geometry.")
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
                    raise NotImplementedError("Have not implemented 2 coated inclusions.")
                elif self.coat_y is not None and self.inc_b_x is None:
                        raise NotImplementedError("Have not implemented 1 coated inclusions.")
                else:
                    raise ValueError("NumBAT doesn't understand you geometry.")
            else:
                raise ValueError("must have at least one nonzero inclusion.")
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
            if self.coat_y is not None:
                msh_template = 'slot_coated'
                self.nb_typ_el = 5
                msh_name = 'slot_c_%(d)s_%(dy)s_%(a)s_%(b)s_%(c)s_%(cc)s_%(ccc)s' % {
                'd': dec_float_str(self.unitcell_x),
                'dy': dec_float_str(self.unitcell_y),
                'a': dec_float_str(self.inc_a_x),
                'b': dec_float_str(self.inc_a_y),
                'c': dec_float_str(self.inc_b_x),
                'cc': dec_float_str(self.slab_a_y),
                'ccc': dec_float_str(self.coat_y)}
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
            else:
                msh_template = 'slot'
                self.nb_typ_el = 4
                msh_name = 'slot_%(d)s_%(dy)s_%(a)s_%(b)s_%(c)s_%(cc)s' % {
                'd': dec_float_str(self.unitcell_x),
                'dy': dec_float_str(self.unitcell_y),
                'a': dec_float_str(self.inc_a_x),
                'b': dec_float_str(self.inc_a_y),
                'c': dec_float_str(self.inc_b_x),
                'cc': dec_float_str(self.slab_a_y)}
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

        elif self.inc_shape in ['onion']:
            msh_template = 'onion'
            self.nb_typ_el = 15
            msh_name = 'onion_%(d)s_%(dy)s_%(a)s_%(b)s_%(c)s_%(d)s_%(e)s_%(f)s_%(g)s' % {
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


        elif self.inc_shape in ['rib_coated']:
                msh_template = 'rib_coated'
                self.nb_typ_el = 4
                msh_name = 'rib_coated%(d)s_%(dy)s_%(a)s_%(b)s_%(cz)s_%(c)s_%(cc)s_%(ccc)s' % {
                'd': dec_float_str(self.unitcell_x),
                'dy': dec_float_str(self.unitcell_y),
                'a': dec_float_str(self.inc_a_x),
                'b': dec_float_str(self.inc_a_y),
                'cz': dec_float_str(self.coat_x),
                'c': dec_float_str(self.coat_y),
                'cc': dec_float_str(self.slab_a_y),
                'ccc': dec_float_str(self.slab_b_y)}
                if not os.path.exists(msh_location + msh_name + '.mail') or self.force_mesh is True:
                    geo_tmp = open(msh_location + '%s_msh_template.geo' % msh_template, "r").read()
                    geo = geo_tmp.replace('d_in_nm = 1000;', "d_in_nm = %f;" % self.unitcell_x)
                    geo = geo.replace('dy_in_nm = 600;', "dy_in_nm = %f;" % self.unitcell_y)
                    geo = geo.replace('ribx = 200;', "ribx = %f;" % self.inc_a_x)
                    geo = geo.replace('riby = 30;', "riby = %f;" % self.inc_a_y)
                    geo = geo.replace('coatx = 20;', "coatx = %f;" % self.coat_x)
                    geo = geo.replace('coaty = 20;', "coaty = %f;" % self.coat_y)
                    geo = geo.replace('slab1 = 100;', "slab1 = %f;" % self.slab_a_y)
                    geo = geo.replace('slab2 = 50;', "slab2 = %f;" % self.slab_b_y)
                    geo = geo.replace('lc = 0;', "lc = %f;" % self.lc)
                    geo = geo.replace('lc2 = lc/1;', "lc2 = lc/%f;" % self.lc2)
                    geo = geo.replace('lc3 = lc/1;', "lc3 = lc/%f;" % self.lc3)

        else:
            raise NotImplementedError("\n Selected inc_shape = '%s' \n " \
            "is not currently implemented. Please make a mesh with gmsh, & \n " \
            "consider contributing this to NumBAT via gitlab." % self.inc_shape)

        self.mesh_file = msh_name + '.mail'
        open(msh_location + msh_name + '.geo', "w").write(geo)
        NumBAT.conv_gmsh(msh_location+msh_name)

        if self.check_msh is True:
            # Automatically show created mesh in gmsh.
            gmsh_cmd = 'gmsh '+ msh_location + msh_name + '.msh'
            os.system(gmsh_cmd)
            gmsh_cmd = 'gmsh '+ msh_location + msh_name + '.geo'
            os.system(gmsh_cmd)


    def calc_EM_modes(self, wl_nm, num_modes, n_eff, Stokes=False, **args):
        """ Run a simulation to find the Struct's EM modes.

            Args:
                wl_nm  (float): Wavelength of EM wave in vacuum.

                num_modes  (int): Number of EM modes to solve for.

                n_eff  (float): Guesstimated effective index of
                    fundamental mode, will be origin of FEM search.

            Returns:
                :Simmo: object
        """
        simmo = Simmo(self, wl_nm, num_modes=num_modes, n_eff=n_eff, Stokes=Stokes)

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
    string = '%8.2f' %dec_float
    string = string.replace('.', '_')
    string = string.replace(" ", "")
    return string
