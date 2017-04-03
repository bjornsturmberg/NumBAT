"""
    mode_calcs.py is a subroutine of NumBAT that contains methods to
    calculate the EM and Acoustic modes of a structure.

    Copyright (C) 2016 Bjorn Sturmberg, Kokou Dossou
"""

import numpy as np
import sys
import os
sys.path.append("../backend/")

import plotting
from fortran import NumBAT


class Simmo(object):
    """ Calculates the modes of a :Struc: object at a wavelength of wl_nm.
    """
    def __init__(self, structure, wl_nm, num_modes=20, n_eff=None, shift_Hz=None, 
                 k_AC=None, EM_sim=None):
        self.structure = structure
        self.wl_m = wl_nm*1e-9
        self.n_eff = n_eff
        self.shift_Hz = shift_Hz
        self.k_AC = k_AC
        self.EM_sim = EM_sim
        self.num_modes = num_modes
        self.mode_pol = None
        self.k_0 = 2 * np.pi / self.wl_m
        # just off normal incidence to avoid degeneracies
        self.k_pll = np.array([1e-16, 1e-16])
        speed_c = 299792458
        self.omega_EM = 2*np.pi*speed_c/self.wl_m # Angular freq in units of Hz


    def calc_EM_modes(self):
        """ Run a Fortran FEM calculation to find the optical modes.

        Returns a :Simmo: object that now has these key values:

        Eig_values: a 1d array of Eigenvalues (propagation constants) in [1/m]

        sol1: the associated Eigenvectors, ie. the fields, stored as
               [field comp, node nu on element, Eig value, el nu]
        """
        self.d_in_m = self.structure.unitcell_x*1e-9
        n_list = []
        n_list_tmp = np.array([self.structure.material_a.n, self.structure.material_b.n, self.structure.material_c.n,
                               self.structure.material_d.n, self.structure.material_e.n, self.structure.material_f.n,
                               self.structure.material_g.n, self.structure.material_h.n, self.structure.material_i.n,
                               self.structure.material_j.n, self.structure.material_k.n, self.structure.material_l.n,
                               self.structure.material_m.n, self.structure.material_n.n, self.structure.material_o.n,
                               self.structure.material_p.n, self.structure.material_q.n, self.structure.material_r.n])
                               # self.structure.bkg_material.n,
                               # self.structure.inc_a_material.n,
                               # self.structure.slab_a_material.n,
                               # self.structure.slab_a_bkg_material.n,
                               # self.structure.slab_b_material.n,
                               # self.structure.slab_b_bkg_material.n,
                               # self.structure.inc_b_material.n,
                               # self.structure.coat_material.n])
        self.el_conv_table_n = {}
        i = 1; j = 1
        for n in n_list_tmp:
            if n != 0:
                n_list.append(n)
                self.el_conv_table_n[i] = j
                j += 1
            i += 1
        self.n_list = np.array(n_list)

        if self.structure.loss is False:
            self.n_list = self.n_list.real

        if self.num_modes < 20:
            self.num_modes = 20
            print("Warning: ARPACK needs >= 20 modes so set num_modes=20.")

        # Parameters that control how FEM routine runs
        self.E_H_field = 1  # Selected formulation (1=E-Field, 2=H-Field)
        i_cond = 2  # Boundary conditions (0=Dirichlet,1=Neumann,2=unitcell_x)
        itermax = 30  # Maximum number of iterations for convergence
        EM_FEM_debug = 0  # Fortran routines will display & save add. info
        # Size of Fortran's complex superarray (scales with mesh)
        # In theory could do some python-based preprocessing
        # on the mesh file to work out RAM requirements
        cmplx_max = 2**28
        real_max = 2**27
        int_max = 2**26

        # Calculate where to center the Eigenmode solver around.
        # (Shift and invert FEM method)
        shift = self.n_eff**2 * self.k_0**2

        if EM_FEM_debug == 1:
            if not os.path.exists("Normed"):
                os.mkdir("Normed")
            if not os.path.exists("Matrices"):
                os.mkdir("Matrices")
            if not os.path.exists("Output"):
                os.mkdir("Output")

        with open("../backend/fortran/msh/"+self.structure.mesh_file) as f:
            self.n_msh_pts, self.n_msh_el = [int(i) for i in f.readline().split()]

        try:
            resm = NumBAT.calc_em_modes(
                self.wl_m, self.num_modes,
                EM_FEM_debug, self.structure.mesh_file, self.n_msh_pts,
                self.n_msh_el, self.structure.nb_typ_el, self.n_list,
                self.k_pll, self.d_in_m, shift, self.E_H_field, i_cond, itermax,
                self.structure.plotting_fields, self.structure.plot_real,
                self.structure.plot_imag, self.structure.plot_abs,
                cmplx_max, real_max, int_max)

            self.Eig_values, self.sol1, self.mode_pol, \
            self.table_nod, self.type_el, self.type_nod, self.x_arr = resm

        except KeyboardInterrupt:
            print("\n\n FEM routine calc_EM_modes",\
            "interrupted by keyboard.\n\n")

        # if not self.structure.plot_field_conc:
        #     self.mode_pol = None

        # if self.structure.plotting_fields != 1:
        #     self.sol1 = None
        #     self.n_list = None
        #     self.E_H_field = None
        #     self.table_nod = None
        #     self.type_el = None
        #     self.x_arr = None
        #     self.n_msh_pts = None
        #     self.n_msh_el = None


### Calc unnormalised power in each EM mode Kokou equiv. of Eq. 8.
        try:
            nnodes = 6
            if self.structure.inc_shape in self.structure.linear_element_shapes:
            ## Integration using analytically evaluated basis function integrals. Fast.
                self.EM_mode_overlap = NumBAT.em_mode_energy_int_v2_ez(
                    self.k_0, self.num_modes, self.n_msh_el, self.n_msh_pts,
                    nnodes, self.table_nod,
                    self.x_arr, self.Eig_values, self.sol1)
            else:
                if self.structure.inc_shape != 'circular':
                    print("Warning: em_mode_energy_int - not sure if mesh contains curvi-linear elements", 
                        "\n using slow quadrature integration by default.\n\n")
            # Integration by quadrature. Slowest.
                self.EM_mode_overlap = NumBAT.em_mode_energy_int_ez(
                    self.k_0, self.num_modes, self.n_msh_el, self.n_msh_pts,
                    nnodes, self.table_nod,
                    self.x_arr, self.Eig_values, self.sol1)
            # Bring Kokou's def into line with CW formulation.
            self.EM_mode_overlap = 2.0*self.EM_mode_overlap

        except KeyboardInterrupt:
            print("\n\n FEM routine EM_mode_energy_int interrupted by keyboard.\n\n")

        ### Not necessary because EM FEM mesh always normalised in area to unity.
        # print area
        # x_tmp = []
        # y_tmp = []
        # for i in np.arange(self.n_msh_pts):
        #     x_tmp.append(self.x_arr[0,i])
        #     y_tmp.append(self.x_arr[1,i])
        # x_min = np.min(x_tmp); x_max=np.max(x_tmp)
        # y_min = np.min(y_tmp); y_max=np.max(y_tmp)
        # area = abs((x_max-x_min)*(y_max-y_min))
        # print area
        # self.EM_mode_overlap = self.EM_mode_overlap*area


    def calc_AC_modes(self):
        """ Run a Fortran FEM calculation to find the acoustic modes.

        Returns a :Simmo: object that now has these key values:

        Eig_values: a 1d array of Eigenvalues (frequencies) in [1/s]

        sol1: the associated Eigenvectors, ie. the fields, stored as
               [field comp, node nu on element, Eig value, el nu]
        """
        self.d_in_m = self.structure.inc_a_x*1e-9

        el_conv_table = {}
        i = 1; j = 1
        for matter in self.structure.acoustic_props_tmp:
            if matter.s != None:
                el_conv_table[i] = j
                j += 1
            i += 1
        final_dict = {}
        for entry in el_conv_table:
            # print entry, self.EM_sim.el_conv_table_n[entry], el_conv_table[entry]
            final_dict[self.EM_sim.el_conv_table_n[entry]] = el_conv_table[entry]
        # print final_dict
        self.typ_el_AC = final_dict

        if self.num_modes < 20:
            self.num_modes = 20
            print("Warning: ARPACK needs >= 20 modes so set num_modes=20.")

        # Parameters that control how FEM routine runs
        i_cond = 1  # Boundary conditions (0=Dirichlet,1=Neumann,2=unitcell_x)
        itermax = 30  # Maximum number of iterations for convergence
        AC_FEM_debug = 0  # Fortran routines will display & save add. info
        ARPACK_tol = 1e-10  # ARPACK accuracy (0.0 for machine precision)
        # Size of Fortran's complex superarray (scales with mesh)
        # In theory could do some python-based preprocessing
        # on the mesh file to work out RAM requirements
        cmplx_max = 2**28
        real_max = 2**27
        int_max = 2**26


        # Calculate where to center the Eigenmode solver around.
        # (Shift and invert FEM method)
        if self.shift_Hz is None:
            # For AC problem shift is a frequency; [shift] = s^-1.
            v_list = []
            for el in range(self.structure.nb_typ_el_AC):
                # Using acoustic velocity of longitudinal mode pg 215 Auld vol 1.
                v_list.append(np.sqrt(self.structure.c_tensor[0,0][el]/self.structure.rho[el]))
                # # Using acoustic velocity of shear mode pg 215 Auld vol 1.
                # v_list.append(np.sqrt(self.structure.c_tensor[3,3][el]/self.structure.rho[el]))
            AC_velocity = np.real(v_list).min()
            shift = np.real(AC_velocity*self.k_AC/(2.*np.pi))
            # print "shift", shift
            shift = 0.9*shift
            # print "shift", shift
        else:
            shift = self.shift_Hz

        # Take existing msh from EM FEM and manipulate mesh to exclude vacuum areas.
        if self.EM_sim:
            suplied_geo_flag = 1
            n_msh_el = self.EM_sim.n_msh_el
            n_msh_pts = self.EM_sim.n_msh_pts
            type_el = self.EM_sim.type_el
            type_nod = self.EM_sim.type_nod
            table_nod = self.EM_sim.table_nod
            x_arr = self.EM_sim.x_arr
            n_el_kept = 0
            n_msh_pts_AC = 0
            type_el_AC = []
            table_nod_AC_tmp = np.zeros(np.shape(table_nod))
            el_convert_tbl = {}
            el_convert_tbl_inv = {}
            node_convert_tbl = {}
            if AC_FEM_debug == 1:
                plotting.plot_msh(x_arr, 'orig')

            for el in range(n_msh_el):
                # print type_el[el]
                if type_el[el] in self.typ_el_AC:
                    # print "in", type_el[el]
                    type_el_AC.append(self.typ_el_AC[type_el[el]])
                    el_convert_tbl[n_el_kept] = el
                    el_convert_tbl_inv[el] = n_el_kept
                    for i in range(6):
                        # Leaves node numbering untouched
                        table_nod_AC_tmp[i][n_el_kept] = table_nod[i][el]
                    n_el_kept += 1
            n_msh_el_AC = n_el_kept
            # Find unique nodes
            node_lst_tmp = []
            for el in range(n_msh_el_AC):
                for i in range(6):
                    node_lst_tmp.append(table_nod_AC_tmp[i][el])
            unique_nodes = list(set(node_lst_tmp))
            n_msh_pts_AC = len(unique_nodes)
            unique_nodes = [int(j) for j in unique_nodes]
            # Mapping unique nodes to start from zero
            for i in range(n_msh_pts_AC):
                node_convert_tbl[unique_nodes[i]] = i
            # Creating finalised table_nod.
            table_nod_AC = []
            for i in range(6):
                el_tbl = []
                for el in range(n_msh_el_AC):
                    # Note table_nod needs to be adjust back to fortran indexing
                    el_tbl.append(node_convert_tbl[table_nod_AC_tmp[i][el]]+1)
                table_nod_AC.append(el_tbl)
            # Find the coordinates of chosen nodes.
            x_arr_AC = np.zeros((2,n_msh_pts_AC))
            for node in unique_nodes:
                # Note x_arr needs to be adjust back to fortran indexing
                x_arr_AC[0,node_convert_tbl[node]] = (x_arr[0,node-1])
                x_arr_AC[1,node_convert_tbl[node]] = (x_arr[1,node-1])

            self.el_convert_tbl = el_convert_tbl
            self.el_convert_tbl_inv = el_convert_tbl_inv
            self.node_convert_tbl = node_convert_tbl


            ### AC FEM uses Neumann B.C.s so type_nod is totally irrelevant!
            # # Find nodes on boundaries of materials
            # node_array = -1*np.ones(n_msh_pts)
            # interface_nodes = []
            # for el in range(n_msh_el):
            #     for i in range(6):
            #         node = table_nod[i][el]
            #         # Check if first time seen this node
            #         if node_array[node - 1] == -1: # adjust to python indexing
            #             node_array[node - 1] = type_el[el]
            #         else:
            #             if node_array[node - 1] != type_el[el]:
            #                 interface_nodes.append(node)
            # interface_nodes = list(set(interface_nodes))
            type_nod_AC = np.zeros(n_msh_pts_AC)

            # import matplotlib
            # matplotlib.use('pdf')
            # import matplotlib.pyplot as plt
            # plt.clf()
            # plt.figure(figsize=(13,13))
            # ax = plt.subplot(1,1,1)
            # for node in unique_nodes:
            #     if node in interface_nodes:
            #         type_nod_AC[node_convert_tbl[node]] = i_cond
            #         plt.plot(x_arr_AC[0,node_convert_tbl[node]], x_arr_AC[1,node_convert_tbl[node]], 'ok')
            # ax.set_aspect('equal')
            # plt.savefig('boundary.pdf', bbox_inches='tight')
            self.n_msh_pts = n_msh_pts_AC
            self.n_msh_el = n_msh_el_AC
        # Default, indicates to use geometry subroutine in FEM routine.
        else:
            suplied_geo_flag = 0
            with open("../backend/fortran/msh/"+self.structure.mesh_file) as f:
                self.n_msh_pts, self.n_msh_el = [int(i) for i in f.readline().split()]
            table_nod_AC = np.zeros((6, self.n_msh_el))
            type_el_AC = np.zeros(self.n_msh_el)
            x_arr_AC = np.zeros((2,self.n_msh_pts))
            type_nod_AC = np.zeros(self.n_msh_pts)

        if AC_FEM_debug == 1:
            print('shift', shift)
            if not os.path.exists("Normed"):
                os.mkdir("Normed")
            if not os.path.exists("Output"):
                os.mkdir("Output")
            if not os.path.exists("Matrices"):
                os.mkdir("Matrices")

        try:
            resm = NumBAT.calc_ac_modes(
                self.k_AC, self.num_modes,
                AC_FEM_debug, self.structure.mesh_file, self.n_msh_pts,
                self.n_msh_el, self.structure.nb_typ_el_AC,
                self.structure.c_tensor, self.structure.rho,
                self.d_in_m, shift, i_cond, itermax, ARPACK_tol,
                self.structure.plotting_fields,
                cmplx_max, real_max, int_max, suplied_geo_flag, type_nod_AC, 
                self.structure.symmetry_flag, table_nod_AC, type_el_AC, x_arr_AC)
            table_nod_out, type_el_out, x_arr_out, \
            self.Eig_values, self.sol1, self.mode_pol = resm

            # FEM Eigenvalue is frequency, rather than angular frequency Omega
            self.Omega_AC = self.Eig_values*2*np.pi

        except KeyboardInterrupt:
            print("\n\n FEM routine calc_AC_modes",\
            "interrupted by keyboard.\n\n")

        if AC_FEM_debug == 1:
            plotting.plot_msh(x_arr_AC, 'in')
            plotting.plot_msh(x_arr_out, 'out')

        # if self.EM_sim is None:
        #     table_nod_out = None
        #     type_el_out = None
        #     x_arr_out = None
        #     self.table_nod = table_nod_AC
        #     self.type_el = type_el_AC
        #     self.x_arr = x_arr_AC
        # else:
        self.table_nod = table_nod_out
        self.type_el = type_el_out
        self.x_arr = x_arr_out

### Calc unnormalised power in each AC mode Eq. 18.
        try:
            nnodes = 6
            # import time
            # start = time.time()
            if self.structure.inc_shape in self.structure.linear_element_shapes:
            # Semi-analytic integration following KD 9/9/16 notes. Fastest!
                self.AC_mode_overlap = NumBAT.ac_mode_energy_int_v4(
                    self.num_modes, self.n_msh_el, self.n_msh_pts,
                    nnodes, self.table_nod, self.type_el, self.x_arr,
                    self.structure.nb_typ_el_AC, self.structure.c_tensor,
                    self.k_AC, self.Omega_AC, self.sol1)
            else:
                if self.structure.inc_shape != 'circular':
                    print("Warning: ac_mode_energy_int - not sure if mesh contains curvi-linear elements", 
                        "\n using slow quadrature integration by default.\n\n")
            # Integration by quadrature. Slowest.
                self.AC_mode_overlap = NumBAT.ac_mode_energy_int(
                    self.num_modes, self.n_msh_el, self.n_msh_pts,
                    nnodes, self.table_nod, self.type_el, self.x_arr,
                    self.structure.nb_typ_el_AC, self.structure.c_tensor_z,
                    self.k_AC, self.Omega_AC, self.sol1, AC_FEM_debug)
        except KeyboardInterrupt:
            print("\n\n FEM routine AC_mode_energy_int interrupted by keyboard.\n\n")
