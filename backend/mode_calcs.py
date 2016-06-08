"""
    mode_calcs.py is a subroutine of NumBAT that contains methods to
    calculate the EM and Acoustic modes of a structure.

    Copyright (C) 2016 Bjorn Sturmberg, Kokou Dossou
"""

import numpy as np
import sys
# from scipy import sqrt
import os
sys.path.append("../backend/")

from fortran import NumBAT

pi = np.pi


class Simmo(object):
    """ Interaction of one :Light: object with one :Struc: object.

        Inherits knowledge of :Struc:, :Light: objects
        Stores the calculated modes of :Struc: for illumination by :Light:
    """
    def __init__(self, structure, wl_nm, q_acoustic=None, num_modes=20):
        self.structure = structure
        self.wl_nm = wl_nm
        self.q_acoustic = q_acoustic
        self.num_modes = num_modes
        self.mode_pol = None
        # just off normal incidence to avoid degeneracies
        self.k_pll = np.array([1e-16, 1e-16])

    def k_pll_norm(self):
        return self.k_pll * self.structure.unitcell_x

    def wl_norm(self):
        """ Return normalised wavelength (wl/unitcell_x). """
        wl = self.wl_nm / self.structure.unitcell_x
        # Avoid Wood Anomalies
        if self.wl_nm % self.structure.unitcell_x == 0:
            wl += 1e-10
        return wl

    def calc_EM_modes(self):
        """ Run a Fortran FEM calculation to find the EM modes \
        of a structure. """
        st = self.structure
        wl = self.wl_nm
        self.n_effs = np.array([st.bkg_material.n(wl), st.inc_a_material.n(wl),
                                st.inc_b_material.n(wl), st.slab_a_material.n(wl),
                                st.slab_a_bkg_material.n(wl), st.slab_b_material.n(wl),
                                st.slab_b_bkg_material.n(wl), st.coating_material.n(wl)])

        self.n_effs = self.n_effs[:self.structure.nb_typ_el]
        if self.structure.loss is False:
            self.n_effs = self.n_effs.real

        if self.num_modes < 20:
            self.num_modes = 20
            print "Warning: ARPACK needs >= 20 modes so set num_modes=20."

        # Parameters that control how FEM routine runs
        self.E_H_field = 1  # Selected formulation (1=E-Field, 2=H-Field)
        i_cond = 2  # Boundary conditions (0=Dirichlet,1=Neumann,2=unitcell_x)
        itermax = 30  # Maximum number of iterations for convergence
        EM_FEM_debug = 0  # Fortran routines will display & save add. info

        # Calculate where to center the Eigenmode solver around.
        # (Shift and invert FEM method)
        max_n = np.real(self.n_effs).max()
        # Take real part so that complex conjugate pair Eigenvalues are
        # equal distance from shift and invert point and therefore both found.
        k_0 = 2 * pi / self.wl_norm()
        shift = 1.1*max_n**2 * k_0**2  \
            - self.k_pll_norm()[0]**2 - self.k_pll_norm()[1]**2

        if EM_FEM_debug == 1:
            print 'shift', shift
            if not os.path.exists("Normed"):
                os.mkdir("Normed")
            if not os.path.exists("Matrices"):
                os.mkdir("Matrices")
            if not os.path.exists("Output"):
                os.mkdir("Output")

        with open("../backend/fortran/msh/"+self.structure.mesh_file) as f:
            self.n_msh_pts, self.n_msh_el = [int(i) for i in f.readline().split()]

        # Size of Fortran's complex superarray (scales with mesh)
        # In theory could do some python-based preprocessing
        # on the mesh file to work out RAM requirements
        cmplx_max = 2**27  # 30
        real_max = 2**23
        int_max = 2**22

        try:
            resm = NumBAT.calc_em_modes(
                self.wl_norm(), self.num_modes,
                EM_FEM_debug, self.structure.mesh_file, self.n_msh_pts,
                self.n_msh_el, self.structure.nb_typ_el, self.n_effs,
                self.k_pll_norm(), shift, self.E_H_field, i_cond, itermax,
                self.structure.plotting_fields, self.structure.plot_real,
                self.structure.plot_imag, self.structure.plot_abs,
                cmplx_max, real_max, int_max)

            self.Eig_value, self.sol1, self.mode_pol, \
            self.table_nod, self.type_el, self.x_arr = resm

            area = self.structure.unitcell_x * self.structure.unitcell_y
            area_norm = area/self.structure.unitcell_x**2

        except KeyboardInterrupt:
            print "\n\n FEM routine calc_EM_modes",\
            "interrupted by keyboard.\n\n"

        # if not self.structure.plot_field_conc:
        #     self.mode_pol = None

        # if self.structure.plotting_fields != 1:
        #     self.sol1 = None
        #     self.n_effs = None
        #     self.E_H_field = None
        #     self.table_nod = None
        #     self.type_el = None
        #     self.x_arr = None
        #     self.n_msh_pts = None
        #     self.n_msh_el = None


    def calc_AC_modes(self):
        """ Run a Fortran FEM calculation to find the EM modes \
        of a structure. """
        st = self.structure
        wl = self.wl_nm
        q_acoustic = self.q_acoustic
        self.d_in_m = self.structure.unitcell_x*1e-9

        if self.num_modes < 20:
            self.num_modes = 20
            print "Warning: ARPACK needs >= 20 modes so set num_modes=20."

        # Parameters that control how FEM routine runs
        i_cond = 1  # Boundary conditions (0=Dirichlet,1=Neumann,2=unitcell_x)
        itermax = 30  # Maximum number of iterations for convergence
        AC_FEM_debug = 0  # Fortran routines will display & save add. info
        ARPACK_tol = 1e-10  # ARPACK accuracy (0.0 for machine precision)

        # Calculate where to center the Eigenmode solver around.
        # (Shift and invert FEM method)
        # For AC problem shift is a frequency - [shift] = s^-1.
        # Using acoustic velocity of longitudinal mode pg 215 Auld vol 1.
        shift1 = np.real(np.sqrt(self.structure.c_tensor[0,0][-1]/self.structure.rho[-1]))
        # Factor 2 from q_acoustic being twice beta.
        shift1 = 0.5*self.q_acoustic*shift1
        # Using acoustic velocity of shear mode pg 215 Auld vol 1.
        shift2 = np.real(np.sqrt(self.structure.c_tensor[3,3][-1]/self.structure.rho[-1]))
        shift2 = 0.5*self.q_acoustic*shift2
        shift = 13.0e9  # used for original test case
        # print repr(shift)
        shift = (shift1 + shift2)/8.
        print 'shift', shift

        if AC_FEM_debug == 1:
            print 'shift', shift
            if not os.path.exists("Normed"):
                os.mkdir("Normed")
            if not os.path.exists("Output"):
                os.mkdir("Output")
            if not os.path.exists("Matrices"):
                os.mkdir("Matrices")

        with open("../backend/fortran/msh/"+self.structure.mesh_file) as f:
            self.n_msh_pts, self.n_msh_el = [int(i) for i in f.readline().split()]

        # Size of Fortran's complex superarray (scales with mesh)
        # In theory could do some python-based preprocessing
        # on the mesh file to work out RAM requirements
        cmplx_max = 2**27  # 30
        real_max = 2**23
        int_max = 2**22

        try:
            resm = NumBAT.calc_ac_modes(
                self.wl_norm(), self.q_acoustic, self.num_modes,
                AC_FEM_debug, self.structure.mesh_file, self.n_msh_pts,
                self.n_msh_el, self.structure.nb_typ_el, 
                self.structure.c_tensor, self.structure.rho,
                self.d_in_m, shift, i_cond, itermax, ARPACK_tol,
                self.structure.plotting_fields,
                cmplx_max, real_max, int_max)
            self.Eig_value, self.sol1, self.mode_pol, \
            self.table_nod, self.type_el, self.x_arr = resm

        except KeyboardInterrupt:
            print "\n\n FEM routine calc_AC_modes",\
            "interrupted by keyboard.\n\n"
