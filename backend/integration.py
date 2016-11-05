"""
    mode_calcs.py is a subroutine of NumBAT that contains methods to
    calculate the EM and Acoustic modes of a structure.

    Copyright (C) 2016  Bjorn Sturmberg, Kokou Dossou
"""

import numpy as np
from scipy import interpolate
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

import plotting
from fortran import NumBAT


def gain_and_qs(sim_EM_wguide, sim_AC_wguide, q_acoustic,
                EM_ival1=0, EM_ival2=0, AC_ival=0):
    """ Calculate interaction integrals and SBS gain.

        Implements Eqs. 33, 41, 45 of
        Wolff et al. PRA 92, 013836 (2015) doi/10.1103/PhysRevA.92.013836
        These are for Q_photoelastic, Q_moving_boundary, and the Acoustic loss
        "alpha" respectively.

        Args:
            sim_EM_wguide  (:Simmo: object): Contains all info on EM modes

            sim_AC_wguide  (:Simmo: object): Contains all info on AC modes

            q_acoustic  (float): Propagation constant of acoustic modes.

        Keyword Args:
            EM_ival1  (int/string): Specify mode number of EM mode 1 (pump mode)
                to calculate interactions for.
                Numbering is python index so runs from 0 to num_EM_modes-1,
                with 0 being fundamental mode (largest prop constant).
                Can also set to 'All' to include all modes.

            EM_ival2  (int/string): Specify mode number of EM mode 2 (stokes mode)
                to calculate interactions for.
                Numbering is python index so runs from 0 to num_EM_modes-1,
                with 0 being fundamental mode (largest prop constant).
                Can also set to 'All' to include all modes.

            AC_ival  (int/string): Specify mode number of AC mode
                to calculate interactions for.
                Numbering is python index so runs from 0 to num_AC_modes-1,
                with 0 being fundamental mode (largest prop constant).
                Can also set to 'All' to include all modes.
    """

    # Notes about internals of fortran integration
    # Calc overlap of basis functions (and PE tensor and epsilon)
    # Then use this multiple times for calc of each mode field values

    # phi is values of Lagrange polynomials (1-6) at that node.
    # grad is value of gradient of Lagrange polynomials (1-6) at that node.
    # i variables refer to E field
    # j variables refer to H field
    # ww weight function
    # coeff numerical integration


    if EM_ival1 == 'All':
        EM_ival1_fortran = -1
    else:
        EM_ival1_fortran = EM_ival1+1  # convert back to Fortran indexing
    if EM_ival2 == 'All':
        EM_ival2_fortran = -1
    else:
        EM_ival2_fortran = EM_ival2+1  # convert back to Fortran indexing
    if AC_ival == 'All':
        AC_ival_fortran = -1
    else:
        AC_ival_fortran = AC_ival+1  # convert back to Fortran indexing

    Fortran_debug = 0
    ncomps = 3
    nnodes = 6
    num_modes_EM = sim_EM_wguide.num_modes
    num_modes_AC = sim_AC_wguide.num_modes
    n_msh_el_AC = sim_AC_wguide.n_msh_el
    trimmed_EM_field = np.zeros((ncomps,nnodes,num_modes_EM,n_msh_el_AC), dtype=complex)
    for el in range(n_msh_el_AC):
        new_el = sim_AC_wguide.el_convert_tbl[el]
        for ival in range(num_modes_EM):
            for n in range(nnodes):
                for x in range(ncomps):
                    trimmed_EM_field[x,n,ival,el] = sim_EM_wguide.sol1[x,n,ival,new_el]
    # sim_EM_wguide.sol1 = trimmed_EM_field
    # sim_EM_wguide.n_msh_el = sim_AC_wguide.n_msh_el
    # sim_EM_wguide.n_msh_pts = sim_AC_wguide.n_msh_pts
    # sim_EM_wguide.type_el = sim_AC_wguide.type_el
    # sim_EM_wguide.table_nod = sim_AC_wguide.table_nod
    # sim_EM_wguide.x_arr = sim_AC_wguide.x_arr
    # plotting.plt_mode_fields(sim_EM_wguide, EM_AC='EM', add_name='trim')

    relevant_eps_effs =[]
    for el_typ in range(sim_EM_wguide.structure.nb_typ_el):
        if el_typ+1 in sim_AC_wguide.structure.typ_el_AC:
            relevant_eps_effs.append(sim_EM_wguide.n_effs[el_typ]**2)

    # Calc alpha (loss) Eq. 45
    try:
        if sim_EM_wguide.structure.inc_shape == 'rectangular':
            alpha = NumBAT.ac_alpha_int_v2(sim_AC_wguide.num_modes,
                sim_AC_wguide.n_msh_el, sim_AC_wguide.n_msh_pts, nnodes,
                sim_AC_wguide.table_nod, sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
                sim_AC_wguide.structure.nb_typ_el_AC, sim_AC_wguide.structure.eta_tensor,
                q_acoustic, sim_AC_wguide.Omega_AC, sim_AC_wguide.sol1,
                sim_AC_wguide.AC_mode_overlap)
        elif sim_EM_wguide.structure.inc_shape == 'circular':
            alpha = NumBAT.ac_alpha_int(sim_AC_wguide.num_modes,
                sim_AC_wguide.n_msh_el, sim_AC_wguide.n_msh_pts, nnodes,
                sim_AC_wguide.table_nod, sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
                sim_AC_wguide.structure.nb_typ_el_AC, sim_AC_wguide.structure.eta_tensor,
                q_acoustic, sim_AC_wguide.Omega_AC, sim_AC_wguide.sol1,
                sim_AC_wguide.AC_mode_overlap, Fortran_debug)
    except KeyboardInterrupt:
        print "\n\n Routine ac_alpha_int interrupted by keyboard.\n\n"
    alpha = np.real(alpha)


    # Calc Q_photoelastic Eq. 33
    try:
        if sim_EM_wguide.structure.inc_shape == 'rectangular':
            Q_PE = NumBAT.photoelastic_int_v2(
                sim_EM_wguide.num_modes, sim_AC_wguide.num_modes, EM_ival1_fortran,
                EM_ival2_fortran, AC_ival_fortran, sim_AC_wguide.n_msh_el,
                sim_AC_wguide.n_msh_pts, nnodes,
                sim_AC_wguide.table_nod, sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
                sim_AC_wguide.structure.nb_typ_el_AC, sim_AC_wguide.structure.p_tensor,
                q_acoustic, trimmed_EM_field, sim_AC_wguide.sol1,
                relevant_eps_effs, sim_EM_wguide.Eig_value, Fortran_debug)
        elif sim_EM_wguide.structure.inc_shape == 'circular':
            Q_PE = NumBAT.photoelastic_int(
                sim_EM_wguide.num_modes, sim_AC_wguide.num_modes, EM_ival1_fortran,
                EM_ival2_fortran, AC_ival_fortran, sim_AC_wguide.n_msh_el,
                sim_AC_wguide.n_msh_pts, nnodes,
                sim_AC_wguide.table_nod, sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
                sim_AC_wguide.structure.nb_typ_el_AC, sim_AC_wguide.structure.p_tensor,
                q_acoustic, trimmed_EM_field, sim_AC_wguide.sol1,
                relevant_eps_effs, sim_EM_wguide.Eig_value, Fortran_debug)
    except KeyboardInterrupt:
        print "\n\n Routine photoelastic_int interrupted by keyboard.\n\n"




    # # Load CW Comsol fields
    # import csv
    # with open('ac_mode_1-5.dat', 'rb') as csvfile:
    #     spamreader = csv.reader(csvfile, delimiter=' ')#, quotechar='|')
    #     for header_rows in range(9):
    #         spamreader.next()
    #     x_coord = []; y_coord = []
    #     u0_x = []; u0_y = []; u0_z = []
    #     u1_x = []; u1_y = []; u1_z = []
    #     u2_x = []; u2_y = []; u2_z = []
    #     u3_x = []; u3_y = []; u3_z = []
    #     u4_x = []; u4_y = []; u4_z = []
    #     for row in spamreader:
    #         row = filter(None, row)
    #         row = [float(x) for x in row]
    #         x_coord.append(row[0])
    #         y_coord.append(row[1])
    #         u0_x.append(row[2] + 1j*row[3])
    #         u0_y.append(row[4] + 1j*row[5])
    #         u0_z.append(row[6] + 1j*row[7])
    #         u1_x.append(row[8] + 1j*row[9])
    #         u1_y.append(row[10] + 1j*row[11])
    #         u1_z.append(row[12] + 1j*row[13])
    #         u2_x.append(row[14] + 1j*row[15])
    #         u2_y.append(row[16] + 1j*row[17])
    #         u2_z.append(row[18] + 1j*row[19])
    #         u3_x.append(row[20] + 1j*row[21])
    #         u3_y.append(row[22] + 1j*row[23])
    #         u3_z.append(row[24] + 1j*row[25])
    #         u4_x.append(row[26] + 1j*row[27])
    #         u4_y.append(row[28] + 1j*row[29])
    #         u4_z.append(row[30] + 1j*row[31])

    # x_coord = np.array(x_coord).reshape(100,100)
    # y_coord = np.array(y_coord).reshape(100,100)
    # u0_x = np.array(u0_x).reshape(100,100)
    # u0_y = np.array(u0_y).reshape(100,100)
    # u0_z = np.array(u0_z).reshape(100,100)
    # u1_x = np.array(u1_x).reshape(100,100)
    # u1_y = np.array(u1_y).reshape(100,100)
    # u1_z = np.array(u1_z).reshape(100,100)
    # u2_x = np.array(u2_x).reshape(100,100)
    # u2_y = np.array(u2_y).reshape(100,100)
    # u2_z = np.array(u2_z).reshape(100,100)
    # u3_x = np.array(u3_x).reshape(100,100)
    # u3_y = np.array(u3_y).reshape(100,100)
    # u3_z = np.array(u3_z).reshape(100,100)
    # u4_x = np.array(u4_x).reshape(100,100)
    # u4_y = np.array(u4_y).reshape(100,100)
    # u4_z = np.array(u4_z).reshape(100,100)

    # x_coord = np.swapaxes(x_coord,0,1)
    # y_coord = np.swapaxes(y_coord,0,1)
    # u0_x = np.swapaxes(u0_x,0,1)
    # u0_y = np.swapaxes(u0_y,0,1)
    # u0_z = np.swapaxes(u0_z,0,1)
    # u1_x = np.swapaxes(u1_x,0,1)
    # u1_y = np.swapaxes(u1_y,0,1)
    # u1_z = np.swapaxes(u1_z,0,1)
    # u2_x = np.swapaxes(u2_x,0,1)
    # u2_y = np.swapaxes(u2_y,0,1)
    # u2_z = np.swapaxes(u2_z,0,1)
    # u3_x = np.swapaxes(u3_x,0,1)
    # u3_y = np.swapaxes(u3_y,0,1)
    # u3_z = np.swapaxes(u3_z,0,1)
    # u4_x = np.swapaxes(u4_x,0,1)
    # u4_y = np.swapaxes(u4_y,0,1)
    # u4_z = np.swapaxes(u4_z,0,1)

    # CW_mat = np.array([[u0_x, u0_y, u0_z],[u1_x, u1_y, u1_z],[u2_x, u2_y, u2_z],[u3_x, u3_y, u3_z],[u4_x, u4_y, u4_z]])


    # with open('opt_mode_closeup.dat', 'rb') as csvfile:
    #     spamreader = csv.reader(csvfile, delimiter=' ')#, quotechar='|')
    #     for header_rows in range(9):
    #         spamreader.next()
    #     E_x_CW = []; E_y_CW = []; E_z_CW = []
    #     for row in spamreader:
    #         row = filter(None, row)
    #         row = [float(x) for x in row]
    #         E_x_CW.append(row[2] + 1j*row[3])
    #         E_y_CW.append(row[4] + 1j*row[5])
    #         E_z_CW.append(row[6] + 1j*row[7])
    # E_x_CW = np.array(E_x_CW).reshape(100,100)
    # E_y_CW = np.array(E_y_CW).reshape(100,100)
    # E_z_CW = np.array(E_z_CW).reshape(100,100)
    # E_x_CW = np.swapaxes(E_x_CW,0,1)
    # E_y_CW = np.swapaxes(E_y_CW,0,1)
    # E_z_CW = np.swapaxes(E_z_CW,0,1)
    # E_mat_CW = np.array([E_x_CW,E_y_CW,E_z_CW])
    # # E_mat_CW = np.array([np.zeros(np.shape(E_z_CW)),E_y_CW,E_z_CW])
    # # E_mat_CW = np.array([np.zeros(np.shape(E_z_CW)),np.zeros(np.shape(E_z_CW)),np.zeros(np.shape(E_z_CW))])
    # # E_mat_CW_2 = np.array([E_x_CW,E_y_CW,E_z_CW])
    # # E_mat_CW_2 = np.array([np.zeros(np.shape(E_z_CW)),E_y_CW,E_z_CW])
    # # E_mat_CW_2 = np.array([np.zeros(np.shape(E_z_CW)),np.zeros(np.shape(E_z_CW)),np.zeros(np.shape(E_z_CW))])
    # # print np.max(abs(E_mat_CW[0]))
    # # print np.max(abs(E_mat_CW[1]))
    # # print np.max(abs(E_mat_CW[2]))



# ### Calc AC mode energy, alpha, Q_photoelastic in python
#     # Note: only picking a single EM mode
#     ival_E=0
#     n_points = 100
#     n_pts_x_CW = 100
#     # field mapping
#     x_tmp = []
#     y_tmp = []
#     for i in np.arange(sim_AC_wguide.n_msh_pts):
#         x_tmp.append(sim_AC_wguide.x_arr[0,i])
#         y_tmp.append(sim_AC_wguide.x_arr[1,i])
#     x_min = np.min(x_tmp); x_max=np.max(x_tmp)
#     y_min = np.min(y_tmp); y_max=np.max(y_tmp)
#     area = abs((x_max-x_min)*(y_max-y_min))
#     n_pts_x = n_points#int(n_points*abs(x_max-x_min)/np.sqrt(area))
#     n_pts_y = n_points#int(n_points*abs(y_max-y_min)/np.sqrt(area))
#     v_x=np.zeros(n_pts_x*n_pts_y)
#     v_y=np.zeros(n_pts_x*n_pts_y)
#     i=0
#     for x in np.linspace(x_min,x_max,n_pts_x):
#         for y in np.linspace(y_max,y_min,n_pts_y):
#             v_x[i] = x
#             v_y[i] = y
#             i+=1
#     v_x = np.array(v_x)
#     v_y = np.array(v_y)

#     # unrolling data for the interpolators
#     table_nod = sim_AC_wguide.table_nod.T
#     x_arr = sim_AC_wguide.x_arr.T

#     alpha_py = np.zeros(len(sim_AC_wguide.Eig_value))
#     F = np.zeros(len(sim_AC_wguide.Eig_value), dtype=np.complex128)
#     F_AC = np.zeros(len(sim_AC_wguide.Eig_value), dtype=np.complex128)
#     F_PE = np.zeros((len(sim_EM_wguide.Eig_value),len(sim_EM_wguide.Eig_value),len(sim_AC_wguide.Eig_value)), dtype=np.complex128)
#     alpha_py_CW = np.zeros(len(sim_AC_wguide.Eig_value))
#     F_CW = np.zeros(len(sim_AC_wguide.Eig_value), dtype=np.complex128)
#     F_AC_CW = np.zeros(len(sim_AC_wguide.Eig_value), dtype=np.complex128)
#     F_PE_CW = np.zeros((len(sim_EM_wguide.Eig_value),len(sim_EM_wguide.Eig_value),len(sim_AC_wguide.Eig_value)), dtype=np.complex128)
#     for ival in range(len(sim_AC_wguide.Eig_value)):
#         # dense triangulation with multiple points
#         v_x6p = np.zeros(6*sim_AC_wguide.n_msh_el)
#         v_y6p = np.zeros(6*sim_AC_wguide.n_msh_el)
#         v_Ex6p = np.zeros(6*sim_AC_wguide.n_msh_el, dtype=np.complex128)
#         v_Ey6p = np.zeros(6*sim_AC_wguide.n_msh_el, dtype=np.complex128)
#         v_Ez6p = np.zeros(6*sim_AC_wguide.n_msh_el, dtype=np.complex128)
#         v_Ex6p_E = np.zeros(6*sim_AC_wguide.n_msh_el, dtype=np.complex128)
#         v_Ey6p_E = np.zeros(6*sim_AC_wguide.n_msh_el, dtype=np.complex128)
#         v_Ez6p_E = np.zeros(6*sim_AC_wguide.n_msh_el, dtype=np.complex128)
#         v_triang6p = []

#         i = 0
#         for i_el in np.arange(sim_AC_wguide.n_msh_el):
#             for i_node in np.arange(6):
#                 # index for the coordinates
#                 i_ex = table_nod[i_el, i_node]-1
#                 # values
#                 v_x6p[i] = x_arr[i_ex, 0]
#                 v_y6p[i] = x_arr[i_ex, 1]
#                 v_Ex6p[i] = sim_AC_wguide.sol1[0,i_node,ival,i_el]
#                 v_Ey6p[i] = sim_AC_wguide.sol1[1,i_node,ival,i_el]
#                 v_Ez6p[i] = sim_AC_wguide.sol1[2,i_node,ival,i_el]
#                 v_Ex6p_E[i] = trimmed_EM_field[0,i_node,ival_E,i_el]
#                 v_Ey6p_E[i] = trimmed_EM_field[1,i_node,ival_E,i_el]
#                 v_Ez6p_E[i] = trimmed_EM_field[2,i_node,ival_E,i_el]
#                 i += 1

#         xy = zip(v_x6p, v_y6p)
#         grid_x, grid_y = np.mgrid[x_min:x_max:n_pts_x*1j, y_min:y_max:n_pts_y*1j]
#         m_ReEx = interpolate.griddata(xy, v_Ex6p.real, (grid_x, grid_y), method='cubic')
#         m_ReEy = interpolate.griddata(xy, v_Ey6p.real, (grid_x, grid_y), method='cubic')
#         m_ReEz = interpolate.griddata(xy, v_Ez6p.real, (grid_x, grid_y), method='cubic')
#         m_ImEx = interpolate.griddata(xy, v_Ex6p.imag, (grid_x, grid_y), method='cubic')
#         m_ImEy = interpolate.griddata(xy, v_Ey6p.imag, (grid_x, grid_y), method='cubic')
#         m_ImEz = interpolate.griddata(xy, v_Ez6p.imag, (grid_x, grid_y), method='cubic')
#         m_Ex = m_ReEx + 1j*m_ImEx
#         m_Ey = m_ReEy + 1j*m_ImEy
#         m_Ez = m_ReEz + 1j*m_ImEz
#         m_Ex = m_Ex.reshape(n_pts_x,n_pts_y)
#         m_Ey = m_Ey.reshape(n_pts_x,n_pts_y)
#         m_Ez = m_Ez.reshape(n_pts_x,n_pts_y)
#         u_mat = np.array([m_Ex, m_Ey, m_Ez])
#         m_ReEx_E = interpolate.griddata(xy, v_Ex6p_E.real, (grid_x, grid_y), method='cubic')
#         m_ReEy_E = interpolate.griddata(xy, v_Ey6p_E.real, (grid_x, grid_y), method='cubic')
#         m_ReEz_E = interpolate.griddata(xy, v_Ez6p_E.real, (grid_x, grid_y), method='cubic')
#         m_ImEx_E = interpolate.griddata(xy, v_Ex6p_E.imag, (grid_x, grid_y), method='cubic')
#         m_ImEy_E = interpolate.griddata(xy, v_Ey6p_E.imag, (grid_x, grid_y), method='cubic')
#         m_ImEz_E = interpolate.griddata(xy, v_Ez6p_E.imag, (grid_x, grid_y), method='cubic')
#         m_Ex_E = m_ReEx_E + 1j*m_ImEx_E
#         m_Ey_E = m_ReEy_E + 1j*m_ImEy_E
#         m_Ez_E = -1j*sim_EM_wguide.Eig_value[ival_E]*(m_ReEz_E + 1j*m_ImEz_E)
#         m_Ex_E = m_Ex_E.reshape(n_pts_x,n_pts_y)
#         m_Ey_E = m_Ey_E.reshape(n_pts_x,n_pts_y)
#         m_Ez_E = m_Ez_E.reshape(n_pts_x,n_pts_y)
#         E_mat = np.array([m_Ex_E, m_Ey_E, m_Ez_E])

#         dx = grid_x[-1,0] - grid_x[-2,0]
#         dy = grid_y[0,-1] - grid_y[0,-2]
#         # print "NumBAT", dx
#         # print "NumBAT", dy
#         del_x_Ex = np.gradient(m_Ex, dx, axis=0)
#         del_y_Ex = np.gradient(m_Ex, dy, axis=1)
#         del_x_Ey = np.gradient(m_Ey, dx, axis=0)
#         del_y_Ey = np.gradient(m_Ey, dy, axis=1)
#         del_x_Ez = np.gradient(m_Ez, dx, axis=0)
#         del_y_Ez = np.gradient(m_Ez, dy, axis=1)
#         del_x_Ex_star = np.gradient(np.conj(m_Ex), dx, axis=0)
#         del_y_Ex_star = np.gradient(np.conj(m_Ex), dy, axis=1)
#         del_x_Ey_star = np.gradient(np.conj(m_Ey), dx, axis=0)
#         del_y_Ey_star = np.gradient(np.conj(m_Ey), dy, axis=1)
#         del_x_Ez_star = np.gradient(np.conj(m_Ez), dx, axis=0)
#         del_y_Ez_star = np.gradient(np.conj(m_Ez), dy, axis=1)
#         del_z_Ex = 1j*q_acoustic*m_Ex
#         del_z_Ey = 1j*q_acoustic*m_Ey
#         del_z_Ez = 1j*q_acoustic*m_Ez
#         del_z_Ex_star = -1j*q_acoustic*np.conj(m_Ex)
#         del_z_Ey_star = -1j*q_acoustic*np.conj(m_Ey)
#         del_z_Ez_star = -1j*q_acoustic*np.conj(m_Ez)

#         del_mat = np.array([[del_x_Ex, del_x_Ey, del_x_Ez], [del_y_Ex, del_y_Ey, del_y_Ez], [del_z_Ex, del_z_Ey, del_z_Ez]])
#         del_mat_star = np.array([[del_x_Ex_star, del_x_Ey_star, del_x_Ez_star], [del_y_Ex_star, del_y_Ey_star, del_y_Ez_star], [del_z_Ex_star, del_z_Ey_star, del_z_Ez_star]])

#         if ival < 5:
#             u_mat_CW = CW_mat[ival]
#             dx_CW = x_coord[-1,0] - x_coord[-2,0]
#             dy_CW = y_coord[0,-1] - y_coord[0,-2]
#             del_x_CWx = np.gradient(CW_mat[ival][0], dx_CW, axis=0)
#             del_y_CWx = np.gradient(CW_mat[ival][0], dy_CW, axis=1)
#             del_x_CWy = np.gradient(CW_mat[ival][1], dx_CW, axis=0)
#             del_y_CWy = np.gradient(CW_mat[ival][1], dy_CW, axis=1)
#             del_x_CWz = np.gradient(CW_mat[ival][2], dx_CW, axis=0)
#             del_y_CWz = np.gradient(CW_mat[ival][2], dy_CW, axis=1)
#             del_z_CWx = 1j*q_acoustic*CW_mat[ival][0]
#             del_z_CWy = 1j*q_acoustic*CW_mat[ival][1]
#             del_z_CWz = 1j*q_acoustic*CW_mat[ival][2]
#             del_x_CWx_star = np.gradient(np.conj(CW_mat[ival][0]), dx_CW, axis=0)
#             del_y_CWx_star = np.gradient(np.conj(CW_mat[ival][0]), dy_CW, axis=1)
#             del_x_CWy_star = np.gradient(np.conj(CW_mat[ival][1]), dx_CW, axis=0)
#             del_y_CWy_star = np.gradient(np.conj(CW_mat[ival][1]), dy_CW, axis=1)
#             del_x_CWz_star = np.gradient(np.conj(CW_mat[ival][2]), dx_CW, axis=0)
#             del_y_CWz_star = np.gradient(np.conj(CW_mat[ival][2]), dy_CW, axis=1)
#             del_z_CWx_star = -1j*q_acoustic*np.conj(CW_mat[ival][0])
#             del_z_CWy_star = -1j*q_acoustic*np.conj(CW_mat[ival][1])
#             del_z_CWz_star = -1j*q_acoustic*np.conj(CW_mat[ival][2])
#             del_mat_CW = np.array([[del_x_CWx, del_x_CWy, del_x_CWz], [del_y_CWx, del_y_CWy, del_y_CWz], [del_z_CWx, del_z_CWy, del_z_CWz]])
#             del_mat_CW_star = np.array([[del_x_CWx_star, del_x_CWy_star, del_x_CWz_star], [del_y_CWx_star, del_y_CWy_star, del_y_CWz_star], [del_z_CWx_star, del_z_CWy_star, del_z_CWz_star]])
#         else:
#             u_mat_CW = np.zeros(len(u_mat))
#             del_mat_CW_star = np.array([[0,0,0],[0,0,0],[0,0,0]])


#         for i in range(3):
#             for k in range(3):
#                 for l in range(3):
#                     integrand_AC = np.conj(u_mat[i])*del_mat[k,l]*sim_AC_wguide.structure.c_tensor_z[i,k,l]
#                     # do a 1-D integral over every row
#                     I = np.zeros( n_pts_x )
#                     for r in range(n_pts_x):
#                         I[r] = np.trapz( np.real(integrand_AC[r,:]), dx=dy )
#                     # then an integral over the result
#                     F_AC[ival] += np.trapz( I, dx=dx )
#                     I = np.zeros( n_pts_x )
#                     for r in range(n_pts_x):
#                         I[r] = np.trapz( np.imag(integrand_AC[r,:]), dx=dy )
#                     F_AC[ival] += 1j*np.trapz( I, dx=dx )

#                 ### CW - start
#                     integrand_AC = np.conj(u_mat_CW[i])*del_mat_CW[k,l]*sim_AC_wguide.structure.c_tensor_z[i,k,l]
#                     # do a 1-D integral over every row
#                     I = np.zeros( n_pts_x_CW )
#                     for r in range(n_pts_x_CW):
#                         I[r] = np.trapz( np.real(integrand_AC[r,:]), dx=dy )
#                     # then an integral over the result
#                     F_AC_CW[ival] += np.trapz( I, dx=dx )
#                     I = np.zeros( n_pts_x_CW )
#                     for r in range(n_pts_x_CW):
#                         I[r] = np.trapz( np.imag(integrand_AC[r,:]), dx=dy )
#                     F_AC_CW[ival] += 1j*np.trapz( I, dx=dx )
#                 ### CW - end

#                     for j in range(3):
#                         integrand = del_mat[i,j]*del_mat_star[k,l]*sim_AC_wguide.structure.eta_tensor[i,j,k,l]
#                         # do a 1-D integral over every row
#                         I = np.zeros( n_pts_x )
#                         for r in range(n_pts_x):
#                             I[r] = np.trapz( np.real(integrand[r,:]), dx=dy )
#                         # then an integral over the result
#                         F[ival] += np.trapz( I, dx=dx )
#                         # # Adding imag comp
#                         I = np.zeros( n_pts_x )
#                         for r in range(n_pts_x):
#                             I[r] = np.trapz( np.imag(integrand[r,:]), dx=dy )
#                         F[ival] += 1j*np.trapz( I, dx=dx )

#                         # integrand_PE = relevant_eps_effs[0]**2 * E_mat[j]*np.conj(E_mat[i])*sim_AC_wguide.structure.p_tensor[i,j,k,l]*del_mat_star[k,l]
#                         # integrand_PE = relevant_eps_effs[0]**2 * E_mat[j]*E_mat[i]*sim_AC_wguide.structure.p_tensor[i,j,k,l]*del_mat_star[k,l]
#                         integrand_PE = relevant_eps_effs[0]**2 * E_mat[j]*E_mat[i]*sim_AC_wguide.structure.p_tensor[i,j,k,l]*del_mat[k,l]
#                         I = np.zeros( n_pts_x )
#                         for r in range(n_pts_x):
#                             I[r] = np.trapz( np.real(integrand_PE[r,:]), dx=dy )
#                         F_PE[ival_E,ival_E,ival] += np.trapz( I, dx=dx )
#                         I = np.zeros( n_pts_x )
#                         for r in range(n_pts_x):
#                             I[r] = np.trapz( np.imag(integrand_PE[r,:]), dx=dy )
#                         F_PE[ival_E,ival_E,ival] += 1j*np.trapz( I, dx=dx )

                # ### CW - start
                #         integrand = del_mat_CW[i,j]*del_mat_CW_star[k,l]*sim_AC_wguide.structure.eta_tensor[i,j,k,l]
                #         # do a 1-D integral over every row
                #         I = np.zeros( n_pts_x_CW )
                #         for r in range(n_pts_x_CW):
                #             I[r] = np.trapz( np.real(integrand[r,:]), dx=dy )
                #         # then an integral over the result
                #         F_CW[ival] += np.trapz( I, dx=dx )
                #         # # Adding imag comp
                #         I = np.zeros( n_pts_x_CW )
                #         for r in range(n_pts_x_CW):
                #             I[r] = np.trapz( np.imag(integrand[r,:]), dx=dy )
                #         F_CW[ival] += 1j*np.trapz( I, dx=dx )

                #         integrand_PE = relevant_eps_effs[0]**2 * E_mat_CW[j]*E_mat_CW[i]*sim_AC_wguide.structure.p_tensor[i,j,k,l]*del_mat_CW[k,l]
                #         # integrand_PE = relevant_eps_effs[0]**2 * E_mat_CW[j]*np.conj(E_mat_CW[i])*sim_AC_wguide.structure.p_tensor[i,j,k,l]*del_mat_CW_star[k,l]
                #         I = np.zeros( n_pts_x_CW )
                #         for r in range(n_pts_x_CW):
                #             I[r] = np.trapz( np.real(integrand_PE[r,:]), dx=dy )
                #         F_PE_CW[ival_E,ival_E,ival] += np.trapz( I, dx=dx )
                #         I = np.zeros( n_pts_x_CW )
                #         for r in range(n_pts_x_CW):
                #             I[r] = np.trapz( np.imag(integrand_PE[r,:]), dx=dy )
                #         F_PE_CW[ival_E,ival_E,ival] += 1j*np.trapz( I, dx=dx )
                # ### CW - end


    # AC_py = -2j*F_AC*sim_AC_wguide.Omega_AC
    # AC_py_CW = -2j*F_AC_CW*sim_AC_wguide.Omega_AC
    # # print "AC", AC_py
    # # print "AC_CW", AC_py_CW
    # # print sim_AC_wguide.AC_mode_overlap
    # # print (AC_py-sim_AC_wguide.AC_mode_overlap)/sim_AC_wguide.AC_mode_overlap
    # # print (AC_py-AC_py_CW)/sim_AC_wguide.AC_mode_overlap

    # alpha_py = F*sim_AC_wguide.Omega_AC**2/AC_py
    # alpha_py_CW = F_CW*sim_AC_wguide.Omega_AC**2/AC_py_CW
    # alpha_py = np.real(alpha_py)
    # alpha_py_CW = np.real(alpha_py_CW)
    # # print "alpha", alpha
    # # print "Q", np.real(sim_AC_wguide.Omega_AC/(2*alpha))
    # # print "alpha_py", alpha_py
    # # print "alpha_py_CW", alpha_py_CW
    # # print "alpha/alpha_py", alpha/alpha_py
    # # print "alpha_py/alpha_py_CW", alpha_py/alpha_py_CW

    # eps_0 = 8.854187817e-12
    # Q_PE_py = F_PE*eps_0
    # Q_PE_py_CW = F_PE_CW*eps_0


    # Calc Q_moving_boundary Eq. 41
    typ_select_in = 1
    typ_select_out = 0
    try:
        # Q_MB = NumBAT.moving_boundary(sim_EM_wguide.num_modes,
        #     sim_AC_wguide.num_modes, EM_ival1_fortran, EM_ival2_fortran,
        #     AC_ival_fortran, sim_EM_wguide.n_msh_el,
        #     sim_EM_wguide.n_msh_pts, nnodes,
        #     sim_EM_wguide.table_nod, sim_EM_wguide.type_el,
        #     sim_EM_wguide.x_arr,
        #     sim_EM_wguide.structure.nb_typ_el, typ_select_in, typ_select_out,
        #     trimmed_EM_field, sim_AC_wguide.sol1,
        #     relevant_eps_effs, Fortran_debug)

        # print sim_EM_wguide.type_el
        # print sim_AC_wguide.type_el

        # sim_AC_wguide.x_arr = 100*sim_AC_wguide.x_arr

        Q_MB = NumBAT.moving_boundary(sim_EM_wguide.num_modes,
            sim_AC_wguide.num_modes, EM_ival1_fortran, EM_ival2_fortran,
            AC_ival_fortran, sim_AC_wguide.n_msh_el,
            sim_AC_wguide.n_msh_pts, nnodes,
            sim_AC_wguide.table_nod, sim_AC_wguide.type_el,
            sim_AC_wguide.x_arr,
            sim_AC_wguide.structure.nb_typ_el_AC, typ_select_in, typ_select_out,
            trimmed_EM_field, sim_AC_wguide.sol1,
            relevant_eps_effs, Fortran_debug)
    except KeyboardInterrupt:
        print "\n\n Routine moving_boundary interrupted by keyboard.\n\n"


    # from collections import Counter

    # n_msh_el = sim_EM_wguide.n_msh_el
    # type_el = sim_EM_wguide.type_el
    # table_nod = sim_EM_wguide.table_nod
    # n_msh_pts = sim_EM_wguide.n_msh_pts
    # x_arr = sim_EM_wguide.x_arr

    # eps_0 = 8.854187817e-12
    # Q_MB = np.zeros((num_modes_EM, num_modes_EM, num_modes_AC), dtype=np.complex128)

    # ### Find nodes that are in elements of various types
    # ### and find elements that have multiple nodes of various types
    # ### ie. are not single vertices on an interface.
    # node_array = -1*np.ones(n_msh_pts)
    # interface_nodes = []
    # edge_el_list = []
    # for el in range(n_msh_el):
    #     el_type = type_el[el]
    #     for i in range(6):
    #         node = table_nod[i][el]
    #         # Check if first time seen this node
    #         if node_array[node - 1] == -1: # adjust to python indexing
    #             node_array[node - 1] = el_type
    #         else:
    #             if node_array[node - 1] != el_type:
    #                 interface_nodes.append(node)
    #                 ## line below is redundant because elements sorted by type w type 1 first
    #                 # if el_type is not bkg_el_type:
    #                 edge_el_list.append(el)
    # interface_nodes = list(set(interface_nodes))
    # edge_els_multi_nodes = [k for (k,v) in Counter(edge_el_list).iteritems() if v > 1]

    # eps_list = [sim_EM_wguide.structure.bkg_material.n(sim_EM_wguide.wl_m),
    #             sim_EM_wguide.structure.inc_a_material.n(sim_EM_wguide.wl_m),
    #             sim_EM_wguide.structure.inc_b_material.n(sim_EM_wguide.wl_m),
    #             sim_EM_wguide.structure.slab_a_material.n(sim_EM_wguide.wl_m),
    #             sim_EM_wguide.structure.slab_a_bkg_material.n(sim_EM_wguide.wl_m),
    #             sim_EM_wguide.structure.slab_b_material.n(sim_EM_wguide.wl_m),
    #             sim_EM_wguide.structure.slab_b_bkg_material.n(sim_EM_wguide.wl_m),
    #             sim_EM_wguide.structure.coating_material.n(sim_EM_wguide.wl_m)]
    # # Line below may be overkill?
    # if sim_EM_wguide.structure.loss is False:
    #     eps_list = np.real(eps_list)

    # test_orient = [1,-1]
    # test1 = [0,0]
    # test2 = [0,0]
    # test3 = [0,0]
    # # x_list = []
    # # y_list = []
    # # x2_list = []
    # # y2_list = []
    # # nx_list = []
    # # ny_list = []
    # for el in edge_els_multi_nodes:
    #     AC_el = sim_AC_wguide.el_convert_tbl_inv[el]
    #     # Below also works
    #     # AC_el = (key for key,tst in sim_AC_wguide.el_convert_tbl.items() if tst==el).next()
    #     # print el
    #     # print AC_el
    #     # print sim_AC_wguide.el_convert_tbl[AC_el]

    #     contour_dr = []
    #     contour_vals = []
    #     first_seg = True
    #     integrand = np.zeros((num_modes_EM, num_modes_EM, num_modes_AC), dtype=np.complex128)
    #     # These are all possible edge line segments.
    #     for [n0,n1] in [[0,3],[3,1],[1,4],[4,2],[2,5],[5,0]]:
    #         node0 = table_nod[n0][el]
    #         node1 = table_nod[n1][el]
    #         if node0 in interface_nodes and node1 in interface_nodes:
    #             # coordinates of line seg. nodes
    #             x1 = x_arr[0,table_nod[n0][el] - 1]
    #             y1 = x_arr[1,table_nod[n0][el] - 1]
    #             x2 = x_arr[0,table_nod[n1][el] - 1]
    #             y2 = x_arr[1,table_nod[n1][el] - 1]
    #             # coordinates of non-vertex nodes, used to test orientation
    #             nodes_all = np.array([0, 1, 2, 3, 4, 5])
    #             nodes_rel = np.array([n0, n1])
    #             nodes_test = np.setdiff1d(nodes_all,nodes_rel)
    #             xt1 = x_arr[0,table_nod[nodes_test[0]][el] - 1]
    #             yt1 = x_arr[1,table_nod[nodes_test[0]][el] - 1]
    #             t1 = np.array([xt1,yt1])
    #             xt2 = x_arr[0,table_nod[nodes_test[1]][el] - 1]
    #             yt2 = x_arr[1,table_nod[nodes_test[1]][el] - 1]
    #             t2 = np.array([xt2,yt2])
    #             xt3 = x_arr[0,table_nod[nodes_test[2]][el] - 1]
    #             yt3 = x_arr[1,table_nod[nodes_test[2]][el] - 1]
    #             t3 = np.array([xt3,yt3])
    #             for i in [0, 1]:
    #                 t = test_orient[i]
    #                 normal_vec = [t*(y2-y1), -1*t*(x2-x1)]
    #                 # start half way along line seg. out along n vector
    #                 test_point = np.array([x1+(x2-x1+normal_vec[0])/2.,
    #                               y1+(y2-y1+normal_vec[1])/2.])
    #                 test1[i] = np.linalg.norm(test_point-t1)
    #                 test2[i] = np.linalg.norm(test_point-t2)
    #                 test3[i] = np.linalg.norm(test_point-t3)
    #             orient = 0
    #             if test1[0] < test1[1]:
    #                 orient += 1
    #             else:
    #                 orient -= 1
    #             if test2[0] < test2[1]:
    #                 orient += 1
    #             else:
    #                 orient -= 1
    #             if test3[0] < test3[1]:
    #                 orient += 1
    #             else:
    #                 orient -= 1
    #             if orient > 0:
    #                 normal_vec = [-1*(y2-y1), (x2-x1)]
    #             else:
    #                 normal_vec = [(y2-y1), -1*(x2-x1)]
    #             n_vec_norm = normal_vec/np.linalg.norm(normal_vec)
    #             # x_list.append(x1)
    #             # y_list.append(y1)
    #             # x2_list.append(x2)
    #             # y2_list.append(y2)
    #             # nx_list.append(n_vec_norm[0])
    #             # ny_list.append(n_vec_norm[1])

    #             # Find el on other side of interface and its epsilon
    #             all_el_w_node0 = np.where(table_nod[:] == node0)
    #             all_el_w_node1 = np.where(table_nod[:] == node1)
    #             all_el_w_node0 = [list(set(all_el_w_node0[0])), list(set(all_el_w_node0[1]))]
    #             all_el_w_node1 = [list(set(all_el_w_node1[0])), list(set(all_el_w_node1[1]))]
    #             all_el_w_node0 = [item for sublist in all_el_w_node0 for item in sublist]
    #             all_el_w_node1 = [item for sublist in all_el_w_node1 for item in sublist]
    #             all_el_w_node0 = [list(set(all_el_w_node0))][0]
    #             all_el_w_node1 = [list(set(all_el_w_node1))][0]
    #             all_el_w_nodes = all_el_w_node0 + all_el_w_node1
    #             all_el_w_node0_and_node1 = [k for (k,v) in Counter(all_el_w_nodes).iteritems() if v > 1]
    #             all_el_w_node0_and_node1.remove(el)
    #             out_side_el = all_el_w_node0_and_node1[0]
    #             # Now finally find actual epsilons
    #             type_el_a = type_el[el]
    #             type_el_b = type_el[out_side_el]
    #             eps_a = eps_list[type_el_a-1] # adjust for fortran indexing
    #             eps_b = eps_list[type_el_b-1] # adjust for fortran indexing

    #             ### Calc integrand on line segment
    #             dr = np.array([x2-x1, y2-y1])
    #             abs_dr = np.sqrt(dr[0]**2 + dr[1]**2)
    #             contour_dr.append(abs_dr)
    #             second_node = False
    #             for n in [n0, n1]:
    #                 if first_seg is True or second_node is True:
    #                     EM_ival1 = EM_ival2 = 0
    #                     for AC_ival in range(num_modes_AC):
    #                         # E-fields
    #                         e1_x = sim_EM_wguide.sol1[0,n,EM_ival1,el]
    #                         e1_y = sim_EM_wguide.sol1[1,n,EM_ival1,el]
    #                         e1_z = sim_EM_wguide.sol1[2,n,EM_ival1,el]
    #                         e2_x = sim_EM_wguide.sol1[0,n,EM_ival2,el]
    #                         e2_y = sim_EM_wguide.sol1[1,n,EM_ival2,el]
    #                         e2_z = sim_EM_wguide.sol1[2,n,EM_ival2,el]
    #                         d1_x = e1_x*eps_a*eps_0
    #                         d1_y = e1_y*eps_a*eps_0
    #                         d1_z = e1_z*eps_a*eps_0
    #                         d2_x = e2_x*eps_a*eps_0
    #                         d2_y = e2_y*eps_a*eps_0
    #                         d2_z = e2_z*eps_a*eps_0
    #                         # Displacement fields
    #                         u_x = sim_AC_wguide.sol1[0,n,AC_ival,AC_el]
    #                         u_y = sim_AC_wguide.sol1[1,n,AC_ival,AC_el]
    #                         u_z = sim_AC_wguide.sol1[2,n,AC_ival,AC_el]

    #                         u_n = np.conj(u_x)*n_vec_norm[0] + np.conj(u_y)*n_vec_norm[1]
    #                         # n_vec_norm[2] = 0 # z-comp!
    #                         n_cross_e1 = np.array([[n_vec_norm[1]*e1_z],
    #                             [-n_vec_norm[0]*e1_z],
    #                             [n_vec_norm[0]*e1_y - n_vec_norm[1]*e1_x]])
    #                         n_cross_e2 = np.array([[n_vec_norm[1]*e2_z],
    #                             [-n_vec_norm[0]*e2_z],
    #                             [n_vec_norm[0]*e2_y - n_vec_norm[1]*e2_x]])
    #                         inter_term1 = (eps_a - eps_b)*eps_0*np.vdot(n_cross_e1,n_cross_e2)
    #                         # print inter_term1
    #                         n_dot_d1 = n_vec_norm[0]*d1_x + n_vec_norm[1]*d1_y
    #                         n_dot_d2 = n_vec_norm[0]*d2_x + n_vec_norm[1]*d2_y
    #                         inter_term2 = (1./eps_b - 1./eps_a)*(1./eps_0)*np.conj(n_dot_d1)*n_dot_d2
    #                         # print inter_term2
    #                         integrand[EM_ival1,EM_ival2,AC_ival] = u_n*(inter_term1 - inter_term2)
    #                     contour_vals.append(integrand)
    #                 second_node = True
    #             first_seg = False
    #     contour_r = [0,contour_dr[0],np.sum(contour_dr)]
    #     for AC_ival in range(num_modes_AC):
    #         reshape_c_vals = []
    #         for i in range(3):
    #             reshape_c_vals.append(contour_vals[i][EM_ival1,EM_ival2,AC_ival])
    #         Q_MB[EM_ival1,EM_ival2,AC_ival] += np.trapz(reshape_c_vals, x=contour_r)

    # print Q_MB[0,0,:]
    # print Q_PE[0,0,:]
    Q = Q_MB
    # Q = Q_PE
    # Q = Q_PE + Q_MB

    # Note: sim_EM_wguide.omega_EM is the optical angular freq in units of Hz
    # Note: sim_AC_wguide.Omega_AC is the acoustic angular freq in units of Hz
    gain = 2*sim_EM_wguide.omega_EM*sim_AC_wguide.Omega_AC*np.real(Q*np.conj(Q))
    # gain_py = 2*sim_EM_wguide.omega_EM*sim_AC_wguide.Omega_AC*np.real(Q_PE_py*np.conj(Q_PE_py))
    # gain_CW = 2*sim_EM_wguide.omega_EM*sim_AC_wguide.Omega_AC*np.real(Q_PE_py_CW*np.conj(Q_PE_py_CW))
    normal_fact = np.zeros((num_modes_EM, num_modes_EM, num_modes_AC), dtype=complex)
    # normal_fact_py = np.zeros((num_modes_EM, num_modes_EM, num_modes_AC), dtype=complex)
    # normal_fact_CW = np.zeros((num_modes_EM, num_modes_EM, num_modes_AC), dtype=complex)
    # P1_CW = 1.954501164316765E-14
    for i in range(num_modes_EM):
        P1 = sim_EM_wguide.EM_mode_overlap[i]
        for j in range(num_modes_EM):
            P2 = sim_EM_wguide.EM_mode_overlap[j]
            for k in range(num_modes_AC):
                P3 = sim_AC_wguide.AC_mode_overlap[k]
                # P3_py = AC_py[k]
                # P3_CW = AC_py_CW[k]
                normal_fact[i, j, k] = P1*P2*P3
                # normal_fact_py[i, j, k] = P1*P2*P3_py
                # normal_fact_CW[i, j, k] = P1_CW*P1_CW*P3_CW
    SBS_gain = np.real(gain/normal_fact)
    # SBS_gain_py = np.real(gain_py/normal_fact_py)
    # SBS_gain_CW = np.real(gain_CW/normal_fact_CW)

    print "MB_gain_ratio", SBS_gain[0,0,2]/alpha[2]/5164.72
    print "MB_gain_ratio", SBS_gain[0,0,4]/alpha[4]/1344.00
    # print "MB_gain_ratio", SBS_gain[0,0,8]/alpha[8]/400.31

    # print "PE_gain_ratio", SBS_gain[0,0,2]/alpha[2]/1152.95
    # print "PE_gain_ratio", SBS_gain[0,0,4]/alpha[4]/6333.18
    # # print "PE_gain_ratio", SBS_gain[0,0,8]/alpha[8]/36.55
    
    # print "SBS_gain_py", SBS_gain_py[0,0,:]/alpha_py
    # print "SBS_gain_CW", SBS_gain_CW[0,0,:]/alpha_py_CW
    # print "gain ratio py", SBS_gain_py[0,0,:]/SBS_gain[0,0,:]
    # print "gain ratio CW", SBS_gain_CW[0,0,:]/SBS_gain[0,0,:]



    # nx_list = np.array(nx_list)
    # ny_list = np.array(ny_list)
    # plt.clf()
    # plt.figure(figsize=(13,13))
    # ax = plt.subplot(1,1,1)
    # # for node in range(np.shape(x_arr)[1]):
    # plt.plot(x_list, y_list, 'o')
    # # plt.plot(x2_list, y2_list, 'o')
    # plt.plot(x_list+nx_list*15e-9, y_list+ny_list*15e-9, 'xr')
    # ax.set_aspect('equal')
    # # ax.set_xlim(0,1000e-9)
    # # ax.set_ylim(-1000e-9,0)
    # plt.savefig('msh_%(add)s.png' %
    #     {'add' : el}, bbox_inches='tight')
    # plt.close()

    return SBS_gain, Q_PE, Q_MB, alpha


# # # Test overlap
# # nel = sim_EM_wguide.n_msh_el
# # type_el = sim_EM_wguide.type_el
# # table_nod = sim_EM_wguide.table_nod
# # x_arr = sim_EM_wguide.x_arr
# # nnodes = 6
# # xel = np.zeros((2,nnodes))
# # nod_el_p = np.zeros(nnodes)
# # xx = [0,0]
# # nquad = 16
# # [wq, xq, yq] = integration.quad_triangle(nquad)


#### Categorise modes by their symmetries #############################################
def symmetries(sim_wguide, n_points=10):
    """ Plot EM mode fields.

        Args:
            sim_wguide : A :Struct: instance that has had calc_modes calculated

        Keyword Args:
            n_points  (int): The number of points across unitcell to \
                interpolate the field onto.
    """

    mode_fields = sim_wguide.sol1

    # field mapping
    x_tmp = []
    y_tmp = []
    for i in np.arange(sim_wguide.n_msh_pts):
        x_tmp.append(sim_wguide.x_arr[0,i])
        y_tmp.append(sim_wguide.x_arr[1,i])
    x_min = np.min(x_tmp); x_max=np.max(x_tmp)
    y_min = np.min(y_tmp); y_max=np.max(y_tmp)
    area = abs((x_max-x_min)*(y_max-y_min))
    n_pts_x = int(n_points*abs(x_max-x_min)/np.sqrt(area))
    n_pts_y = int(n_points*abs(y_max-y_min)/np.sqrt(area))
    v_x=np.zeros(n_pts_x*n_pts_y)
    v_y=np.zeros(n_pts_x*n_pts_y)
    i=0
    for x in np.linspace(x_min,x_max,n_pts_x):
        for y in np.linspace(y_min,y_max,n_pts_y):
            v_x[i] = x
            v_y[i] = y
            i+=1
    v_x = np.array(v_x)
    v_y = np.array(v_y)

    # unrolling data for the interpolators
    table_nod = sim_wguide.table_nod.T
    x_arr = sim_wguide.x_arr.T

    sym_list = []

    for ival in range(len(sim_wguide.Eig_value)):
        # dense triangulation with multiple points
        v_x6p = np.zeros(6*sim_wguide.n_msh_el)
        v_y6p = np.zeros(6*sim_wguide.n_msh_el)
        v_Ex6p = np.zeros(6*sim_wguide.n_msh_el, dtype=np.complex128)
        v_Ey6p = np.zeros(6*sim_wguide.n_msh_el, dtype=np.complex128)
        # v_Ez6p = np.zeros(6*sim_wguide.n_msh_el, dtype=np.complex128)
        v_triang6p = []

        i = 0
        for i_el in np.arange(sim_wguide.n_msh_el):

            # triangles
            idx = np.arange(6*i_el, 6*(i_el+1))
            triangles = [[idx[0], idx[3], idx[5]],
                         [idx[1], idx[4], idx[3]],
                         [idx[2], idx[5], idx[4]],
                         [idx[3], idx[4], idx[5]]]
            v_triang6p.extend(triangles)

            for i_node in np.arange(6):
                # index for the coordinates
                i_ex = table_nod[i_el, i_node]-1
                # values
                v_x6p[i] = x_arr[i_ex, 0]
                v_y6p[i] = x_arr[i_ex, 1]
                v_Ex6p[i] = mode_fields[0,i_node,ival,i_el]
                v_Ey6p[i] = mode_fields[1,i_node,ival,i_el]
                # v_Ez6p[i] = mode_fields[2,i_node,ival,i_el]
                i += 1

        # dense triangulation with unique points
        v_triang1p = []
        for i_el in np.arange(sim_wguide.n_msh_el):
            # triangles
            triangles = [[table_nod[i_el,0]-1,table_nod[i_el,3]-1,table_nod[i_el,5]-1],
                         [table_nod[i_el,1]-1,table_nod[i_el,4]-1,table_nod[i_el,3]-1],
                         [table_nod[i_el,2]-1,table_nod[i_el,5]-1,table_nod[i_el,4]-1],
                         [table_nod[i_el,3]-1,table_nod[i_el,4]-1,table_nod[i_el,5]-1]]
            v_triang1p.extend(triangles)

        # triangulations
        triang6p = matplotlib.tri.Triangulation(v_x6p,v_y6p,v_triang6p)
        triang1p = matplotlib.tri.Triangulation(x_arr[:,0],x_arr[:,1],v_triang1p)

        # building interpolators: triang1p for the finder, triang6p for the values
        finder = matplotlib.tri.TrapezoidMapTriFinder(triang1p)
        ReEx = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ex6p.real,trifinder=finder)
        ImEx = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ex6p.imag,trifinder=finder)
        ReEy = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ey6p.real,trifinder=finder)
        ImEy = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ey6p.imag,trifinder=finder)
        # ReEz = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ez6p.real,trifinder=finder)
        # ImEz = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ez6p.imag,trifinder=finder)

        ### plotting
        # interpolated fields
        m_ReEx = ReEx(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ReEy = ReEy(v_x,v_y).reshape(n_pts_x,n_pts_y)
        # m_ReEz = ReEz(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ImEx = ImEx(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ImEy = ImEy(v_x,v_y).reshape(n_pts_x,n_pts_y)
        # m_ImEz = ImEz(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_Ex = m_ReEx + 1j*m_ImEx
        m_Ey = m_ReEy + 1j*m_ImEy
        # m_Ez = m_ReEz + 1j*m_ImEz

        m_Ex_ymirror = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        m_Ex_xmirror = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        m_Ex_rotated = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        m_Ey_ymirror = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        m_Ey_xmirror = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        m_Ey_rotated = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        # m_Ez_ymirror = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        # m_Ez_xmirror = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        # m_Ez_rotated = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        Ex_sigma_y = 0
        Ey_sigma_y = 0
        # Ez_sigma_y = 0
        Ex_sigma_x = 0
        Ey_sigma_x = 0
        # Ez_sigma_x = 0
        Ex_C_2 = 0
        Ey_C_2 = 0
        # Ez_C_2 = 0
        # max_E = max(np.max(np.abs(m_Ex)), np.max(np.abs(m_Ey)), np.max(np.abs(m_Ez)))

        for ix in range(n_pts_x):
            for iy in range(n_pts_y):
                m_Ex_ymirror[ix,iy] = (m_Ex[ix,n_pts_y-iy-1])
                m_Ey_ymirror[ix,iy] = -1*(m_Ey[ix,n_pts_y-iy-1])
                # m_Ez_ymirror[ix,iy] = m_Ez[ix,n_pts_y-iy-1]
                m_Ex_xmirror[ix,iy] = -1*(m_Ex[n_pts_x-ix-1,iy])
                m_Ey_xmirror[ix,iy] = (m_Ey[n_pts_x-ix-1,iy])
                # m_Ez_xmirror[ix,iy] = m_Ez[n_pts_x-ix-1,iy]
                m_Ex_rotated[ix,iy] = -1*(m_Ex[n_pts_x-ix-1,n_pts_y-iy-1])
                m_Ey_rotated[ix,iy] = -1*(m_Ey[n_pts_x-ix-1,n_pts_y-iy-1])
                # m_Ez_rotated[ix,iy] = m_Ez[n_pts_x-ix-1,iy]
                # Ex_sigma_y += np.abs(m_Ex[ix,iy] - m_Ex_ymirror[ix,iy])/max_E
                # Ey_sigma_y += np.abs(m_Ey[ix,iy] - m_Ey_ymirror[ix,iy])/max_E
                # Ez_sigma_y += np.abs(m_Ez[ix,iy] - m_Ez_ymirror[ix,iy])/max_E
                # Ex_sigma_x += np.abs(m_Ex[ix,iy] - m_Ex_xmirror[ix,iy])/max_E
                # Ey_sigma_x += np.abs(m_Ey[ix,iy] - m_Ey_xmirror[ix,iy])/max_E
                # Ez_sigma_x += np.abs(m_Ez[ix,iy] - m_Ez_xmirror[ix,iy])/max_E
                # Ex_C_2 += np.abs(m_Ex[ix,iy] - m_Ex_rotated[ix,iy])/max_E
                # Ey_C_2 += np.abs(m_Ey[ix,iy] - m_Ey_rotated[ix,iy])/max_E
                # Ez_C_2 += np.abs(m_Ez[ix,iy] - m_Ez_rotated[ix,iy])/max_E

        Ex_sigma_y = np.sum(np.abs(m_Ex - m_Ex_ymirror))
        Ey_sigma_y = np.sum(np.abs(m_Ey - m_Ey_ymirror))
        # Ez_sigma_y = np.sum(np.abs(m_Ez - m_Ez_ymirror))
        Ex_sigma_x = np.sum(np.abs(m_Ex - m_Ex_xmirror))
        Ey_sigma_x = np.sum(np.abs(m_Ey - m_Ey_xmirror))
        # Ez_sigma_x = np.sum(np.abs(m_Ez - m_Ez_xmirror))
        Ex_C_2 = np.sum(np.abs(m_Ex - m_Ex_rotated))
        Ey_C_2 = np.sum(np.abs(m_Ey - m_Ey_rotated))
        # Ez_C_2 = np.sum(np.abs(m_Ez - m_Ez_rotated))
        # sigma_y = (Ex_sigma_y + Ey_sigma_y + Ez_sigma_y)/(n_pts_x*n_pts_y)
        # sigma_x = (Ex_sigma_x + Ey_sigma_x + Ez_sigma_x)/(n_pts_x*n_pts_y)
        # C_2 = (Ex_C_2 + Ey_C_2 + Ez_C_2)/(n_pts_x*n_pts_y)
        sigma_y = (Ex_sigma_y + Ey_sigma_y)/(n_pts_x*n_pts_y)
        sigma_x = (Ex_sigma_x + Ey_sigma_x)/(n_pts_x*n_pts_y)
        C_2 = (Ex_C_2 + Ey_C_2)/(n_pts_x*n_pts_y)

        if abs(C_2) > 0.2:
            C_2_print = -1
        else:
            C_2_print = 1
        if abs(sigma_y) > 0.1:
            sigma_y_print = -1
        else:
            sigma_y_print = 1
        if abs(sigma_x) > 0.1:
            sigma_x_print = -1
        else:
            sigma_x_print = 1
        # print '------'
        # print ival
        # print 'C_2', C_2_print
        # print 'C_2', C_2
        # print 'sigma_x', sigma_x_print
        # print 'sigma_x', sigma_x
        # print 'sigma_y', sigma_y_print
        # print 'sigma_y', sigma_y
        sym_list.append([C_2_print, sigma_y_print, sigma_x_print])

        # v_plots = [np.real(m_Ex_ymirror),np.real(m_Ey_ymirror),np.real(m_Ez_ymirror),
        #     np.imag(m_Ex_ymirror),np.imag(m_Ey_ymirror),np.imag(m_Ez_ymirror)]
        # # field plots
        # plt.clf()
        # plt.figure(figsize=(13,13))
        # for i_p,plot in enumerate(v_plots):
        #     ax = plt.subplot(3,3,i_p+1)
        #     im = plt.imshow(plot.T,cmap='inferno');
        #     # colorbar
        #     divider = make_axes_locatable(ax)
        #     cax = divider.append_axes("right", size="5%", pad=0.1)
        #     cbar = plt.colorbar(im, cax=cax)
        # plt.savefig('fields/field_%(i)i-ymirror.pdf' %
        #         {'i' : ival}, bbox_inches='tight')
        # plt.close()
        # v_plots = [np.real(m_Ex_xmirror),np.real(m_Ey_xmirror),np.real(m_Ez_xmirror),
        #     np.imag(m_Ex_xmirror),np.imag(m_Ey_xmirror),np.imag(m_Ez_xmirror)]
        # # field plots
        # plt.clf()
        # plt.figure(figsize=(13,13))
        # for i_p,plot in enumerate(v_plots):
        #     ax = plt.subplot(3,3,i_p+1)
        #     im = plt.imshow(plot.T,cmap='inferno');
        #     # colorbar
        #     divider = make_axes_locatable(ax)
        #     cax = divider.append_axes("right", size="5%", pad=0.1)
        #     cbar = plt.colorbar(im, cax=cax)
        # plt.savefig('fields/field_%(i)i-xmirror.pdf' %
        #         {'i' : ival}, bbox_inches='tight')
        # plt.close()
        # v_plots = [np.real(m_Ex_rotated),np.real(m_Ey_rotated),np.real(m_Ez_rotated),
        #     np.imag(m_Ex_rotated),np.imag(m_Ey_rotated),np.imag(m_Ez_rotated)]
        # # field plots
        # plt.clf()
        # plt.figure(figsize=(13,13))
        # for i_p,plot in enumerate(v_plots):
        #     ax = plt.subplot(3,3,i_p+1)
        #     im = plt.imshow(plot.T,cmap='inferno');
        #     # colorbar
        #     divider = make_axes_locatable(ax)
        #     cax = divider.append_axes("right", size="5%", pad=0.1)
        #     cbar = plt.colorbar(im, cax=cax)
        # plt.savefig('fields/field_%(i)i-rotated.pdf' %
        #         {'i' : ival}, bbox_inches='tight')
        # plt.close()
        # v_plots = [m_ReEx,m_ReEy,m_ReEz,
        #     m_ImEx,m_ImEy,m_ImEz]
        # # field plots
        # plt.clf()
        # plt.figure(figsize=(13,13))
        # for i_p,plot in enumerate(v_plots):
        #     ax = plt.subplot(3,3,i_p+1)
        #     im = plt.imshow(plot.T,cmap='inferno');
        #     # colorbar
        #     divider = make_axes_locatable(ax)
        #     cax = divider.append_axes("right", size="5%", pad=0.1)
        #     cbar = plt.colorbar(im, cax=cax)
        # plt.savefig('fields/field_%(i)i.pdf' %
        #         {'i' : ival}, bbox_inches='tight')
        # plt.close()

    return sym_list




def quad_triangle(nquad):
    """ Implementation of quad_triangle

                   Integration numerique
                   ---------------------
    Evalue les integrales elementaires des composantes convectives
    sur chaque triangle. on utilise ici la methode de hammer a
    seize points de gauss qui integre exactement des polynomes du
    huitieme degre.

    Google Translate:
    Evaluates integrals elementary convective components
    on each triangle. here we use the method of a hammer
    sixteen points of gauss that integrates exactly the polynomials
    eighth degree.

    Reference
    J. N. Lyness and D. Jespersen
    "Moderate Degree Symmetric Quadrature Rules for the Triangle"
    J. Inst. Math. Appl., 1975, 15(1), pp. 19-32
    "J. Inst. Math. Appl." is now Continued as "IMA J. Appl. Math."
    J. Inst. Math. Appl. = Journal of the Institute of Mathematics and its Applications
    IMA J. Appl. Math.   = IMA Journal of Applied Mathematics
    """

    if nquad is not 16:
        raise ValueError, 'integration.quad_triangle nquad must == 16.'
    wq = np.zeros(nquad)
    xq = np.zeros(nquad)
    yq = np.zeros(nquad)

    poidbar = 1.443156076777862e-1 / 2.0
    poid1   = 2.852749028018549e-1 / 6.0
    poid2   = 9.737549286959440e-2 / 6.0
    poid3   = 3.096521116041552e-1 / 6.0
    poid4   = 1.633818850466092e-1 / 12.0

    coorbar = 1.0 / 3.0

    coor1grp1 = 4.592925882927229e-1
    coor2grp1 = 8.141482341455413e-2
    coor1grp2 = 5.054722831703103e-2
    coor2grp2 = 8.989055433659379e-1
    coor1grp3 = 1.705693077517601e-1
    coor2grp3 = 6.588613844964797e-1
    coor1grp4 = 7.284923929554041e-1
    coor2grp4 = 2.63112829634638689e-1
    coor3grp4 = 8.394777409957211e-3

    i = 0
    xq[i] = coorbar
    yq[i] = coorbar
    wq[i] = poidbar
    i = 1
    xq[i] = coor1grp1
    yq[i] = coor1grp1
    wq[i] = poid1
    i = 2
    xq[i] = coor1grp1
    yq[i] = coor2grp1
    wq[i] = poid1
    i = 3
    xq[i] = coor2grp1
    yq[i] = coor1grp1
    wq[i] = poid1
    i = 4
    xq[i] = coor1grp2
    yq[i] = coor1grp2
    wq[i] = poid2
    i = 5
    xq[i] = coor1grp2
    yq[i] = coor2grp2
    wq[i] = poid2
    i = 6
    xq[i] = coor2grp2
    yq[i] = coor1grp2
    wq[i] = poid2
    i = 7
    xq[i] = coor1grp3
    yq[i] = coor1grp3
    wq[i] = poid3
    i = 8
    xq[i] = coor1grp3
    yq[i] = coor2grp3
    wq[i] = poid3
    i = 9
    xq[i] = coor2grp3
    yq[i] = coor1grp3
    wq[i] = poid3
    i = 10
    xq[i] = coor1grp4
    yq[i] = coor2grp4
    wq[i] = poid4
    i = 11
    xq[i] = coor2grp4
    yq[i] = coor1grp4
    wq[i] = poid4
    i = 12
    xq[i] = coor2grp4
    yq[i] = coor3grp4
    wq[i] = poid4
    i = 13
    xq[i] = coor3grp4
    yq[i] = coor2grp4
    wq[i] = poid4
    i = 14
    xq[i] = coor3grp4
    yq[i] = coor1grp4
    wq[i] = poid4
    i = 15
    xq[i] = coor1grp4
    yq[i] = coor3grp4
    wq[i] = poid4

    return [wq, xq, yq]
