"""
    mode_calcs.py is a subroutine of NumBAT that contains methods to
    calculate the EM and Acoustic modes of a structure.

    Copyright (C) 2016  Bjorn Sturmberg, Kokou Dossou
"""

import numpy as np
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable

# import materials
# import objects
# import mode_calcs
import plotting
from fortran import NumBAT


def gain_and_qs(sim_EM_wguide, sim_AC_wguide, q_acoustic,
                EM_ival1=0, EM_ival2=0, AC_ival=0):
    """ Calculate interaction integrals and SBS gain.

    Calc Qs of a range of selected modes
    Pass in the q_acoustic as this is beta of AC modes

    By default considers the interactions between all modes,
    can also specify specific modes.

    """

### Notes about internals of fortran integration
# Calc overlap of basis functions (and PE tensor and epsilon)
# Then use this multiple times for calc of each mode field values

# phi is values of Lagrange polynomials (1-6) at that node.
# grad is value of gradient of Lagrange polynomials (1-6) at that node.


    if EM_ival1 == 'All':
        EM_ival1_fortran = -1
        # P1 = sim_EM_wguide.EM_mode_overlap
    else:
        EM_ival1_fortran = EM_ival1+1  # convert back to fortran indexing
        # P1 = np.zeros(num_modes_EM, dtype=complex)
        # P1[EM_ival1] = sim_EM_wguide.EM_mode_overlap[EM_ival1]

    if EM_ival2 == 'All':
        EM_ival2_fortran = -1
        # P2 = sim_EM_wguide.EM_mode_overlap
    else:
        EM_ival2_fortran = EM_ival2+1  # convert back to fortran indexing
        # P2 = np.zeros(num_modes_EM, dtype=complex)
        # P2[EM_ival2] = sim_EM_wguide.EM_mode_overlap[EM_ival2]

    if AC_ival == 'All':
        AC_ival_fortran = -1
        # Note: sim_AC_wguide.Omega_AC if the acoustic angular freq in units of Hz
        AC_freq_Omega = sim_AC_wguide.Omega_AC
        # P3 = sim_AC_wguide.AC_mode_overlap
    else:
        AC_freq_Omega = sim_AC_wguide.Omega_AC[AC_ival]
        AC_ival_fortran = AC_ival+1  # convert back to fortran indexing
        # P3 = np.zeros(num_modes_AC, dtype=complex)
        # P3[AC_ival] = sim_AC_wguide.AC_mode_overlap[AC_ival]

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

    # sim_AC_wguide.AC_mode_overlap[0] = sim_AC_wguide.AC_mode_overlap[0]*2.0541115841
    # sim_AC_wguide.AC_mode_overlap[1] = sim_AC_wguide.AC_mode_overlap[1]*1.62055717426
    # sim_AC_wguide.AC_mode_overlap[2] = sim_AC_wguide.AC_mode_overlap[2]*1.97449348579
    # sim_AC_wguide.AC_mode_overlap[3] = sim_AC_wguide.AC_mode_overlap[3]*1.7819220338
    # sim_AC_wguide.AC_mode_overlap[4] = sim_AC_wguide.AC_mode_overlap[4]*1.05417483114
    # sim_AC_wguide.AC_mode_overlap[5] = sim_AC_wguide.AC_mode_overlap[5]*1.59968666271
    # sim_AC_wguide.AC_mode_overlap[6] = sim_AC_wguide.AC_mode_overlap[6]*1.06616883215
    # sim_AC_wguide.AC_mode_overlap[7] = sim_AC_wguide.AC_mode_overlap[7]*0.917209648528
    # sim_AC_wguide.AC_mode_overlap[8] = sim_AC_wguide.AC_mode_overlap[8]*1.01503248635

### Calc alpha (loss) Eq. 45
    try:
        if sim_EM_wguide.structure.inc_shape == 'rectangular':
            alpha, basis_overlap_alpha = NumBAT.ac_alpha_int_v2(sim_AC_wguide.num_modes,
                sim_AC_wguide.n_msh_el, sim_AC_wguide.n_msh_pts, nnodes,
                sim_AC_wguide.table_nod, sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
                sim_AC_wguide.structure.nb_typ_el_AC, sim_AC_wguide.structure.eta_tensor,
                q_acoustic, sim_AC_wguide.Omega_AC, sim_AC_wguide.sol1,
                sim_AC_wguide.AC_mode_overlap)
        elif sim_EM_wguide.structure.inc_shape == 'circular':
            alpha, basis_overlap_alpha = NumBAT.ac_alpha_int(sim_AC_wguide.num_modes,
                sim_AC_wguide.n_msh_el, sim_AC_wguide.n_msh_pts, nnodes,
                sim_AC_wguide.table_nod, sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
                sim_AC_wguide.structure.nb_typ_el_AC, sim_AC_wguide.structure.eta_tensor,
                q_acoustic, sim_AC_wguide.Omega_AC, sim_AC_wguide.sol1,
                sim_AC_wguide.AC_mode_overlap, Fortran_debug)
    except KeyboardInterrupt:
        print "\n\n Routine ac_alpha_int interrupted by keyboard.\n\n"
    alpha = np.real(alpha)

    # alpha[0] = alpha[0]/2.0541115841
    # alpha[1] = alpha[1]/1.62055717426
    # alpha[2] = alpha[2]/1.97449348579
    # alpha[3] = alpha[3]/1.7819220338
    # alpha[4] = alpha[4]/1.05417483114
    # alpha[5] = alpha[5]/1.59968666271
    # alpha[6] = alpha[6]/1.06616883215
    # alpha[7] = alpha[7]/0.917209648528
    # alpha[8] = alpha[8]/1.01503248635


### Calc Q_photoelastic Eq. 33
    try:
        #TODO: removes basis_overlaps
        if sim_EM_wguide.structure.inc_shape == 'rectangular':
            Q_PE, basis_overlap_PE, field_overlap_PE = NumBAT.photoelastic_int_v2(
                sim_EM_wguide.num_modes, sim_AC_wguide.num_modes, EM_ival1_fortran,
                EM_ival2_fortran, AC_ival_fortran, sim_AC_wguide.n_msh_el,
                sim_AC_wguide.n_msh_pts, nnodes,
                sim_AC_wguide.table_nod, sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
                sim_AC_wguide.structure.nb_typ_el_AC, sim_AC_wguide.structure.p_tensor,
                q_acoustic, trimmed_EM_field, sim_AC_wguide.sol1,
                relevant_eps_effs, sim_EM_wguide.Eig_value, Fortran_debug)
        elif sim_EM_wguide.structure.inc_shape == 'circular':
            Q_PE, basis_overlap_PE = NumBAT.photoelastic_int(
                sim_EM_wguide.num_modes, sim_AC_wguide.num_modes, EM_ival1_fortran,
                EM_ival2_fortran, AC_ival_fortran, sim_AC_wguide.n_msh_el,
                sim_AC_wguide.n_msh_pts, nnodes,
                sim_AC_wguide.table_nod, sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
                sim_AC_wguide.structure.nb_typ_el_AC, sim_AC_wguide.structure.p_tensor,
                q_acoustic, trimmed_EM_field, sim_AC_wguide.sol1,
                relevant_eps_effs, sim_EM_wguide.Eig_value, Fortran_debug)
    except KeyboardInterrupt:
        print "\n\n Routine photoelastic_int interrupted by keyboard.\n\n"

    # n_points = 100
    # x_tmp = []
    # y_tmp = []
    # for i in np.arange(sim_AC_wguide.n_msh_pts):
    #     x_tmp.append(sim_AC_wguide.x_arr[0,i])
    #     y_tmp.append(sim_AC_wguide.x_arr[1,i])
    # x_min = np.min(x_tmp); x_max=np.max(x_tmp)
    # y_min = np.min(y_tmp); y_max=np.max(y_tmp)
    # area = abs((x_max-x_min)*(y_max-y_min))
    # n_pts_x = int(n_points*abs(x_max-x_min)/np.sqrt(area))
    # n_pts_y = int(n_points*abs(y_max-y_min)/np.sqrt(area))
    # v_x=np.zeros(n_pts_x*n_pts_y)
    # v_y=np.zeros(n_pts_x*n_pts_y)
    # i=0
    # for x in np.linspace(x_min,x_max,n_pts_x):
    #     # for y in np.linspace(y_min,y_max,n_pts_y):
    #     for y in np.linspace(y_max,y_min,n_pts_y):
    #         v_x[i] = x
    #         v_y[i] = y
    #         i+=1
    # v_x = np.array(v_x)
    # v_y = np.array(v_y)

    # # unrolling data for the interpolators
    # table_nod = sim_AC_wguide.table_nod.T
    # x_arr = sim_AC_wguide.x_arr.T

    # for ival in [0]:
    #     # dense triangulation with multiple points
    #     v_x6p = np.zeros(6*sim_AC_wguide.n_msh_el)
    #     v_y6p = np.zeros(6*sim_AC_wguide.n_msh_el)
    #     v_Ex6p = np.zeros(6*sim_AC_wguide.n_msh_el, dtype=np.complex128)
    #     v_Ey6p = np.zeros(6*sim_AC_wguide.n_msh_el, dtype=np.complex128)
    #     v_Ez6p = np.zeros(6*sim_AC_wguide.n_msh_el, dtype=np.complex128)
    #     v_triang6p = []

    #     i = 0
    #     for i_el in np.arange(sim_AC_wguide.n_msh_el):
    #         # triangles
    #         idx = np.arange(6*i_el, 6*(i_el+1))
    #         triangles = [[idx[0], idx[3], idx[5]],
    #                      [idx[1], idx[4], idx[3]],
    #                      [idx[2], idx[5], idx[4]],
    #                      [idx[3], idx[4], idx[5]]]
    #         v_triang6p.extend(triangles)

    #         for i_node in np.arange(6):
    #             # index for the coordinates
    #             i_ex = table_nod[i_el, i_node]-1
    #             # values
    #             v_x6p[i] = x_arr[i_ex, 0]
    #             v_y6p[i] = x_arr[i_ex, 1]
    #             v_Ex6p[i] = sim_EM_wguide.sol1[0,i_node,ival,i_el]
    #             v_Ey6p[i] = sim_EM_wguide.sol1[1,i_node,ival,i_el]
    #             v_Ez6p[i] = -1j*sim_EM_wguide.Eig_value[ival]*sim_EM_wguide.sol1[2,i_node,ival,i_el]
    #             i += 1

    #     v_E6p = np.sqrt(np.abs(v_Ex6p)**2 +
    #                     np.abs(v_Ey6p)**2 +
    #                     np.abs(v_Ez6p)**2)

    # eta = np.zeros(len(sim_AC_wguide.Eig_value), dtype=np.complex128)
    # for ival in range(len(sim_AC_wguide.Eig_value)):
    #     # dense triangulation with multiple points
    #     AC_v_x6p = np.zeros(6*sim_AC_wguide.n_msh_el)
    #     AC_v_y6p = np.zeros(6*sim_AC_wguide.n_msh_el)
    #     AC_v_Ex6p = np.zeros(6*sim_AC_wguide.n_msh_el, dtype=np.complex128)
    #     AC_v_Ey6p = np.zeros(6*sim_AC_wguide.n_msh_el, dtype=np.complex128)
    #     AC_v_Ez6p = np.zeros(6*sim_AC_wguide.n_msh_el, dtype=np.complex128)
    #     AC_v_triang6p = []

    #     i = 0
    #     for i_el in np.arange(sim_AC_wguide.n_msh_el):
    #         # triangles
    #         idx = np.arange(6*i_el, 6*(i_el+1))
    #         triangles = [[idx[0], idx[3], idx[5]],
    #                      [idx[1], idx[4], idx[3]],
    #                      [idx[2], idx[5], idx[4]],
    #                      [idx[3], idx[4], idx[5]]]
    #         AC_v_triang6p.extend(triangles)

    #         for i_node in np.arange(6):
    #             # index for the coordinates
    #             i_ex = table_nod[i_el, i_node]-1
    #             # values
    #             AC_v_x6p[i] = x_arr[i_ex, 0]
    #             AC_v_y6p[i] = x_arr[i_ex, 1]
    #             AC_v_Ex6p[i] = sim_AC_wguide.sol1[0,i_node,ival,i_el]
    #             AC_v_Ey6p[i] = sim_AC_wguide.sol1[1,i_node,ival,i_el]
    #             AC_v_Ez6p[i] = sim_AC_wguide.sol1[2,i_node,ival,i_el]
    #             i += 1

    #     AC_v_E6p = np.sqrt(np.abs(AC_v_Ex6p)**2 +
    #                     np.abs(AC_v_Ey6p)**2 +
    #                     np.abs(AC_v_Ez6p)**2)


    #     # Scalar calc from Poulton Josa 2013

    #     top = np.sum(v_E6p**2*AC_v_E6p)**2
    #     bot = np.sum(v_E6p**2)**2 + np.sum(AC_v_E6p)**2
    #     eta[ival] = top/bot
    #     # eta[ival] = np.sum(v_E6p**2*AC_v_E6p)**2
    #     # Power normalise
    #     speed_c = 299792458
    #     omega_p = 2.0*np.pi/(1550*1e-9)
    #     n = np.sqrt(12.25)
    #     s = 2330
    #     p_12 = 0.017
    #     V_L = 8440
    #     Q = 1000
    #     prefactor = omega_p**2*n**7*p_12**2*Q/(speed_c**3*s*V_L*sim_AC_wguide.Omega_AC)

    # print eta*(314.7e-9*0.9*314.7e-9)
    # print prefactor*eta



































    # print Q_PE[0,0,2]
    # print Q_PE2[0,0,2]
    # Q_PE=0
    Q_MB = 0.0 # Haven't implemented Moving Boundary integral (but nor did Rakich)
    Q = Q_PE + Q_MB
    # Note: sim_EM_wguide.omega_EM is the optical angular freq in units of Hz
    gain = 2*sim_EM_wguide.omega_EM*AC_freq_Omega*np.real(Q*np.conj(Q))
    normal_fact = np.zeros((num_modes_EM, num_modes_EM, num_modes_AC), dtype=complex)
    for i in range(num_modes_EM):
        P1 = sim_EM_wguide.EM_mode_overlap[i]
        for j in range(num_modes_EM):
            P2 = sim_EM_wguide.EM_mode_overlap[j]
            for k in range(num_modes_AC):
                P3 = sim_AC_wguide.AC_mode_overlap[k]
                normal_fact[i, j, k] = P1*P2*P3
    SBS_gain = np.real(gain/normal_fact)

    return SBS_gain, Q_PE, Q_MB, alpha
















# ### Calc Q_moving_boundary Eq. 41
# from collections import Counter

# eps_0 = 8.854187817e-12

# n_msh_el = sim_EM_wguide.n_msh_el
# type_el = sim_EM_wguide.type_el
# table_nod = sim_EM_wguide.table_nod
# n_msh_pts = sim_EM_wguide.n_msh_pts
# x_arr = sim_EM_wguide.x_arr

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

# Q_MB = 0
# eps_list = [sim_EM_wguide.structure.bkg_material.n(sim_EM_wguide.wl_nm), sim_EM_wguide.structure.inc_a_material.n(sim_EM_wguide.wl_nm), sim_EM_wguide.structure.inc_b_material.n(sim_EM_wguide.wl_nm), sim_EM_wguide.structure.slab_a_material.n(sim_EM_wguide.wl_nm), sim_EM_wguide.structure.slab_a_bkg_material.n(sim_EM_wguide.wl_nm), sim_EM_wguide.structure.slab_b_material.n(sim_EM_wguide.wl_nm),sim_EM_wguide.structure.slab_b_bkg_material.n(sim_EM_wguide.wl_nm), sim_EM_wguide.structure.coating_material.n(sim_EM_wguide.wl_nm)]
# # Line below may be overkill?
# if sim_EM_wguide.structure.loss is False:
#     eps_list = np.real(eps_list)

# test_orient = [1,-1]
# test1 = [0,0]
# test2 = [0,0]
# test3 = [0,0]
# # for el in edge_els_multi_nodes:
# for el in [edge_els_multi_nodes[0]]:
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
#             xt1 = x_arr[0,table_nod[3][el] - 1]
#             yt1 = x_arr[1,table_nod[3][el] - 1]
#             t1 = np.array([xt1,yt1])
#             xt2 = x_arr[0,table_nod[4][el] - 1]
#             yt2 = x_arr[1,table_nod[4][el] - 1]
#             t2 = np.array([xt2,yt2])
#             xt3 = x_arr[0,table_nod[5][el] - 1]
#             yt3 = x_arr[1,table_nod[5][el] - 1]
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
#             elif test1[0] > test1[1]:
#                 orient -= 1
#             if test2[0] < test2[1]:
#                 orient += 1
#             elif test2[0] > test2[1]:
#                 orient -= 1
#             if test3[0] < test3[1]:
#                 orient += 1
#             elif test3[0] > test3[1]:
#                 orient -= 1
#             if orient > 0:
#                 normal_vec = [-1*(y2-y1), (x2-x1)]
#             elif orient < 0:
#                 normal_vec = [(y2-y1), -1*(x2-x1)]
#             else:
#                 raise Warning, \
#                 'Cannot find orientation of normal vector'
#             n_vec_norm = normal_vec/np.linalg.norm(normal_vec)

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
#             dr = np.array([x2-x1, y2-y1, 0.0])
#             for n in [n0, n1]:
#                 # E-fields
#                 e1_x = sim_EM_wguide.sol1[0,n,EM_ival1,el]
#                 e1_y = sim_EM_wguide.sol1[1,n,EM_ival1,el]
#                 e1_z = sim_EM_wguide.sol1[2,n,EM_ival1,el]
#                 e2_x = sim_EM_wguide.sol1[0,n,EM_ival2,el]
#                 e2_y = sim_EM_wguide.sol1[1,n,EM_ival2,el]
#                 e2_z = sim_EM_wguide.sol1[2,n,EM_ival2,el]
#                 # Displacement fields
#                 u = sim_EM_wguide.sol1 #ToDo: replace with actual u field
#                 u_x = u[0,n,AC_ival,el]
#                 u_y = u[1,n,AC_ival,el]
#                 u_z = u[2,n,AC_ival,el]

#                 u_n = u_x*n_vec_norm[0] + u_y*n_vec_norm[1]
#                 # n_vec_norm[2] = 0 # z-comp!
#                 # n_cross_e1 = np.array([n_vec_norm[1]*e1_z - n_vec_norm[2]*e1_y],
#                 #     [-n_vec_norm[0]*e1_z + n_vec_norm[2]*e1_x],
#                 #     [n_vec_norm[0]*e1_y - n_vec_norm[1]*e1_x])
#                 n_cross_e1 = np.array([[n_vec_norm[1]*e1_z],
#                     [-n_vec_norm[0]*e1_z],
#                     [n_vec_norm[0]*e1_y - n_vec_norm[1]*e1_x]])
#                 n_cross_e2 = np.array([[n_vec_norm[1]*e2_z],
#                     [-n_vec_norm[0]*e2_z],
#                     [n_vec_norm[0]*e2_y - n_vec_norm[1]*e2_x]])
#                 inter_term1 = (eps_a - eps_b)*eps_0*np.conj(n_cross_e1)*n_cross_e2
#                 # ToDo: what is d field?
#                 n_cross_d_1 = n_cross_e1
#                 n_cross_d_2 = n_cross_e2
#                 inter_term2 = (1./eps_a - 1./eps_b)*(1./eps_0)*np.conj(n_cross_d_1)*n_cross_d_2
#                 integrand = u_n*(inter_term1 - inter_term2)
#                 Q_MB += integrand/2.0
#                 print Q_MB



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
