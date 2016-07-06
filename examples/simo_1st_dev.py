import time
import datetime
import numpy as np
import sys
from multiprocessing import Pool
sys.path.append("../backend/")

import materials
import objects
import mode_calcs
import integration
import plotting
from fortran import NumBAT


### Geometric parameters
wl_nm = 1550
speed_c = 299792458
opt_freq_GHz = speed_c/wl_nm # putting in wl in nm gives you GHz
unitcell_x = 2.5*1550
inc_a_x = 314.7
unitcell_y = unitcell_x
inc_a_y = 0.9*inc_a_x
inc_shape = 'rectangular'
# inc_shape = 'circular'

### Optical parameters
eps = 12.25
num_EM_modes = 20
num_AC_modes = 20

### Acoustic parameters
# Inclusion a
s = 2330  # kg/m3
c_11 = 165.7e9; c_12 = 63.9e9; c_44 = 79.6e9  # Pa
p_11 = -0.044; p_12 = 0.017; p_44 = -0.051
eta_11 = 5.9e-3 ; eta_12 = 5.16e-3  # Pa s
eta_44 = 620e-6   # Pa s
inc_a_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]

wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        bkg_material=materials.Material(1.0 + 0.0j),
                        inc_a_material=materials.Material(np.sqrt(eps)),
                        loss=False, inc_a_AC=inc_a_AC_props,
                        lc_bkg=0.2, lc2=20.0, lc3=20.0)#,
                        # lc_bkg=0.2, lc2=70.0, lc3=100.0)#,
                        # make_mesh_now=False, plotting_fields=False,
                        # mesh_file='rect_acoustic_3.mail')


### Calculate Electromagnetic Modes
# sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes)
# np.savez('wguide_data', sim_EM_wguide=sim_EM_wguide)
npzfile = np.load('wguide_data.npz')
sim_EM_wguide = npzfile['sim_EM_wguide'].tolist()
# print 'k_z of EM wave \n', sim_EM_wguide.Eig_value
# plotting.plt_mode_fields(sim_EM_wguide, xlim=0.4, ylim=0.4, EM_AC='EM')


### Calculate Acoustic Modes
# Acoustic k has to push optical mode from -ve lightline to +ve, hence factor 2.
# Backward SBS
q_acoustic = 2*sim_EM_wguide.Eig_value[0]
# # Forward (intramode) SBS
# q_acoustic = 0.0
sim_AC_wguide = wguide.calc_AC_modes(wl_nm, q_acoustic, num_AC_modes, EM_sim=sim_EM_wguide)
np.savez('wguide_data_AC', sim_AC_wguide=sim_AC_wguide)
# npzfile = np.load('wguide_data_AC.npz')
# sim_AC_wguide = npzfile['sim_AC_wguide'].tolist()
# print 'Omega of AC wave \n', sim_AC_wguide.Eig_value*1e-9 # GHz
# prop_AC_modes = np.array([np.real(x) for x in sim_AC_wguide.Eig_value if abs(np.real(x)) > abs(np.imag(x))])
# prop_AC_modes = np.array([x for x in prop_AC_modes if np.real(x) > 0.0])
# print 'Omega of AC wave \n', prop_AC_modes*1e-9/(2*np.pi*8.54e3/inc_a_x) # GHz
# plotting.plt_mode_fields(sim_AC_wguide, EM_AC='AC')


### Calculate interaction integrals
# SBS_gain, Q_PE, Q_MB = integration.gain_and_qs(sim_EM_wguide, sim_AC_wguide)
# integration.gain_and_qs(sim_EM_wguide, sim_AC_wguide)


def gain_and_qs(sim_EM_wguide, sim_AC_wguide, q_acoustic,
                EM_mode1=0, EM_mode2=0, AC_mode=0):
    """ Calculate interaction integrals and SBS gain.
    """
    ncomps = 3
    nnodes = 6
    num_EM_modes = len(sim_EM_wguide.Eig_value)
    n_msh_el_AC = sim_AC_wguide.n_msh_el
    trimmed_EM_field = np.zeros((ncomps,nnodes,num_EM_modes,n_msh_el_AC), dtype=complex)
    for el in range(n_msh_el_AC):
        new_el = sim_AC_wguide.el_convert_tbl[el]
        for ival in range(num_EM_modes):
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

### Calc Q_photoelastic Eq. 33
    try:
        Fortran_debug = 0
        Q_PE = NumBAT.photoelastic_int(
            sim_EM_wguide.num_modes, sim_AC_wguide.num_modes, EM_mode1,
            EM_mode2, AC_mode, sim_AC_wguide.n_msh_el, sim_AC_wguide.n_msh_pts, nnodes,
            sim_AC_wguide.table_nod, sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
            sim_AC_wguide.structure.nb_typ_el_AC, sim_AC_wguide.structure.p_tensor,
            q_acoustic, trimmed_EM_field, sim_AC_wguide.sol1,
            relevant_eps_effs, Fortran_debug)

        # EM_mode_overlap2 = NumBAT.em_mode_energy_int_v2(
        #     sim_EM_wguide.wl_norm(), sim_EM_wguide.num_modes, sim_AC_wguide.n_msh_el, 
        #     sim_AC_wguide.n_msh_pts, nnodes, sim_AC_wguide.table_nod,
        #     sim_AC_wguide.x_arr, sim_EM_wguide.Eig_value, trimmed_EM_field)
        # # EM_mode_overlap2 = NumBAT.em_mode_energy_int_v2(
        # #     sim_EM_wguide.wl_norm(), sim_EM_wguide.num_modes, sim_EM_wguide.n_msh_el, 
        # #     sim_EM_wguide.n_msh_pts, nnodes, sim_EM_wguide.table_nod,
        # #     sim_EM_wguide.x_arr, sim_EM_wguide.Eig_value, sim_EM_wguide.sol1)

        # print sim_EM_wguide.EM_mode_overlap
        # print EM_mode_overlap2

        alpha_numbat = NumBAT.ac_alpha_int(sim_AC_wguide.num_modes, sim_AC_wguide.n_msh_el, 
            sim_AC_wguide.n_msh_pts, nnodes, sim_AC_wguide.table_nod, 
            sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
            sim_AC_wguide.structure.nb_typ_el_AC, sim_AC_wguide.structure.eta_tensor,
            q_acoustic, sim_AC_wguide.sol1, sim_AC_wguide.sol1, 
            sim_AC_wguide.AC_mode_overlap, Fortran_debug)

    except KeyboardInterrupt:
        print "\n\n Routine photoelastic_int",\
        "interrupted by keyboard.\n\n"

# Christians values for alpha of first 3 modes
    print '---------'
    print alpha_numbat[2]
    alpha = 1/98.70e-6
    print alpha
    print alpha_numbat[2]/alpha
    print '---------'
    print alpha_numbat[0]
    alpha_0 = 1/186.52e-6
    print alpha_0
    print alpha_numbat[0]/alpha_0
    print '---------'
    print alpha_numbat[1]
    alpha_1 = 1/142.79e-6
    print alpha_1
    print alpha_numbat[1]/alpha_1
    print '---------'
    # print '---------'
    # print 98.70/186.52
    # print alpha_numbat[0]/alpha_numbat[2]
    # print 98.70/142.79
    # print alpha_numbat[1]/alpha_numbat[2]
    # print '---------'




    # eps_0 = 8.854187817e-12
    # Q_MB = 0.0 # Haven't implemented Moving Boundary integral (but nor did Rakich)
    # Q = Q_PE + Q_MB
    # print "Q", Q
    # gain = 2*opt_freq_GHz*1e9*sim_AC_wguide.Eig_value[AC_mode]*np.real(Q*np.conj(Q))
    # normal_fact = sim_EM_wguide.EM_mode_overlap[EM_mode1]#*eps_0*unitcell_x*1e-9*unitcell_y*1e-9
    # normal_fact = normal_fact*sim_EM_wguide.EM_mode_overlap[EM_mode2]#*eps_0*unitcell_x*1e-9*unitcell_y*1e-9
    # normal_fact = normal_fact*sim_AC_wguide.AC_mode_overlap[AC_mode]#*inc_a_x*1e-9*inc_a_y*1e-9

    # print "gain", gain
    # print "EM mode 1 power", sim_EM_wguide.EM_mode_overlap[EM_mode1]#*eps_0*unitcell_x*1e-9*unitcell_y*1e-9
    # print "EM mode 2 power", sim_EM_wguide.EM_mode_overlap[EM_mode2]#*eps_0*unitcell_x*1e-9*unitcell_y*1e-9
    # print "AC mode power", sim_AC_wguide.AC_mode_overlap[AC_mode]#*inc_a_x*1e-9*inc_a_y*1e-9

    # gain = gain/normal_fact
    # alpha = 1/98.70e-6
    # print alpha_numbat[2]/alpha
    # gain = gain/alpha
    # print "gain normalised", gain
    # print "gain error factor", gain/310.




    # sim_EM_wguide.EM_mode_overlap_v2
    # print sim_AC_wguide.el_convert_tbl
    # sim_AC_wguide.node_convert_tbl

### Calc Q_moving_boundary Eq. 41


gain_and_qs(sim_EM_wguide, sim_AC_wguide, q_acoustic, AC_mode=2)







# ### Calc Q_moving_boundary Eq. 41
# from collections import Counter

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
# EM_ival_1 = 0
# EM_ival_2 = 0
# Ac_ival = 0
# eps_0 = 1.0 # ToDo: update this value
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
#                 e1_x = sim_EM_wguide.sol1[0,n,EM_ival_1,el]
#                 e1_y = sim_EM_wguide.sol1[1,n,EM_ival_1,el]
#                 e1_z = sim_EM_wguide.sol1[2,n,EM_ival_1,el]
#                 e2_x = sim_EM_wguide.sol1[0,n,EM_ival_2,el]
#                 e2_y = sim_EM_wguide.sol1[1,n,EM_ival_2,el]
#                 e2_z = sim_EM_wguide.sol1[2,n,EM_ival_2,el]
#                 # Displacement fields
#                 u = sim_EM_wguide.sol1 #ToDo: replace with actual u field
#                 u_x = u[0,n,Ac_ival,el]
#                 u_y = u[1,n,Ac_ival,el]
#                 u_z = u[2,n,Ac_ival,el]

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
# # integrand = 0.0
# # for ival in [0]:
# # # for ival in range(len(sim_EM_wguide.k_z)):
# #     NumBAT.EM_mode_energy_int()


