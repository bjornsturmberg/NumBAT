"""
    mode_calcs.py is a subroutine of NumBAT that contains methods to
    calculate the EM and Acoustic modes of a structure.

    Copyright (C) 2016  Bjorn Sturmberg, Kokou Dossou 
"""

import numpy as np

# import materials
# import objects
# import mode_calcs
# import plotting
from fortran import NumBAT


def gain_and_qs(sim_EM_wguide, sim_AC_wguide, q_acoustic,
                EM_ival1=0, EM_ival2=0, AC_ival=0):
    """ Calculate interaction integrals and SBS gain.
    """

# Calc Qs of a range of selected modes
# Pass in the q_acoustic as this is beta of AC mode (which gives z derivative)

### Notes about internals of fortran integration
# Calc overlap of basis functions (and PE tensor and epsilon) 
# Then use this multiple times for calc of each mode field values

# phi is values of Lagrange polynomials (1-6) at that node.
# grad is value of gradient of Lagrange polynomials (1-6) at that node.

    ncomps = 3
    nnodes = 6
    speed_c = 299792458
    opt_freq_GHz = speed_c/sim_EM_wguide.wl_nm # putting in wl in nm gives you GHz
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

    Fortran_debug = 0
### Calc Q_photoelastic Eq. 33
    try:
        Q_PE = NumBAT.photoelastic_int(
            sim_EM_wguide.num_modes, sim_AC_wguide.num_modes, EM_ival1,
            EM_ival2, AC_ival, sim_AC_wguide.n_msh_el, sim_AC_wguide.n_msh_pts, nnodes,
            sim_AC_wguide.table_nod, sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
            sim_AC_wguide.structure.nb_typ_el_AC, sim_AC_wguide.structure.p_tensor,
            q_acoustic, trimmed_EM_field, sim_AC_wguide.sol1,
            relevant_eps_effs, Fortran_debug)
    except KeyboardInterrupt:
        print "\n\n Routine photoelastic_int interrupted by keyboard.\n\n"

### Calc loss alphar Eq. 45
    try:
        alpha_numbat = NumBAT.ac_alpha_int(sim_AC_wguide.num_modes, sim_AC_wguide.n_msh_el, 
            sim_AC_wguide.n_msh_pts, nnodes, sim_AC_wguide.table_nod, 
            sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
            sim_AC_wguide.structure.nb_typ_el_AC, sim_AC_wguide.structure.eta_tensor,
            q_acoustic, sim_AC_wguide.sol1, sim_AC_wguide.sol1, 
            sim_AC_wguide.AC_mode_overlap, Fortran_debug)
    except KeyboardInterrupt:
        print "\n\n Routine ac_alpha_int interrupted by keyboard.\n\n"


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

    eps_0 = 8.854187817e-12
    Q_MB = 0.0 # Haven't implemented Moving Boundary integral (but nor did Rakich)
    Q = Q_PE + Q_MB
    gain = 2*opt_freq_GHz*1e9*sim_AC_wguide.Eig_value[AC_ival]*np.real(Q*np.conj(Q))
    P1 = sim_EM_wguide.EM_mode_overlap[EM_ival1]#*eps_0*unitcell_x*1e-9*unitcell_y*1e-9
    P2 = sim_EM_wguide.EM_mode_overlap[EM_ival2]#*eps_0*unitcell_x*1e-9*unitcell_y*1e-9
    P3 = sim_AC_wguide.AC_mode_overlap[AC_ival]#*inc_a_x*1e-9*inc_a_y*1e-9
    normal_fact = P1*P2*P3

    print "Q", Q
    print "gain", gain
    print "EM mode 1 power", P1
    print "EM mode 2 power", P2
    print "AC mode power", P3

    gain = gain/normal_fact
    alpha = 1/98.70e-6
    SBS_gain = gain/alpha

    return SBS_gain, Q_PE, Q_MB


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
# EM_ival1 = 0
# EM_ival2 = 0
# AC_ival = 0
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