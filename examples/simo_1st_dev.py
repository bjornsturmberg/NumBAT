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


unitcell_x = 2.5*1550
inc_a_x = 314.7
unitcell_y = unitcell_x
inc_a_y = 0.9*inc_a_x
inc_shape = 'rectangular'

### Optical parameters
eps = 12.25

wl_nm = 1550
num_EM_modes = 30

### Acoustic parameters
s = 2330  # kg/m3

c_11 = 165.7  # GPa
c_12 = 63.9  # GPa
c_44 = 79.6  # GPa

p_11 = -0.044
p_12 = 0.017
p_44 = -0.051

eta_11 = 5.9  # m Pa s
eta_12 = 5.16  # m Pa s
eta_44 = 620  # mu Pa s

# wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
#                         bkg_material=materials.Material(1.0 + 0.0j),
#                         inc_a_material=materials.Material(np.sqrt(eps)),
#                         lc_bkg=0.2, lc2=1.0, lc3=1.0, check_msh=False)

# sim_wguide = wguide.calc_modes(wl_nm, num_EM_modes)
# np.savez('wguide_data', sim_wguide=sim_wguide)
npzfile = np.load('wguide_data.npz')
sim_wguide = npzfile['sim_wguide'].tolist()

# betas = sim_wguide.k_z
# print 'k_z of EM wave \n', betas
# plotting.plot_EM_modes(sim_wguide)

# # Test overlap
# nel = sim_wguide.n_msh_el
# type_el = sim_wguide.type_el
# table_nod = sim_wguide.table_nod
x_arr = sim_wguide.x_arr
# nnodes = 6
# xel = np.zeros((2,nnodes))
# nod_el_p = np.zeros(nnodes)
# xx = [0,0]
# nquad = 16
# [wq, xq, yq] = integration.quad_triangle(nquad)

# integrand = 0.0
# print np.shape(sim_wguide.sol1)

# for ival in [0]:
# # for ival in range(len(sim_wguide.k_z)):
#     NumBAT.EM_mode_energy_int()

n_msh_el = sim_wguide.n_msh_el
type_el = sim_wguide.type_el
table_nod = sim_wguide.table_nod
n_msh_pts = sim_wguide.n_msh_pts
print 'n_msh_el', n_msh_el
print 'n_msh_pts', sim_wguide.n_msh_pts
print 'table_nod', np.shape(table_nod)
# print table_nod[0]
print 'x_arr', np.shape(x_arr)
print 'type_el', np.shape(type_el)

# find triangles associated with node
# find triangles type
# if types not all equal -> node on interface
# for node in range(1,sim_wguide.n_msh_pts+1):
#     print np

# bkg_el_type = 1
interface_nodes = []
edge_el_list = []
node_array = -1*np.ones(n_msh_pts)
### Find nodes that are in elements of various types
### and find elements that have multiple nodes of various types
### ie. are not single verticies on an interface.
for el in range(n_msh_el):
    el_type = type_el[el]
    for i in range(6):
        node = table_nod[i][el]
        # Check if first time seen this node
        if node_array[node - 1] == -1: # adjust to python indexing
            node_array[node - 1] = el_type
        else:
            if node_array[node - 1] != el_type:
                interface_nodes.append(node)
                ## line below is redundant because elements sorted by type w type 1 first
                # if el_type is not bkg_el_type:
                edge_el_list.append(el)

interface_nodes = list(set(interface_nodes))
from collections import Counter
edge_els_multi_nodes = [k for (k,v) in Counter(edge_el_list).iteritems() if v > 1]

test_orient = [1,-1]
test1 = [0,0]
test2 = [0,0]
test3 = [0,0]
for el in edge_els_multi_nodes:
    # These are all possible edge line segments.
    for [n1,n2] in [[0,3],[3,1],[1,4],[4,2],[2,5],[5,0]]:
        node0 = table_nod[n1][el]
        node1 = table_nod[n2][el]
        if node0 in interface_nodes and node1 in interface_nodes:
            # coordinates of line seg. nodes
            x1 = x_arr[0,table_nod[n1][el] - 1]
            y1 = x_arr[1,table_nod[n1][el] - 1]
            x2 = x_arr[0,table_nod[n2][el] - 1]
            y2 = x_arr[1,table_nod[n2][el] - 1]
            # coordinates of non-vertex nodes, used to test orientation
            xt1 = x_arr[0,table_nod[3][el] - 1]
            yt1 = x_arr[1,table_nod[3][el] - 1]
            t1 = np.array([xt1,yt1])
            xt2 = x_arr[0,table_nod[4][el] - 1]
            yt2 = x_arr[1,table_nod[4][el] - 1]
            t2 = np.array([xt2,yt2])
            xt3 = x_arr[0,table_nod[5][el] - 1]
            yt3 = x_arr[1,table_nod[5][el] - 1]
            t3 = np.array([xt3,yt3])
            for i in [0, 1]:
                t = test_orient[i]
                normal_vec = [t*(y2-y1), -1*t*(x2-x1)]
                # start half way along line seg. out along n vector
                test_point = np.array([x1+(x2-x1+normal_vec[0])/2.,
                              y1+(y2-y1+normal_vec[1])/2.])
                test1[i] = np.linalg.norm(test_point-t1)
                test2[i] = np.linalg.norm(test_point-t2)
                test3[i] = np.linalg.norm(test_point-t3)
            orient = 0
            if test1[0] < test1[1]:
                orient += 1
            elif test1[0] > test1[1]:
                orient -= 1
            if test2[0] < test2[1]:
                orient += 1
            elif test2[0] > test2[1]:
                orient -= 1
            if test3[0] < test3[1]:
                orient += 1
            elif test3[0] > test3[1]:
                orient -= 1
            if orient > 0:
                normal_vec = [-1*(y2-y1), (x2-x1)]
            elif orient < 0:
                normal_vec = [(y2-y1), -1*(x2-x1)]
            else:
                raise Warning, \
                'Cannot find orientation of normal vector'
            normal_vec_norm = normal_vec/np.linalg.norm(normal_vec)
            # print normal_vec_norm
            # print x1
            # print y1
            # print x2
            # print y2


# import matplotlib
# matplotlib.use('pdf')
# import matplotlib.pyplot as plt
# # plt.figure(figsize=(13,13))
# # count = 0
# interface_nodes2 = []
# for el in [0,1,2,3,4,5,6,7,8]:#edge_els_multi_nodes:
#     plt.clf()
#     # below is incorrect - need to include all possible combinations of the elements nodes
#     for i in range(0,6):
#         print table_nod[i][el] - 1
#         x = x_arr[0,table_nod[i][el] - 1]
#         y = x_arr[1,table_nod[i][el] - 1]
#         print 'x1 = ', x_arr[0,table_nod[i][el] - 1]
#         print 'y1 = ', x_arr[1,table_nod[i][el] - 1]
#         plt.plot(x, y, 'o')
#         plt.text(x+0.001, y+0.001, str(i))
#     plt.savefig('triangle_%i.png' %el)




### Calc Q_photoelastic Eq. 33
### Calc Q_deformation_pol Eq. 36
### Calc Q_moving_boundary Eq. 41
