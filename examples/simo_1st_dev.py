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
#                         lc_bkg=0.09, lc2=3.0, lc3=3.0, check_msh=False)

# sim_wguide = wguide.calc_modes(wl_nm, num_EM_modes)
# np.savez('wguide_data', sim_wguide=sim_wguide)
# print repr(wguide_data)
npzfile = np.load('sim_wguide.npz')
sim_wguide = npzfile['sim_wguide'].tolist()

# betas = sim_wguide.k_z
# print 'k_z of EM wave \n', betas
# plotting.plot_EM_modes(sim_wguide)

# Test overlap
nel = sim_wguide.n_msh_el
type_el = sim_wguide.type_el
table_nod = sim_wguide.table_nod
x_arr = sim_wguide.x_arr
nnodes = 6
xel = np.zeros((2,nnodes))
nod_el_p = np.zeros(nnodes)
xx = [0,0]
nquad = 16
[wq, xq, yq] = integration.quad_triangle(nquad)

integrand = 0.0
print np.shape(sim_wguide.sol1)

for ival in [0]:
# for ival in range(len(sim_wguide.k_z)):
    # loop over all elements 
    for iel in range(nel):
        typ_e = type_el[iel]
        for inode in range(nnodes):
            tmp_1 = table_nod[inode,iel] - 1 # adjust for python index from 0
            nod_el_p[inode] = tmp_1
            xel[0,inode] = x_arr[0,tmp_1]
            xel[1,inode] = x_arr[1,tmp_1]

        for iq in range(nquad):
            xx[0] = xq[iq]
            xx[1] = yq[iq]
            ww = wq[iq]

            x_g, det_jacobian = NumBAT.jacobian_p1_2d(xx, xel, nnodes)
            coeff_1 = ww * abs(det_jacobian)

            for wave_comp in [0,1,2]:
                point_value = coeff_1*sim_wguide.sol1[wave_comp,inode,ival,iel]
                integrand += point_value
