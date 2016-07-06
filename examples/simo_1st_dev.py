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
eta_11 = 5.9e-3 ; eta_12 = 5.16e-3 ; eta_44 = 620e-6  # Pa s
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
# Backward SBS
# Acoustic k has to push optical mode from +ve lightline to -ve, hence factor 2.
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
SBS_gain, Q_PE, Q_MB = integration.gain_and_qs(sim_EM_wguide, 
                           sim_AC_wguide, q_acoustic, AC_ival=2)

