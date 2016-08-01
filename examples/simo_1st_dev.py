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


speed_c = 299792458
### Geometric parameters 
## All spacial variables given in nm!
# opt_freq_GHz = 11
# wl_nm = 2*np.pi*speed_c/(opt_freq_GHz)
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
                        lc_bkg=0.1, lc2=50.0, lc3=50.0)#,
                        # make_mesh_now=False, plotting_fields=False,
                        # mesh_file='rect_acoustic_3.mail')


# opt_freq_GHz = np.linspace(11, 25, 50)
# wl_nms = 2*np.pi*speed_c/opt_freq_GHz
# gain_array = []

# for wl_nm in wl_nms:
### Calculate Electromagnetic Modes
sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes)
np.savez('wguide_data', sim_EM_wguide=sim_EM_wguide)
# npzfile = np.load('wguide_data.npz')
# sim_EM_wguide = npzfile['sim_EM_wguide'].tolist()
print 'k_z of EM wave \n', sim_EM_wguide.Eig_value
# plotting.plt_mode_fields(sim_EM_wguide, xlim=0.4, ylim=0.4, EM_AC='EM')


### Calculate Acoustic Modes
# Backward SBS
# Acoustic k has to push optical mode from +ve lightline to -ve, hence factor 2.
q_acoustic = 2*sim_EM_wguide.Eig_value[0]
# # Forward (intramode) SBS
# q_acoustic = 0.0
sim_AC_wguide = wguide.calc_AC_modes(wl_nm, q_acoustic, num_AC_modes, EM_sim=sim_EM_wguide)
# sim_AC_wguide = wguide.calc_AC_modes(wl_nm, q_acoustic, num_AC_modes, EM_sim=None)
# np.savez('wguide_data_AC', sim_AC_wguide=sim_AC_wguide)
# npzfile = np.load('wguide_data_AC.npz')
# sim_AC_wguide = npzfile['sim_AC_wguide'].tolist()
print 'Omega of AC wave \n', sim_AC_wguide.Eig_value/(2*np.pi)*1e-9 # GHz
# prop_AC_modes = np.array([np.real(x) for x in sim_AC_wguide.Eig_value if abs(np.real(x)) > abs(np.imag(x))])
# prop_AC_modes = np.array([x for x in prop_AC_modes if np.real(x) > 0.0])
# print 'Omega of AC wave \n', prop_AC_modes*1e-9/(2*np.pi*8.54e3/inc_a_x) # GHz
# plotting.plt_mode_fields(sim_AC_wguide, EM_AC='AC')


### Calculate interaction integrals
SBS_gain, Q_PE, Q_MB, alpha, P1, P3 = integration.gain_and_qs(sim_EM_wguide, 
                           sim_AC_wguide, q_acoustic, 
                           EM_ival1=0, EM_ival2=0, AC_ival=2)

print "num_EM_modes", num_EM_modes
print "num_AC_modes", num_AC_modes
print "lc_bkg", wguide.lc
print "lc_bkg", wguide.lc2
print "lc_bkg", wguide.lc3
print "alpha[2]", alpha[2]
alpha_2 = 1/98.70e-6
print "CW_alpha/alpha[2]", alpha_2/alpha[2]
print "EM power", P1
print "AC power", P3
print "Q_PE", Q_PE
print "Gain", SBS_gain

    # gain_array.append(np.real(SBS_gain))


# import matplotlib
# matplotlib.use('pdf')
# import matplotlib.pyplot as plt
# plt.figure(figsize=(13,13))
# plt.clf()
# # plt.plot(gain_array,opt_freq_GHz,'r',linewidth=3)
# plt.plot(opt_freq_GHz,gain_array,'r',linewidth=3)
# plt.savefig('gain.pdf')
# plt.close()