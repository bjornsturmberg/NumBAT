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
unitcell_y = unitcell_x
inc_a_x = 300
inc_a_y = 280
# print inc_a_x/wl_nm
# print 2*np.pi/wl_nm
# print (2*np.pi/wl_nm)/(2*np.pi*speed_c/inc_a_x)
# inc_a_x = .35*wl_nm
# print inc_a_x
inc_shape = 'rectangular'


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
                        # lc_bkg=0.2, lc2=20.0, lc3=20.0)#,
                        lc_bkg=0.1, lc2=30.0, lc3=20.0)


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

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
# marks = ['ok','ob','or','og','oc','om','^k','^b','^r','^g','^c','^m','sk','sb','sr','sg','sc','sm']
plt.clf()
plt.figure(figsize=(10,6))
ax = plt.subplot(1,1,1)
for q_ac in np.linspace(0.0,q_acoustic,20):
    sim_AC_wguide = wguide.calc_AC_modes(wl_nm, q_ac, num_AC_modes, EM_sim=sim_EM_wguide)
    prop_AC_modes = np.array([np.real(x) for x in sim_AC_wguide.Eig_value if abs(np.real(x)) > abs(np.imag(x))])
    # prop_AC_modes = np.array([np.real(x) for x in sim_AC_wguide.Eig_value if abs(np.imag(x)) < 1e-0])
    # prop_AC_modes = np.array([np.real(x) for x in sim_AC_wguide.Eig_value])
    # prop_AC_modes = np.array([x for x in prop_AC_modes if np.real(x) > 0.0])
    # print 'Omega of AC wave \n', sim_AC_wguide.Eig_value*1e-9 # GHz
    # plotting.plt_mode_fields(sim_AC_wguide, EM_AC='AC')
    sym_list = integration.symmetries(sim_AC_wguide)

    # prop_AC_modes_allowed = integration.allowed_symmetries(sim_AC_wguide)
    for i in range(len(prop_AC_modes)):
        Om = prop_AC_modes[i]*1e-9
        # plt.plot(q_ac/q_acoustic, Om, marks[i%len(marks)])
        plt.plot(q_ac/q_acoustic, Om, 'ok')
        if sym_list[i][0] == 1 and sym_list[i][1] == 1 and sym_list[i][2] == 1:
            plt.plot(np.real(q_ac/q_acoustic), Om, 'or')
        if sym_list[i][0] == -1 and sym_list[i][1] == 1 and sym_list[i][2] == -1:
            plt.plot(np.real(q_ac/q_acoustic), Om, 'vc')
        if sym_list[i][0] == 1 and sym_list[i][1] == -1 and sym_list[i][2] == -1:
            plt.plot(np.real(q_ac/q_acoustic), Om, 'sb')
        if sym_list[i][0] == -1 and sym_list[i][1] == -1 and sym_list[i][2] == 1:
            plt.plot(np.real(q_ac/q_acoustic), Om, '^g')
    ax.set_ylim(0,20)
    ax.set_xlim(0,1)
plt.xlabel(r'Axial wavevector (normalised)', fontsize=16)
plt.ylabel(r'Frequency (GHz)', fontsize=16)
plt.savefig('disp.pdf', bbox_inches='tight')
plt.close()
