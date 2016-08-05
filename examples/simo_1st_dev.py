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

# start = time.time()

speed_c = 299792458
### Geometric parameters 
## All spacial variables given in nm!
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

EM_ival1=0
EM_ival2=0
AC_ival='All'

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
                        lc_bkg=0.1, lc2=20.0, lc3=20.0)#,
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
# # Acoustic k has to push optical mode from +ve lightline to -ve, hence factor 2.
q_acoustic = 2*sim_EM_wguide.Eig_value[0]
# # Forward (intramode) SBS
# q_acoustic = 0.0
# sim_AC_wguide = wguide.calc_AC_modes(wl_nm, q_acoustic, num_AC_modes, EM_sim=sim_EM_wguide)
# np.savez('wguide_data_AC', sim_AC_wguide=sim_AC_wguide)
npzfile = np.load('wguide_data_AC.npz')
sim_AC_wguide = npzfile['sim_AC_wguide'].tolist()
# print 'Res freq of AC wave (GHz) \n', sim_AC_wguide.Eig_value*1e-9
# prop_AC_modes = np.array([np.real(x) for x in sim_AC_wguide.Eig_value if abs(np.real(x)) > abs(np.imag(x))])
# prop_AC_modes = np.array([x for x in prop_AC_modes if np.real(x) > 0.0])
# print 'Omega of AC wave \n', prop_AC_modes*1e-9/(2*np.pi*8.54e3/inc_a_x) # GHz
# plotting.plt_mode_fields(sim_AC_wguide, EM_AC='AC')


### Calculate interaction integrals
SBS_gain, Q_PE, Q_MB, alpha, P1, P3 = integration.gain_and_qs(sim_EM_wguide, 
                           sim_AC_wguide, q_acoustic, 
                           EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival)
# np.savez('wguide_data_AC_gain', SBS_gain=SBS_gain, alpha=alpha)
# npzfile = np.load('wguide_data_AC_gain.npz')
# SBS_gain = npzfile['SBS_gain']
# alpha = npzfile['alpha']

### TMP TEST
q_acoustic = 1j*q_acoustic
### Calculate interaction integrals
SBS_gain, Q_PE, Q_MB, alpha2, P1, P3 = integration.gain_and_qs(sim_EM_wguide, 
                           sim_AC_wguide, q_acoustic, 
                           EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival)

print alpha
print alpha2
print alpha - alpha2


# # # trim1 = 5
# # # trim = 13
# # syms = integration.symmetries(sim_AC_wguide)
# # # for i in range(trim-trim1):
# # for i in range(num_AC_modes):
# #    # i = i+trim1
# #    sym = syms[i]
# #    print 'Res freq', sim_AC_wguide.Eig_value[i]*1e-9
# #    print 'SES gain', SBS_gain[0,0,i]/alpha[i]
# #    print 'alpha', alpha[i]
# #    if sym[0] == 1 and sym[1] == 1 and sym[2] == 1:
# #       print 'A'
# #    if sym[0] == -1 and sym[1] == 1 and sym[2] == -1:
# #       print 'B1'
# #    if sym[0] == 1 and sym[1] == -1 and sym[2] == -1:
# #       print 'B2'
# #    if sym[0] == -1 and sym[1] == -1 and sym[2] == 1:
# #       print 'B3'

# # # print "num_EM_modes", num_EM_modes
# # # print "num_AC_modes", num_AC_modes
# # # print "lc_bkg", wguide.lc
# # # print "lc_bkg", wguide.lc2
# # # print "lc_bkg", wguide.lc3

# print "alpha[1]", alpha[1]
# alpha_1 = 1/142.79e-6
# print alpha_1
# print "CW_alpha/alpha[1]", alpha_1/alpha[1]

# print "alpha[2]", alpha[2]
# alpha_2 = 1/98.70e-6
# print alpha_2
# print "CW_alpha/alpha[2]", alpha_2/alpha[2]

# print "alpha[3]", alpha[3]
# alpha_3 = 1/203.98e-6
# print alpha_3
# print "CW_alpha/alpha[3]", alpha_3/alpha[3]

# print "alpha[4]", alpha[4]
# alpha_4 = 1/27.75e-6
# print alpha_4
# print "CW_alpha/alpha[4]", alpha_4/alpha[4]

# print "alpha[5]", alpha[5]
# alpha_5 = 1/66.69e-6
# print alpha_5
# print "CW_alpha/alpha[5]", alpha_5/alpha[5]

# print "alpha[8]", alpha[8]
# alpha_8 = 1/43.90e-6
# print alpha_8
# print "CW_alpha/alpha[8]", alpha_8/alpha[8]


# # # print "EM power", P1
# # # print "AC power", P3
# # # # print "Q_PE", Q_PE
# # print "Gain", SBS_gain[0,0,2]/alpha[2]
# # print "Gain", SBS_gain[0,0,2]/alpha_2
# # print "Gain", SBS_gain[0,0,4]/alpha[4]
# # print "Gain", SBS_gain[0,0,4]/alpha_4
# # print "Gain", SBS_gain[0,0,8]/alpha[8]
# # print "Gain", SBS_gain[0,0,8]/alpha_8
# # # # print "Gain", SBS_gain[0,0,:]/alpha
# # # print "Gain", SBS_gain[0,0,AC_ival]/alpha[AC_ival]


# # import matplotlib
# # matplotlib.use('pdf')
# # import matplotlib.pyplot as plt
# # plt.figure(figsize=(13,13))
# # plt.clf()
# # AC_detuning_range = np.linspace(-2e9, 2e9, 1e3)
# # # Line width of resonances should be v_g * alpha, but we don't have convenient access to v_g
# # # speed_in_Si = 9620 # m/s 
# # # LW = speed_in_Si*alpha
# # phase_v = sim_AC_wguide.Eig_value/q_acoustic # phase velocity as approximation to group velocity
# # LW = phase_v*alpha
# # # print LW
# # for AC_i in range(num_AC_modes):
# # # for AC_i in range(trim-trim1):
# #    # AC_i = AC_i + trim1
# #    gain_list = SBS_gain[EM_ival1,EM_ival2,AC_i]*LW[AC_i]/(LW[AC_i]**2 + AC_detuning_range**2)
# #    # gain_list = SBS_gain[0,0,AC_i]*alpha[AC_i]/(alpha[AC_i]**2 + AC_detuning_range**2)
# #    freq_list_GHz = np.real(sim_AC_wguide.Eig_value[AC_i] + AC_detuning_range)*1e-9
# #    plt.plot(freq_list_GHz, np.real(gain_list),linewidth=3)
# # plt.xlim(10,25)
# # plt.savefig('gain.pdf')
# # plt.close()


# # AC_detuning_range = np.linspace(-1e9, 1e9, 1e3)
# # gain_list = SBS_gain*alpha[AC_ival]/(alpha[AC_ival]**2 + AC_detuning_range**2)
# # freq_list_GHz = np.real(sim_AC_wguide.Eig_value[AC_ival] + AC_detuning_range)*1e-9
# # plt.plot(freq_list_GHz, np.real(gain_list),'r',linewidth=3)
# # plt.savefig('gain.pdf')
# # plt.close()