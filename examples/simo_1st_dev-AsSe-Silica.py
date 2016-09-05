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
wl_nm = 1550
unitcell_x = 2.5*1550
unitcell_y = unitcell_x
inc_shape = 'rectangular'
# inc_shape = 'circular'
inc_a_x = 500
inc_a_y = inc_a_x
# inc_a_x = 1500
# inc_a_y = 680

### Optical parameters
n_b = 1.44
n_i = 2.83
num_EM_modes = 20
num_AC_modes = 20#400

EM_ival1=0
EM_ival2=0
AC_ival='All'

### Acoustic parameters
def isotropic_stiffness(E, v):
   """
   Calculate the stiffness matrix components of isotropic 
   materials, given the two free parameters:
   E: Youngs_modulus
   v: Poisson_ratio

   Ref: http://www.efunda.com/formulae/solid_mechanics/mat_mechanics/hooke_isotropic.cfm
   """
   c_11 = E*(1-v)/((1+v)*(1-2*v))
   c_12 = c_11
   c_44 = E*(1-2*v)/((1+v)*(1-2*v))

   return c_11, c_12, c_44

# Background - Silca
s = 2200  # kg/m3
E = 7.3e10
v = 0.17
c_11, c_12, c_44 = isotropic_stiffness(E, v)
p_11 = 0.121; p_12 = 0.270; p_44 = -0.075
eta_11 = 1.6e-3 ; eta_12 = 1.29e-3 ; eta_44 = 0.16e-3  # Pa s
bkg_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]
# Inclusion a - Chalc
s = 4460  # kg/m3
E = 1.73e10
v = 0.29
c_11, c_12, c_44 = isotropic_stiffness(E, v)
p_11 = 0.314; p_12 = 0.266; p_44 = 0.024
eta_11 = 9e-3 ; eta_12 = 7.5e-3 ; eta_44 = 0.75e-3  # Pa s
inc_a_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]

wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        bkg_material=materials.Material(n_b),
                        inc_a_material=materials.Material(n_i),
                        loss=False, 
                        bkg_AC=bkg_AC_props,
                        inc_a_AC=inc_a_AC_props,
                        lc_bkg=0.1, lc2=40.0, lc3=20.0)


### Calculate Electromagnetic Modes
n_eff = 2.5
shift_Hz = n_eff**2 * (2*np.pi/(wl_nm*1e-9))**2
# sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes, shift_Hz=shift_Hz)
# np.savez('wguide_data-chalc', sim_EM_wguide=sim_EM_wguide)
npzfile = np.load('wguide_data-chalc.npz')
sim_EM_wguide = npzfile['sim_EM_wguide'].tolist()
# print 'k_z of EM wave \n', sim_EM_wguide.Eig_value
# plotting.plt_mode_fields(sim_EM_wguide, xlim=0.4, ylim=0.4, EM_AC='EM')


### Calculate Acoustic Modes
# Backward SBS
# Acoustic k has to push optical mode from +ve lightline to -ve, hence factor 2.
q_acoustic = 2*sim_EM_wguide.Eig_value[0]
# # Forward (intramode) SBS
# q_acoustic = 0.0
# As2S3_bulk_velocity = 2595m/s - NOT AsSe!
# @1550 gives 7.6 GHz as freq in normal waveguides
shift_Hz = 7.6e9
# shift_Hz = 6.0e9
# shift_Hz = shift_Hz/3.
print 'shift_Hz', shift_Hz
sim_AC_wguide = wguide.calc_AC_modes(wl_nm, q_acoustic, num_AC_modes, 
   EM_sim=sim_EM_wguide)#, shift_Hz=shift_Hz )
# np.savez('wguide_data_AC', sim_AC_wguide=sim_AC_wguide)
# npzfile = np.load('wguide_data_AC.npz')
# sim_AC_wguide = npzfile['sim_AC_wguide'].tolist()
print 'Res freq of AC wave (GHz) \n', sim_AC_wguide.Eig_value*1e-9
# plotting.plt_mode_fields(sim_AC_wguide, EM_AC='AC', n_points=300, xlim=0.6, ylim=0.3)


# # import time
# # start = time.time()
# ### Calculate interaction integrals
# SBS_gain, Q_PE, Q_MB, alpha, P1, P3 = integration.gain_and_qs(sim_EM_wguide, 
#                            sim_AC_wguide, q_acoustic, 
#                            EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival)
# # elapsed = (time.time() - start)
# # print 'TIME', elapsed
# np.savez('wguide_data_AC_gain', SBS_gain=SBS_gain, alpha=alpha)
# # npzfile = np.load('wguide_data_AC_gain.npz')
# # SBS_gain = npzfile['SBS_gain']
# # alpha = npzfile['alpha']

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

# # print "Gain", SBS_gain[0,0,AC_ival]/alpha[AC_ival]


# import matplotlib
# matplotlib.use('pdf')
# import matplotlib.pyplot as plt
# plt.figure(figsize=(13,13))
# plt.clf()
# AC_detuning_range = np.linspace(-2e9, 2e9, 1e3)
# # Line width of resonances should be v_g * alpha, but we don't have convenient access to v_g
# # speed_in_Si = 9620 # m/s 
# # LW = speed_in_Si*alpha
# phase_v = sim_AC_wguide.Eig_value/q_acoustic # phase velocity as approximation to group velocity
# LW = phase_v*alpha
# # print LW
# for AC_i in range(num_AC_modes):
# # for AC_i in range(trim-trim1):
#    # AC_i = AC_i + trim1
#    gain_list = SBS_gain[EM_ival1,EM_ival2,AC_i]*LW[AC_i]/(LW[AC_i]**2 + AC_detuning_range**2)
#    # gain_list = SBS_gain[0,0,AC_i]*alpha[AC_i]/(alpha[AC_i]**2 + AC_detuning_range**2)
#    freq_list_GHz = np.real(sim_AC_wguide.Eig_value[AC_i] + AC_detuning_range)*1e-9
#    plt.plot(freq_list_GHz, np.real(gain_list),linewidth=3)
# plt.xlim(10,25)
# plt.savefig('gain.pdf')
# plt.close()
