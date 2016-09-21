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
# inc_a_x = 300
# inc_a_y = 280

### Optical parameters
eps = 12.25
num_EM_modes = 20
num_AC_modes = 20

EM_ival1=0
EM_ival2=EM_ival1
AC_ival='All'
# AC_ival=2

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
   c_12 = E*v/((1+v)*(1-2*v))
   c_44 = E*(1-2*v)/((1+v)*(1-2*v))

   return c_11, c_12, c_44

### Acoustic parameters
# Inclusion a
s = 2330  # kg/m3
c_11 = 165.7e9; c_12 = 63.9e9; c_44 = 79.6e9  # Pa
p_11 = -0.094; p_12 = 0.017; p_44 = -0.051
eta_11 = 5.9e-3 ; eta_12 = 5.16e-3 ; eta_44 = 0.620e-3  # Pa
# E = 170e9
# v = 0.28
# c_11, c_12, c_44 = isotropic_stiffness(E, v)
# p_11 = 0.09; p_12 = -0.017; p_44 = -0.051
# # print c_11, c_12, c_44
inc_a_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]

wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        bkg_material=materials.Material(1.0 + 0.0j),
                        inc_a_material=materials.Material(np.sqrt(eps)),
                        loss=False, inc_a_AC=inc_a_AC_props,
                        lc_bkg=0.1, lc2=40.0, lc3=20.0)
                        # make_mesh_now=True)#, plotting_fields=True, plot_imag=1)#,
                        # mesh_file='rect_acoustic_3.mail')

### Calculate Electromagnetic Modes
# sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes)
# np.savez('wguide_data', sim_EM_wguide=sim_EM_wguide)
npzfile = np.load('wguide_data.npz')
sim_EM_wguide = npzfile['sim_EM_wguide'].tolist()
# print 'k_z of EM wave \n', sim_EM_wguide.Eig_value
# plotting.plt_mode_fields(sim_EM_wguide, xlim=0.4, ylim=0.4, EM_AC='EM')#,
    # n_points=1000, quiver_steps=10)
# plotting.plt_mode_fields(sim_EM_wguide, xlim=0.45, ylim=0.45, EM_AC='EM')


### Calculate Acoustic Modes
# Backward SBS
# Acoustic k has to push optical mode from +ve lightline to -ve, hence factor 2.
q_acoustic = 2*np.real(sim_EM_wguide.Eig_value[0])
# print q_acoustic*inc_a_x*1e-9/np.pi
# Forward (intramode) SBS
# q_acoustic = 0.0
# sim_AC_wguide = wguide.calc_AC_modes(wl_nm, q_acoustic,
#     num_AC_modes, EM_sim=sim_EM_wguide, shift_Hz=12e9)# shift_Hz=18e9)
# np.savez('wguide_data_AC', sim_AC_wguide=sim_AC_wguide)
npzfile = np.load('wguide_data_AC.npz')
sim_AC_wguide = npzfile['sim_AC_wguide'].tolist()
print 'Res freq of AC wave (GHz) \n', np.real(sim_AC_wguide.Eig_value)*1e-9
# plotting.plt_mode_fields(sim_AC_wguide, EM_AC='AC')#, add_name='-check')


# # Try to test with a simple field we know the answer to
# sim_AC_wguide.Omega_AC = np.ones(num_AC_modes)
# sim_AC_wguide.AC_mode_overlap = np.ones(num_AC_modes)
# grad = 1.0
# for el in range(sim_AC_wguide.n_msh_el):
#     for ival in range(num_AC_modes):
#         for n in range(6):
#             local_nod = sim_AC_wguide.table_nod[n][el]-1
#             # sim_AC_wguide.sol1[0,n,ival,el] = grad*sim_AC_wguide.x_arr[0,local_nod]
#             # sim_AC_wguide.sol1[1,n,ival,el] = grad*sim_AC_wguide.x_arr[0,local_nod]
#             # sim_AC_wguide.sol1[2,n,ival,el] = grad*sim_AC_wguide.x_arr[0,local_nod]
#             sim_AC_wguide.sol1[0,n,ival,el] = 0.0
#             sim_AC_wguide.sol1[1,n,ival,el] = 0.0
#             sim_AC_wguide.sol1[2,n,ival,el] = 1.0
# plotting.plt_mode_fields(sim_AC_wguide, EM_AC='AC')


### Calculate interaction integrals
SBS_gain, Q_PE, Q_MB, alpha = integration.gain_and_qs(
    sim_EM_wguide, sim_AC_wguide, q_acoustic,
    EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival)
# np.savez('wguide_data_AC_gain', SBS_gain=SBS_gain, alpha=alpha)
# npzfile = np.load('wguide_data_AC_gain.npz')
# SBS_gain = npzfile['SBS_gain']
# alpha = npzfile['alpha']




# Q_Rakich = 1000
# alpha_Rakich = sim_AC_wguide.Omega_AC/(2*Q_Rakich)
# # print SBS_gain[0,0,:]#/alpha

# print 310.25*(1./98.70e-6)/alpha_Rakich[2]
# print 2464.98*(1./27.75e-6)/alpha_Rakich[4]
# print 36.55*(1./43.90e-6)/alpha_Rakich[8]

# print 1307.0*2*(1./98.70e-6)/(2*np.pi*12.34e9)
# print "Q", sim_AC_wguide.Omega_AC/(2*alpha)

# print 'lc_bkg = ', wguide.lc
# print 'lc_2 = ', wguide.lc2
# print 'lc_3 = ', wguide.lc3

# print 'alpha 2 / CW alpha', alpha[2]/(1./98.70e-6)
# print 'alpha 4 / CW alpha', alpha[4]/(1./27.75e-6)
# print 'alpha 8 / CW alpha', alpha[8]/(1./43.90e-6)

# print 'SBS_gain 2', SBS_gain[0,0,2]/alpha[2]
# print 'SBS_gain 4', SBS_gain[0,0,4]/alpha[4]
# print 'SBS_gain 8', SBS_gain[0,0,8]/alpha[8]
# print 'SBS_gain 2 / CW gain', SBS_gain[0,0,2]/alpha[2]/1141#310.25
# print 'SBS_gain 4 / CW gain', SBS_gain[0,0,4]/alpha[4]/6000#2464.98
# print 'SBS_gain 8 / CW gain', SBS_gain[0,0,8]/alpha[8]/36.55
# print 'SBS_gain 2 / CW gain (using CW alpha)', SBS_gain[0,0,2]/(1./98.70e-6)/310.25
# print 'SBS_gain 4 / CW gain (using CW alpha)', SBS_gain[0,0,4]/(1./27.75e-6)/2464.98
# print 'SBS_gain 8 / CW gain (using CW alpha)', SBS_gain[0,0,8]/(1./43.90e-6)/36.55

# SBS_gain[0,0,2] = (1./98.70e-6)*310.25
# SBS_gain[0,0,4] = (1./27.75e-6)*2464.98
# SBS_gain[0,0,8] = (1./43.90e-6)*36.55

# alpha = sim_AC_wguide.Omega_AC/1000
# print alpha[2]*310.25
# print alpha[4]*2464.98
# print alpha[8]*36.55

import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
plt.figure(figsize=(13,13))
plt.clf()
tune_range = 2e4
AC_detuning_range = np.append(np.linspace(-2, 0, tune_range), np.linspace(0, 2, tune_range)[1:])*1e9
# print min(abs(AC_detuning_range))
# Line width of resonances should be v_g * alpha, but we don't have convenient access to v_g
# speed_in_Si = 9620 # m/s
# LW = speed_in_Si*alpha
phase_v = sim_AC_wguide.Eig_value/q_acoustic # phase velocity as approximation to group velocity
LW = phase_v*alpha
# LW = np.ones(len(alpha))
for AC_i in range(num_AC_modes):
   gain_list = SBS_gain[EM_ival1,EM_ival2,AC_i]*LW[AC_i]/(LW[AC_i]**2 + AC_detuning_range**2)
   # gain_list = SBS_gain[0,0,AC_i]*alpha[AC_i]/(alpha[AC_i]**2 + AC_detuning_range**2)
   freq_list_GHz = np.real(sim_AC_wguide.Eig_value[AC_i] + AC_detuning_range)*1e-9
   plt.plot(freq_list_GHz, np.real(gain_list),linewidth=3)
plt.xlim(10,25)
plt.savefig('gain.pdf')
plt.close()
