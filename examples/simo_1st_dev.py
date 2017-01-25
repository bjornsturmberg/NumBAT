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
unitcell_y = unitcell_x
inc_a_x = 314.7
inc_a_y = 0.9*inc_a_x
inc_shape = 'rectangular'
# inc_shape = 'circular'


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
# p_11 = -0.09; p_12 = -0.017; p_44 = -0.051
# # print c_11, c_12, c_44
inc_a_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]

wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        bkg_material=materials.Material(1.0 + 0.0j),
                        inc_a_material=materials.Material(np.sqrt(eps)),
                        loss=False, inc_a_AC=inc_a_AC_props,
                        lc_bkg=0.1, lc2=20.0, lc3=20.0)
                        # make_mesh_now=True)#, plotting_fields=True, plot_imag=1)#,
                        # mesh_file='rect_acoustic_3.mail')

### Calculate Electromagnetic Modes
sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes)
# np.savez('wguide_data', sim_EM_wguide=sim_EM_wguide)
# npzfile = np.load('wguide_data.npz')
# sim_EM_wguide = npzfile['sim_EM_wguide'].tolist()
# print 'k_z of EM wave \n', sim_EM_wguide.Eig_value
# plotting.plt_mode_fields(sim_EM_wguide, xlim=0.4, ylim=0.4, EM_AC='EM')#,
    # n_points=1000, quiver_steps=10)
# plotting.plt_mode_fields(sim_EM_wguide, xlim=0.45, ylim=0.45, EM_AC='EM')


### Calculate Acoustic Modes
# Backward SBS
# Acoustic k has to push optical mode from +ve lightline to -ve, hence factor 2.
k_AC = 2*np.real(sim_EM_wguide.Eig_value[0])
# print k_AC*inc_a_x*1e-9/np.pi
# Forward (intramode) SBS
# k_AC = 0.0
sim_AC_wguide = wguide.calc_AC_modes(wl_nm, k_AC,
    num_AC_modes, EM_sim=sim_EM_wguide)
# np.savez('wguide_data_AC', sim_AC_wguide=sim_AC_wguide)
# npzfile = np.load('wguide_data_AC.npz')
# sim_AC_wguide = npzfile['sim_AC_wguide'].tolist()
print 'Res freq of AC wave (GHz) \n', np.real(sim_AC_wguide.Eig_value)*1e-9
# plotting.plt_mode_fields(sim_AC_wguide, EM_AC='AC')#, add_name='-check')


### Calculate interaction integrals
SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha = integration.gain_and_qs(
    sim_EM_wguide, sim_AC_wguide, k_AC,
    EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival, fixed_Q=None)
# np.savez('wguide_data_AC_gain', SBS_gain=SBS_gain, alpha=alpha)
# npzfile = np.load('wguide_data_AC_gain.npz')
# SBS_gain = npzfile['SBS_gain']
# alpha = npzfile['alpha']


print "SBS_gain", SBS_gain[0,0,:]/alpha
print "SBS_gain_MB", SBS_gain_MB[0,0,:]/alpha
print "SBS_gain_PE", SBS_gain_PE[0,0,:]/alpha


import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt


fixed_q = 1000.
alpha = (np.pi*k_AC/fixed_q)*np.ones(sim_AC_wguide.num_modes)


# freqs_gains = []
# Calculate the EM effective index of the waveguide.
n_eff_sim = round(np.real(sim_EM_wguide.Eig_value[0]*((wl_nm*1e-9)/(2.*np.pi))), 4)
print "n_eff", n_eff_sim

# Construct the SBS gain spectrum, built up
# from Lorentzian peaks of the individual modes.
tune_steps = 5e4
tune_range = 10 # GHz
# Construct an odd range of frequencies that is guaranteed to include the central resonance freq.
detuning_range = np.append(np.linspace(-1*tune_range, 0, tune_steps),
                   np.linspace(0, tune_range, tune_steps)[1:])*1e9 # GHz
k_AC = 2*np.real(sim_EM_wguide.Eig_value[0])
phase_v = sim_AC_wguide.Eig_value/k_AC
line_width = phase_v*alpha
interp_grid_points = 10000
interp_grid = np.linspace(10, 25, interp_grid_points)
interp_values = np.zeros(interp_grid_points)

plt.figure(figsize=(13,13))
plt.clf()
for AC_i in range(len(alpha)):
    gain_list = np.real(SBS_gain[EM_ival1,EM_ival2,AC_i]/alpha[AC_i]
                 *line_width[AC_i]**2/(line_width[AC_i]**2 + detuning_range**2))
    freq_list_GHz = np.real(sim_AC_wguide.Eig_value[AC_i] + detuning_range)*1e-9
    plt.plot(freq_list_GHz, gain_list)#,linewidth=2)
    # set up an interpolation for summing all the gain peaks
    interp_spectrum = np.interp(interp_grid, freq_list_GHz, gain_list)
    interp_values += interp_spectrum
# freqs_gains.append(zip(interp_grid, interp_values))
# plt.plot(interp_grid, interp_values, 'k', linewidth=3)
plt.xlim(10,25)
plt.xlabel('Frequency (GHz)')#, fontsize=16)
plt.ylabel('Gain 1/(Wm)')#, fontsize=16)
plt.savefig('gain_spectra-mode_comps.pdf')
plt.close()


interp_values = np.zeros(interp_grid_points)
interp_values_PE = np.zeros(interp_grid_points)
interp_values_MB = np.zeros(interp_grid_points)
plt.figure(figsize=(13,13))
plt.clf()
for AC_i in range(len(alpha)):
    gain_list = np.real(SBS_gain[EM_ival1,EM_ival2,AC_i]/alpha[AC_i]
                 *line_width[AC_i]**2/(line_width[AC_i]**2 + detuning_range**2))
    freq_list_GHz = np.real(sim_AC_wguide.Eig_value[AC_i] + detuning_range)*1e-9
    interp_spectrum = np.interp(interp_grid, freq_list_GHz, gain_list)
    interp_values += interp_spectrum

    gain_list_PE = np.real(SBS_gain_PE[EM_ival1,EM_ival2,AC_i]/alpha[AC_i]
                 *line_width[AC_i]**2/(line_width[AC_i]**2 + detuning_range**2))
    interp_spectrum_PE = np.interp(interp_grid, freq_list_GHz, gain_list_PE)
    interp_values_PE += interp_spectrum_PE

    gain_list_MB = np.real(SBS_gain_MB[EM_ival1,EM_ival2,AC_i]/alpha[AC_i]
                 *line_width[AC_i]**2/(line_width[AC_i]**2 + detuning_range**2))
    interp_spectrum_MB = np.interp(interp_grid, freq_list_GHz, gain_list_MB)
    interp_values_MB += interp_spectrum_MB

# freqs_gains.append(zip(interp_grid, interp_values))
plt.plot(interp_grid, interp_values_PE, 'r', linewidth=3)
plt.plot(interp_grid, interp_values_MB, 'g', linewidth=3)
plt.plot(interp_grid, interp_values, 'b', linewidth=3)
# plt.xlim(10,25)
plt.xlim(12,14)
plt.xlabel('Frequency (GHz)')#, fontsize=16)
plt.ylabel('Gain 1/(Wm)')#, fontsize=16)
plt.savefig('gain_spectra-MB_PE_comps.pdf')
plt.close()