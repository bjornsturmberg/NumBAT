""" Calculate the backward SBS gain spectra of a
    chalcogenide (As2S3) waveguide surrounded in silica.
"""

import time
import datetime
import numpy as np
import sys
sys.path.append("../backend/")
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

import materials
import objects
import mode_calcs
import integration
import plotting
from fortran import NumBAT


# Geometric Parameters - all in nm.
wl_nm = 1550
unitcell_x = 4.5*wl_nm
unitcell_y = unitcell_x
inc_a_x = 850
inc_a_y = 680
inc_shape = 'rectangular'

# Optical parameters
n_b = 1.44
n_i = 2.83
num_EM_modes = 20
# There are lots of leaky acoustic modes!
num_AC_modes = 250
EM_ival1=0
EM_ival2=EM_ival1
AC_ival='All'

# Acoustic Parameters
# Background - silca
# Use the isotropic parameters for silca.
s = 2200  # kg/m3
E = 7.3e10
v = 0.17
c_11, c_12, c_44 = materials.isotropic_stiffness(E, v)
p_11 = 0.121; p_12 = 0.270; p_44 = -0.075
eta_11 = 1.6e-3 ; eta_12 = 1.29e-3 ; eta_44 = 0.16e-3  # Pa s
bkg_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]
# Inclusion a - As2S3
s = 3210  # kg/m3
c_11 = 2.104e10; c_12 = 8.363e9; c_44 =6.337e9 # Pa
p_11 = 0.25; p_12 = 0.24; p_44 = 0.005
eta_11 = 9e-3 ; eta_12 = 7.5e-3 ; eta_44 = 0.75e-3  # Pa s
inc_a_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]

# Use all specified parameters to create a waveguide object.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        bkg_material=materials.Material(n_b),
                        inc_a_material=materials.Material(n_i),
                        loss=False, 
                        bkg_AC=bkg_AC_props,
                        inc_a_AC=inc_a_AC_props,plotting_fields=True,
                        lc_bkg=0.5, lc2=30.0, lc3=30.0)


# Calculate Electromagnetic Modes
# Provide an estimate of the effective index of the fundamental optical mode.
n_eff = 2.2
shift_Hz = n_eff**2 * (2*np.pi/(wl_nm*1e-9))**2
sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes, shift_Hz=shift_Hz)
# np.savez('wguide_data-chalc', sim_EM_wguide=sim_EM_wguide)
# npzfile = np.load('wguide_data-chalc.npz')
# sim_EM_wguide = npzfile['sim_EM_wguide'].tolist()
plotting.plt_mode_fields(sim_EM_wguide, xlim=0.4, ylim=0.4, EM_AC='EM')


# Print the wavevectors of EM modes.
print 'k_z of EM wave \n', (sim_EM_wguide.Eig_value, 4)
n_eff_sim = ((sim_EM_wguide.Eig_value[0]*((wl_nm*1e-9)/(2.*np.pi))), 4)
print 'n_eff of fund. EM mode \n', n_eff_sim
n_eff_sim = ((sim_EM_wguide.Eig_value[2]*((wl_nm*1e-9)/(2.*np.pi))), 4)
print 'n_eff of fund. EM mode \n', n_eff_sim


# # Choose acoustic wavenumber to solve for
# # Backward SBS
# # AC mode couples EM modes on +ve to -ve lightline, hence factor 2.
# k_AC = 2*np.real(sim_EM_wguide.Eig_value[0])


# # Calculate Acoustic Modes
# # As2S3_bulk_velocity = 2595m/s @1550 gives 7.6 GHz as freq in reg waveguides
# shift_Hz = 5.8e9
# print shift_Hz
# sim_AC_wguide = wguide.calc_AC_modes(wl_nm, k_AC, num_AC_modes, 
#    EM_sim=sim_EM_wguide)#, shift_Hz=shift_Hz )
# # np.savez('wguide_data_AC', sim_AC_wguide=sim_AC_wguide)

# # The previous two lines can be commented out and the following
# # two uncommented to provide precisely the same objects for the
# # remainder of the simulation.
# # npzfile = np.load('wguide_data_AC.npz')
# # sim_AC_wguide = npzfile['sim_AC_wguide'].tolist()

# # Print the frequencies of AC modes.
# print 'Res freq of AC wave (GHz) \n', np.real(sim_AC_wguide.Eig_value)*1e-9



# # Calculate interaction integrals and SBS gain
# # Don't trust the eta tensor parameters so choose to use a fixed Q factor for all modes
# fixed_Q = 500
# SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha = integration.gain_and_qs(
#     sim_EM_wguide, sim_AC_wguide, k_AC,
#     EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival, fixed_Q=fixed_Q)
# np.savez('wguide_data_AC_gain', SBS_gain=SBS_gain, alpha=alpha)

# # The previous two lines can be commented out and the following
# # three uncommented to provide precisely the same objects for the
# # remainder of the simulation.
# # npzfile = np.load('wguide_data_AC_gain.npz')
# # SBS_gain = npzfile['SBS_gain']
# # alpha = npzfile['alpha']


# # Construct the SBS gain spectrum, built up
# # from Lorentzian peaks of the individual modes.
# plt.figure(figsize=(13,13))
# plt.clf()
# tune_steps = 5e4
# tune_range = 10 # GHz
# # Construct an odd range of frequencies that is guaranteed to include the central resonance freq.
# detuning_range = np.append(np.linspace(-1*tune_range, 0, tune_steps),
#                    np.linspace(0, tune_range, tune_steps)[1:])*1e9 # GHz
# # Line width of resonances should be v_g * alpha,
# # but we don't have convenient access to v_g, therefore
# # phase velocity as approximation to group velocity
# phase_v = sim_AC_wguide.Eig_value/k_AC
# line_width = phase_v*alpha
# interp_grid_points = 10000
# interp_grid = np.linspace(10, 25, interp_grid_points)
# interp_values = np.zeros(interp_grid_points)
# for AC_i in range(num_AC_modes):
#     gain_list = np.real(SBS_gain[EM_ival1,EM_ival2,AC_i]/alpha[AC_i]
#                  *line_width[AC_i]**2/(line_width[AC_i]**2 + detuning_range**2))
#     freq_list_GHz = np.real(sim_AC_wguide.Eig_value[AC_i] + detuning_range)*1e-9
#     plt.plot(freq_list_GHz, gain_list,linewidth=3)
#     # set up an interpolation for summing all the gain peaks
#     interp_spectrum = np.interp(interp_grid, freq_list_GHz, gain_list)
#     interp_values += interp_spectrum

# plt.plot(interp_grid, interp_values, 'k', linewidth=4)
# plt.xlim(10,25)
# plt.xlabel('Frequency (GHz)', fontsize=16)
# plt.ylabel('Gain 1/(Wm)', fontsize=16)
# plt.savefig('gain_spectra-chalc.pdf')
# plt.close()

