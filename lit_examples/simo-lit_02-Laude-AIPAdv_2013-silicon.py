""" Calculate the backward SBS gain spectra of a
    silicon waveguide surrounded in air.

    Show how to save simulation objects (eg EM mode calcs)
    to expedite the process of altering later parts of
    simulations.
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
unitcell_x = 7*wl_nm
unitcell_y = unitcell_x
inc_a_x = 1500
inc_a_y = 1000
inc_shape = 'rectangular'

# Optical Parameters
n_inc_a = 3.47
num_EM_modes = 20
num_AC_modes = 800
EM_ival1=0
EM_ival2=EM_ival1
AC_ival='All'

# Silicon - Laude AIP Advances 2013
s = 2331  # kg/m3
c_11 = 166e9; c_12 = 64e9; c_44 = 79e9
p_11 = -0.1; p_12 = -0.01; p_44 = -0.051
eta_11 = 1.6e-3 ; eta_12 = 1.29e-3 ; eta_44 = 0.16e-3  # Pa s
inc_a_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]

# Use all specified parameters to create a waveguide object.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        bkg_material=materials.Material(1.0 + 0.0j),
                        inc_a_material=materials.Material(n_inc_a),
                        loss=False, inc_a_AC=inc_a_AC_props,
                        lc_bkg=3, lc2=2000.0, lc3=20.0)

# Expected effective index of fundamental guided mode.
n_eff=3.4

# Calculate Electromagnetic Modes
# sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes, n_eff=n_eff)
# np.savez('wguide_data-Laude_Si', sim_EM_wguide=sim_EM_wguide)
npzfile = np.load('wguide_data-Laude_Si.npz')
sim_EM_wguide = npzfile['sim_EM_wguide'].tolist()
# plotting.plt_mode_fields(sim_EM_wguide, xlim=0.4, ylim=0.4, EM_AC='EM', add_name='Laude_Si')
# plotting.plt_mode_fields(sim_EM_wguide, EM_AC='EM', add_name='NW')

# Print the wavevectors of EM modes.
print 'k_z of EM modes \n', np.round(np.real(sim_EM_wguide.Eig_values), 4)

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_wguide.Eig_values*((wl_nm*1e-9)/(2.*np.pi)))
print "n_eff = ", np.round(n_eff_sim, 4)

k_AC = 2*np.real(sim_EM_wguide.Eig_values[0])

shift_Hz = 31e9

# Calculate Acoustic Modes
sim_AC_wguide = wguide.calc_AC_modes(wl_nm, num_AC_modes, k_AC=k_AC,
    EM_sim=sim_EM_wguide, shift_Hz=shift_Hz)
# np.savez('wguide_data_AC-Laude_Si', sim_AC_wguide=sim_AC_wguide)
# npzfile = np.load('wguide_data_AC-Laude_Si.npz')
# sim_AC_wguide = npzfile['sim_AC_wguide'].tolist()
# plotting.plt_mode_fields(sim_AC_wguide, EM_AC='AC', add_name='Laude_Si')

# Print the frequencies of AC modes.
print 'Freq of AC modes (GHz) \n', np.round(np.real(sim_AC_wguide.Eig_values)*1e-9, 4)

# Do not calculate the acoustic loss from our fields, but instead set a 
# predetirmined Q factor. (Useful for instance when replicating others results).
set_q_factor = 15000.

# Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# as well as just for PE, and just for MB. Also calculate acoustic loss alpha.
SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha = integration.gain_and_qs(
    sim_EM_wguide, sim_AC_wguide, k_AC,
    EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival)#, fixed_Q=set_q_factor)
# np.savez('wguide_data_AC-Laude_Si_gain', SBS_gain=SBS_gain, alpha=alpha)
# npzfile = np.load('wguide_data_AC-Laude_Si_gain.npz')
# SBS_gain = npzfile['SBS_gain']
# alpha = npzfile['alpha']


# Construct the SBS gain spectrum, built up
# from Lorentzian peaks of the individual modes.
tune_steps = 5e4
tune_range = 10 # GHz
# Construct an odd range of frequencies that is guaranteed to include 
# the central resonance frequency.
detuning_range = np.append(np.linspace(-1*tune_range, 0, tune_steps),
                   np.linspace(0, tune_range, tune_steps)[1:])*1e9 # GHz
# Line width of resonances should be v_g * alpha,
# but we don't have convenient access to v_g, therefore
# phase velocity as approximation to group velocity
phase_v = sim_AC_wguide.Eig_values/k_AC
linewidth = phase_v*alpha

freq_min = 20  # GHz
freq_max = 45  # GHz
interp_grid_points = 10000
interp_grid = np.linspace(freq_min, freq_max, interp_grid_points)
interp_values = np.zeros(interp_grid_points)

freq_threshold = 4  # GHz

plt.figure()
plt.clf()
for AC_i in range(len(alpha)):
    if np.real(sim_AC_wguide.Eig_values[AC_i])*1e-9 > freq_threshold:
        # print SBS_gain[EM_ival1,EM_ival2,AC_i]
        # print alpha[AC_i]
        gain_list = np.real(SBS_gain[EM_ival1,EM_ival2,AC_i]/alpha[AC_i]
                     *linewidth[AC_i]**2/(linewidth[AC_i]**2 + detuning_range**2))
        freq_list_GHz = np.real(sim_AC_wguide.Eig_values[AC_i] + detuning_range)*1e-9
        plt.plot(freq_list_GHz, gain_list)
        # set up an interpolation for summing all the gain peaks
        interp_spectrum = np.interp(interp_grid, freq_list_GHz, gain_list)
        interp_values += interp_spectrum
plt.plot(interp_grid, interp_values, 'k', linewidth=3, label="Total")
plt.legend(loc=0)
plt.xlim(freq_min,freq_max)
plt.xlabel('Frequency (GHz)')
plt.ylabel('Gain (1/Wm)')
plt.savefig('gain_spectra-mode_comps.pdf')
plt.close()


plt.figure()
plt.clf()
for AC_i in range(len(alpha)):
    if np.real(sim_AC_wguide.Eig_values[AC_i])*1e-9 > freq_threshold:
        gain_list = np.real(SBS_gain[EM_ival1,EM_ival2,AC_i]/alpha[AC_i]
                     *linewidth[AC_i]**2/(linewidth[AC_i]**2 + detuning_range**2))
        freq_list_GHz = np.real(sim_AC_wguide.Eig_values[AC_i] + detuning_range)*1e-9
        plt.plot(freq_list_GHz, 20*np.log10(abs(gain_list)))
        # set up an interpolation for summing all the gain peaks
        interp_spectrum = np.interp(interp_grid, freq_list_GHz, gain_list)
        interp_values += interp_spectrum
plt.plot(interp_grid, 20*np.log10(abs(interp_values)), 'k', linewidth=3, label="Total")
plt.legend(loc=0)
plt.xlim(freq_min,freq_max)
plt.ylim(-80,50)
plt.xlabel('Frequency (GHz)')
plt.ylabel('Gain (dB)')
plt.savefig('gain_spectra-mode_comps-dB.pdf')
plt.close()


freq_min_zoom = freq_min
freq_max_zoom = freq_max
interp_values = np.zeros(interp_grid_points)
interp_values_PE = np.zeros(interp_grid_points)
interp_values_MB = np.zeros(interp_grid_points)
plt.figure()
plt.clf()
for AC_i in range(len(alpha)):
    gain_list = np.real(SBS_gain[EM_ival1,EM_ival2,AC_i]/alpha[AC_i]
                 *linewidth[AC_i]**2/(linewidth[AC_i]**2 + detuning_range**2))
    freq_list_GHz = np.real(sim_AC_wguide.Eig_values[AC_i] + detuning_range)*1e-9
    interp_spectrum = np.interp(interp_grid, freq_list_GHz, gain_list)
    interp_values += interp_spectrum

    gain_list_PE = np.real(SBS_gain_PE[EM_ival1,EM_ival2,AC_i]/alpha[AC_i]
                 *linewidth[AC_i]**2/(linewidth[AC_i]**2 + detuning_range**2))
    interp_spectrum_PE = np.interp(interp_grid, freq_list_GHz, gain_list_PE)
    interp_values_PE += interp_spectrum_PE

    gain_list_MB = np.real(SBS_gain_MB[EM_ival1,EM_ival2,AC_i]/alpha[AC_i]
                 *linewidth[AC_i]**2/(linewidth[AC_i]**2 + detuning_range**2))
    interp_spectrum_MB = np.interp(interp_grid, freq_list_GHz, gain_list_MB)
    interp_values_MB += interp_spectrum_MB
plt.plot(interp_grid, interp_values, 'k', linewidth=3, label="Total")
plt.plot(interp_grid, interp_values_PE, 'r', linewidth=3, label="PE")
plt.plot(interp_grid, interp_values_MB, 'g', linewidth=3, label="MB")
plt.legend(loc=0)
plt.xlim(freq_min_zoom,freq_max_zoom)
plt.xlabel('Frequency (GHz)')
plt.ylabel('Gain (1/Wm)')
plt.savefig('gain_spectra-MB_PE_comps.pdf')
plt.close()


plt.figure()
plt.clf()
plt.plot(interp_grid, 20*np.log10(abs(interp_values)), 'k', linewidth=3, label="Total")
plt.plot(interp_grid, 20*np.log10(abs(interp_values_PE)), 'r', linewidth=3, label="PE")
plt.plot(interp_grid, 20*np.log10(abs(interp_values_MB)), 'g', linewidth=3, label="MB")
plt.legend(loc=0)
plt.xlim(freq_min_zoom,freq_max_zoom)
plt.xlabel('Frequency (GHz)')
plt.ylabel('Gain (dB)')
plt.savefig('gain_spectra-MB_PE_comps-dB.pdf')
plt.close()