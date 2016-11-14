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
unitcell_x = 2.5*wl_nm
unitcell_y = unitcell_x
inc_a_x = 314.7
inc_a_y = 0.9*inc_a_x
inc_shape = 'rectangular'

# Optical Parameters
eps = 12.25
num_EM_modes = 20
num_AC_modes = 20
EM_ival1=0
EM_ival2=EM_ival1
AC_ival='All'

# Acoustic Parameters
s = 2330  # kg/m3
c_11 = 165.7e9; c_12 = 63.9e9; c_44 = 79.6e9  # Pa
p_11 = -0.094; p_12 = 0.017; p_44 = -0.051
eta_11 = 5.9e-3 ; eta_12 = 5.16e-3 ; eta_44 = 0.620e-3  # Pa
inc_a_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]

# Use all specified parameters to create a waveguide object.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        bkg_material=materials.Material(1.0 + 0.0j),
                        inc_a_material=materials.Material(np.sqrt(eps)),
                        loss=False, inc_a_AC=inc_a_AC_props,
                        lc_bkg=0.1, lc2=40.0, lc3=20.0)


# Calculate Electromagnetic Modes
sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes)
# Save calculated :Simmo: object for EM calculation.
np.savez('wguide_data', sim_EM_wguide=sim_EM_wguide)

# The previous two lines can be commented out and the following
# two uncommented to provide precisely the same objects for the
# remainder of the simulation.
# npzfile = np.load('wguide_data.npz')
# sim_EM_wguide = npzfile['sim_EM_wguide'].tolist()

# Print the wavevectors of EM modes.
print 'k_z of EM wave \n', sim_EM_wguide.Eig_value


# Choose acoustic wavenumber to solve for
# Backward SBS
# AC mode couples EM modes on +ve to -ve lightline, hence factor 2.
k_AC = 2*np.real(sim_EM_wguide.Eig_value[0])


# Calculate Acoustic Modes
sim_AC_wguide = wguide.calc_AC_modes(wl_nm, k_AC,
    num_AC_modes, EM_sim=sim_EM_wguide)
# np.savez('wguide_data_AC', sim_AC_wguide=sim_AC_wguide)

# The previous two lines can be commented out and the following
# two uncommented to provide precisely the same objects for the
# remainder of the simulation.
# npzfile = np.load('wguide_data_AC.npz')
# sim_AC_wguide = npzfile['sim_AC_wguide'].tolist()

# Print the frequencies of AC modes.
print 'Res freq of AC wave (GHz) \n', np.real(sim_AC_wguide.Eig_value)*1e-9



# Calculate interaction integrals and SBS gain
# - for PE and MB effects combined, as well as
# just for PE, and just for MB. Also calculate acoustic loss alpha.
SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha = integration.gain_and_qs(
    sim_EM_wguide, sim_AC_wguide, k_AC,
    EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival)
np.savez('wguide_data_AC_gain', SBS_gain=SBS_gain, alpha=alpha)

# The previous two lines can be commented out and the following
# three uncommented to provide precisely the same objects for the
# remainder of the simulation.
# npzfile = np.load('wguide_data_AC_gain.npz')
# SBS_gain = npzfile['SBS_gain']
# alpha = npzfile['alpha']


# Construct the SBS gain spectrum, built up
# from Lorentzian peaks of the individual modes.
plt.figure(figsize=(13,13))
plt.clf()
tune_steps = 5e4
tune_range = 10 # GHz
# Construct an odd range of frequencies that is guaranteed to include the central resonance freq.
detuning_range = np.append(np.linspace(-1*tune_range, 0, tune_steps),
                   np.linspace(0, tune_range, tune_steps)[1:])*1e9 # GHz
# Line width of resonances should be v_g * alpha,
# but we don't have convenient access to v_g, therefore
# phase velocity as approximation to group velocity
phase_v = sim_AC_wguide.Eig_value/k_AC
line_width = phase_v*alpha
interp_grid_points = 10000
interp_grid = np.linspace(10, 25, interp_grid_points)
interp_values = np.zeros(interp_grid_points)
for AC_i in range(num_AC_modes):
    gain_list = np.real(SBS_gain[EM_ival1,EM_ival2,AC_i]/alpha[AC_i]
                 *line_width[AC_i]**2/(line_width[AC_i]**2 + detuning_range**2))
    freq_list_GHz = np.real(sim_AC_wguide.Eig_value[AC_i] + detuning_range)*1e-9
    plt.plot(freq_list_GHz, gain_list,linewidth=3)
    # set up an interpolation for summing all the gain peaks
    interp_spectrum = np.interp(interp_grid, freq_list_GHz, gain_list)
    interp_values += interp_spectrum

plt.plot(interp_grid, interp_values, 'k', linewidth=4)
plt.xlim(10,25)
plt.xlabel('Frequency (GHz)', fontsize=16)
plt.ylabel('Gain 1/(Wm)', fontsize=16)
plt.savefig('gain_spectra.pdf')
plt.close()

