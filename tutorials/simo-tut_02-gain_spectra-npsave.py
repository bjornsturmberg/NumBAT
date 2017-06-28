""" Calculate the backward SBS gain spectra of a
    silicon waveguide surrounded in air.

    Show how to save simulation objects (eg. EM mode calcs)
    to expedite the process of altering later parts of
    simulations.
"""

import time
import datetime
import numpy as np
import sys
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt

sys.path.append("../backend/")
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

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 20
EM_ival_pump = 0
EM_ival_Stokes = EM_ival_pump
AC_ival = 'All'

# Use of a more refined mesh to produce field plots.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        material_bkg=materials.Air,
                        material_a=materials.Si,
                        lc_bkg=3, lc2=2000.0, lc3=1000.0)


# Expected effective index of fundamental guided mode.
n_eff = wguide.material_a.n-0.1

# Calculate Electromagnetic modes.
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff)
# Save calculated :Simmo: object for EM calculation.
np.savez('wguide_data', sim_EM_pump=sim_EM_pump)

# Once npz files have been saved from one simulation run,
# the previous three lines can be commented and the following
# two line uncommented. This provides precisely the same objects
# for the remainder of the simulation.
# npzfile = np.load('wguide_data.npz')
# sim_EM_pump = npzfile['sim_EM_pump'].tolist()

sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)
np.savez('wguide_data2', sim_EM_Stokes=sim_EM_Stokes)
# npzfile = np.load('wguide_data2.npz')
# sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()

# Print the wavevectors of EM modes.
print('k_z of EM modes \n', np.round(np.real(sim_EM_pump.Eig_values), 4))

# Plot the E fields of the EM modes fields - specified with EM_AC='EM_E'.
# Zoom in on the central region (of big unitcell) with xlim_, ylim_ args.
plotting.plt_mode_fields(sim_EM_pump, xlim_min=0.4, xlim_max=0.4, 
                         ylim_min=0.4, ylim_max=0.4, EM_AC='EM_E')
# Plot the H fields of the EM modes - specified with EM_AC='EM_H'.
plotting.plt_mode_fields(sim_EM_pump, xlim_min=0.4, xlim_max=0.4, 
                         ylim_min=0.4, ylim_max=0.4, EM_AC='EM_H')

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.Eig_values[0]*((wl_nm*1e-9)/(2.*np.pi)))
print("n_eff", np.round(n_eff_sim, 4))
# Acoustic wavevector
k_AC = np.real(sim_EM_pump.Eig_values[0] - sim_EM_Stokes.Eig_values[0])

# Calculate Acoustic modes.
sim_AC = wguide.calc_AC_modes(num_modes_AC, k_AC, EM_sim=sim_EM_pump)
# Save calculated :Simmo: object for AC calculation.
np.savez('wguide_data_AC', sim_AC=sim_AC)

# npzfile = np.load('wguide_data_AC.npz')
# sim_AC = npzfile['sim_AC'].tolist()

# Print the frequencies of AC modes.
print('Freq of AC modes (GHz) \n', np.round(np.real(sim_AC.Eig_values)*1e-9, 4))

# Plot the AC modes fields, important to specify this with EM_AC='AC'.
# The AC modes are calculated on a subset of the full unitcell,
# which excludes vacuum regions, so no need to restrict area plotted.
# We want to get pdf files so set pdf_png='pdf' 
# (default is png as these are easier to flick through).
plotting.plt_mode_fields(sim_AC, EM_AC='AC', pdf_png='pdf')

# Do not calculate the acoustic loss from our fields, instead set a Q factor.
set_q_factor = 1000.

# Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha, Q_factors = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_q_factor)
# Save the gain calculation results
np.savez('wguide_data_AC_gain', SBS_gain=SBS_gain, SBS_gain_PE=SBS_gain_PE, 
            SBS_gain_MB=SBS_gain_MB, alpha=alpha)

# Once npz files have been saved from one simulation run,
# the previous six lines can be commented and the following
# five line uncommented. This provides precisely the same objects
# for the remainder of the simulation.
# npzfile = np.load('wguide_data_AC_gain.npz')
# SBS_gain = npzfile['SBS_gain']
# SBS_gain_PE = npzfile['SBS_gain_PE']
# SBS_gain_MB = npzfile['SBS_gain_MB']
# alpha = npzfile['alpha']

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = np.real(sim_AC.Eig_values[0])*1e-9 - 2  # GHz
freq_max = np.real(sim_AC.Eig_values[-1])*1e-9 + 2  # GHz
plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha, k_AC,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max)
# Zoomed in version
freq_min = 11  # GHz
freq_max = 15  # GHz
plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha, k_AC,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max, add_name='_zoom')
