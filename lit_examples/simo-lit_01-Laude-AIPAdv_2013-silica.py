""" Replicating the results of
    Generation of phonons from electrostriction in small-core optical waveguides
    Laude et al.
    http://dx.doi.org/10.1063/1.4801936
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
num_EM_modes = 20
num_AC_modes = 120
EM_ival1 = 0
EM_ival2 = EM_ival1
AC_ival = 'All'

# Material parameters as in paper 
# Silicon
n = 1.44
s = 2203  # kg/m3
c_11 = 78e9; c_12 = 16e9; c_44 = 31e9
p_11 = 0.12; p_12 = 0.270; p_44 = -0.073
eta_11 = 1.6e-3 ; eta_12 = 1.29e-3 ; eta_44 = 0.16e-3  # Pa s
SiO2_props = [n, s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]

# Use all specified parameters to create a waveguide object.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        bkg_material=materials.Air,
                        inc_a_material=materials.Material(SiO2_props),
                        lc_bkg=3, lc2=2000.0, lc3=20.0)

# Expected effective index of fundamental guided mode.
n_eff = 1.3

# Calculate Electromagnetic modes.
sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes, n_eff=n_eff)

# Print the wavevectors of EM modes.
print('k_z of EM modes \n', np.round(np.real(sim_EM_wguide.Eig_values), 4))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_wguide.Eig_values*((wl_nm*1e-9)/(2.*np.pi)))
print("n_eff = ", np.round(n_eff_sim, 4))

k_AC = 2*np.real(sim_EM_wguide.Eig_values[0])

shift_Hz = 8e9

# Calculate Acoustic modes.
sim_AC_wguide = wguide.calc_AC_modes(wl_nm, num_AC_modes, k_AC=k_AC,
    EM_sim=sim_EM_wguide, shift_Hz=shift_Hz)

# Print the frequencies of AC modes.
print('Freq of AC modes (GHz) \n', np.round(np.real(sim_AC_wguide.Eig_values)*1e-9, 4))

# Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha = integration.gain_and_qs(
    sim_EM_wguide, sim_AC_wguide, k_AC,
    EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 4  # GHz
freq_max = 13  # GHz
plotting.gain_specta(sim_AC_wguide, SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha, k_AC,
    EM_ival1, EM_ival2, AC_ival, freq_min=freq_min, freq_max=freq_max)
