""" 
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
unitcell_x = 4*wl_nm
unitcell_y = unitcell_x
inc_a_x = 150
inc_a_y = 190
inc_shape = 'slot'
inc_b_x = 250
# Assumed in current mesh template
# inc_b_y = inc_a_y
slab_a_y = wl_nm

# Optical Parameters
n_As2S3 = 2.44
n_Si = 3.5
n_SiO2 = 1.44
num_EM_modes = 20
num_AC_modes = 40
EM_ival1=0
EM_ival2=EM_ival1
AC_ival='All'


# As2S3
s = 3210  # kg/m3
c_11 = 2.10e10; c_12 = 8.36e9; c_44 =6.34e9 # Pa
p_11 = 0.25; p_12 = 0.24; p_44 = 0.005
eta_11 = 9e-3 ; eta_12 = 7.5e-3 ; eta_44 = 0.75e-3  # Pa s
inc_a_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]
#  Silicon
s = 2330  # kg/m3
c_11 = 165.7e9; c_12 = 63.9e9; c_44 = 79.6e9  # Pa
p_11 = -0.09; p_12 = 0.017; p_44 = -0.051
eta_11 = 5.9e-3 ; eta_12 = 5.16e-3 ; eta_44 = 0.620e-3  # Pa
inc_b_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]
# Silca
s = 2203  # kg/m3
c_11 = 78.5e9; c_12 = 16.1e9; c_44 = 31.2e9
p_11 = 0.121; p_12 = 0.270; p_44 = -0.075
eta_11 = 1.6e-3 ; eta_12 = 1.29e-3 ; eta_44 = 0.16e-3  # Pa s
slab_a_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]


# Use all specified parameters to create a waveguide object.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        inc_b_x =inc_b_x, slab_a_y=slab_a_y,
                        bkg_material=materials.Material(1.0),
                        inc_a_material=materials.Material(n_As2S3),
                        inc_b_material=materials.Material(n_Si),
                        slab_a_material=materials.Material(n_SiO2),
                        loss=False, inc_a_AC=inc_a_AC_props,
                        inc_b_AC=inc_b_AC_props, slab_a_AC=slab_a_AC_props,
                        lc_bkg=3, lc2=2000.0, lc3=1000.0)
# In this case lc3 is meshing around ribs encasing the slot (the Si)

# Expected effective index of fundamental guided mode.
n_eff=2.8

# Calculate Electromagnetic Modes
sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes, n_eff=n_eff)
# np.savez('wguide_data', sim_EM_wguide=sim_EM_wguide)
# npzfile = np.load('wguide_data.npz')
# sim_EM_wguide = npzfile['sim_EM_wguide'].tolist()

# plotting.plt_mode_fields(sim_EM_wguide, xlim_min=0.4, xlim_max=0.4, 
#                           ylim_min=0.1, ylim_max=0.8, EM_AC='EM', add_name='slot')

# Print the wavevectors of EM modes.
print 'k_z of EM modes \n', np.round(np.real(sim_EM_wguide.Eig_values), 4)

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_wguide.Eig_values*((wl_nm*1e-9)/(2.*np.pi)))
print "n_eff = ", np.round(n_eff_sim, 4)

k_AC = 2*np.real(sim_EM_wguide.Eig_values[0])

shift_Hz = 4e9

# Calculate Acoustic Modes
sim_AC_wguide = wguide.calc_AC_modes(wl_nm, num_AC_modes, k_AC=k_AC,
    EM_sim=sim_EM_wguide, shift_Hz=shift_Hz)
# np.savez('wguide_data_AC', sim_AC_wguide=sim_AC_wguide)
# npzfile = np.load('wguide_data_AC.npz')
# sim_AC_wguide = npzfile['sim_AC_wguide'].tolist()

# plotting.plt_mode_fields(sim_AC_wguide, xlim_min=0.4, xlim_max=0.4, 
#                           ylim_min=0.7, ylim_max=0.0, EM_AC='AC', add_name='slot')

# Print the frequencies of AC modes.
print 'Freq of AC modes (GHz) \n', np.round(np.real(sim_AC_wguide.Eig_values)*1e-9, 4)

# Do not calculate the acoustic loss from our fields, but instead set a 
# predetirmined Q factor. (Useful for instance when replicating others results).
set_q_factor = 1000.

# Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# as well as just for PE, and just for MB. Also calculate acoustic loss alpha.
SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha = integration.gain_and_qs(
    sim_EM_wguide, sim_AC_wguide, k_AC,
    EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival, fixed_Q=set_q_factor)
# np.savez('wguide_data_AC_gain', SBS_gain=SBS_gain, SBS_gain_PE=SBS_gain_PE, SBS_gain_MB=SBS_gain_MB, alpha=alpha)
# npzfile = np.load('wguide_data_AC_gain.npz')
# SBS_gain = npzfile['SBS_gain']
# SBS_gain_PE = npzfile['SBS_gain_PE']
# SBS_gain_MB = npzfile['SBS_gain_MB']
# alpha = npzfile['alpha']

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5  # GHz
freq_max = 10  # GHz
plotting.gain_specta(sim_AC_wguide, SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha, k_AC,
    EM_ival1, EM_ival2, AC_ival, freq_min=freq_min, freq_max=freq_max, add_name='_Sayyed')