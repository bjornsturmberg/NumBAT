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
unitcell_x = 6.5*wl_nm
unitcell_y = unitcell_x
inc_a_x = 1900
inc_a_y = 680
inc_shape = 'rectangular'

# Optical parameters
n_b = 1.44
n_i = 2.44
num_EM_modes = 20
# There are lots of leaky acoustic modes!
num_AC_modes = 350
EM_ival1=0
EM_ival2=EM_ival1
AC_ival='All'

# Acoustic Parameters
# Background - silca
# # Use the isotropic parameters for silca.
# s = 2200  # kg/m3
# E = 7.3e10
# v = 0.17
# c_11, c_12, c_44 = materials.isotropic_stiffness(E, v)
# p_11 = 0.121; p_12 = 0.270; p_44 = -0.075

# Silca - Laude AIP Advances 2013
s = 2203  # kg/m3
c_11 = 78e9; c_12 = 16e9; c_44 = 31e9
p_11 = 0.12; p_12 = 0.270; p_44 = -0.073

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
                        inc_a_AC=inc_a_AC_props,plotting_fields=False,
                        lc_bkg=3, lc2=2500.0, lc3=10.0)

# Expected effective index of fundamental guided mode.
n_eff = 2.

# Calculate Electromagnetic Modes
# sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes, n_eff=n_eff)
# np.savez('wguide_data-chalc', sim_EM_wguide=sim_EM_wguide)
npzfile = np.load('wguide_data-chalc.npz')
sim_EM_wguide = npzfile['sim_EM_wguide'].tolist()
# plotting.plt_mode_fields(sim_EM_wguide, xlim_min=0.4, xlim_max=0.4, 
#                           ylim_min=0.4, ylim_max=0.4, EM_AC='EM', add_name='As2S3')

# Print the wavevectors of EM modes.
print 'k_z of EM wave \n', np.round(np.real(sim_EM_wguide.Eig_values), 4)
n_eff_sim = np.real((sim_EM_wguide.Eig_values[0]*((wl_nm*1e-9)/(2.*np.pi))))
print 'n_eff of fund. EM mode \n', np.round(n_eff_sim, 4)
n_eff_sim = np.real((sim_EM_wguide.Eig_values[2]*((wl_nm*1e-9)/(2.*np.pi))))
print 'n_eff of 3rd EM mode \n', np.round(n_eff_sim, 4)

# Choose acoustic wavenumber to solve for
# Backward SBS
# AC mode couples EM modes on +ve to -ve lightline, hence factor 2.
k_AC = 2*np.real(sim_EM_wguide.Eig_values[0])

# Calculate Acoustic Modes
# sim_AC_wguide = wguide.calc_AC_modes(wl_nm, num_AC_modes, k_AC=k_AC,
#     EM_sim=sim_EM_wguide)
# np.savez('wguide_data_AC', sim_AC_wguide=sim_AC_wguide)
npzfile = np.load('wguide_data_AC.npz')
sim_AC_wguide = npzfile['sim_AC_wguide'].tolist()
# plotting.plt_mode_fields(sim_AC_wguide, xlim_min=0.4, xlim_max=0.4, 
#                           ylim_min=0.4, ylim_max=0.4, EM_AC='AC', add_name='As2S3')

# Print the frequencies of AC modes.
print 'Res freq of AC wave (GHz) \n', np.round(np.real(sim_AC_wguide.Eig_values*1e-9), 4)

# Experimentaly obtained Q factor.
set_q_factor = 1600.

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
    EM_ival1, EM_ival2, AC_ival, freq_min=freq_min, freq_max=freq_max, add_name='_As2S3')
# Zoomed in version
freq_min = 7.5  # GHz
freq_max = 8.6  # GHz
plotting.gain_specta(sim_AC_wguide, SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha, k_AC,
    EM_ival1, EM_ival2, AC_ival, freq_min=freq_min, freq_max=freq_max, add_name='_As2S3_zoom')
