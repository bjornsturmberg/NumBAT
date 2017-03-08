""" Calculate the backward SBS gain spectra of a
    silicon waveguide on a pedestal surrounded in air.
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
inc_a_x = 450
inc_a_y = 230
inc_shape = 'rectangular'
slab_a_x = 15
slab_a_y = 300
slab_b_x = unitcell_x-400
slab_b_y = 800

# Optical Parameters
n_silicon = 3.5
n_silica = 1.44
num_EM_modes = 20
num_AC_modes = 80
EM_ival1=0
EM_ival2=EM_ival1
AC_ival='All'

# Acoustic Parameters
# Silicon
s = 2330  # kg/m3
c_11 = 166e9; c_12 = 64e9; c_44 = 79e9  # Pa
p_11 = -0.09; p_12 = 0.017; p_44 = -0.051
eta_11 = 5.9e-3 ; eta_12 = 5.16e-3 ; eta_44 = 0.620e-3  # Pa
inc_a_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]


# Silica also
s = 2203  # kg/m3
c_11 = 78e9; c_12 = 16e9; c_44 = 31e9
p_11 = 0.12; p_12 = 0.270; p_44 = -0.073
eta_11 = 1.6e-3 ; eta_12 = 1.29e-3 ; eta_44 = 0.16e-3  # Pa s
slab_a_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]

# Silica
s = 2203  # kg/m3
c_11 = 78e9; c_12 = 16e9; c_44 = 31e9
p_11 = 0.12; p_12 = 0.270; p_44 = -0.073
eta_11 = 1.6e-3 ; eta_12 = 1.29e-3 ; eta_44 = 0.16e-3  # Pa s
slab_b_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]


# Use all specified parameters to create a waveguide object.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        slab_a_x=slab_a_x, slab_a_y=slab_a_y,
                        slab_b_x=slab_b_x, slab_b_y=slab_b_y,
                        bkg_material=materials.Material(1.0 + 0.0j),
                        inc_a_material=materials.Material(n_silicon),
                        slab_a_material=materials.Material(n_silica),
                        slab_a_bkg_material=materials.Material(1.0 + 0.0j),
                        slab_b_material=materials.Material(n_silica),
                        slab_b_bkg_material=materials.Material(1.0 + 0.0j),
                        loss=False, 
                        inc_a_AC=inc_a_AC_props,
                        slab_a_AC = slab_a_AC_props,
                        slab_b_AC = slab_b_AC_props,
                        lc_bkg=2, lc2=4000.0, lc3=20.0)


# Expected effective index of fundamental guided mode.
n_eff = np.real(n_silicon)-0.1

# Calculate Electromagnetic Modes
# sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes, n_eff)
# np.savez('wguide_data', sim_EM_wguide=sim_EM_wguide)
npzfile = np.load('wguide_data.npz')
sim_EM_wguide = npzfile['sim_EM_wguide'].tolist()

# plotting.plt_mode_fields(sim_EM_wguide, 
#                          xlim_min=0.35, xlim_max=0.35, ylim_min=0.1, ylim_max=0.55, 
#                          EM_AC='EM', add_name='slab', pdf_png='pdf')

# Print the wavevectors of EM modes.
print 'k_z of EM modes \n', np.round(np.real(sim_EM_wguide.Eig_values), 4)

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_wguide.Eig_values[0]*((wl_nm*1e-9)/(2.*np.pi)))

k_AC = 2*np.real(sim_EM_wguide.Eig_values[0])
shift_Hz=10e9


# Calculate Acoustic Modes
sim_AC_wguide = wguide.calc_AC_modes(wl_nm, num_AC_modes, k_AC,
    EM_sim=sim_EM_wguide, shift_Hz=shift_Hz)
np.savez('wguide_data_AC', sim_AC_wguide=sim_AC_wguide)
# npzfile = np.load('wguide_data_AC.npz')
# sim_AC_wguide = npzfile['sim_AC_wguide'].tolist()

plotting.plt_mode_fields(sim_AC_wguide, EM_AC='AC', add_name='slab', pdf_png='png')

# # Print the frequencies of AC modes.
# print 'Freq of AC modes (GHz) \n', np.round(np.real(sim_AC_wguide.Eig_values)*1e-9, 4)


# # Do not calculate the acoustic loss from our fields, but instead set a 
# # predetirmined Q factor. (Useful for instance when replicating others results).
# set_q_factor = 400.

# # Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# # as well as just for PE, and just for MB. Also calculate acoustic loss alpha.
# SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha = integration.gain_and_qs(
#     sim_EM_wguide, sim_AC_wguide, k_AC,
#     EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival, fixed_Q=set_q_factor)
# np.savez('wguide_data_AC_gain', SBS_gain=SBS_gain, SBS_gain_PE=SBS_gain_PE, SBS_gain_MB=SBS_gain_MB, alpha=alpha)
# # npzfile = np.load('wguide_data_AC_gain.npz')
# # SBS_gain = npzfile['SBS_gain']
# # SBS_gain_PE = npzfile['SBS_gain_PE']
# # SBS_gain_MB = npzfile['SBS_gain_MB']
# # alpha = npzfile['alpha']

# print SBS_gain[EM_ival1,EM_ival2,:]/alpha

# # Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
# freq_min = 0  # GHz
# freq_max = 20  # GHz
# plotting.gain_specta(sim_AC_wguide, SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha, k_AC,
#     EM_ival1, EM_ival2, AC_ival, freq_min=freq_min, freq_max=freq_max)
# # Zoomed in version
# freq_min = 10  # GHz
# freq_max = 14  # GHz
# plotting.gain_specta(sim_AC_wguide, SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha, k_AC,
#     EM_ival1, EM_ival2, AC_ival, freq_min=freq_min, freq_max=freq_max, add_name='_zoom')
