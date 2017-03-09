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

num_EM_modes = 20
# There are lots of leaky acoustic modes!
num_AC_modes = 350
EM_ival1 = 0
EM_ival2 = EM_ival1
AC_ival = 'All'

# Use all specified parameters to create a waveguide object.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        bkg_material=materials.SiO2,
                        inc_a_material=materials.As2S3,
                        lc_bkg=3, lc2=2000.0, lc3=10.0)

# Expected effective index of fundamental guided mode.
n_eff = 2.

# Calculate Electromagnetic Modes
sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes, n_eff=n_eff)
# np.savez('wguide_data-chalc', sim_EM_wguide=sim_EM_wguide)
# npzfile = np.load('wguide_data-chalc.npz')
# sim_EM_wguide = npzfile['sim_EM_wguide'].tolist()
# plotting.plt_mode_fields(sim_EM_wguide, xlim_min=0.4, xlim_max=0.4, 
#                           ylim_min=0.4, ylim_max=0.4, EM_AC='EM', add_name='As2S3')

# Print the wavevectors of EM modes.
print 'k_z of EM wave \n', np.round(np.real(sim_EM_wguide.Eig_values), 4)
n_eff_sim = np.real((sim_EM_wguide.Eig_values[0]*((wl_nm*1e-9)/(2.*np.pi))))
print 'n_eff of fund. EM mode \n', np.round(n_eff_sim, 4)

k_AC = 2*np.real(sim_EM_wguide.Eig_values[0])

# Calculate Acoustic Modes
sim_AC_wguide = wguide.calc_AC_modes(wl_nm, num_AC_modes, k_AC=k_AC,
    EM_sim=sim_EM_wguide)
# np.savez('wguide_data_AC', sim_AC_wguide=sim_AC_wguide)
# npzfile = np.load('wguide_data_AC.npz')
# sim_AC_wguide = npzfile['sim_AC_wguide'].tolist()
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
