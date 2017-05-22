""" Calculate the backward SBS gain spectra of a Si
    slot waveguide containing As2S3 surrounded by Air,
    sitting on a SiO2 substrate.
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
# Current mesh template assume inc_b_y = inc_a_y
slab_a_y = wl_nm

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 40
EM_ival_pump = 0
EM_ival_Stokes = EM_ival_pump
AC_ival = 'All'

wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        inc_b_x =inc_b_x, slab_a_y=slab_a_y,
                        material_a=materials.Air,
                        material_b=materials.As2S3_exp,
                        material_c=materials.SiO2,
                        material_d=materials.Si,
                        lc_bkg=3, lc2=2000.0, lc3=1000.0)
# In this case lc3 is meshing around ribs encasing the slot (the Si)

# Expected effective index of fundamental guided mode.
n_eff = wguide.material_b.n-0.1

# Calculate Electromagnetic modes.
sim_EM_pump = wguide.calc_EM_modes(wl_nm, num_modes_EM_pump, n_eff=n_eff)
# np.savez('wguide_data', sim_EM_pump=sim_EM_pump)
# npzfile = np.load('wguide_data.npz')
# sim_EM_pump = npzfile['sim_EM_pump'].tolist()

sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)
# np.savez('wguide_data2', sim_EM_Stokes=sim_EM_Stokes)
# npzfile = np.load('wguide_data2.npz')
# sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()

# plotting.plt_mode_fields(sim_EM_pump, xlim_min=0.4, xlim_max=0.4, 
#                           ylim_min=0.1, ylim_max=0.8, EM_AC='EM', add_name='slot')

# Print the wavevectors of EM modes.
print('k_z of EM modes \n', np.round(np.real(sim_EM_pump.Eig_values), 4))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.Eig_values*((wl_nm*1e-9)/(2.*np.pi)))
print("n_eff = ", np.round(n_eff_sim, 4))

k_AC = np.real(sim_EM_pump.Eig_values[0] - sim_EM_Stokes.Eig_values[0])

# Specify the expected acoustic frequency (slightly low balled).
shift_Hz = 4e9

# Calculate Acoustic modes.
sim_AC_wguide = wguide.calc_AC_modes(wl_nm, num_modes_AC, k_AC=k_AC,
    EM_sim=sim_EM_pump, shift_Hz=shift_Hz)
# np.savez('wguide_data_AC', sim_AC_wguide=sim_AC_wguide)
# npzfile = np.load('wguide_data_AC.npz')
# sim_AC_wguide = npzfile['sim_AC_wguide'].tolist()

# plotting.plt_mode_fields(sim_AC_wguide, xlim_min=0.4, xlim_max=0.4, 
#                           ylim_min=0.7, ylim_max=0.0, EM_AC='AC', add_name='slot')

# Print the frequencies of AC modes.
print('Freq of AC modes (GHz) \n', np.round(np.real(sim_AC_wguide.Eig_values)*1e-9, 4))

set_q_factor = 1000.

SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha, Q_factors = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC_wguide, k_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_q_factor)
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
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max, add_name='_slot')