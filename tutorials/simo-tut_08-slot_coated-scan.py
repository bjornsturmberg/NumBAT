""" Calculate the backward SBS gain spectra of a Si
    slot waveguide containing As2S3 surrounded in SiO2.

    This time include a capping layer of SiO2 and 
    investigate the effect of this layer's thickness.
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

num_EM_modes = 20
num_AC_modes = 40
EM_ival1 = 0
EM_ival2 = EM_ival1
AC_ival = 'All'


coat_y_list = np.linspace(50,200,4)
for coat_y in coat_y_list:
    wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                            inc_b_x =inc_b_x, slab_a_y=slab_a_y, coat_y=coat_y,
                            bkg_material=materials.Air,
                            inc_a_material=materials.As2S3,
                            inc_b_material=materials.Si,
                            slab_a_material=materials.SiO2,
                            coat_material=materials.SiO2,
                            lc_bkg=3, lc2=1500.0, lc3=700.0)

    # Expected effective index of fundamental guided mode.
    n_eff = 2.8

    # Calculate Electromagnetic modes.
    sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes, n_eff=n_eff)

    k_AC = 2*np.real(sim_EM_wguide.Eig_values[0])

    shift_Hz = 4e9

    # Calculate Acoustic modes.
    sim_AC_wguide = wguide.calc_AC_modes(wl_nm, num_AC_modes, k_AC=k_AC,
        EM_sim=sim_EM_wguide, shift_Hz=shift_Hz)

    # plotting.plt_mode_fields(sim_AC_wguide, xlim_min=0.4, xlim_max=0.4, 
    #                           ylim_min=0.7, ylim_max=0.0, EM_AC='AC', add_name='_%i' %int(coat_y))

    set_q_factor = 1000.

    SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha = integration.gain_and_qs(
        sim_EM_wguide, sim_AC_wguide, k_AC,
        EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival, fixed_Q=set_q_factor)

    # Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
    freq_min = 5  # GHz
    freq_max = 15  # GHz
    plotting.gain_specta(sim_AC_wguide, SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha, k_AC,
        EM_ival1, EM_ival2, AC_ival, freq_min=freq_min, freq_max=freq_max, add_name='_%i' %int(coat_y))