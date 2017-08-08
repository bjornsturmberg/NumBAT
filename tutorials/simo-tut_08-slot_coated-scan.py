""" Calculate the backward SBS gain spectra of a Si
    slot waveguide containing As2S3 surrounded in SiO2.

    This time include a capping layer of SiO2 and 
    investigate the effect of this layer's thickness.
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

start = time.time()

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


coat_y_list = np.linspace(50,200,4)
for coat_y in coat_y_list:
    wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                            inc_b_x =inc_b_x, slab_a_y=slab_a_y, coat_y=coat_y,
                            material_bkg=materials.Vacuum,
                            material_a=materials.As2S3_2017_Morrison,
                            material_b=materials.SiO2_2013_Laude,
                            material_c=materials.Si_2016_Smith,
                            material_d=materials.SiO2_2013_Laude,
                            lc_bkg=3, lc2=3000.0, lc3=2000.0)

    # Expected effective index of fundamental guided mode.
    n_eff = wguide.material_a.n-0.1

    # Calculate Electromagnetic modes.
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)
    sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

    k_AC = np.real(sim_EM_pump.Eig_values[0] - sim_EM_Stokes.Eig_values[0])

    shift_Hz = 4e9

    # Calculate Acoustic modes.
    sim_AC = wguide.calc_AC_modes(num_modes_AC, k_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)

    # plotting.plt_mode_fields(sim_AC, xlim_min=0.4, xlim_max=0.4, 
    #                           ylim_min=0.7, ylim_max=0.0, EM_AC='AC', 
    #                           prefix_str='tut_08-', suffix_str='_%i' %int(coat_y))

    set_q_factor = 1000.

    SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
        sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC,
        EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_q_factor)

    # Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
    freq_min = np.real(sim_AC.Eig_values[0])*1e-9 - 2  # GHz
    freq_max = np.real(sim_AC.Eig_values[-1])*1e-9 + 2  # GHz
    
    plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
        EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max, 
        prefix_str='tut_08-', suffix_str='_%i' %int(coat_y))
    

end = time.time()
print("\n Simulation time (sec.)", (end - start))
