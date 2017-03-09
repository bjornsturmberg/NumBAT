""" Replicating the results of
    Brillouin scattering self-cancellation
    Florez et al.
    http://dx.doi.org/10.1038/ncomms11759
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
unitcell_x = 5*wl_nm
unitcell_y = unitcell_x
inc_a_x = 550
inc_a_y = inc_a_x
inc_shape = 'circular'

num_EM_modes = 20
num_AC_modes = 40
EM_ival1 = 0
EM_ival2 = EM_ival1
AC_ival = 'All'

# Use all specified parameters to create a waveguide object.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        bkg_material=materials.Air,
                        inc_a_material=materials.SiO2,
                        lc_bkg=3, lc2=2000.0, lc3=20.0)

# Expected effective index of fundamental guided mode.
n_eff = 1.4

# Calculate Electromagnetic Modes
sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes, n_eff=n_eff)

plotting.plt_mode_fields(sim_EM_wguide, xlim_min=0.4, xlim_max=0.4, 
                         ylim_min=0.4, ylim_max=0.4, EM_AC='EM', pdf_png='pdf')

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

plotting.plt_mode_fields(sim_AC_wguide, EM_AC='AC', add_name='NW')

# Print the frequencies of AC modes.
print 'Freq of AC modes (GHz) \n', np.round(np.real(sim_AC_wguide.Eig_values)*1e-9, 4)

set_q_factor = 1000.

# Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha = integration.gain_and_qs(
    sim_EM_wguide, sim_AC_wguide, k_AC,
    EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival, fixed_Q=set_q_factor)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 0  # GHz
freq_max = 12  # GHz
plotting.gain_specta(sim_AC_wguide, SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha, k_AC,
    EM_ival1, EM_ival2, AC_ival, freq_min=freq_min, freq_max=freq_max)