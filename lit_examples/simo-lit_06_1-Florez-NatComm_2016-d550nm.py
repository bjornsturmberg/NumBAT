""" Replicating the results of
    Brillouin scattering self-cancellation
    Florez et al.
    http://dx.doi.org/10.1038/ncomms11759
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
from plotting import FieldDecorator
from fortran import NumBAT

# use this class to add or alter features to the final plots
class EMDecorator(FieldDecorator):
  def __init__(self):
    super().__init__()

    #title_font=24
    #self._multi_sizes= {'title':title_font-2, 'subplot_title':title_font-5, 'cbar_tick':title_font-10, 'ax_tick':title_font-10, 'ax_label':title_font-10 }
    #self._single_sizes= {'ax_label':80, 'subplot_title':80, 'cbar_tick':60, 'ax_tick':70}
    ##self._single_sizes= {'ax_label':60, 'subplot_title':60, 'cbar_tick':40, 'ax_tick':40}
    #self._is_single=True

  def extra_axes_commands(self, ax):
    circle1 = plt.Circle((0, 0), 0.275, color='black', fill=False)
    ax.add_artist(circle1)


class ACDecorator(FieldDecorator):
  def __init__(self):
    super().__init__()
    #title_font=24
    #self._multi_sizes= {'title':title_font-2, 'subplot_title':title_font-5, 'cbar_tick':title_font-10, 'ax_tick':title_font-10, 'ax_label':title_font-10 }
#
#    self._single_sizes= {'ax_label':30, 'subplot_title':40, 'cbar_tick':20, 'ax_tick':30, 'title_pad':20}
#    self._single_sizes= {'ax_label':80, 'subplot_title':80, 'cbar_tick':60, 'ax_tick':70, 'title_pad':25}
#    self._is_single=True

  def extra_axes_commands(self, ax):
    circle1 = plt.Circle((0, 0), 0.275, color='black', fill=False)
    ax.add_artist(circle1)

emdecorate=EMDecorator()
acdecorate=ACDecorator()


start = time.time()

# Geometric Parameters - all in nm.
wl_nm = 1550
unitcell_x = 2*wl_nm
unitcell_y = unitcell_x
inc_a_x = 550 # Diameter
inc_a_y = inc_a_x
inc_shape = 'circular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 40
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix_str = 'fig10-'

# Use all specified parameters to create a waveguide object.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        material_bkg=materials.Vacuum,
                        material_a=materials.SiO2_2013_Laude,
                        lc_bkg=.25, lc2=170.0, lc3=85.0)

# Expected effective index of fundamental guided mode.
n_eff = 1.4

doem=True
doac=True
new_calcs=False

if doem:
  # Calculate Electromagnetic Modes
  if new_calcs:
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)
    np.savez('wguide_data_florez', sim_EM_pump=sim_EM_pump)
  else:
    npzfile = np.load('wguide_data_florez.npz', allow_pickle=True)
    sim_EM_pump = npzfile['sim_EM_pump'].tolist()
  sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)
  
  plotting.plt_mode_fields(sim_EM_pump, xlim_min=0.3, xlim_max=0.3, ivals=[0],
                           ylim_min=0.3, ylim_max=0.3, EM_AC='EM_E', 
                           prefix_str=prefix_str, pdf_png='png', ticks=True,
                           decorator=emdecorate, quiver_steps=20)
  
  # Print the wavevectors of EM modes.
  print('k_z of EM modes \n', np.round(np.real(sim_EM_pump.Eig_values), 4))
  
  # Calculate the EM effective index of the waveguide.
  n_eff_sim = np.real(sim_EM_pump.Eig_values*((wl_nm*1e-9)/(2.*np.pi)))
  print("n_eff = ", np.round(n_eff_sim, 4))

if not doac: sys.exit(0)

k_AC = np.real(sim_EM_pump.Eig_values[EM_ival_pump] - sim_EM_Stokes.Eig_values[EM_ival_Stokes])

shift_Hz = 4e9

# Calculate Acoustic Modes
if new_calcs:
  sim_AC = wguide.calc_AC_modes(num_modes_AC, k_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)
  np.savez('wguide_data_florez_AC', sim_AC=sim_AC)
else:
  npzfile = np.load('wguide_data_florez_AC.npz', allow_pickle=True)
  sim_AC = npzfile['sim_AC'].tolist()

plotting.plt_mode_fields(sim_AC, EM_AC='AC', prefix_str=prefix_str, suffix_str='',
ticks=True, ivals=[4,5], comps=('ut','uabs'), 
xlim_min=-.1, ylim_min=-.1, xlim_max=-.1, ylim_max=-.1, 
decorator=acdecorate, quiver_steps=20, pdf_png='png',colorbar=True)

# Print the frequencies of AC modes.
print('Freq of AC modes (GHz) \n', np.round(np.real(sim_AC.Eig_values)*1e-9, 4))

set_q_factor = 1000.

# Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_q_factor)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5  # GHz
freq_max = 12  # GHz
plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    prefix_str=prefix_str, pdf_png='png')

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 5.86  # GHz
freq_max = 5.9  # GHz
plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    prefix_str=prefix_str, suffix_str='-5', pdf_png='png')

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 6.28  # GHz
freq_max = 6.32  # GHz
plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    prefix_str=prefix_str, suffix_str='-6', pdf_png='png')

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 8.09  # GHz
freq_max = 8.13  # GHz
plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    prefix_str=prefix_str, suffix_str='-8', pdf_png='png')

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 11.65  # GHz
freq_max = 11.69  # GHz
plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    prefix_str=prefix_str, suffix_str='-11', pdf_png='png')

end = time.time()
print("\n Simulation time (sec.)", (end - start))
