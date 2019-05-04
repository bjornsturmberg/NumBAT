""" Replicating the results of
    Compact Brillouin devices through hybrid 
    integration on Silicon
    Morrison et al.
    https://doi.org/10.1364/OPTICA.4.000847
"""

import time
import datetime
import numpy as np
import sys

sys.path.append("../backend/")
import materials
import objects
import mode_calcs
import integration
import plotting
from fortran import NumBAT

# Naming conventions
# AC: acoustic
# EM: electromagnetic
# k_AC: acoustic wavenumber

start = time.time()

# Geometric Parameters - all in nm.
wl_nm = 1550 # Wavelength of EM wave in vacuum.
# Unit cell must be large to ensure fields are zero at boundary.
unitcell_x = 6*wl_nm
unitcell_y = 0.75*unitcell_x
# Waveguide widths.
inc_a_x = 1900
inc_a_y = 680
# Shape of the waveguide.
inc_shape = 'rib_coated'

slab_a_x = 4000
slab_a_y = 1000

coat_x = 200 
coat_y = 1000


# Number of electromagnetic modes to solve for.
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
# Number of acoustic modes to solve for.
num_modes_AC = 100
# The EM pump mode(s) for which to calculate interaction with AC modes.
# Can specify a mode number (zero has lowest propagation constant) or 'All'.
EM_ival_pump = 0
# The EM Stokes mode(s) for which to calculate interaction with AC modes.
EM_ival_Stokes = 0
# The AC mode(s) for which to calculate interaction with EM modes.
AC_ival = 'All'

prefix_str = 'lit_09-'

# Use specified parameters to create a waveguide object.
# Note use of rough mesh for demonstration purposes.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        slab_a_x=slab_a_x, slab_a_y=slab_a_y, coat_x=coat_x, coat_y=coat_y,
                        material_bkg=materials.Vacuum,
                        material_a=materials.As2S3_2017_Morrison, # waveguide
                        material_b=materials.SiO2_2016_Smith,     # slab
                        material_c=materials.SiO2_2016_Smith,     # coating
                        lc_bkg=1, lc2=800.0, lc3=400.0)

# Expected effective index of fundamental guided mode.
n_eff = wguide.material_a.n-0.1

# Calculate the Electromagnetic modes of the pump field.
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff)
# # np.savez('wguide_data', sim_EM_pump=sim_EM_pump)
# npzfile = np.load('wguide_data.npz')
# sim_EM_pump = npzfile['sim_EM_pump'].tolist()

plotting.plt_mode_fields(sim_EM_pump, xlim_min=0.4, xlim_max=0.4, ivals=[0], 
                         ylim_min=0.3, ylim_max=0.3, EM_AC='EM_E', num_ticks=3,
                         prefix_str=prefix_str, pdf_png='png')

# Calculate the Electromagnetic modes of the Stokes field.
sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)
# np.savez('wguide_data2', sim_EM_Stokes=sim_EM_Stokes)
# npzfile = np.load('wguide_data2.npz')
# sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()

# Print the wavevectors of EM modes.
print('\n k_z of EM modes \n', np.round(np.real(sim_EM_pump.Eig_values),4))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.Eig_values[0]*((wl_nm*1e-9)/(2.*np.pi)))
print("\n n_eff = ", np.round(n_eff_sim, 4))

k_AC = np.real(sim_EM_pump.Eig_values[EM_ival_pump] - sim_EM_Stokes.Eig_values[EM_ival_Stokes])
print('\n AC wavenumber (1/m) = ', np.round(k_AC, 4))


k_AC= 2.*9173922.1698

# Calculate Acoustic modes.
shift_Hz = 7.5*1e9 # select the lowest frequency to start FEM search from.
sim_AC = wguide.calc_AC_modes(num_modes_AC, k_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)
# # np.savez('wguide_data_AC', sim_AC=sim_AC)
# npzfile = np.load('wguide_data_AC.npz')
# sim_AC = npzfile['sim_AC'].tolist()

plotting.plt_mode_fields(sim_AC, EM_AC='AC', prefix_str=prefix_str,
     num_ticks=3, xlim_min=0.1, xlim_max=0.1, pdf_png='png')

# Print the frequencies of AC modes.
print('\n Freq of AC modes (GHz) \n', np.round(np.real(sim_AC.Eig_values)*1e-9, 4))

set_Q_factor = 190 # set the mechanic Q manually

# Calculate interaction integrals and SBS gain for PE and MB effects combined,
# as well as just for PE, and just for MB. Also calculate acoustic loss alpha.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_Q_factor)
# Print the Backward SBS gain of the AC modes.
print("\n SBS_gain PE contribution \n", SBS_gain_PE[EM_ival_pump,EM_ival_Stokes,:])
print("SBS_gain MB contribution \n", SBS_gain_MB[EM_ival_pump,EM_ival_Stokes,:])
print("SBS_gain total \n", SBS_gain[EM_ival_pump,EM_ival_Stokes,:])
# Mask negligible gain values to improve clarity of print out.
threshold = -1e-3
masked_PE = np.ma.masked_inside(SBS_gain_PE[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
masked_MB = np.ma.masked_inside(SBS_gain_MB[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
masked = np.ma.masked_inside(SBS_gain[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
print("\n SBS_gain PE contribution \n", masked_PE)
print("SBS_gain MB contribution \n", masked_MB)
print("SBS_gain total \n", masked)


# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 7.2  # GHz
freq_max = 8.1  # GHz
plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    prefix_str=prefix_str, pdf_png='png')

end = time.time()
print("\n Simulation time (sec.)", (end - start))