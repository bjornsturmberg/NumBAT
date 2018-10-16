""" Replicating the results of
    On-chip inter-modal Brillouin scattering
    Kittlaus et al.
    http://dx.doi.org/10.1038/ncomms15819
"""

import time
import datetime
import numpy as np
import sys
import copy

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
unitcell_x = 7*wl_nm
unitcell_y = 0.7*unitcell_x
# Waveguide widths.
inc_a_x = 1500
inc_a_y = 80
# Shape of the waveguide.
# Use double coated geometry to control meshing around rib waveguide.
inc_shape = 'rib_double_coated'

slab_a_x = 2850
slab_a_y = 135

# areas included purely
slab_b_y = 100
coat_x = 50 
coat_y = 100
coat2_x = 100
coat2_y = 200
lc_bkg = 4  # background
lc2 = 8000  # edge of rib
lc3 = 3000   # edge of slab_a 
lc4 = 50    # edge of coat
lc5 = 20    # edge of slab_b
lc6 = 4     # edge of coat2

# Number of electromagnetic modes to solve for.
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
# Number of acoustic modes to solve for.
num_modes_AC = 35
# The EM pump mode(s) for which to calculate interaction with AC modes.
# Can specify a mode number (zero has lowest propagation constant) or 'All'.
EM_ival_pump = 0
# The EM Stokes mode(s) for which to calculate interaction with AC modes.
EM_ival_Stokes = 1 # INTERMODE SBS TE0 to TE1
# The AC mode(s) for which to calculate interaction with EM modes.
AC_ival = 'All'

# Si_110 = copy.deepcopy(materials.Si_2015_Van_Laer)
Si_110 = copy.deepcopy(materials.Si_2016_Smith)
Si_110.rotate_axis(np.pi/4,'z-axis', save_rotated_tensors=True)


print("c_11", Si_110.c_11)
print("c_12", Si_110.c_12)
print("c_13", Si_110.c_13)
print("c_14", Si_110.c_14)
print("c_15", Si_110.c_15)
print("c_16", Si_110.c_16)
print("c_21", Si_110.c_21)
print("c_22", Si_110.c_22)
print("c_23", Si_110.c_23)
print("c_24", Si_110.c_24)
print("c_25", Si_110.c_25)
print("c_26", Si_110.c_26)
print("c_31", Si_110.c_31)
print("c_32", Si_110.c_32)
print("c_33", Si_110.c_33)
print("c_34", Si_110.c_34)
print("c_35", Si_110.c_35)
print("c_36", Si_110.c_36)
print("c_41", Si_110.c_41)
print("c_42", Si_110.c_42)
print("c_43", Si_110.c_43)
print("c_44", Si_110.c_44)
print("c_45", Si_110.c_45)
print("c_46", Si_110.c_46)
print("c_51", Si_110.c_51)
print("c_52", Si_110.c_52)
print("c_53", Si_110.c_53)
print("c_54", Si_110.c_54)
print("c_55", Si_110.c_55)
print("c_56", Si_110.c_56)
print("c_61", Si_110.c_61)
print("c_62", Si_110.c_62)
print("c_63", Si_110.c_63)
print("c_64", Si_110.c_64)
print("c_65", Si_110.c_65)
print("c_66", Si_110.c_66)

prefix_str = 'lit_08-'

# Use specified parameters to create a waveguide object.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        slab_a_x=slab_a_x, slab_a_y=slab_a_y, slab_b_y=slab_b_y, 
                        coat_x=coat_x, coat_y=coat_y, coat2_x=coat2_x, coat2_y=coat2_y,
                        material_bkg=materials.Vacuum,
                        material_a=Si_110, #plt_mesh=True,
                        material_b=Si_110, material_c=materials.Vacuum,
                        material_d=materials.Vacuum, material_e=materials.Vacuum,
                        symmetry_flag=False,
                        lc_bkg=lc_bkg, lc2=lc2, lc3=lc3,
                        lc4=lc4, lc5=lc5, lc6=lc6)
# Expected effective index of fundamental guided mode.
n_eff = wguide.material_a.n-0.1

# Calculate Electromagnetic Modes
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)
# np.savez('wguide_data', sim_EM_pump=sim_EM_pump)
# npzfile = np.load('wguide_data.npz')
# sim_EM_pump = npzfile['sim_EM_pump'].tolist()

sim_EM_Stokes = mode_calcs.fwd_Stokes_modes(sim_EM_pump)
# np.savez('wguide_data2', sim_EM_Stokes=sim_EM_Stokes)
# npzfile = np.load('wguide_data2.npz')
# sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()

# plotting.plt_mode_fields(sim_EM_pump, xlim_min=0.35, xlim_max=0.35, ivals=[0,1], 
#                          ylim_min=0.3, ylim_max=0.3, EM_AC='EM_E', num_ticks=3,
#                          prefix_str=prefix_str, pdf_png='png')

# Print the wavevectors of EM modes.
print('k_z of EM modes \n', np.round(np.real(sim_EM_pump.Eig_values), 4))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.Eig_values*((wl_nm*1e-9)/(2.*np.pi)))
print("n_eff = ", np.round(n_eff_sim, 4))

k_AC = np.real(sim_EM_pump.Eig_values[EM_ival_pump] - sim_EM_Stokes.Eig_values[EM_ival_Stokes])
print('Intermode q_AC (Hz) \n', k_AC)

shift_Hz = 2e9

# Calculate Acoustic Modes
sim_AC = wguide.calc_AC_modes(num_modes_AC, k_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)
# np.savez('wguide_data_AC', sim_AC=sim_AC)
# npzfile = np.load('wguide_data_AC.npz')
# sim_AC = npzfile['sim_AC'].tolist()

# Print the frequencies of AC modes.
print('Freq of AC modes (GHz) \n', np.round(np.real(sim_AC.Eig_values)*1e-9, 4))

# plotting.plt_mode_fields(sim_AC, EM_AC='AC', prefix_str=prefix_str,
#      num_ticks=3, xlim_min=0.1, xlim_max=0.1, pdf_png='png')

set_q_factor = 460.

# Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_q_factor)

# Mask negligible gain values to improve clarity of print out.
threshold = 1e-3
masked_PE = np.ma.masked_inside(SBS_gain_PE[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
masked_MB = np.ma.masked_inside(SBS_gain_MB[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
masked = np.ma.masked_inside(SBS_gain[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)

print("\n SBS_gain PE contribution \n", masked_PE)
print("SBS_gain MB contribution \n", masked_MB)
print("SBS_gain total \n", masked)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 0.5  # GHz
freq_max = 9.5  # GHz
plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    prefix_str=prefix_str, suffix_str='', pdf_png='pdf')

end = time.time()
print("\n Simulation time (sec.)", (end - start))