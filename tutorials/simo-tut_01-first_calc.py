""" Calculate the backward SBS gain for modes in a
    silicon waveguide surrounded in air.
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
unitcell_x = 2.5*wl_nm
unitcell_y = unitcell_x
# Waveguide widths.
inc_a_x = 314.7
inc_a_y = 0.9*inc_a_x
# Shape of the waveguide.
inc_shape = 'rectangular'

# Number of electromagnetic modes to solve for.
num_EM_modes = 20
# Number of acoustic modes to solve for.
num_AC_modes = 20
# The first EM mode(s) for which to calculate interaction with AC modes.
# Can specify a mode number (zero has lowest propagation constant) or 'All'.
EM_ival1 = 0
# The second EM mode(s) for which to calculate interaction with AC modes.
EM_ival2 = EM_ival1
# The AC mode(s) for which to calculate interaction with EM modes.
AC_ival = 'All'

# Use specified parameters to create a waveguide object.
# Note use of rough mesh for demonstration purposes.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        material_a=materials.Air,
                        material_b=materials.Si,
                        lc_bkg=2, lc2=200.0, lc3=5.0)

# Expected effective index of fundamental guided mode.
n_eff = wguide.material_b.n-0.1

# Calculate Electromagnetic modes.
sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes, n_eff)
# Print the wavevectors of EM modes.
print('\n k_z of EM modes \n', np.round(np.real(sim_EM_wguide.Eig_values),4))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_wguide.Eig_values[0]*((wl_nm*1e-9)/(2.*np.pi)))
print("\n n_eff = ", np.round(n_eff_sim, 4))

# Choose acoustic wavenumber to solve for.
# Backward SBS - AC mode couples EM modes on +ve to -ve lightline, hence factor 2.
k_AC = 2*np.real(sim_EM_wguide.Eig_values[0])
print('\n AC wavenumber (1/m) = ', np.round(k_AC, 4))
# Forward (intramode) SBS - EM modes on same lightline.
# k_AC = 0.0

# Calculate Acoustic modes.
sim_AC_wguide = wguide.calc_AC_modes(wl_nm, num_AC_modes, 
    k_AC=k_AC, EM_sim=sim_EM_wguide)
# Print the frequencies of AC modes.
print('\n Freq of AC modes (GHz) \n', np.round(np.real(sim_AC_wguide.Eig_values)*1e-9, 4))

# Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# as well as just for PE, and just for MB. Also calculate acoustic loss alpha.
SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha = integration.gain_and_qs(
    sim_EM_wguide, sim_AC_wguide, k_AC,
    EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival)
# Print the Backward SBS gain of the AC modes.
print("\n SBS_gain PE contribution \n", SBS_gain_PE[EM_ival1,EM_ival2,:]/alpha)
print("SBS_gain MB contribution \n", SBS_gain_MB[EM_ival1,EM_ival2,:]/alpha)
print("SBS_gain total \n", SBS_gain[EM_ival1,EM_ival2,:]/alpha)
# Mask negligible gain values to improve clarity of print out.
threshold = 1e-3
masked_PE = np.ma.masked_inside(SBS_gain_PE[EM_ival1,EM_ival2,:]/alpha, 0, threshold)
masked_MB = np.ma.masked_inside(SBS_gain_MB[EM_ival1,EM_ival2,:]/alpha, 0, threshold)
masked = np.ma.masked_inside(SBS_gain[EM_ival1,EM_ival2,:]/alpha, 0, threshold)
print("\n SBS_gain PE contribution \n", masked_PE)
print("SBS_gain MB contribution \n", masked_MB)
print("SBS_gain total \n", masked)

end = time.time()
print("\n Simulation time (sec.)", (end - start))
