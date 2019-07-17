"""
    test_case_0_rect_silicon.py is a simulation example for NumBAT.

    Copyright (C) 2015  Bjorn Sturmberg, Kokou Dossou.

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

"""
Test simulation of a simple rectangular waveguide made of silicon.
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

from numpy.testing import assert_allclose as assert_ac
from numpy.testing import assert_equal

start = time.time()

casefile_name = 'case_0'

# Geometric Parameters - all in nm.
wl_nm = 1550 # Wavelength of EM wave in vacuum.
# Unit cell must be large to ensure fields are zero at boundary.
unitcell_x = 2.5*wl_nm
unitcell_y = unitcell_x
# Waveguide width (x direction).
inc_a_x = 314.7
# Waveguide height (y direction).
inc_a_y = 0.9*inc_a_x
# Shape of the waveguide could also be 'circular'.
inc_shape = 'rectangular'

# Number of electromagnetic modes to solve for.
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
# Number of acoustic modes to solve for.
num_AC_modes = 20
# The EM pump mode(s) for which to calculate interaction with AC modes.
# Can specify a mode number (zero has lowest propagation constant) or 'All'.
EM_ival_pump = 0
# The EM Stokes mode(s) for which to calculate interaction with AC modes.
EM_ival_Stokes = 0
# The AC mode(s) for which to calculate interaction with EM modes.
AC_ival='All'

# Use all specified parameters to create a waveguide object.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        material_bkg=materials.Vacuum,
                        material_a=materials.Si_2016_Smith,
                        lc_bkg=1, lc2=1000.0, lc3=400.0)

# Expected effective index of fundamental guided mode.
n_eff = wguide.material_a.n-0.1

# Calculate Electromagnetic Modes
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)

sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

k_AC = np.real(sim_EM_pump.Eig_values[0] - sim_EM_Stokes.Eig_values[0])

# Calculate Acoustic Modes
sim_AC_wguide = wguide.calc_AC_modes(num_AC_modes, k_AC=k_AC, EM_sim=sim_EM_pump)

# Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# as well as just for PE, and just for MB. Also calculate acoustic loss alpha.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC_wguide, k_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)
# Mask negligible gain values to improve clarity of print out.
threshold = 1e-3
threshold_indices = abs(SBS_gain_PE) < threshold
SBS_gain_PE[threshold_indices] = 0
threshold_indices = abs(SBS_gain_MB) < threshold
SBS_gain_MB[threshold_indices] = 0
threshold_indices = abs(SBS_gain) < threshold
SBS_gain[threshold_indices] = 0
masked_PE = SBS_gain_PE[EM_ival_pump,EM_ival_Stokes,:]
masked_MB = SBS_gain_MB[EM_ival_pump,EM_ival_Stokes,:]
masked = SBS_gain[EM_ival_pump,EM_ival_Stokes,:]

test_list1 = list(zip(sim_EM_pump.Eig_values, sim_AC_wguide.Eig_values))
test_list2 = list(zip(masked_PE, masked_MB, masked))

# # SAVE DATA AS REFERENCE
# # Only run this after changing what is simulated - this
# # generates a new set of reference answers to check against
# # in the future
# np.savez_compressed("ref/%s.npz" % casefile_name, 
#         test_list1 = test_list1, test_list2 = test_list2)
# assert False, "Reference results saved successfully, \n tests would pass trivially so we'll skip them."

def test_list_matches_saved(casefile_name = casefile_name):
    rtol = 1e-6
    atol = 1e-6
    ref = np.load("ref/%s.npz" % casefile_name)
    for case, rcase in zip(test_list1, ref['test_list1']):
        yield assert_ac, case, rcase, rtol, atol
    for case, rcase in zip(test_list2, ref['test_list2']):
        yield assert_ac, case, rcase, rtol, atol
