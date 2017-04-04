"""
    test_case_0_rect_silicon.py is a simulation example for EMUstack.

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
num_EM_modes = 20
# Number of acoustic modes to solve for.
num_AC_modes = 20
# The first EM mode(s) for which to calculate interaction with AC modes.
# Can specify a mode number (zero has lowest propagation constant) or 'All'.
EM_ival1=0
# The second EM mode(s) for which to calculate interaction with AC modes.
EM_ival2=EM_ival1
# The AC mode(s) for which to calculate interaction with EM modes.
AC_ival='All'

# Acoustic Parameters
# Silicon
n = 3.48
# Density
s = 2329  # kg/m3
# Stiffness tensor components.
c_11 = 165.7e9; c_12 = 63.9e9; c_44 = 79.6e9  # Pa
# Photoelastic tensor components
p_11 = -0.094; p_12 = 0.017; p_44 = -0.051
# Acoustic loss tensor components.
eta_11 = 5.9e-3 ; eta_12 = 5.16e-3 ; eta_44 = 0.620e-3  # Pa
# Put acoustic parameters together for convenience.
Si_props = [n, s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]

start = time.time()

# Use all specified parameters to create a waveguide object.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        material_a=materials.Air,
                        material_b=materials.Material(Si_props),
                        lc_bkg=2, lc2=2000.0, lc3=20.0)

# Expected effective index of fundamental guided mode.
n_eff = wguide.material_b.n-0.1

# Calculate Electromagnetic Modes
sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes, n_eff=n_eff)

k_AC = 2*np.real(sim_EM_wguide.Eig_values[0])

# Calculate Acoustic Modes
sim_AC_wguide = wguide.calc_AC_modes(wl_nm, num_AC_modes, 
    k_AC=k_AC, EM_sim=sim_EM_wguide)

# Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# as well as just for PE, and just for MB. Also calculate acoustic loss alpha.
SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha = integration.gain_and_qs(
    sim_EM_wguide, sim_AC_wguide, k_AC,
    EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival)
# Mask negligible gain values to improve clarity of print out.
threshold = 1e-3
masked_PE = np.ma.masked_inside(SBS_gain_PE[EM_ival1,EM_ival2,:]/alpha, 0, threshold)
masked_MB = np.ma.masked_inside(SBS_gain_MB[EM_ival1,EM_ival2,:]/alpha, 0, threshold)
masked = np.ma.masked_inside(SBS_gain[EM_ival1,EM_ival2,:]/alpha, 0, threshold)

# print("\n SBS_gain PE contribution \n", masked_PE)
# print("SBS_gain MB contribution \n", masked_MB)
# print("SBS_gain total \n", masked)

test_list1 = list(zip(sim_EM_wguide.Eig_values, sim_AC_wguide.Eig_values))
test_list2 = list(zip(masked_PE, masked_MB, masked))

# SAVE DATA AS REFERENCE
# Only run this after changing what is simulated - this
# generates a new set of reference answers to check against
# in the future
# np.savez_compressed("ref/%s.npz" % casefile_name, 
#         test_list1 = test_list1, test_list2 = test_list2)
# assert False, "Reference results saved successfully, \
# but tests will now pass trivially so let's not run them now."

def test_list_matches_saved(casefile_name = casefile_name):
    rtol = 1e-6
    atol = 1e-6
    ref = np.load("ref/%s.npz" % casefile_name)
# print(repr(ref))
# print(repr(ref['test_list']))
# print(repr(ref['test_list'][0]))
    for case, rcase in zip(test_list1, ref['test_list1']):
        yield assert_ac, case, rcase, rtol, atol
    for case, rcase in zip(test_list2, ref['test_list2']):
        yield assert_ac, case, rcase, rtol, atol
