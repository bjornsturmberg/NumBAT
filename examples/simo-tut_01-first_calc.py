""" Calculate the backward SBS gain for modes in a
    silicon waveguide surrounded in air.
"""

import time
import datetime
import numpy as np
import sys
from multiprocessing import Pool
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

# Geometric Parameters - all in nm.
wl_nm = 1550 # Wavelength of EM wave in vacuum.
# Unit cell must be large to ensure fields are zero at boundary.
unitcell_x = 2.5*1550
unitcell_y = unitcell_x
# Waveguide width (x direction).
inc_a_x = 314.7
# Waveguide height (y direction).
inc_a_y = 0.9*inc_a_x
# Shape of the waveguide could also be 'circular'.
inc_shape = 'rectangular'

# Optical Parameters
# Permittivity
eps = 12.25
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
# Density
s = 2330  # kg/m3
# Stiffness tensor components.
c_11 = 165.7e9; c_12 = 63.9e9; c_44 = 79.6e9  # Pa
# Photoelastic tensor components
p_11 = -0.094; p_12 = 0.017; p_44 = -0.051
# Acoustic loss tensor components.
eta_11 = 5.9e-3 ; eta_12 = 5.16e-3 ; eta_44 = 0.620e-3  # Pa
# Put acoustic parameters together for convenience.
inc_a_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]

# Use all specified parameters to create a waveguide object.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        bkg_material=materials.Material(1.0 + 0.0j),
                        inc_a_material=materials.Material(np.sqrt(eps)),
                        loss=False, inc_a_AC=inc_a_AC_props,
                        lc_bkg=0.1, lc2=40.0, lc3=20.0)


# Calculate Electromagnetic Modes
sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes)
# Print the wavevectors of EM modes.
print 'k_z of EM wave \n', sim_EM_wguide.Eig_value
# Plot the EM modes fields, important to specify this with EM_AC='EM'.
# Zoom in on the central region (of big unitcell) with xlim, ylim args.
plotting.plt_mode_fields(sim_EM_wguide, xlim=0.4, ylim=0.4, EM_AC='EM')


# Choose acoustic wavenumber to solve for
# Backward SBS
# AC mode couples EM modes on +ve to -ve lightline, hence factor 2.
k_AC = 2*np.real(sim_EM_wguide.Eig_value[0])
print 'AC wavenumber (1/m) \n', k_AC
# Forward (intramode) SBS
# EM modes on same lightline.
# k_AC = 0.0


# Calculate Acoustic Modes
sim_AC_wguide = wguide.calc_AC_modes(wl_nm, k_AC,
    num_AC_modes, EM_sim=sim_EM_wguide, shift_Hz=12e9)
# Print the frequencies of AC modes.
print 'Res freq of AC wave (GHz) \n', np.real(sim_AC_wguide.Eig_value)*1e-9
# Plot the AC modes fields, important to specify this with EM_AC='AC'.
# The AC modes are calculated on a subset of the full unitcell,
# which excludes vacuum regions, so no need to restrict area plotted.
plotting.plt_mode_fields(sim_AC_wguide, EM_AC='AC')


# Calculate interaction integrals and SBS gain
# - for PE and MB effects combined, as well as
# just for PE, and just for MB. Also calculate acoustic loss alpha.
SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha = integration.gain_and_qs(
    sim_EM_wguide, sim_AC_wguide, k_AC,
    EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival)
# Print the Backward SBS gain of the AC modes.
print "SBS_gain \n", SBS_gain[EM_ival1,EM_ival2,:]/alpha
