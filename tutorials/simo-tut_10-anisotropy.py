""" Calculate the backward SBS gain spectra of a
    silicon waveguide surrounded in air.

    Show how to save simulation objects (eg. EM mode calcs)
    to expedite the process of altering later parts of
    simulations.
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


# Geometric Parameters - all in nm.
wl_nm = 1550
unitcell_x = 2.5*wl_nm
unitcell_y = unitcell_x
inc_a_x = 314.7
inc_a_y = 0.9*inc_a_x
inc_shape = 'rectangular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 20
EM_ival_pump = 0
EM_ival_Stokes = EM_ival_pump
AC_ival = 'All'

# Anisotropic Acoustic Parameters
n = 3.48
# Density
s = 2329  # kg/m3
# Stiffness tensor components in Pa
c_11 = 165.6e9; c_12 = 63.9e9; c_13 = 63.9e9; c_14 = 0.0; c_15 = 0.0; c_16 = 0.0 
c_21 = 63.9e9; c_22 = 165.6e9; c_23 = 63.9e9; c_24 = 0.0; c_25 = 0.0; c_26 = 0.0 
c_31 = 63.9e9; c_32 = 63.9e9; c_33 = 165.6e9; c_34 = 0.0; c_35 = 0.0; c_36 = 0.0 
c_41 = 0.0; c_42 = 0.0; c_43 = 0.0; c_44 = 79.5e9; c_45 = 0.0; c_46 = 0.0 
c_51 = 0.0; c_52 = 0.0; c_53 = 0.0; c_54 = 0.0; c_55 = 79.5e9; c_56 = 0.0 
c_61 = 0.0; c_62 = 0.0; c_63 = 0.0; c_64 = 0.0; c_65 = 0.0; c_66 = 79.5e9 
# Photoelastic tensor components
p_11 = -0.094; p_12 = 0.017; p_13 = 0.017; p_14 = 0.0; p_15 = 0.0; p_16 = 0.0 
p_21 = 0.017; p_22 = -0.094; p_23 = 0.017; p_24 = 0.0; p_25 = 0.0; p_26 = 0.0
p_31 = 0.017; p_32 = 0.017; p_33 = -0.094; p_34 = 0.0; p_35 = 0.0; p_36 = 0.0 
p_41 = 0.0; p_42 = 0.0; p_43 = 0.0; p_44 = -0.051; p_45 = 0.0; p_46 = 0.0
p_51 = 0.0; p_52 = 0.0; p_53 = 0.0; p_54 = 0.0; p_55 = -0.051; p_56 = 0.0 
p_61 = 0.0; p_62 = 0.0; p_63 = 0.0; p_64 = 0.0; p_65 = 0.0; p_66 = -0.051
# Acoustic loss tensor components in Pa
eta_11 = 5.9e-3; eta_12 = 5.16e-3; eta_13 = 5.16e-3; eta_14 = 0.0; eta_15 = 0.0; eta_16 = 0.0
eta_21 = 5.16e-3; eta_22 = 5.9e-3; eta_23 = 5.16e-3; eta_24 = 0.0; eta_25 = 0.0; eta_26 = 0.0 
eta_31 = 5.16e-3; eta_32 = 5.16e-3; eta_33 = 5.9e-3; eta_34 = 0.0; eta_35 = 0.0; eta_36 = 0.0
eta_41 = 0.0; eta_42 = 0.0; eta_43 = 0.0; eta_44 = 0.620e-3; eta_45 = 0.0; eta_46 = 0.0 
eta_51 = 0.0; eta_52 = 0.0; eta_53 = 0.0; eta_54 = 0.0; eta_55 = 0.620e-3; eta_56 = 0.0
eta_61 = 0.0; eta_62 = 0.0; eta_63 = 0.0; eta_64 = 0.0; eta_65 = 0.0; eta_66 = 0.620e-3

# Put acoustic parameters together for convenience.
test_props = [n, s, 
c_11, c_12, c_13, c_14, c_15, c_16, c_21, c_22, c_23, c_24, c_25, c_26, 
c_31, c_32, c_33, c_34, c_35, c_36, c_41, c_42, c_43, c_44, c_45, c_46, 
c_51, c_52, c_53, c_54, c_55, c_56, c_61, c_62, c_63, c_64, c_65, c_66, 
p_11, p_12, p_13, p_14, p_15, p_16, p_21, p_22, p_23, p_24, p_25, p_26, 
p_31, p_32, p_33, p_34, p_35, p_36, p_41, p_42, p_43, p_44, p_45, p_46, 
p_51, p_52, p_53, p_54, p_55, p_56, p_61, p_62, p_63, p_64, p_65, p_66, 
eta_11, eta_12, eta_13, eta_14, eta_15, eta_16, eta_21, eta_22, eta_23, eta_24, eta_25, eta_26, 
eta_31, eta_32, eta_33, eta_34, eta_35, eta_36, eta_41, eta_42, eta_43, eta_44, eta_45, eta_46, 
eta_51, eta_52, eta_53, eta_54, eta_55, eta_56, eta_61, eta_62, eta_63, eta_64, eta_65, eta_66]

# Use of a more refined mesh to produce field plots.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        material_a=materials.Air,
                        # material_b=materials.Si,
                        material_b=materials.Material(test_props),
                        symmetry_flag=False,
                        lc_bkg=3, lc2=2000.0, lc3=1000.0)


# Expected effective index of fundamental guided mode.
n_eff = wguide.material_b.n-0.1

# Calculate Electromagnetic modes.
sim_EM_pump = wguide.calc_EM_modes(wl_nm, num_modes_EM_pump, n_eff)
# np.savez('wguide_data', sim_EM_pump=sim_EM_pump)
# npzfile = np.load('wguide_data.npz')
# sim_EM_pump = npzfile['sim_EM_pump'].tolist()
sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)
# np.savez('wguide_data2', sim_EM_Stokes=sim_EM_Stokes)
# npzfile = np.load('wguide_data2.npz')
# sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()

# Print the wavevectors of EM modes.
print('k_z of EM modes \n', np.round(np.real(sim_EM_pump.Eig_values), 4))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.Eig_values[0]*((wl_nm*1e-9)/(2.*np.pi)))
print("n_eff", np.round(n_eff_sim, 4))

# Choose acoustic wavenumber to solve for backward SBS
k_AC = np.real(sim_EM_pump.Eig_values[0] - sim_EM_Stokes.Eig_values[0])

# Calculate Acoustic modes.
sim_AC = wguide.calc_AC_modes(wl_nm, num_modes_AC, k_AC, EM_sim=sim_EM_pump)
# np.savez('wguide_data_AC', sim_AC=sim_AC)
# npzfile = np.load('wguide_data_AC.npz')
# sim_AC = npzfile['sim_AC'].tolist()

# Print the frequencies of AC modes.
print('Freq of AC modes (GHz) \n', np.round(np.real(sim_AC.Eig_values)*1e-9, 4))

plotting.plt_mode_fields(sim_AC, EM_AC='AC', pdf_png='png')

# set_q_factor = 1000.

# Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha, Q_factors = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)#, fixed_Q=set_q_factor)

# Print the Backward SBS gain of the AC modes.
print("\n SBS_gain PE contribution \n", SBS_gain_PE[EM_ival_Stokes,EM_ival_pump,:]/alpha)
print("SBS_gain MB contribution \n", SBS_gain_MB[EM_ival_Stokes,EM_ival_pump,:]/alpha)
print("SBS_gain total \n", SBS_gain[EM_ival_Stokes,EM_ival_pump,:]/alpha)
# Mask negligible gain values to improve clarity of print out.
threshold = 1e-3
masked_PE = np.ma.masked_inside(SBS_gain_PE[EM_ival_Stokes,EM_ival_pump,:]/alpha, 0, threshold)
masked_MB = np.ma.masked_inside(SBS_gain_MB[EM_ival_Stokes,EM_ival_pump,:]/alpha, 0, threshold)
masked = np.ma.masked_inside(SBS_gain[EM_ival_Stokes,EM_ival_pump,:]/alpha, 0, threshold)
print("\n SBS_gain PE contribution \n", masked_PE)
print("SBS_gain MB contribution \n", masked_MB)
print("SBS_gain total \n", masked)