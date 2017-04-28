""" Replicating the results of
    Interaction between light and highly confined hypersound in a silicon photonic nanowire
    Van Laer et al.
    http://dx.doi.org/10.1038/nphoton.2015.11
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
unitcell_x = 2.5*wl_nm
unitcell_y = unitcell_x
inc_a_x = 450
inc_a_y = 230
inc_shape = 'rectangular'
slab_a_x = 15
slab_a_y = 300
slab_b_x = unitcell_x-400
slab_b_y = 800

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 80
EM_ival_pump=0
EM_ival_Stokes=EM_ival_pump
AC_ival='All'

# Material parameters as in paper 
# Silicon
n = 3.5
s = 2330  # kg/m3
c_11 = 166e9; c_12 = 64e9; c_44 = 79e9  # Pa
p_11 = -0.09; p_12 = 0.017; p_44 = -0.051
eta_11 = 5.9e-3 ; eta_12 = 5.16e-3 ; eta_44 = 0.620e-3  # Pa
Si_props = [n, s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]

# Silica
n = 1.44
s = 2203  # kg/m3
c_11 = 78e9; c_12 = 16e9; c_44 = 31e9
p_11 = 0.12; p_12 = 0.270; p_44 = -0.073
eta_11 = 1.6e-3 ; eta_12 = 1.29e-3 ; eta_44 = 0.16e-3  # Pa s
SiO2_props = [n, s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]

# Use all specified parameters to create a waveguide object.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        slab_a_x=slab_a_x, slab_a_y=slab_a_y,
                        slab_b_x=slab_b_x, slab_b_y=slab_b_y,
                        material_a=materials.Air,
                        material_b=materials.Material(Si_props),
                        material_c=materials.Material(SiO2_props),
                        material_d=materials.Air,
                        material_e=materials.Material(SiO2_props),
                        material_f=materials.Air,
                        lc_bkg=2, lc2=4000.0, lc3=20.0)

# Expected effective index of fundamental guided mode.
n_eff = wguide.material_b.n-0.1

# Calculate Electromagnetic Modes
sim_EM_pump = wguide.calc_EM_modes(wl_nm, num_modes_EM_pump, n_eff)
sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

# plotting.plt_mode_fields(sim_EM_pump, 
#                          xlim_min=0.35, xlim_max=0.35, ylim_min=0.1, ylim_max=0.55, 
#                          EM_AC='EM', add_name='slab', pdf_png='pdf')

# Print the wavevectors of EM modes.
print('k_z of EM modes \n', np.round(np.real(sim_EM_pump.Eig_values), 4))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.Eig_values[0]*((wl_nm*1e-9)/(2.*np.pi)))

k_AC = np.real(sim_EM_pump.Eig_values[0] - sim_EM_Stokes.Eig_values[0])

shift_Hz = 10e9

# Calculate Acoustic Modes
sim_AC_wguide = wguide.calc_AC_modes(wl_nm, num_modes_AC, k_AC,
    EM_sim=sim_EM_pump, shift_Hz=shift_Hz)

plotting.plt_mode_fields(sim_AC_wguide, EM_AC='AC', add_name='slab', pdf_png='png')
