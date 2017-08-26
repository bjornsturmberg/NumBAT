""" Replicating the results of
    Interaction between light and highly confined 
    hypersound in a silicon photonic nanowire
    Van Laer et al.
    http://dx.doi.org/10.1038/nphoton.2015.11
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


start = time.time()

# Geometric Parameters - all in nm.
wl_nm = 1550
unitcell_x = 2.5*wl_nm
unitcell_y = unitcell_x
inc_a_x = 450
inc_a_y = 230
inc_shape = 'pedestal'
slab_a_x = 15
slab_a_y = 300
slab_b_x = unitcell_x-400
slab_b_y = 800

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 80
EM_ival_pump = 0
EM_ival_Stokes = EM_ival_pump
AC_ival = 'All'

prefix_str = 'lit_04-'

# Use all specified parameters to create a waveguide object.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        slab_a_x=slab_a_x, slab_a_y=slab_a_y,
                        slab_b_x=slab_b_x, slab_b_y=slab_b_y,
                        material_bkg=materials.Vacuum,            # background
                        material_a=materials.Si_2015_Van_Laer,    # rib
                        material_b=materials.SiO2_2015_Van_Laer,  # slab
                        material_c=materials.SiO2_2015_Van_Laer,  # pillar
                        lc_bkg=2, lc2=4000.0, lc3=1000.0)

# Expected effective index of fundamental guided mode.
n_eff = wguide.material_a.n-0.1

# Calculate Electromagnetic Modes
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff)
sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

# plotting.plt_mode_fields(sim_EM_pump, 
#                          xlim_min=0.35, xlim_max=0.35, ylim_min=0.1, ylim_max=0.55, 
#                          EM_AC='EM_E', pdf_png='pdf', prefix_str=prefix_str, suffix_str='slab')

# Print the wavevectors of EM modes.
print('k_z of EM modes \n', np.round(np.real(sim_EM_pump.Eig_values), 4))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.Eig_values[0]*((wl_nm*1e-9)/(2.*np.pi)))

k_AC = 5

shift_Hz = 8e9

# Calculate Acoustic Modes
sim_AC = wguide.calc_AC_modes(num_modes_AC, k_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)

plotting.plt_mode_fields(sim_AC, EM_AC='AC', prefix_str=prefix_str, suffix_str='slab', pdf_png='png')

end = time.time()
print("\n Simulation time (sec.)", (end - start))

