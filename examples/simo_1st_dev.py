import time
import datetime
import numpy as np
import sys
from multiprocessing import Pool
sys.path.append("../backend/")

import materials
import objects
import mode_calcs
import plotting


unitcell_x = 2.5*1550
inc_a_x = 314.7
unitcell_y = unitcell_x
inc_a_y = 0.9*inc_a_x
inc_shape = 'rectangular'
# slab_a_x=unitcell_x
# slab_a_y=11

wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        bkg_material=materials.Material(1.0 + 0.0j),
                        inc_a_material=materials.Material(np.sqrt(12.25)),
                        lc_bkg=0.09, lc2=6.0, lc3=6.0, check_msh=False)

wl_nm = 1550
num_modes = 30

sim_wguide = wguide.calc_modes(wl_nm, num_modes)

betas = sim_wguide.k_z
# print 'k_z of EM wave \n', betas

plotting.plot_EM_modes(sim_wguide)
