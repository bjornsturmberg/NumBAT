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


unitcell_x = 200
inc_a_x = 70
unitcell_y = 60
inc_a_y = 20
inc_shape = 'rectangular'
slab_a_x=unitcell_x
slab_a_y=11

wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        bkg_material=materials.Material(1.0 + 0.0j),
                        inc_a_material=materials.Material(2.5 + 0.0j))

wl_nm = 900
num_modes = 30

sim_wguide = wguide.calc_modes(wl_nm, num_modes)

betas = sim_wguide.k_z
# print 'k_z of EM wave \n', betas

plotting.plot_EM_modes(sim_wguide)
