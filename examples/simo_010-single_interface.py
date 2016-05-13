import time
import datetime
import numpy as np
import sys
from multiprocessing import Pool
sys.path.append("../backend/")

import objects


unitcell_x = 200
inc_a_x = 70
unitcell_y = 60
inc_a_y = 20
inc_shape = 'rectangular'
slab_a_x=unitcell_x
slab_a_y=11
# slab_b_x=None, slab_b_y=None,

objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,slab_a_x,slab_a_x)