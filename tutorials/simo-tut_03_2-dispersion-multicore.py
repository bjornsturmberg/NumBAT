""" Calculate a dispersion diagram of the acoustic modes
    from k_AC ~ 0 (forward SBS) to k_AC = 2*k_EM (backward SBS).
    Use python's (embarrassing parallel) multiprocessing package.
"""

import time
import datetime
import numpy as np
import sys
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from multiprocessing import Pool

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
unitcell_x = 3.0*wl_nm
unitcell_y = unitcell_x
inc_a_x = 800.
inc_a_y = 220.
inc_shape = 'rectangular'
# Choose modes to include.
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 60
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix_str = 'tut_03_2-'

# Note that this mesh is quite fine, may not be required if purely using dispersive sims
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        material_bkg=materials.Vacuum,
                        material_a=materials.Si_2016_Smith,
                        lc_bkg=1, lc2=600.0, lc3=300.0)

# Estimated effective index of fundamental guided mode.
n_eff = wguide.material_a.n-0.1

# Calculate Electromagnetic modes.
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff)
sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

# Will scan from forward to backward SBS so need to know k_AC of backward SBS.
k_AC = np.real(sim_EM_pump.Eig_values[EM_ival_pump] - sim_EM_Stokes.Eig_values[EM_ival_Stokes])

# Rather than calculating with a loop we can use pool to do a multi core sim
def ac_mode_freqs(k_ac):
    print('Commencing mode calculation for k_ac = %f'% k_ac)

    # Calculate the modes, grab the output frequencies only and convert to GHz
    sim_AC = wguide.calc_AC_modes(num_modes_AC, k_ac, EM_sim=sim_EM_pump)
    prop_AC_modes = np.array([np.real(x) for x in sim_AC.Eig_values if abs(np.real(x)) > abs(np.imag(x))])
    mode_freqs = prop_AC_modes*1.e-9
    # Clear memory
    sim_AC = None

    print('Completed mode calculation for width a_x = %f'% k_ac)

    # Return the frequencies and simulated k_ac value in a list
    return mode_freqs


# Now we utilise multi-core calculations to perform parallel simulations and speed up the simulation
test_name = 'dispersion_multicore'
nu_ks = 5  # start with a low number of k_ac values to get an idea
acoustic_ks = np.linspace(5., k_AC*1.1, nu_ks)

num_cores = 5  # should be appropriate for individual machine/vm, and memory!
pool = Pool(num_cores)
pooled_mode_freqs = pool.map(ac_mode_freqs, acoustic_ks)
# Note pool.map() doesn't pass errors back from fortran routines very well.
# It's good practise to run the extrema of your simulation range through map()
# before launcing full multicore simulation.

# We will pack the above values into a single array for plotting purposes, initialise first
freq_arr = np.empty((nu_ks, num_modes_AC))
for i_w, sim_freqs in enumerate(pooled_mode_freqs):
    # Set the value to the values in the frequency array
    freq_arr[i_w] = sim_freqs

# Now that we have packed will save to a numpy file for better plotting and reference
file_name = 'freq_array_200'
np.save(file_name, freq_arr)
np.save(file_name+'_qs', acoustic_ks)  # and the q values

# Also plot a figure for reference
plot_range = num_modes_AC
plt.clf()
plt.figure(figsize=(10,6))
ax = plt.subplot(1,1,1)
for idx in range(plot_range):
    # slicing in the row direction for plotting purposes
    freq_slice = freq_arr[:, idx]
    plt.plot(acoustic_ks/k_AC, freq_slice, 'r')

# Set the limits and plot axis labels
ax.set_ylim(0,35)
ax.set_xlim(0,1.1)
plt.xlabel(r'Axial wavevector (normalised)')
plt.ylabel(r'Frequency (GHz)')
plt.savefig(prefix_str+test_name+'.pdf', bbox_inches='tight')
plt.savefig(prefix_str+test_name+'.png', bbox_inches='tight')
plt.close()

# Output the normalisation k value for reference
print("The 2kp is: %f" % k_AC)

end = time.time()
print("\n Simulation time (sec.)", (end - start))

