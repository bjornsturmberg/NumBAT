""" Calculate the backward SBS gain spectra of a Si
    slot waveguide containing As2S3 on a SiO2 slab.

    This time include a capping layer of SiO2 and 
    investigate the effect of this layer's thickness.
"""

import time
import datetime
import numpy as np
import sys
from multiprocessing import Pool
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
unitcell_x = 4*wl_nm
unitcell_y = 0.3*unitcell_x
inc_shape = 'slot_coated'
inc_a_x = 150
inc_a_y = 190
inc_b_x = 250
# Current mesh template assume inc_b_y = inc_a_y
slab_a_x = 1000
slab_a_y = 100

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 40
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix_str = 'tut_08-'

# Function to return ac freqs for given coating thickness
def ac_mode_freqs(coat_y):
    print('Commencing mode calculation for coat_y = %f'% coat_y)

    wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                            slab_a_x=slab_a_x, slab_a_y=slab_a_y, inc_b_x=inc_b_x,
                            coat_y=coat_y,
                            material_bkg=materials.Vacuum,            # background
                            material_a=materials.As2S3_2017_Morrison, # slot
                            material_b=materials.SiO2_2013_Laude,     # slab
                            material_c=materials.Si_2016_Smith,       # walls of slot
                            material_d=materials.SiO2_2013_Laude,     # coating
                            lc_bkg=1, lc2=400.0, lc3=200.0)

    # Expected effective index of fundamental guided mode.
    n_eff = wguide.material_a.n-0.1

    # Calculate Electromagnetic modes.
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff)
    sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

    k_AC = np.real(sim_EM_pump.Eig_values[EM_ival_pump] - sim_EM_Stokes.Eig_values[EM_ival_Stokes])

    shift_Hz = 4e9

    # Calculate Acoustic modes.
    sim_AC = wguide.calc_AC_modes(num_modes_AC, k_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz)

    # plotting.plt_mode_fields(sim_AC, xlim_min=0.4, xlim_max=0.4, 
    #                           ylim_min=0.7, ylim_max=0.0, EM_AC='AC', 
    #                           prefix_str=prefix_str, suffix_str='_%i' %int(coat_y))

    set_q_factor = 1000.

    SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
        sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC,
        EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_q_factor)

    # Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
    freq_min = 4 # np.real(sim_AC.Eig_values[0])*1e-9 - 2  # GHz
    freq_max = 14 # np.real(sim_AC.Eig_values[-1])*1e-9 + 2  # GHz
    
    plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
        EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max, 
        prefix_str=prefix_str, suffix_str='_%i' %int(coat_y))

    # Convert to GHz
    mode_freqs = sim_AC.Eig_values*1.e-9
    # Clear memory
    wguide = sim_EM_pump = sim_EM_Stokes = sim_AC = None
    SBS_gain = SBS_gain_PE = SBS_gain_MB = linewidth_Hz = Q_factors = alpha = None

    print('Completed mode calculation for coating coat_y = %f'% coat_y)

    # Return the frequencies and simulated k_ac value in a list
    return mode_freqs


nu_coats = 5
coat_min = 5
coat_max = 200
coat_y_list = np.linspace(coat_min,coat_max,nu_coats)

num_cores = 5  # should be appropriate for individual machine/vm, and memory!
pool = Pool(num_cores)
pooled_mode_freqs = pool.map(ac_mode_freqs, coat_y_list)
# Note pool.map() doesn't pass errors back from fortran routines very well.
# It's good practise to run the extrema of your simulation range through map()
# before launcing full multicore simulation.

# We will pack the above values into a single array for plotting purposes, initialise first
freq_arr = np.empty((nu_coats, num_modes_AC))
for i_w, sim_freqs in enumerate(pooled_mode_freqs):
    # Set the value to the values in the frequency array
    freq_arr[i_w] = sim_freqs

# Also plot a figure for reference
plot_range = num_modes_AC
plt.clf()
plt.figure(figsize=(10,6))
ax = plt.subplot(1,1,1)
for idx in range(plot_range):
    # slicing in the row direction for plotting purposes
    freq_slice = freq_arr[:, idx]
    plt.plot(coat_y_list, freq_slice, 'g')

# Set the limits and plot axis labels
ax.set_xlim(coat_min,coat_max)
plt.xlabel(r'Coating Thickness (nm)')
plt.ylabel(r'Frequency (GHz)')
plt.savefig(prefix_str+'freq_changes.pdf', bbox_inches='tight')
plt.savefig(prefix_str+'freq_changes.png', bbox_inches='tight')
plt.close()


end = time.time()
print("\n Simulation time (sec.)", (end - start))

