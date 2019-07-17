""" Calculate the backward SBS gain spectra as a function of
    waveguide width, for silicon waveguides surrounded in air.

    Also shows how to use python multiprocessing library.
"""

import time
import datetime
import numpy as np
import sys
from multiprocessing import Pool
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter

sys.path.append("../backend/")
import materials
import objects
import mode_calcs
import integration
import plotting
from fortran import NumBAT


start = time.time()

# Select the number of CPUs to use in simulation.
num_cores = 6

# Geometric Parameters - all in nm.
wl_nm = 1550
inc_shape = 'rectangular'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 20
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix_str = 'tut_04-'

# Width previous simo's done for, with known meshing params
known_geo = 315.

def modes_n_gain(wguide):
    print ('Commencing mode calculation for width a_x = %f' % wguide.inc_a_x)
    # Expected effective index of fundamental guided mode.
    n_eff = (wguide.material_a.n-0.1) * wguide.inc_a_x/known_geo
    # Calculate Electromagnetic modes.
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff)
    sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)
    k_AC = np.real(sim_EM_pump.Eig_values[EM_ival_pump] - sim_EM_Stokes.Eig_values[EM_ival_Stokes])
    # Calculate Acoustic modes.
    sim_AC = wguide.calc_AC_modes(num_modes_AC, k_AC, EM_sim=sim_EM_pump)
    # Calculate interaction integrals and SBS gain.
    SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
        sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC,
        EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)
    ## Clear memory
    #sim_EM_pump = sim_EM_Stokes = sim_AC = None

    print ('Completed mode calculation for width a_x = %f' % wguide.inc_a_x)
    return [sim_EM_pump, sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC]


nu_widths = 6
waveguide_widths = np.linspace(300,350,nu_widths)
geo_objects_list = []
# Scale meshing to new structures.
for width in waveguide_widths:
    msh_ratio = (width/known_geo)
    unitcell_x = 2.5*wl_nm*msh_ratio
    unitcell_y = unitcell_x
    inc_a_x = width
    inc_a_y = 0.9*inc_a_x

    wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,
                            inc_a_y,inc_shape,
                            material_bkg=materials.Vacuum,
                            material_a=materials.Si_2016_Smith,
                            lc_bkg=1, lc2=600.0, lc3=300.0)
    geo_objects_list.append(wguide)


new_calcs=True
if new_calcs:
  # Run widths in parallel across num_cores CPUs using multiprocessing package.
  pool = Pool(num_cores)
  
  # Note pool.map() doesn't pass errors back from fortran routines very well.
  # It's good practise to run the extrema of your simulation range through map()
  # before launcing full multicore simulation.

  width_objs = pool.map(modes_n_gain, geo_objects_list)
  np.savez('Simo_results', width_objs=width_objs)
  
else:
  npzfile = np.load('Simo_results.npz', allow_pickle=True)
  width_objs = npzfile['width_objs'].tolist()
  
n_effs = []
freqs_gains = []
interp_grid_points = 10000
int_min = 10
int_max = 26
interp_grid = np.linspace(int_min, int_max, interp_grid_points)
for i_w, width_obj in enumerate(width_objs):
    interp_values = np.zeros(interp_grid_points)
    sim_EM = width_obj[0]
    sim_AC = width_obj[1]
    SBS_gain = width_obj[2]
    SBS_gain_PE = width_obj[3]
    SBS_gain_MB = width_obj[4]
    linewidth_Hz = width_obj[5]
    k_AC = width_obj[6]
    # Calculate the EM effective index of the waveguide (k_AC = 2*k_EM).
    n_eff_sim = np.round(np.real((k_AC/2.)*((wl_nm*1e-9)/(2.*np.pi))), 4)
    n_effs.append(n_eff_sim)

    print(sim_AC)
    # Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
    freq_min = np.real(sim_AC.Eig_values[0])*1e-9 - 5  # GHz
    freq_max = np.real(sim_AC.Eig_values[-1])*1e-9 + 5  # GHz
    plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
        EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max, 
        prefix_str=prefix_str, suffix_str='_scan%i' % i_w)

    # Repeat calc to collect data for waterfall plot.
    tune_steps = 5e4
    tune_range = 10 # GHz
    detuning_range = np.append(np.linspace(-1*tune_range, 0, tune_steps),
                       np.linspace(0, tune_range, tune_steps)[1:])*1e9 # GHz
    # Linewidth of Lorentzian is half the FWHM style linewidth.
    linewidth = linewidth_Hz/2
    for AC_i in range(len(linewidth_Hz)):
        gain_list = np.real(SBS_gain[EM_ival_Stokes,EM_ival_pump,AC_i]
                     *linewidth[AC_i]**2/(linewidth[AC_i]**2 + detuning_range**2))
        freq_list_GHz = np.real(sim_AC.Eig_values[AC_i] + detuning_range)*1e-9
        interp_spectrum = np.interp(interp_grid, freq_list_GHz, gain_list)
        interp_values += interp_spectrum
    freqs_gains.append(list(zip(interp_grid, abs(interp_values))))

print('Widths', waveguide_widths)
print('n_effs', n_effs)

# Plot a 'waterfall' plot.
fig = plt.figure()
ax = fig.gca(projection='3d')
poly = PolyCollection(freqs_gains)
poly.set_alpha(0.7)
ax.add_collection3d(poly, zs=waveguide_widths, zdir='y')
ax.set_xlabel('Frequency (GHz)', fontsize=14)
ax.set_xlim3d(int_min,int_max)
ax.set_ylabel('Width (nm)', fontsize=14)
ax.set_ylim3d(waveguide_widths[0], waveguide_widths[-1])
ax.set_zlabel('|Gain| (1/Wm)', fontsize=14)
ax.set_zlim3d(0,1500)
# We change the fontsize of minor ticks label 
plt.tick_params(axis='both', which='major', labelsize=12, pad=-2)
plt.savefig(prefix_str+'gain_spectra-waterfall.pdf')
plt.savefig(prefix_str+'gain_spectra-waterfall.png')
plt.close()

end = time.time()
print("\n Simulation time (sec.)", (end - start))

