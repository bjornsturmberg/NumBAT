""" Calculate the backward SBS gain spectra as a function of
    waveguide width, for silicon waveguides surrounded in air.

    Also shows how to use python multiprocessing library.
"""

import time
import datetime
import numpy as np
import sys
from multiprocessing import Pool
sys.path.append("../backend/")
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PolyCollection
from matplotlib.colors import colorConverter

import materials
import objects
import mode_calcs
import integration
import plotting
from fortran import NumBAT


# Select the number of CPUs to use in simulation.
num_cores = 6

# Geometric Parameters - all in nm.
wl_nm = 1550
inc_shape = 'rectangular'

# Optical Parameters
eps = 12.25
num_EM_modes = 20
num_AC_modes = 20
EM_ival1=0
EM_ival2=EM_ival1
AC_ival='All'

# Acoustic Parameters
s = 2330  # kg/m3
c_11 = 165.7e9; c_12 = 63.9e9; c_44 = 79.6e9  # Pa
p_11 = -0.094; p_12 = 0.017; p_44 = -0.051
eta_11 = 5.9e-3 ; eta_12 = 5.16e-3 ; eta_44 = 0.620e-3  # Pa
inc_a_AC_props = [s, c_11, c_12, c_44, p_11, p_12, p_44,
                  eta_11, eta_12, eta_44]


def modes_n_gain(wguide):
    # Calculate Electromagnetic Modes
    sim_EM_wguide = wguide.calc_EM_modes(wl_nm, num_EM_modes)
    # Backward SBS
    k_AC = 2*np.real(sim_EM_wguide.Eig_value[0])
    # Calculate Acoustic Modes
    sim_AC_wguide = wguide.calc_AC_modes(wl_nm, k_AC,
        num_AC_modes, EM_sim=sim_EM_wguide)
    # Calculate interaction integrals and SBS gain
    SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha = integration.gain_and_qs(
        sim_EM_wguide, sim_AC_wguide, k_AC,
        EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival)

    return [sim_EM_wguide, sim_AC_wguide, SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha, k_AC]


nu_widths = 12
waveguide_widths = np.linspace(300,400,nu_widths)
geo_objects_list = []
# Width previous simo's done for, with known meshing params
known_geo = 315.
# Scale meshing to new structures
for width in waveguide_widths:
    msh_ratio = (width/known_geo)
    # Geometric Parameters - all in nm.
    unitcell_x = 2.5*wl_nm*msh_ratio
    unitcell_y = unitcell_x
    inc_a_x = width
    inc_a_y = 0.9*inc_a_x

    lc_bkg = 0.1#/msh_ratio

    # Use all specified parameters to create a waveguide object.
    wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,
                            inc_a_y,inc_shape,
                            bkg_material=materials.Material(1.0 + 0.0j),
                            inc_a_material=materials.Material(np.sqrt(eps)),
                            loss=False, inc_a_AC=inc_a_AC_props,
                            lc_bkg=lc_bkg, lc2=40.0, lc3=10.0)
    geo_objects_list.append(wguide)


# Run widths in parallel across num_cores CPUs using multiprocessing package.
pool = Pool(num_cores)
width_objs = pool.map(modes_n_gain, geo_objects_list)
# np.savez('Simo_results', width_objs=width_objs)
# npzfile = np.load('Simo_results.npz')
# width_objs = npzfile['width_objs'].tolist()


n_effs = []
freqs_gains = []
interp_grid_points = 10000
interp_grid = np.linspace(10, 25, interp_grid_points)
for i_w, width_obj in enumerate(width_objs):
    interp_values = np.zeros(interp_grid_points)
    sim_EM = width_obj[0]
    sim_AC = width_obj[1]
    SBS_gain = width_obj[2]
    alpha = width_obj[5]
    k_AC = width_obj[6]
    # Calculate the EM effective index of the waveguide.
    n_eff_sim = round(np.real(k_AC*((wl_nm*1e-9)/(2.*np.pi))), 4)
    n_effs.append(n_eff_sim)

    # Construct the SBS gain spectrum, built up
    # from Lorentzian peaks of the individual modes.
    tune_steps = 5e4
    tune_range = 10 # GHz
    # Construct an odd range of frequencies that is guaranteed to include 
    # the central resonance frequency
    detuning_range = np.append(np.linspace(-1*tune_range, 0, tune_steps),
                       np.linspace(0, tune_range, tune_steps)[1:])*1e9 # GHz
    phase_v = sim_AC.Eig_value/k_AC
    line_width = phase_v*alpha
    
    plt.figure(figsize=(13,13))
    plt.clf()
    for AC_i in range(len(alpha)):
        gain_list = np.real(SBS_gain[EM_ival1,EM_ival2,AC_i]/alpha[AC_i]
                     *line_width[AC_i]**2/(line_width[AC_i]**2 + detuning_range**2))
        freq_list_GHz = np.real(sim_AC.Eig_value[AC_i] + detuning_range)*1e-9
        plt.plot(freq_list_GHz, gain_list,linewidth=3)
        # set up an interpolation for summing all the gain peaks
        interp_spectrum = np.interp(interp_grid, freq_list_GHz, gain_list)
        interp_values += interp_spectrum
    freqs_gains.append(zip(interp_grid, interp_values))
    plt.plot(interp_grid, interp_values, 'k', linewidth=4)
    plt.xlim(10,25)
    plt.xlabel('Frequency (GHz)')
    plt.ylabel('Gain 1/(Wm)')
    plt.savefig('gain_spectra_%i.pdf' % i_w)
    plt.close()

print 'Widths', waveguide_widths
print 'n_effs', n_effs

# Plot a 'waterfall' plot.
fig = plt.figure()
ax = fig.gca(projection='3d')
poly = PolyCollection(freqs_gains)
poly.set_alpha(0.7)
ax.add_collection3d(poly, zs=waveguide_widths, zdir='y')
ax.set_xlabel('Frequency (GHz)', fontsize=14)
ax.set_xlim3d(10,25)
ax.set_ylabel('Width (nm)', fontsize=14)
ax.set_ylim3d(waveguide_widths[0], waveguide_widths[-1])
ax.set_zlabel('Gain 1/(Wm)', fontsize=14)
ax.set_zlim3d(0,1500)
# We change the fontsize of minor ticks label 
plt.tick_params(axis='both', which='major', labelsize=12, pad=-2)
plt.savefig('gain_spectra_waterfall.pdf')
plt.close()
