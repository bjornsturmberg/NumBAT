""" Calculate the convergence as a function of FEM mesh
    for backward SBS gain spectra of a
    silicon waveguide surrounded in air.

    Again utilises python multiprocessing.
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


# Select the number of CPUs to use in simulation.
num_cores = 6

# Geometric Parameters - all in nm.
wl_nm = 1550
unitcell_x = 2.5*wl_nm
unitcell_y = unitcell_x
inc_a_x = 314.7
inc_a_y = 0.9*inc_a_x
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
        num_AC_modes, EM_sim=sim_EM_wguide, shift_Hz=12e9)
    # Calculate interaction integrals and SBS gain
    SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha = integration.gain_and_qs(
        sim_EM_wguide, sim_AC_wguide, k_AC,
        EM_ival1=EM_ival1, EM_ival2=EM_ival2, AC_ival=AC_ival)

    return [sim_EM_wguide, sim_AC_wguide, SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha]

nu_lcs = 15
lc_bkg_list = np.linspace(10,0.4,nu_lcs)
lc_list = np.linspace(1,50,nu_lcs)
geo_objects_list = []
for i_lc, lc_ref in enumerate(lc_list):
    lc_bkg = lc_bkg_list[i_lc]
    wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,
                            inc_a_y,inc_shape,
                            bkg_material=materials.Material(1.0 + 0.0j),
                            inc_a_material=materials.Material(np.sqrt(eps)),
                            loss=False, inc_a_AC=inc_a_AC_props,
                            lc_bkg=lc_bkg, lc2=lc_ref, lc3=lc_ref, force_mesh=True)
    geo_objects_list.append(wguide)


# Do not run in parallel, otherwise there are confusions reading the msh files!
lc_objs = map(modes_n_gain, geo_objects_list)
np.savez('Simo_results', lc_objs=lc_objs)
# npzfile = np.load('Simo_results.npz')
# lc_objs = npzfile['lc_objs'].tolist()


rel_modes = [2,4,8]
rel_mode_freq_EM = np.zeros(nu_lcs,dtype=complex)
rel_mode_freq_AC = np.zeros((nu_lcs,len(rel_modes)),dtype=complex)
rel_mode_gain = np.zeros((nu_lcs,len(rel_modes)),dtype=complex)
# rel_mode_alpha = np.zeros((nu_lcs,len(rel_modes)),dtype=complex)
for i_lc, lc_obj in enumerate(lc_objs):
    rel_mode_freq_EM[i_lc] = lc_obj[0].Eig_value[0]
    for i_m, rel_mode in enumerate(rel_modes):
        rel_mode_freq_AC[i_lc,i_m] = lc_obj[1].Eig_value[rel_mode]
        rel_mode_gain[i_lc,i_m] = lc_obj[2][EM_ival1,EM_ival2,rel_mode]/lc_obj[5][rel_mode]
        # rel_mode_alpha[i_lc,i_m] = lc_obj[5][rel_mode]



xlabel = "Mesh Refinement Factor"
fig = plt.figure()#figsize=(13,13))
plt.clf()
ax1 = fig.add_subplot(1,1,1)
ax2 = ax1.twinx()
ax2.yaxis.tick_left()
ax2.yaxis.set_label_position("left")
EM_plot_Mk = rel_mode_freq_EM*1e-6
error0 = np.abs((np.array(EM_plot_Mk[0:-1])-EM_plot_Mk[-1])/EM_plot_Mk[-1])
ax2.plot(lc_list[0:-1], error0, 'b-v',label=r'Error')
ax1.plot(lc_list, np.real(EM_plot_Mk), 'r-.o',label=r'EM k$_z$')
ax1.yaxis.tick_right()
ax1.spines['right'].set_color('red')
ax1.yaxis.label.set_color('red')
ax1.yaxis.set_label_position("right")
ax1.tick_params(axis='y', colors='red')
# import matplotlib.ticker as mtick
# ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.5f'))
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels, loc='center left')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc='center right')
ax1.set_xlabel(xlabel)
ax1.set_ylabel(r"EM k$_z$ ($\times 10^6$ 1/m)")
ax2.set_ylabel(r"Relative Error EM k$_z$")
ax2.set_yscale('log', nonposx='clip')
plt.savefig('convergence-freq_EM.pdf')#, bbox_inches='tight')
plt.close()


fig = plt.figure()#figsize=(13,13))
plt.clf()
ax1 = fig.add_subplot(1,1,1)
ax2 = ax1.twinx()
ax2.yaxis.tick_left()
ax2.yaxis.set_label_position("left")
for i_m, rel_mode in enumerate(rel_modes):
    rel_mode_freq_AC_plot_GHz = rel_mode_freq_AC[:,i_m]*1e-9
    error0 = np.abs((np.array(rel_mode_freq_AC_plot_GHz[0:-1])-rel_mode_freq_AC_plot_GHz[-1])/rel_mode_freq_AC_plot_GHz[-1])
    ax2.plot(lc_list[0:-1], error0, '-v',label='Error mode #%i'%rel_mode)
    ax1.plot(lc_list, np.real(rel_mode_freq_AC_plot_GHz), '-.o',label=r'AC Freq mode #%i'%rel_mode)
ax1.yaxis.tick_right()
ax1.spines['right'].set_color('red')
ax1.yaxis.label.set_color('red')
ax1.yaxis.set_label_position("right")
ax1.tick_params(axis='y', colors='red')
# import matplotlib.ticker as mtick
# ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.5f'))
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels, loc='center left')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc='center right')
ax1.set_xlabel(xlabel)
ax1.set_ylabel(r"AC Freq (GHz)")
ax2.set_ylabel(r"Relative Error AC Freq")
ax2.set_yscale('log', nonposx='clip')
plt.savefig('convergence-freq_AC.pdf')#, bbox_inches='tight')
plt.close()



fig = plt.figure()#figsize=(13,13))
plt.clf()
ax1 = fig.add_subplot(1,1,1)
ax2 = ax1.twinx()
ax2.yaxis.tick_left()
ax2.yaxis.set_label_position("left")
for i_m, rel_mode in enumerate(rel_modes):
    rel_mode_freq_AC_plot = rel_mode_gain[:,i_m]
    error0 = np.abs((np.array(rel_mode_freq_AC_plot[0:-1])-rel_mode_freq_AC_plot[-1])/rel_mode_freq_AC_plot[-1])
    ax2.plot(lc_list[0:-1], error0, '-v',label=r'Error mode #%i'%rel_mode)
    ax1.plot(lc_list, np.real(rel_mode_freq_AC_plot), '-.o',label=r'Gain mode #%i'%rel_mode)
ax1.yaxis.tick_right()
ax1.spines['right'].set_color('red')
ax1.yaxis.label.set_color('red')
ax1.yaxis.set_label_position("right")
ax1.tick_params(axis='y', colors='red')
# import matplotlib.ticker as mtick
# ax1.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.5f'))
handles, labels = ax2.get_legend_handles_labels()
ax2.legend(handles, labels, loc='center left')
handles, labels = ax1.get_legend_handles_labels()
ax1.legend(handles, labels, loc='center right')
ax1.set_xlabel(xlabel)
ax1.set_ylabel(r"Gain")
ax2.set_ylabel(r"Relative Error Gain")
ax2.set_yscale('log', nonposx='clip')
plt.savefig('convergence-Gain.pdf')#, bbox_inches='tight')
plt.close()
