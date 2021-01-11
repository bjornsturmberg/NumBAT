""" Example showing how the 'onion' geometry template can be used
    to simulate a circular Si waveguide clad in SiO2.
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
from math import *


start = time.time()

# Geometric Parameters - all in nm.
wl_nm = 1550
unitcell_x = 200*wl_nm
unitcell_y = unitcell_x
inc_a_x = 2*4070
inc_b_x = 2*127000
inc_shape = 'onion2'

num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 20
EM_ival_pump = 0
EM_ival_Stokes = 0
AC_ival = 'All'

prefix_str = 'tut_12a-'

# Use of a more refined mesh to produce field plots.
wguide = objects.Struct(unitcell_x,inc_a_x,inc_shape=inc_shape,
                        inc_b_x=inc_b_x,
                        unitcell_y=unitcell_y,
                        material_bkg=materials.materials_dict["Vacuum"],
                        material_a=materials.materials_dict["SiO2GeO2_smf28"],
                        material_b=materials.materials_dict["SiO2_smf28"],
                        lc_bkg=.5, lc_refine_1=10.0, lc_refine_2=200, plt_mesh=True)


# Expected effective index of fundamental guided mode.
n_eff = wguide.material_a.n-0.001

new_calcs=True

micron=1e-6
lam_hi=1.65
lam_lo=1.45
klo=2*pi/(lam_hi)
khi=2*pi/(lam_lo)
ksteps=2
kvec=np.linspace(klo, khi, ksteps)

kzmat=np.zeros([len(kvec), num_modes_EM_pump], dtype=float)
neffmat=np.zeros([len(kvec), num_modes_EM_pump], dtype=float)
for ik, tk in enumerate(kvec):

  t_wl_nm=2*pi/tk*1000
  print(t_wl_nm)

  # Calculate Electromagnetic modes.
  if new_calcs:
    sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, t_wl_nm, n_eff)
    np.savez('wguide_data_{0}'.format(ik), sim_EM_pump=sim_EM_pump)
  
    sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)
    np.savez('wguide_data2_{0}'.format(ik), sim_EM_Stokes=sim_EM_Stokes)
    kzmat[ik]= np.real(sim_EM_pump.Eig_values)
    neffmat[ik]= kzmat[ik]/tk*micron
    print (ik, kzmat[ik], neffmat[ik])
    plotting.plt_mode_fields(sim_EM_pump, xlim_min=0.3, xlim_max=0.3, ylim_min=0.3,
                         ylim_max=0.3, ivals=[0], contours=True, EM_AC='EM_E', 
                         prefix_str=prefix_str+'{0}'.format(ik), ticks=True, quiver_steps=20, comps=['Et'])

    plotting.plt_mode_fields(sim_EM_pump, xlim_min=0.3, xlim_max=0.3, ylim_min=0.3,
                         ylim_max=0.3, ivals=[0], contours=True, EM_AC='EM_H', 
                         prefix_str=prefix_str+'{0}'.format(ik), ticks=True, quiver_steps=20, comps=['Ht'])

  else:
    npzfile = np.load('wguide_data{0}.npz'.format(ik), allow_pickle=True)
    sim_EM_pump = npzfile['sim_EM_pump'].tolist()
    npzfile = np.load('wguide_data2{0}.npz'.format(ik), allow_pickle=True)
    sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()
  
#print (kvec)
#print (kzmat)
#print (neffmat)
out=np.zeros([len(kvec), num_modes_EM_pump+2], dtype=float)
out[:,0]=kvec
out[:,1]=2*pi/kvec
out[:,2:]=neffmat
np.savetxt('neffdata', out)
sys.exit(0)
# Print the wavevectors of EM modes.
print('k_z of EM modes \n', np.round(np.real(sim_EM_pump.Eig_values), 4))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.Eig_values[0]*((wl_nm*1e-9)/(2.*np.pi)))
print("n_eff", np.round(n_eff_sim, 4))

# # Plot the E fields of the EM modes fields - specified with EM_AC='EM_E'.
# # Zoom in on the central region (of big unitcell) with xlim_, ylim_ args.
# # Only plot fields of fundamental (ival = 0) mode.
plotting.plt_mode_fields(sim_EM_pump, xlim_min=0.3, xlim_max=0.3, ylim_min=0.3,
                         ylim_max=0.3, ivals=[0], contours=True, EM_AC='EM_E', 
                         prefix_str=prefix_str, ticks=True, quiver_steps=20, comps=['Et'])

plotting.plt_mode_fields(sim_EM_pump, xlim_min=0.3, xlim_max=0.3, ylim_min=0.3,
                         ylim_max=0.3, ivals=[0], contours=True, EM_AC='EM_H', 
                         prefix_str=prefix_str, ticks=True, quiver_steps=20, comps=['Ht'])

sys.exit(0)
# Acoustic wavevector
k_AC = np.real(sim_EM_pump.Eig_values[EM_ival_pump] - sim_EM_Stokes.Eig_values[EM_ival_Stokes])

# Calculate Acoustic modes.
if new_calcs:
  sim_AC = wguide.calc_AC_modes(num_modes_AC, k_AC, EM_sim=sim_EM_pump)
  np.savez('wguide_data_AC', sim_AC=sim_AC)
else:
  npzfile = np.load('wguide_data_AC.npz', allow_pickle=True)
  sim_AC = npzfile['sim_AC'].tolist()

# Print the frequencies of AC modes.
print('Freq of AC modes (GHz) \n', np.round(np.real(sim_AC.Eig_values)*1e-9, 4))

plotting.plt_mode_fields(sim_AC, EM_AC='AC', pdf_png='png', contours=False, 
                         prefix_str=prefix_str, ticks=True, ivals=[0], quiver_steps=20)

# Calculate the acoustic loss from our fields.
# Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC, EM_ival_pump=EM_ival_pump, 
    EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = np.real(sim_AC.Eig_values[0])*1e-9 - 2  # GHz
freq_max = np.real(sim_AC.Eig_values[-1])*1e-9 + 2  # GHz
plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max, 
    prefix_str=prefix_str)

end = time.time()
print("\n Simulation time (sec.)", (end - start))

