""" Replicating the results of
    On-chip inter-modal Brillouin scattering
    Kittlaus et al.
    http://dx.doi.org/10.1038/ncomms15819
"""

import time
import datetime
import numpy as np
import sys
import copy
from matplotlib.ticker import AutoMinorLocator

sys.path.append("../backend/")
import materials
import objects
import mode_calcs
import integration
import plotting
from plotting import FieldDecorator
from fortran import NumBAT

class EMDecorator(FieldDecorator):
  def __init__(self):
    super().__init__()

  def extra_axes_commands(self, ax):

    ax.tick_params(axis='x',color='gray', which='both')
    ax.tick_params(axis='y',color='gray', which='both')
    if self._is_single:
      ax.tick_params(axis='x',length=20)
      ax.tick_params(axis='y',length=20)
      ax.tick_params(axis='x',width=2)
      ax.tick_params(axis='y',width=2)

      ax.tick_params(axis='x',length=10,which='minor')
      ax.tick_params(axis='y',length=10,which='minor')
      ax.tick_params(axis='x',width=2, which='minor')
      ax.tick_params(axis='y',width=2, which='minor')

      ax.set_xticks(np.arange(-1.00,1.01,.1), minor=True)
      ax.set_yticks([-.3,0,.3])
      ax.set_yticks([-.5,-.4,-.2,-.1,.1,.2,.4,.5], minor=True)

    ax.set_aspect('equal')

emdecorate=EMDecorator()

class ACDecorator(plotting.FieldDecorator):
  def __init__(self):
    super().__init__()

    self.set_singleplot_fontsize('ax_label',30) 
    self.set_singleplot_fontsize('subplot_title',30)
    self.set_singleplot_fontsize('cbar_tick',20)
    self.set_singleplot_fontsize('ax_tick',30)

    self.set_singleplot_axes_property('cbar_pad','-30%')  # compensate for call to set_aspect() below
    self.set_multiplot_axes_property('cbar_pad','-30%')   # compensate for call to set_aspect() below
    self.set_singleplot_axes_property('cbar_size','2%')   # compensate for call to set_aspect() below
    self.set_multiplot_axes_property('cbar_size','2%')    # compensate for call to set_aspect() below

  def extra_axes_commands(self, ax):
    ax.set_aspect(3)
    pass

acdecorate=ACDecorator()
# Naming conventions
# AC: acoustic
# EM: electromagnetic
# k_AC: acoustic wavenumber

start = time.time()

# Geometric Parameters - all in nm.
wl_nm = 1550 # Wavelength of EM wave in vacuum.
# Unit cell must be large to ensure fields are zero at boundary.
unitcell_x = 7*wl_nm
unitcell_y = 0.7*unitcell_x
# Waveguide widths.
inc_a_x = 1500
inc_a_y = 80
# Shape of the waveguide.
# Use double coated geometry to control meshing around rib waveguide.
inc_shape = 'rib_double_coated'

slab_a_x = 2850
slab_a_y = 135

# areas included purely
slab_b_y = 500
coat_x = 50 
coat_y = 100
coat2_x = 100
coat2_y = 200

#original old gmsh set that works
lc_bkg = 4  # background
lc2 = 8000  # edge of rib
lc3 = 3000   # edge of slab_a 
lc4 = 50    # edge of coat
lc5 = 20    # edge of slab_b
lc6 = 4     # edge of coat2

##working set
lc_bkg = .5  # background
lc2 = 1000  # edge of rib
lc3 = 400   # edge of slab_a 
lc4 = 10    # edge of coat
lc5 = 5   # edge of slab_b
lc6 = 1   # edge of coat2

#scaled working set: doesn't work
lc_bkg = 1  # background
lc2 = 2000  # edge of rib
lc3 = 600   # edge of slab_a 
lc4 = 10    # edge of coat
lc5 = 4   # edge of slab_b
lc6 = 1  # edge of coat2

##scaled working set: doesn't work
#lc_bkg = .1  # background
#lc2 = 200  # edge of rib
#lc3 = 75   # edge of slab_a 
#lc4 = 50    # edge of coat
#c5  = 50  # edge of slab_b
#lc6 = 50  # edge of coat2





# Number of electromagnetic modes to solve for.
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
# Number of acoustic modes to solve for.
num_modes_AC = 35
# The EM pump mode(s) for which to calculate interaction with AC modes.
# Can specify a mode number (zero has lowest propagation constant) or 'All'.
EM_ival_pump = 0
# The EM Stokes mode(s) for which to calculate interaction with AC modes.
EM_ival_Stokes = 1 # INTERMODE SBS TE0 to TE1
# The AC mode(s) for which to calculate interaction with EM modes.
AC_ival = 'All'

# Si_110 = copy.deepcopy(materials.Si_2015_Van_Laer)
Si_110 = copy.deepcopy(materials.Si_2016_Smith)
Si_110.rotate_axis(np.pi/4,'z-axis', save_rotated_tensors=True)

prefix_str = 'fig16-'

# Use specified parameters to create a waveguide object.
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        slab_a_x=slab_a_x, slab_a_y=slab_a_y, slab_b_y=slab_b_y, 
                        coat_x=coat_x, coat_y=coat_y, coat2_x=coat2_x, coat2_y=coat2_y,
                        material_bkg=materials.Vacuum,
                        material_a=Si_110, #plt_mesh=True,
                        material_b=Si_110, material_c=materials.Vacuum,
                        material_d=materials.Vacuum, material_e=materials.Vacuum,
                        symmetry_flag=False,
                        lc_bkg=lc_bkg, lc2=lc2, lc3=lc3,
                        lc4=lc4, lc5=lc5, lc6=lc6)
# Expected effective index of fundamental guided mode.
n_eff = wguide.material_a.n-0.1

# Calculate Electromagnetic Modes
print("starting EM pump modes")
#sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff=n_eff, debug=True)
#np.savez('wguide_data', sim_EM_pump=sim_EM_pump)
npzfile = np.load('wguide_data.npz', allow_pickle=True)
sim_EM_pump = npzfile['sim_EM_pump'].tolist()

print("starting EM Stokes modes")
sim_EM_Stokes = mode_calcs.fwd_Stokes_modes(sim_EM_pump)
# np.savez('wguide_data2', sim_EM_Stokes=sim_EM_Stokes)
# npzfile = np.load('wguide_data2.npz', allow_pickle=True)
# sim_EM_Stokes = npzfile['sim_EM_Stokes'].tolist()

print("starting EM field plotting ")
plotting.plt_mode_fields(sim_EM_pump, xlim_min=0.4 , xlim_max=0.4 , ivals=[0,1], 
                         ylim_min=0.435 , ylim_max=0.435 , EM_AC='EM_E', num_ticks=3,
                         prefix_str=prefix_str, pdf_png='png', ticks=True,
                          decorator=emdecorate, quiver_steps=20, 
                          comps=('Ex','Ey', 'Ez','Eabs','Et'), n_points=2000, colorbar=True)

plotting.plt_mode_fields(sim_EM_pump, xlim_min=0.4 , xlim_max=0.4 , ivals=[0,1], 
                         ylim_min=0.435 , ylim_max=0.435 , EM_AC='EM_H', num_ticks=3,
                         prefix_str=prefix_str, pdf_png='png', ticks=True,
                          decorator=emdecorate, quiver_steps=20, 
                          comps=('Hx','Hy', 'Hz','Habs','Ht'), n_points=2000, colorbar=True)

# Print the wavevectors of EM modes.
print('k_z of EM modes \n', np.round(np.real(sim_EM_pump.Eig_values), 4))

# Calculate the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.Eig_values*((wl_nm*1e-9)/(2.*np.pi)))
print("n_eff = ", np.round(n_eff_sim, 4))

k_AC = np.real(sim_EM_pump.Eig_values[EM_ival_pump] - sim_EM_Stokes.Eig_values[EM_ival_Stokes])
print('Intermode q_AC (Hz) \n', k_AC)

shift_Hz = 2e9

# Calculate Acoustic Modes
print("starting acoustic modes")
#sim_AC = wguide.calc_AC_modes(num_modes_AC, k_AC, EM_sim=sim_EM_pump, shift_Hz=shift_Hz, debug=True)
#np.savez('wguide_data_AC', sim_AC=sim_AC)
npzfile = np.load('wguide_data_AC.npz', allow_pickle=True)
sim_AC = npzfile['sim_AC'].tolist()

# Print the frequencies of AC modes.
print('Freq of AC modes (GHz) \n', np.round(np.real(sim_AC.Eig_values)*1e-9, 4))

print("plotting acoustic modes")
plotting.plt_mode_fields(sim_AC, EM_AC='AC', prefix_str=prefix_str, ivals=(7, 13, 23,), num_ticks=3, 
     xlim_min=-.05, xlim_max=-0.05, ylim_min=-.1, ylim_max=-0.1, quiver_steps=20,
     pdf_png='png',ticks=True, comps=('ux','ut','uz','uabs'), decorator=acdecorate, colorbar=True)


set_q_factor = 460.

# Calculate interaction integrals and SBS gain for PE and MB effects combined, 
# as well as just for PE, and just for MB.
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC,
    EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival, fixed_Q=set_q_factor)

# Mask negligible gain values to improve clarity of print out.
threshold = 1e-3
masked_PE = np.ma.masked_inside(SBS_gain_PE[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
masked_MB = np.ma.masked_inside(SBS_gain_MB[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
masked = np.ma.masked_inside(SBS_gain[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)

print("\n SBS_gain PE contribution \n", masked_PE)
print("SBS_gain MB contribution \n", masked_MB)
print("SBS_gain total \n", masked)

# Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
freq_min = 0.5  # GHz
freq_max = 9.5  # GHz
plotting.gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
    EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min=freq_min, freq_max=freq_max,
    prefix_str=prefix_str, suffix_str='', pdf_png='png')

end = time.time()
print("\n Simulation time (sec.)", (end - start))
