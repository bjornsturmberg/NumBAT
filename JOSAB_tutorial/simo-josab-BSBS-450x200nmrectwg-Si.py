"""
Script to evaluate backward Brillouin scattering in a rectangular Si waveguide
"""

# Import the necessary packages
import time
import datetime
import numpy as np
import sys
import math
sys.path.append("../backend/")
import materials
import objects
import mode_calcs
import integration
import plotting
from fortran import NumBAT
from matplotlib import pyplot as plt
# Naming conventions
# AC: acoustic
# EM: electromagnetic
# k_AC: acoustic wavevector

start = time.time()

# Specify Geometric Parameters - all in [nm].
wl_nm = 1550 # Wavelength of EM wave in vacuum.
inc_a_x = 450 # Waveguide widths.
inc_a_y = 200
# Unit cell dimensions must be sufficiently large to ensure fields are zero at outermost boundary.
unitcell_x = 3.01*wl_nm
unitcell_y = unitcell_x #be careful to ensure not whole integer multiples
inc_shape = 'rectangular' # Shape of the waveguide.

# Specify number of electromagnetic modes and acoustic modes involved in the
# calculation for BSBS
num_modes_EM_pump = 20
num_modes_EM_Stokes = num_modes_EM_pump
num_modes_AC = 100
# The mode number for the optical field. Typically 0 for BSBS.
EM_ival_pump = 0
# The EM Stokes mode number for which to calculate interaction with AC modes. Typically 0 for BSBS.
EM_ival_Stokes = EM_ival_pump
# The AC mode(s) for which to calculate interaction with EM modes.
AC_ival = 'All'

# Output files are generated in a folder with the following prefix
prefix_str = 'bsbs-josab-450x200nmSi'

# Use all specified parameters to create a waveguide object
wguide = objects.Struct(unitcell_x,inc_a_x,unitcell_y,inc_a_y,inc_shape,
                        material_bkg=materials.materials_dict["Vacuum"],
                        material_a=materials.materials_dict["Si_2021_Poulton"],
                        lc_bkg=0.05, # mesh coarseness in background, larger lc_bkg = coarser along horizontal outer edge
                        lc_refine_1=20.0, # mesh refinement factor near the interface of waveguide, larger = finer along horizontal interface
                        lc_refine_2=30.0, # mesh refinement factor near the origin/centre of waveguide
                        plt_mesh=False, # creates png file of geometry and mesh in backend/fortran/msh/
                        check_mesh=False) # note requires x-windows configuration to work

# Print information on material data in terminal
print('\nUsing %s material data from' % wguide.material_a.chemical)
print('Author:', wguide.material_a.author)
print('Year:', wguide.material_a.date)
print('Ref:', wguide.material_a.doi)

# Initial guess for the EM effective index of the waveguide
n_eff = wguide.material_a.n-0.1

# Calculate the Electromagnetic modes of the pump field.
sim_EM_pump = wguide.calc_EM_modes(num_modes_EM_pump, wl_nm, n_eff)

# Print the wavevectors of EM modes.
print('\n k_z of EM modes \n', np.round(np.real(sim_EM_pump.Eig_values),4))

# A computation interruption if needed
#sys.exit("We interrupt your regularly scheduled computation to bring you ... for now")

# Calculate the Electromagnetic modes of the Stokes field.
sim_EM_Stokes = mode_calcs.bkwd_Stokes_modes(sim_EM_pump)

# Generating images for the EM modes involved in the calculation
print("Plotting EM fields ")
plotting.plt_mode_fields(sim_EM_pump,
                         ivals=[EM_ival_pump],
                         EM_AC='EM_E', num_ticks=3,xlim_min=0.4, xlim_max=0.4, ylim_min=0.4, ylim_max=0.4,
                         prefix_str=prefix_str, pdf_png='png', ticks=True, quiver_steps=10,
                         comps=['Et','Eabs'], n_points=1000, colorbar=True)

# Calculating the EM effective index of the waveguide.
n_eff_sim = np.real(sim_EM_pump.Eig_values[EM_ival_pump]*((wl_nm*1e-9)/(2.*np.pi)))
print("\n Fundamental optical mode ")
print(" n_eff = ", np.round(n_eff_sim, 4))

# Calculating the acoustic wavevector
k_AC = np.real(sim_EM_pump.Eig_values[EM_ival_pump] - sim_EM_Stokes.Eig_values[EM_ival_Stokes])
print('\n AC wavenumber (1/m) = ', np.round(k_AC, 4))

# Calculating Acoustic modes, using the mesh from the EM calculation.
sim_AC = wguide.calc_AC_modes(num_modes_AC, k_AC, EM_sim=sim_EM_pump)

# Print the frequencies of AC modes.
AC_freqs_Hz = np.round(np.real(sim_AC.Eig_values)*1e-9, 4)
print('\n Freq of AC modes (GHz) \n', AC_freqs_Hz)

# Calculate total SBS gain, photoelastic and moving boundary contributions, as
# well as other important quantities
SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha = integration.gain_and_qs(
    sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC, EM_ival_pump=EM_ival_pump,
    EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)

# Mask negligible gain values to improve clarity of print out.
threshold = -1e-3
masked_PE = np.ma.masked_inside(SBS_gain_PE[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
masked_MB = np.ma.masked_inside(SBS_gain_MB[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)
masked = np.ma.masked_inside(SBS_gain[EM_ival_pump,EM_ival_Stokes,:], 0, threshold)

# Display these in terminal
print("\n Displaying results with negligible components masked out")
print("SBS_gain [1/(Wm)] PE contribution \n", masked_PE)
print("SBS_gain [1/(Wm)] MB contribution \n", masked_MB)
print("SBS_gain [1/(Wm)] total \n", masked)
#determining the location of the maximum gain
maxGainloc=np.argmax(abs(masked.data)) ;

print("Plotting acoustic mode corresponding to maximum")
plotting.plt_mode_fields(sim_AC, EM_AC='AC', prefix_str=prefix_str, ivals=[maxGainloc],
                         num_ticks=3, quiver_steps=40, pdf_png='png',ticks=True, comps=['ut','uabs'], colorbar=True)

# Displaying results for the maximum found in the selection
print("-----------------")
print("Displaying results for maximum gain value found:")
maxGainloc=np.argmax(abs(masked.data)) ;
print("Greatest SBS_gain  [1/(Wm)] total \n", masked.data[maxGainloc])
print("displaying corresponding acoustic mode number (i.e., AC_field_#) for reference \n",maxGainloc )
print("EM Pump Power [Watts] \n", sim_EM_pump.EM_mode_power[EM_ival_pump] )
print("EM Stokes Power [Watts] \n", sim_EM_Stokes.EM_mode_power[EM_ival_Stokes] )
print("EM angular frequency [THz] \n", sim_EM_pump.omega_EM/1e12 )
print("AC Energy Density [J*m^{-1}] \n", sim_AC.AC_mode_energy_elastic[maxGainloc] )
print("AC loss alpha [1/s] \n", alpha[maxGainloc] )
print("AC frequency [GHz] \n", sim_AC.Omega_AC[maxGainloc]/(1e9*2*math.pi) )
print("AC linewidth [MHz] \n", linewidth_Hz[maxGainloc]/1e6)

#since the overlap is not returned directly we'll have to deduce it
absQtot2 = (alpha[maxGainloc]*sim_EM_pump.EM_mode_power[EM_ival_pump]*sim_EM_Stokes.EM_mode_power[EM_ival_Stokes]*sim_AC.AC_mode_energy_elastic[maxGainloc]*masked.data[maxGainloc])/(2*sim_EM_pump.omega_EM*sim_AC.Omega_AC[maxGainloc]);
absQtot = pow(absQtot2,1/2)
print("Total coupling |Qtot| [W*m^{-1}*s] \n", absQtot )

end = time.time()
print("\n Simulation time (sec.)", (end - start))
