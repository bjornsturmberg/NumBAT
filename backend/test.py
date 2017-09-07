
import os
import numpy as np
# from scipy import sqrt
# import subprocess
# from matplotlib.mlab import griddata
# from scipy import interpolate
# import matplotlib
# matplotlib.use('pdf')
# import matplotlib.pyplot as plt
# import matplotlib.gridspec as gridspec
# from mpl_toolkits.axes_grid1 import make_axes_locatable
# from matplotlib import ticker

# from fortran import NumBAT

#### Short utility functions ##################################################
def zeros_int_str(zero_int):
    r""" Construct the SBS gain spectrum, built from Lorentzian peaks of the individual modes.
            

        Args:
            sim_AC : An AC ``Struct`` instance that has had calc_modes calculated

            SBS_gain  (array): Totlat SBS gain of modes.

            SBS_gain_PE  (array): Moving Bountary gain of modes.

            SBS_gain_MB  (array): Photoelastic gain of modes.

            linewidth_Hz  (array): Linewidth of each mode [Hz].

            k_AC  (float): Acoustic wavevector.

            EM_ival_pump  (int or 'All'): Which EM pump mode(s) to consider.

            EM_ival_Stokes  (int or 'All'): Which EM Stokes mode(s) to consider.

            AC_ival  (int or 'All'):  Which AC mode(s) to consider.

            freq_min  (float): Minimum of frequency range.

            freq_max  (float): Maximum of frequency range.

        Keyword Args:
            num_interp_pts  (int): Number of frequency points to interpolate to.

            dB  (bool): Save a set of spectra in dB units.
"""
    string = '%4.0f' % zero_int
    fmt_string = string.replace(' ','0')
    return fmt_string
