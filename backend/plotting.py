# plotting.py is a subroutine of NumBAT that contains numerous plotting
# routines.

# Copyright (C) 2017  Bjorn Sturmberg.

# NumBAT is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program. If not, see <http://www.gnu.org/licenses/>.


import os
import sys
import numpy as np
from scipy import sqrt
import subprocess
from matplotlib.mlab import griddata
from scipy import interpolate
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib import ticker

from fortran import NumBAT

try: 
    plt.style.use('NumBATstyle')
except (ValueError, IOError, AttributeError): print("Preferred matplotlib style file not found.")
colors = [color['color'] for color in list(plt.rcParams['axes.prop_cycle'])]


# font = {'family' : 'normal',
#         'weight' : 'bold',
#         'size'   : 18}
# matplotlib.rc('font', **font)
linesstrength = 2.5

class FieldDecorator(object):
  def __init__(self):
    base_font=24
    self._multiplot_fontsizes= {'title':base_font-2, 'subplot_title':base_font-5, 'ax_label':base_font-10, 
        'ax_tick':base_font-10, 'cbar_tick':base_font-10}
    self._multiplot_axesprops={'linewidth':'.75', 'edgecolor':'gray', 'title_pad':5 , 'cbar_size':'5%', 'cbar_pad': '2%'}

    self._singleplot_fontsizes= {'title':base_font-2, 'ax_label':50, 'subplot_title':50,'cbar_tick':30, 'ax_tick':40 }
    self._singleplot_axesprops={'linewidth':'.75', 'edgecolor':'gray', 'title_pad':20, 'cbar_size':'5%', 'cbar_pad': '2%'}

    self._is_single=True

  def _fontsizes(self):
    if self._is_single: return self._singleplot_fontsizes
    else: return self._multiplot_fontsizes

  def _get_axes_prop(self):
    if self._is_single: return self._singleplot_axesprops
    else: return self._multiplot_axesprops

  def _set_for_single(self): 
    self._is_single=True
  def _set_for_multi(self):  
    self._is_single=False

  def get_font_size(self, lab):
    fs=10
    try:
      fs=self._fontsizes()[lab]
    except:
      print('Warning: unknown fontsize label "{0}" in FieldDecorator::get_font_size()'.format(lab))
    return fs

  def get_axes_property(self, lab):
    prop=''
    try:
      prop=self._get_axes_prop()[lab]
    except:
      print('Warning: unknown axes property label "{0}" in FieldDecorator::get_axes_property()'.format(lab))
      print (self._get_axes_prop())
      sys.exit(1)
    return prop

  def is_single_plot(self): return self._is_single

  def set_singleplot_axes_property(self, lab, sz): 
    self._singleplot_axesprops[lab]=sz
  def set_multiplot_axes_property(self, lab, sz): 
    self._multiplot_axesprops[lab]=sz

  def set_singleplot_fontsize(self, lab, sz): 
    self._singleplot_fontsizes[lab]=sz
  def set_multiplot_fontsize(self, lab, sz): 
    self._multiplot_fontsizes[lab]=sz

  def extra_axes_commands(self, ax):
    pass



#### Natural constants ########################################################
ASTM15_tot_I   = 900.084            # Integral ASTM 1.5 solar irradiance W/m**2
Plancks_h      = 6.62607015e-34     # Planck's constant in Js (exact)
speed_c        = 299792458          # Speed of light in vacuum in m/s (exact)
charge_e       = 1.602176634e-19    # Charge of an electron in C (exact)
###############################################################################


#### Short utility functions ##################################################
def zeros_int_str(zero_int):
    """ Convert integer into string with '0' in place of ' '. """
    # if zero_int == 0:
    #     fmt_string = '0000'
    # else:
    #     string = '%4.0f' % zero_int
    #     fmt_string = string.replace(' ','0')
    string = '%4.0f' % zero_int
    fmt_string = string.replace(' ','0')
    return fmt_string

###############################################################################


def gain_spectra(sim_AC, SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, k_AC,
                EM_ival_pump, EM_ival_Stokes, AC_ival, freq_min, freq_max, num_interp_pts=3000,
                save_fig=True, dB=False, dB_peak_amp=10, mode_comps=False, semilogy=False,
                pdf_png='png', save_txt=False, prefix_str='', suffix_str=''):
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

            dB_peak_amp  (float): Set the peak amplitude of highest gain mode in dB.

            mode_comps  (bool): Plot decomposition of spectra into individual modes.

            semilogy  (bool): PLot y-axis on log scale.

            save_fig  (bool): Save figure at all.

            pdf_png  (str): Save figures as 'png' or 'pdf'.

            save_txt  (bool): Save spectra data to txt file.

            prefix_str  (str): String to be appended to start of file name.

            suffix_str  (str): String to be appended to end of file name.
    """
    tune_steps = 5e4
    tune_range = 10 # GHz
    # Construct an odd range of freqs guaranteed to include central resonance frequency.
    detuning_range = np.append(np.linspace(-1*tune_range, 0, tune_steps),
                       np.linspace(0, tune_range, tune_steps)[1:])*1e9 # GHz
    interp_grid = np.linspace(freq_min, freq_max, num_interp_pts)

    # Linewidth of Lorentzian is half the FWHM style linewidth.
    linewidth = linewidth_Hz/2

    # Plot decomposition of spectra into individual modes.
    interp_values = np.zeros(num_interp_pts)
    if mode_comps and save_fig:
        plt.figure()
        plt.clf()
    if AC_ival == 'All':
        for AC_i in range(len(linewidth)):
            gain_list = np.real(SBS_gain[EM_ival_pump,EM_ival_Stokes,AC_i]
                         *linewidth[AC_i]**2/(linewidth[AC_i]**2 + detuning_range**2))
            freq_list_GHz = np.real(sim_AC.Eig_values[AC_i] + detuning_range)*1e-9
            if mode_comps:
                plt.plot(freq_list_GHz, np.abs(gain_list), linewidth=2)
                if save_txt:
                    save_array = (freq_list_GHz, gain_list)
                    np.savetxt('gain_spectra-mode_comps%(add)s-%(mode)i.csv' 
                                % {'add' : suffix_str, 'mode' : AC_i}, 
                                save_array, delimiter=',')
            # set up an interpolation for summing all the gain peaks
            interp_spectrum = np.interp(interp_grid, freq_list_GHz, gain_list)
            interp_values += interp_spectrum
    else: raise NotImplementedError("Spectrum plotting for limited AC modes not implemented.")
    return_interp_values = interp_values
    if mode_comps:
        plt.plot(interp_grid, np.abs(interp_values), 'b', linewidth=3, label="Total")
        plt.legend(loc=0)
        if freq_min and freq_max:
            plt.xlim(freq_min,freq_max)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('|Gain| (1/Wm)')
        if save_txt:
            save_array = (interp_grid, interp_values)
            np.savetxt('gain_spectra-mode_comps%(add)s-Total.csv' 
                        % {'add' : suffix_str}, save_array, delimiter=',')

    if mode_comps and save_fig:
        if pdf_png=='png':
            plt.savefig('%(pre)sgain_spectra-mode_comps%(add)s.png' % {'pre' : prefix_str, 'add' : suffix_str})
        elif pdf_png=='pdf':
            plt.savefig('%(pre)sgain_spectra-mode_comps%(add)s.pdf' % {'pre' : prefix_str, 'add' : suffix_str})
        plt.close()


    interp_values = np.zeros(num_interp_pts)
    interp_values_PE = np.zeros(num_interp_pts)
    interp_values_MB = np.zeros(num_interp_pts)
    plt.figure()
    plt.clf()
    if AC_ival == 'All':
        for AC_i in range(len(linewidth)):
            gain_list = np.real(SBS_gain[EM_ival_pump,EM_ival_Stokes,AC_i]
                         *linewidth[AC_i]**2/(linewidth[AC_i]**2 + detuning_range**2))
            freq_list_GHz = np.real(sim_AC.Eig_values[AC_i] + detuning_range)*1e-9
            interp_spectrum = np.interp(interp_grid, freq_list_GHz, gain_list)
            interp_values += interp_spectrum

            gain_list_PE = np.real(SBS_gain_PE[EM_ival_pump,EM_ival_Stokes,AC_i]
                         *linewidth[AC_i]**2/(linewidth[AC_i]**2 + detuning_range**2))
            interp_spectrum_PE = np.interp(interp_grid, freq_list_GHz, gain_list_PE)
            interp_values_PE += interp_spectrum_PE

            gain_list_MB = np.real(SBS_gain_MB[EM_ival_pump,EM_ival_Stokes,AC_i]
                         *linewidth[AC_i]**2/(linewidth[AC_i]**2 + detuning_range**2))
            interp_spectrum_MB = np.interp(interp_grid, freq_list_GHz, gain_list_MB)
            interp_values_MB += interp_spectrum_MB
    else: raise NotImplementedError("Spectrum plotting for limited AC modes not implemented.")
    if save_fig:
        plt.plot(interp_grid, np.abs(interp_values_PE), 'r', linewidth=3, label="PE")
        plt.plot(interp_grid, np.abs(interp_values_MB), 'g', linewidth=3, label="MB")
        plt.plot(interp_grid, np.abs(interp_values), 'b', linewidth=3, label="Total")
        plt.legend(loc=0)
        if freq_min and freq_max:
            plt.xlim(freq_min,freq_max)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('|Gain| (1/Wm)')

        if pdf_png=='png':
            plt.savefig('%(pre)sgain_spectra-MB_PE_comps%(add)s.png' % {'pre' : prefix_str, 'add' : suffix_str})
        elif pdf_png=='pdf':
            plt.savefig('%(pre)sgain_spectra-MB_PE_comps%(add)s.pdf' % {'pre' : prefix_str, 'add' : suffix_str})
    plt.close()

    if save_txt:
        save_array = (interp_grid, interp_values)
        np.savetxt('%(pre)sgain_spectra-MB_PE_comps%(add)s-Total.csv' 
                    % {'pre' : prefix_str, 'add' : suffix_str}, 
                    save_array, delimiter=',')
        save_array = (interp_grid, interp_values_PE)
        np.savetxt('%(pre)sgain_spectra-MB_PE_comps%(add)s-PE.csv' 
                    % {'pre' : prefix_str, 'add' : suffix_str}, 
                    save_array, delimiter=',')
        save_array = (interp_grid, interp_values_MB)
        np.savetxt('%(pre)sgain_spectra-MB_PE_comps%(add)s-MB.csv' 
                    % {'pre' : prefix_str, 'add' : suffix_str}, 
                    save_array, delimiter=',')

    if dB:
        plt.figure()
        plt.clf()

        max_G = np.max(interp_values)
        dB_const = dB_peak_amp/(4.34*max_G)
        plt.plot(interp_grid, np.abs(10*np.log10(np.exp(abs(interp_values)*dB_const))), 'b', linewidth=3, label="Total")
        plt.legend(loc=0)
        if freq_min and freq_max:
            plt.xlim(freq_min,freq_max)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('Amplification (dB)')

        if pdf_png=='png':
            plt.savefig('%(pre)sgain_spectra-dB%(add)s.png' % {'pre' : prefix_str, 'add' : suffix_str})
        elif pdf_png=='pdf':
            plt.savefig('%(pre)sgain_spectra-dB%(add)s.pdf' % {'pre' : prefix_str, 'add' : suffix_str})
        plt.close()

        if save_txt:
            save_array = (interp_grid, 10*np.log10(np.exp(abs(interp_values)*dB_const)))
            np.savetxt('%(pre)sgain_spectra-dB%(add)s.csv' 
                        % {'pre' : prefix_str, 'add' : suffix_str}, 
                        save_array, delimiter=',')

    if semilogy and save_fig:
        plt.figure()
        plt.clf()
        plt.semilogy(interp_grid, abs(interp_values_PE), 'r', linewidth=3, label="PE")
        plt.semilogy(interp_grid, abs(interp_values_MB), 'g', linewidth=3, label="MB")
        plt.semilogy(interp_grid, abs(interp_values), 'b', linewidth=2, label="Total")
        plt.legend(loc=0)
        if freq_min and freq_max:
            plt.xlim(freq_min,freq_max)
        plt.xlabel('Frequency (GHz)')
        plt.ylabel('|Gain| (1/Wm)')

        if pdf_png=='png':
            plt.savefig('%(pre)sgain_spectra-MB_PE_comps-logy%(add)s.png' % {'pre' : prefix_str, 'add' : suffix_str})
        elif pdf_png=='pdf':
            plt.savefig('%(pre)sgain_spectra-MB_PE_comps-logy%(add)s.pdf' % {'pre' : prefix_str, 'add' : suffix_str})
        plt.close()

    return return_interp_values

def plot_component_axes(ax, v_x, v_y, v_XX, v_YY, plot, v_label, plps):
  decorator = plps['decorator']

  plot_threshold = 1e-4 # set negligible components to explicitly zero

  if plps['ticks']: #set tick range
    xm=v_x[-1]+v_x[0]
    ym=v_y[-1]+v_y[0]
    extents=np.array([v_x[0]-xm/2, v_x[-1]-xm/2, v_y[0]-ym/2, v_y[-1]-ym/2])*1e6  # Convert to length units in microns

  else:
    extents=None

  if np.max(np.abs(plot[~np.isnan(plot)])) < plot_threshold: # if the data is all noise, just plot zeros
      # im = plt.imshow(plot.T,cmap='viridis');
      im = plt.imshow(np.zeros(np.shape(plot.T)), origin='lower', extent=extents);
  else:
      interp=None
      #interp='bilinear';
      im = plt.imshow(plot.T, origin='lower', extent=extents, interpolation=interp)

  # limits
  axes = plt.gca()
  xmin, xmax = axes.get_xlim()
  ymin, ymax = axes.get_ylim()
  width_x = xmax-xmin
  width_y = ymax-ymin
  
  if plps['xlim_min'] != None:
      ax.set_xlim(xmin+plps['xlim_min']*width_x,xmax-plps['xlim_max']*width_x)
  if plps['ylim_min'] != None:
      ax.set_ylim(ymin+plps['ylim_min']*width_y,ymax-plps['ylim_max']*width_y)

  if plps['ticks']:
    plt.xticks()
    plt.yticks()
    ax.tick_params(labelsize=decorator.get_font_size('ax_tick'))
    ax.xaxis.set_tick_params(width=1.0)
    ax.yaxis.set_tick_params(width=1.0)
    ax.set_xlabel('$x$ [$\mathrm{\mu}$m]', size=decorator.get_font_size('ax_label'))
    ax.set_ylabel('$y$ [$\mathrm{\mu}$m]', size=decorator.get_font_size('ax_label'))
    ax.grid(False)
  else:
    plt.xticks([])
    plt.yticks([])

  plt.rc('axes',titlepad=decorator.get_axes_property('title_pad'))
  plt.title(v_label,fontsize=decorator.get_font_size('subplot_title'))


  # colorbar
  docbar=plps['colorbar']
  if docbar:
    divider = make_axes_locatable(ax)
    #cax = divider.append_axes("right", size="2%", pad=-6)
    cax = divider.append_axes("right", size=decorator.get_axes_property('cbar_size'), pad=decorator.get_axes_property('cbar_pad'))
    cbar = plt.colorbar(im, cax=cax)
    if plps['num_ticks']:
        cbarticks = np.linspace(np.min(plot), np.max(plot), num=plps['num_ticks'])                
    elif plps['ylim_min'] != None and plps['ylim_min'] >0 : # this is very strange logic. what is going on here? probably meant to be saying something different
        if plps['xlim_min']/plps['ylim_min'] > 3:
            cbarticks = np.linspace(np.min(plot), np.max(plot), num=3)
        if plps['xlim_min']/plps['ylim_min'] > 1.5:
            cbarticks = np.linspace(np.min(plot), np.max(plot), num=5)
        else:
            cbarticks = np.linspace(np.min(plot), np.max(plot), num=7)
    else:
        cbarticks = np.linspace(np.min(plot), np.max(plot), num=7)
    cbar.set_ticks(cbarticks)
    cbarlabels = ['%.2f' %t for t in cbarticks]
    cbar.set_ticklabels(cbarlabels)
    cbar.ax.tick_params(labelsize=decorator.get_font_size('cbar_tick'))

  if plps['contours']:
    if plps['contour_lst']:
        if docbar: cbarticks = plps['contour_lst']
    if np.max(np.abs(plot[~np.isnan(plot)])) > plot_threshold:
        CS2 = ax.contour(v_XX, v_YY, plot.T, levels=cbarticks, colors=colors[::-1], linewidths=(1.5,))
        if docbar: cbar.add_lines(CS2)

  if decorator!=None: decorator.extra_axes_commands(ax)

def plot_supertitle(plps, sim_wguide, ival):
  if plps['EM_AC']=='EM_E' or plps['EM_AC']=='EM_H':
     kz=sim_wguide.Eig_values[ival]
     n_eff = kz * sim_wguide.wl_m / (2*np.pi)
     if np.imag(kz) < 0:
         k_str = r'$k_z = %(re_k)f %(im_k)f i$'% {'re_k' : np.real(kz), 'im_k' : np.imag(kz)}
         n_str = r'$n_{eff} = %(re_k)f %(im_k)f i$'%  {'re_k' : np.real(n_eff), 'im_k' : np.imag(n_eff)}
     else:
         k_str = r'$k_z = %(re_k)f + %(im_k)f i$'% {'re_k' : np.real(kz), 'im_k' : np.imag(kz)}
         n_str = r'$n_{eff} = %(re_k)f + %(im_k)f i$'%  {'re_k' : np.real(n_eff), 'im_k' : np.imag(n_eff)}
  else:
     n_str = ''
     kz=sim_wguide.Eig_values[ival] * 1e-9
     if np.imag(sim_wguide.Eig_values[ival]) < 0:
         k_str = r'$\Omega/2\pi = %(re_k)f %(im_k)f i$ GHz'% {'re_k' : np.real(kz), 'im_k' : np.imag(kz)} 
     else:
         k_str = r'$\Omega/2\pi = %(re_k)f + %(im_k)f i$ GHz'% {'re_k' : np.real(kz), 'im_k' : np.imag(kz)} 
  fulltitle='Mode #' + str(ival) + '   ' + k_str + '   ' + n_str+"\n"
  return fulltitle

def plot_filename(plps, ival, label=None):
  filestart='%(pre)sfields/%(s)s_field_%(i)i%(add)s' % {'pre' : plps['prefix_str'], 
          's' : plps['EM_AC'], 'i' : ival, 'add' : plps['suffix_str']}
  if label!=None: filestart+='_'+label
 
  if plps['pdf_png']=='png': fig_filename=filestart+'.png'
  elif plps['pdf_png']=='pdf': fig_filename=filestart+'.pdf'
  else: raise ValueError("pdf_png must be either 'png' or 'pdf'.")

  return fig_filename

def plot_component_quiver(ax, v_x_q, v_y_q, vq_plots, plps):
  decorator = plps['decorator']

  m_ReEx_q = vq_plots[0]
  m_ReEy_q = vq_plots[1]
  m_ImEx_q = vq_plots[2]
  m_ImEy_q = vq_plots[3]

  # convert to microns
  v_x_q_um=v_x_q*1e6
  v_y_q_um=v_y_q*1e6

  # centre at zero
  xm=v_x_q_um[-1,0]+v_x_q_um[0,0]
  ym=v_y_q_um[0,-1]+v_y_q_um[0,0]
  v_x_q_um-=xm/2
  v_y_q_um-=ym/2

  #TODO: for some reason, the elastic fields make much nicer quiver plots than the E/H fields. why?
  m_x_q = m_ReEx_q+m_ImEx_q
  m_y_q = m_ReEy_q+m_ImEy_q
  plt.quiver(v_x_q_um, v_y_q_um, m_x_q, m_y_q, 
      np.sqrt(np.real((m_ReEx_q+1j*m_ImEx_q)*(m_ReEx_q-1j*m_ImEx_q)
      +(m_ReEy_q+1j*m_ImEy_q)*(m_ReEy_q-1j*m_ImEy_q))),  #colour the arrows based on this array
      linewidths=(0.2,), edgecolors=('k'), pivot='mid', headlength=5) # length of the arrows

  ax.set_aspect('equal')
  axes = plt.gca()
  ax.set_xlim(v_x_q_um[0,0],v_x_q_um[-1,0]) # this step is needed because quiver doesn't seem
 
  ax.set_ylim(v_y_q_um[0,0],v_y_q_um[0,-1]) # to use its input x and y vectors to set range limits
                                            # clean this up so as to avoid seemingly circular calls following
  xmin, xmax = axes.get_xlim()
  ymin, ymax = axes.get_ylim()
  width_x = xmax-xmin
  width_y = ymax-ymin

  if plps['ticks']:
    plt.xticks()
    plt.yticks()
    ax.tick_params(labelsize=decorator.get_font_size('ax_tick'))
    ax.xaxis.set_tick_params(width=1.0)
    ax.yaxis.set_tick_params(width=1.0)
    ax.set_xlabel('$x$ [$\mathrm{\mu}$m]', size=decorator.get_font_size('ax_label'))
    ax.set_ylabel('$y$ [$\mathrm{\mu}$m]', size=decorator.get_font_size('ax_label'))
    ax.grid(False)
  else:
    plt.xticks([])
    plt.yticks([])

  if plps['xlim_min'] != None:
      ax.set_xlim(xmin+plps['xlim_min']*width_x,xmax-plps['xlim_max']*width_x)
  if plps['ylim_min'] != None:
      ax.set_ylim(ymin+plps['ylim_min']*width_y,ymax-plps['ylim_max']*width_y)

  plt.rc('axes',titlepad=decorator.get_axes_property('title_pad'))

  if plps['EM_AC']=='EM_E':   plt.title('$(E_x, E_y)$',fontsize=decorator.get_font_size('subplot_title'))
  elif plps['EM_AC']=='EM_H': plt.title('$(H_x, H_y)$',fontsize=decorator.get_font_size('subplot_title'))
  elif plps['EM_AC']=='AC':   plt.title('$(u_x, u_y)$',fontsize=decorator.get_font_size('subplot_title'))

  decorator = plps['decorator']
  if decorator!=None: decorator.extra_axes_commands(ax)

def plot_all_components(v_x, v_y, v_x_q, v_y_q, v_XX, v_YY, v_plots, vq_plots, v_labels, plps, sim_wguide, ival):
  decorator = plps['decorator']
  plt.clf()
  fig = plt.figure(figsize=(15,15))
  plt.rc('axes', linewidth=decorator.get_axes_property('linewidth'))
  plt.rc('axes', edgecolor=decorator.get_axes_property('edgecolor'))

        # field plots
  for i_p,plot in enumerate(v_plots):
     ax = plt.subplot(3,3,i_p+1)
     plot_component_axes(ax, v_x, v_y, v_XX, v_YY, plot, v_labels[i_p], plps)  # the scalar plots

  ax = plt.subplot(3,3,i_p+2)
  plot_component_quiver(ax, v_x_q, v_y_q, vq_plots, plps)  # the transverse vector plot

   
  fulltitle=plot_supertitle(plps, sim_wguide, ival)
  plt.suptitle(fulltitle, fontsize=decorator.get_font_size('title'))

  # TODO: add get_axes_property() options to play with spacing of these
  # plt.tight_layout(pad=2.5, w_pad=0.5, h_pad=1.0)
  fig.set_tight_layout(True)

  fig.subplots_adjust(hspace=.05, wspace=0.05)

  figfile=plot_filename(plps, ival)
  save_figure(plt, figfile)

  plt.close()


  if plps['EM_AC']=='AC' and plps['stress_fields']:
     ### Interpolate onto rectangular Cartesian grid
     xy = list(zip(v_x6p, v_y6p))
     grid_x, grid_y = np.mgrid[x_min:x_max:n_pts_x*1j, y_min:y_max:n_pts_y*1j]
     m_ReEx = interpolate.griddata(xy, v_Ex6p.real, (grid_x, grid_y), method='linear')
     m_ReEy = interpolate.griddata(xy, v_Ey6p.real, (grid_x, grid_y), method='linear')
     m_ReEz = interpolate.griddata(xy, v_Ez6p.real, (grid_x, grid_y), method='linear')
     m_ImEx = interpolate.griddata(xy, v_Ex6p.imag, (grid_x, grid_y), method='linear')
     m_ImEy = interpolate.griddata(xy, v_Ey6p.imag, (grid_x, grid_y), method='linear')
     m_ImEz = interpolate.griddata(xy, v_Ez6p.imag, (grid_x, grid_y), method='linear')
     m_AbsE = interpolate.griddata(xy, v_E6p.real, (grid_x, grid_y), method='linear')
     dx = grid_x[-1,0] - grid_x[-2,0]
     dy = grid_y[0,-1] - grid_y[0,-2]
     m_Ex = m_ReEx + 1j*m_ImEx
     m_Ey = m_ReEy + 1j*m_ImEy
     m_Ez = m_ReEz + 1j*m_ImEz
     m_Ex = m_Ex.reshape(n_pts_x,n_pts_y)
     m_Ey = m_Ey.reshape(n_pts_x,n_pts_y)
     m_Ez = m_Ez.reshape(n_pts_x,n_pts_y)
     m_AbsE = m_AbsE.reshape(n_pts_x,n_pts_y)

     m_ReEx = np.real(m_Ex)
     m_ReEy = np.real(m_Ey)
     m_ReEz = np.real(m_Ez)
     m_ImEx = np.imag(m_Ex)
     m_ImEy = np.imag(m_Ey)
     m_ImEz = np.imag(m_Ez)

     del_x_Ex = np.gradient(m_Ex, dx, axis=0)
     del_y_Ex = np.gradient(m_Ex, dy, axis=1)
     del_x_Ey = np.gradient(m_Ey, dx, axis=0)
     del_y_Ey = np.gradient(m_Ey, dy, axis=1)
     del_x_Ez = np.gradient(m_Ez, dx, axis=0)
     del_y_Ez = np.gradient(m_Ez, dy, axis=1)
     del_z_Ex = 1j*sim_wguide.k_AC*m_Ex
     del_z_Ey = 1j*sim_wguide.k_AC*m_Ey
     del_z_Ez = 1j*sim_wguide.k_AC*m_Ez

     # Flip y order as imshow has origin at top left
     del_mat = np.array([del_x_Ex[:,::-1].real, del_x_Ey[:,::-1].real, del_x_Ez[:,::-1].real, del_x_Ex[:,::-1].imag, del_x_Ey[:,::-1].imag, del_x_Ez[:,::-1].imag, del_y_Ex[:,::-1].real, del_y_Ey[:,::-1].real, del_y_Ez[:,::-1].real, del_y_Ex[:,::-1].imag, del_y_Ey[:,::-1].imag, del_y_Ez[:,::-1].imag, del_z_Ex[:,::-1].real, del_z_Ey[:,::-1].real, del_z_Ez[:,::-1].real, del_z_Ex[:,::-1].imag, del_z_Ey[:,::-1].imag, del_z_Ez[:,::-1].imag])
     v_labels = ["Re($S_{xx}$)","Re($S_{xy}$)","Re($S_{xz}$)","Im($S_{xx}$)","Im($S_{xy}$)","Im($S_{xz}$)","Re($S_{yx}$)","Re($S_{yy}$)","Re($S_{yz}$)","Im($S_{yx}$)","Im($S_{yy}$)","Im($S_{yz}$)","Re($S_{zx}$)","Re($S_{zy}$)","Re($S_{zz}$)","Im($S_{zx}$)","Im($S_{zy}$)","Im($S_{zz}$)"]

     # stress field plots
     plt.clf()
     fig = plt.figure(figsize=(15,30))
     for i_p,plot in enumerate(del_mat):
         ax = plt.subplot(6,3,i_p+1)
         im = plt.imshow(plot.T);
         # no ticks
         plt.xticks([])
         plt.yticks([])
         # limits
         if xlim_min != None:
             ax.set_xlim(xlim_min*n_points,(1-xlim_max)*n_points)
         if ylim_min != None:
             ax.set_ylim((1-ylim_min)*n_points,ylim_max*n_points)
         # titles
         plt.title(v_labels[i_p],fontsize=decorator.get_font_size('subplot_title'))
         # colorbar
         divider = make_axes_locatable(ax)
         cax = divider.append_axes("right", size="5%", pad=0.1)
         cbar = plt.colorbar(im, cax=cax, format='%.2e')
         if num_ticks:
             cbarticks = np.linspace(np.min(plot), np.max(plot), num=num_ticks)                
         elif ylim_min != None:
             if xlim_min/ylim_min > 3:
                 cbarlabels = np.linspace(np.min(plot), np.max(plot), num=3)
             if xlim_min/ylim_min > 1.5:
                 cbarlabels = np.linspace(np.min(plot), np.max(plot), num=5)
             else:
                 cbarlabels = np.linspace(np.min(plot), np.max(plot), num=7)
         else:
             cbarlabels = np.linspace(np.min(plot), np.max(plot), num=7)
         cbar.set_ticks(cbarlabels)
         cbarlabels = ['%.2f' %t for t in cbarlabels]
         cbar.set_ticklabels(cbarlabels)
         if contours:
             if contour_lst:
                 cbarticks = contour_lst
             if np.max(np.abs(plot[~np.isnan(plot)])) > plot_threshold:
                 CS2 = ax.contour(v_XX, v_YY, plot.T, levels=cbarticks, colors=colors[::-1], linewidths=(1.5,))
             cbar.add_lines(CS2)
         cbar.ax.tick_params(labelsize=decorator.get_font_size('cbar_tick'))
     fig.set_tight_layout(True)
     n_str = ''
     if np.imag(sim_wguide.Eig_values[ival]) < 0:
         k_str = r'$\Omega/2\pi = %(re_k)f %(im_k)f i$ GHz'% \
             {'re_k' : np.real(sim_wguide.Eig_values[ival]*1e-9),
             'im_k' : np.imag(sim_wguide.Eig_values[ival]*1e-9)}
     else:
         k_str = r'$\Omega/2\pi = %(re_k)f + %(im_k)f i$ GHz'% \
             {'re_k' : np.real(sim_wguide.Eig_values[ival]*1e-9),
             'im_k' : np.imag(sim_wguide.Eig_values[ival]*1e-9)}
     plt.suptitle('Mode #' + str(ival) + '   ' + k_str + '   ' + n_str, fontsize=decorator.get_font_size('title'))

     if pdf_png=='png':
         plt.savefig('%(pre)sfields/%(s)s_S_field_%(i)i%(add)s.png' %
             {'pre' : prefix_str, 's' : EM_AC, 'i' : ival, 'add' : suffix_str})
     elif pdf_png=='pdf':
         plt.savefig('%(pre)sfields/%(s)s_S_field_%(i)i%(add)s.pdf' %
             {'pre' : prefix_str, 's' : EM_AC, 'i' : ival, 'add' : suffix_str}, bbox_inches='tight')
     plt.close()


def save_figure(plt, figfile):
  if figfile[-3:-1]=='png': plt.savefig(figfile)
  else: plt.savefig(figfile, bbox_inches='tight')


def plot_component(v_x, v_y, v_XX, v_YY, plot, label, plps, sim_wguide, ival, comp):
  plt.clf()
  fig = plt.figure(figsize=(15,15))

  decorator = plps['decorator']
  plt.rc('axes', linewidth=decorator.get_axes_property('linewidth'))
  plt.rc('axes', edgecolor=decorator.get_axes_property('edgecolor'))

  ax = plt.subplot(111)
  if comp in ('Et', 'Ht', 'ut'):
    plot_component_quiver(ax, v_x, v_y, plot, plps)
  else: 
    plot_component_axes(ax, v_x, v_y, v_XX, v_YY, plot, label, plps)

  figfile=plot_filename(plps, ival, comp)
  save_figure(plt, figfile)
  plt.close()


#### Standard plotting of spectra #############################################
def plt_mode_fields(sim_wguide, ivals=None, n_points=501, quiver_steps=50, 
                  xlim_min=None, xlim_max=None, ylim_min=None, ylim_max=None,
                  EM_AC='EM_E', num_ticks=None, colorbar=True, contours=False, contour_lst=None,
                  stress_fields=False, pdf_png='png', 
                  prefix_str='', suffix_str='', ticks=False, comps=None, decorator=None):
    """ Plot E or H fields of EM mode, or the AC modes displacement fields.

        Args:
            sim_wguide : A ``Struct`` instance that has had calc_modes calculated

        Keyword Args:
            ivals  (list): mode numbers of modes you wish to plot

            n_points  (int): The number of points across unitcell to
                interpolate the field onto

            xlim_min  (float): Limit plotted xrange to xlim_min:(1-xlim_max) of unitcell

            xlim_max  (float): Limit plotted xrange to xlim_min:(1-xlim_max) of unitcell

            ylim_min  (float): Limit plotted yrange to ylim_min:(1-ylim_max) of unitcell

            ylim_max  (float): Limit plotted yrange to ylim_min:(1-ylim_max) of unitcell

            EM_AC  (str): Either 'EM' or 'AC' modes

            num_ticks  (int): Number of tick marks

            contours  (bool): Controls contours being overlaid on fields

            contour_lst  (list): Specify contour values

            stress_fields  (bool): Calculate acoustic stress fields

            pdf_png  (str): File type to save, either 'png' or 'pdf' 

            prefix_str  (str): Add a string to start of file name
            
            suffix_str  (str): Add a string to end of file name.
    """

    if EM_AC is not 'EM_E' and EM_AC is not 'EM_H' and EM_AC is not 'AC':
        raise ValueError("EM_AC must be either 'AC', 'EM_E' or 'EM_H'.")

    # Calculate the magnetic field from the electric field
    if EM_AC == 'EM_H':
        nnodes = 6
        sim_wguide.sol1_H = NumBAT.h_mode_field_ez(sim_wguide.k_0, sim_wguide.num_modes, 
            sim_wguide.n_msh_el, sim_wguide.n_msh_pts, nnodes, sim_wguide.table_nod, 
            sim_wguide.x_arr, sim_wguide.Eig_values, sim_wguide.sol1)


    plt.clf()

    # field mapping
    x_tmp = []
    y_tmp = []
    for i in np.arange(sim_wguide.n_msh_pts):
        x_tmp.append(sim_wguide.x_arr[0,i])
        y_tmp.append(sim_wguide.x_arr[1,i])
    x_min = np.min(x_tmp); x_max=np.max(x_tmp)
    y_min = np.min(y_tmp); y_max=np.max(y_tmp)
    area = abs((x_max-x_min)*(y_max-y_min))
    n_pts_x = int(n_points*abs(x_max-x_min)/np.sqrt(area))
    n_pts_y = int(n_points*abs(y_max-y_min)/np.sqrt(area))
    v_x = np.zeros(n_pts_x*n_pts_y)
    v_y = np.zeros(n_pts_x*n_pts_y)
    v_XX=None
    v_YY=None
    if contours:
        v_XX, v_YY = np.meshgrid(range(n_pts_x), range(n_pts_y))
    i=0
    for x in np.linspace(x_min,x_max,n_pts_x):
        for y in np.linspace(y_min,y_max,n_pts_y):
            v_x[i] = x
            v_y[i] = y
            i+=1
    v_x = np.array(v_x)
    v_y = np.array(v_y)

    # unrolling data for the interpolators
    table_nod = sim_wguide.table_nod.T
    x_arr = sim_wguide.x_arr.T

    if ivals:
        ival_range = ivals
    else:
        ival_range = range(len(sim_wguide.Eig_values))

    for ival in ival_range:
        # dense triangulation with multiple points
        v_x6p = np.zeros(6*sim_wguide.n_msh_el)
        v_y6p = np.zeros(6*sim_wguide.n_msh_el)
        v_Ex6p = np.zeros(6*sim_wguide.n_msh_el, dtype=np.complex128)
        v_Ey6p = np.zeros(6*sim_wguide.n_msh_el, dtype=np.complex128)
        v_Ez6p = np.zeros(6*sim_wguide.n_msh_el, dtype=np.complex128)
        v_triang6p = []

        i = 0
        for i_el in np.arange(sim_wguide.n_msh_el):
            # triangles
            idx = np.arange(6*i_el, 6*(i_el+1))
            triangles = [[idx[0], idx[3], idx[5]],
                         [idx[1], idx[4], idx[3]],
                         [idx[2], idx[5], idx[4]],
                         [idx[3], idx[4], idx[5]]]
            v_triang6p.extend(triangles)

            for i_node in np.arange(6):
                # index for the coordinates
                i_ex = table_nod[i_el, i_node]-1
                # values
                v_x6p[i] = x_arr[i_ex, 0]
                v_y6p[i] = x_arr[i_ex, 1]

                if EM_AC == 'EM_E' or EM_AC == 'AC' :
                    v_Ex6p[i] = sim_wguide.sol1[0,i_node,ival,i_el]
                    v_Ey6p[i] = sim_wguide.sol1[1,i_node,ival,i_el]
                    v_Ez6p[i] = sim_wguide.sol1[2,i_node,ival,i_el]
                if EM_AC == 'EM_H':
                    v_Ex6p[i] = sim_wguide.sol1_H[0,i_node,ival,i_el]
                    v_Ey6p[i] = sim_wguide.sol1_H[1,i_node,ival,i_el]
                    v_Ez6p[i] = sim_wguide.sol1_H[2,i_node,ival,i_el]
                i += 1

        v_E6p = np.sqrt(np.abs(v_Ex6p)**2 + np.abs(v_Ey6p)**2 + np.abs(v_Ez6p)**2)

        ### Interpolate onto triangular grid - honest to FEM elements
        # dense triangulation with unique points
        v_triang1p = []
        for i_el in np.arange(sim_wguide.n_msh_el):
            # triangles
            triangles = [[table_nod[i_el,0]-1,table_nod[i_el,3]-1,table_nod[i_el,5]-1],
                         [table_nod[i_el,1]-1,table_nod[i_el,4]-1,table_nod[i_el,3]-1],
                         [table_nod[i_el,2]-1,table_nod[i_el,5]-1,table_nod[i_el,4]-1],
                         [table_nod[i_el,3]-1,table_nod[i_el,4]-1,table_nod[i_el,5]-1]]
            v_triang1p.extend(triangles)

        # triangulations
        triang6p = matplotlib.tri.Triangulation(v_x6p,v_y6p,v_triang6p)
        triang1p = matplotlib.tri.Triangulation(x_arr[:,0],x_arr[:,1],v_triang1p)
        
        #fields are called Ex, Ey, etc regardless of whether we are plotting E, H, or u/S

        # building interpolators: triang1p for the finder, triang6p for the values
        #TODO: could be more efficient only interpolating the fields which are ultimately to be used?
        finder = matplotlib.tri.TrapezoidMapTriFinder(triang1p)
        ReEx = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ex6p.real,trifinder=finder)
        ImEx = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ex6p.imag,trifinder=finder)
        ReEy = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ey6p.real,trifinder=finder)
        ImEy = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ey6p.imag,trifinder=finder)
        ReEz = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ez6p.real,trifinder=finder)
        ImEz = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ez6p.imag,trifinder=finder)
        AbsE = matplotlib.tri.LinearTriInterpolator(triang6p,v_E6p,trifinder=finder)
        # interpolated fields
        m_ReEx = ReEx(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ReEy = ReEy(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ReEz = ReEz(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ImEx = ImEx(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ImEy = ImEy(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ImEz = ImEz(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_AbsE = AbsE(v_x,v_y).reshape(n_pts_x,n_pts_y)

        v_x_q = v_x.reshape(n_pts_x,n_pts_y)
        v_y_q = v_y.reshape(n_pts_x,n_pts_y)
        v_x_q = v_x_q[0::quiver_steps,0::quiver_steps]
        v_y_q = v_y_q[0::quiver_steps,0::quiver_steps]
        m_ReEx_q = m_ReEx[0::quiver_steps,0::quiver_steps]
        m_ReEy_q = m_ReEy[0::quiver_steps,0::quiver_steps]
        m_ImEx_q = m_ImEx[0::quiver_steps,0::quiver_steps]
        m_ImEy_q = m_ImEy[0::quiver_steps,0::quiver_steps]

        ### No longer needed as imshow is using origin=lower
        # Flip y order as imshow has origin at top left
        #v_plots = [m_ReEx[:,::-1],m_ReEy[:,::-1],m_ReEz[:,::-1],m_ImEx[:,::-1],m_ImEy[:,::-1],m_ImEz[:,::-1],m_AbsE[:,::-1]]

        v_plots = [m_ReEx, m_ReEy, m_ReEz, m_ImEx, m_ImEy, m_ImEz, m_AbsE]
        vq_plots = [m_ReEx_q, m_ReEy_q, m_ImEx_q, m_ImEy_q]

        if EM_AC=='EM_E':
            v_labels = [r"Re($E_x$)",r"Re($E_y$)",r"Re($E_z$)",r"Im($E_x$)",r"Im($E_y$)",r"Im($E_z$)",r"$|E|$"]
        elif EM_AC == 'EM_H':
            v_labels = [r"Re($H_x$)",r"Re($H_y$)",r"Re($H_z$)",r"Im($H_x$)",r"Im($H_y$)",r"Im($H_z$)",r"$|H|$"]
        else:
            v_labels = [r"Re($u_x$)",r"Re($u_y$)",r"Re($u_z$)",r"Im($u_x$)",r"Im($u_y$)",r"Im($u_z$)",r"$|u|$"]



        if decorator == None: decorator = FieldDecorator()

        plot_params={ 'xlim_min': xlim_min, 'xlim_max': xlim_max, 'ylim_min': ylim_min, 
                     'ylim_max': ylim_max, 'ticks': ticks, 'num_ticks':num_ticks,
                      'colorbar':colorbar, 'contours':contours, 'contour_lst':contour_lst, 'EM_AC':EM_AC,
                      'prefix_str': prefix_str, 'suffix_str': suffix_str, 'pdf_png': pdf_png, 'stress_fields':stress_fields, 'decorator': decorator }

        if not os.path.exists("%sfields" % prefix_str): os.mkdir("%sfields" % prefix_str)

        decorator._set_for_multi()
        plot_all_components(v_x, v_y, v_x_q, v_y_q, v_XX, v_YY, v_plots, vq_plots, v_labels, plot_params, sim_wguide, ival)

        decorator._set_for_single()
        if comps!=None:
          for comp in comps:
            if   comp=='Ex'   and EM_AC=='EM_E': plot_component(v_x, v_y, v_XX, v_YY, v_plots[0], v_labels[0], plot_params, sim_wguide, ival, comp)
            elif comp=='Ey'   and EM_AC=='EM_E': plot_component(v_x, v_y, v_XX, v_YY, v_plots[1], v_labels[1], plot_params, sim_wguide, ival, comp)
            elif comp=='Ez'   and EM_AC=='EM_E': plot_component(v_x, v_y, v_XX, v_YY, v_plots[5], v_labels[5], plot_params, sim_wguide, ival, comp)
            elif comp=='Eabs' and EM_AC=='EM_E': plot_component(v_x, v_y, v_XX, v_YY, v_plots[6], v_labels[6], plot_params, sim_wguide, ival, comp)
            elif comp=='Et'   and EM_AC=='EM_E': plot_component(v_x_q, v_y_q, v_XX, v_YY, vq_plots, '$(E_x,E_y)$', plot_params, sim_wguide, ival, comp)
            elif comp=='Hx'   and EM_AC=='EM_H': plot_component(v_x, v_y, v_XX, v_YY, v_plots[0], v_labels[0], plot_params, sim_wguide, ival, comp)
            elif comp=='Hy'   and EM_AC=='EM_H': plot_component(v_x, v_y, v_XX, v_YY, v_plots[1], v_labels[1], plot_params, sim_wguide, ival, comp)
            elif comp=='Hz'   and EM_AC=='EM_H': plot_component(v_x, v_y, v_XX, v_YY, v_plots[5], v_labels[5], plot_params, sim_wguide, ival, comp)
            elif comp=='Habs' and EM_AC=='EM_H': plot_component(v_x, v_y, v_XX, v_YY, v_plots[6], v_labels[6], plot_params, sim_wguide, ival, comp)
            elif comp=='Ht'   and EM_AC=='EM_H': plot_component(v_x_q, v_y_q, v_XX, v_YY, vq_plots, '$(H_x,H_y)$', plot_params, sim_wguide, ival, comp)
            elif comp=='ux'   and EM_AC=='AC':   plot_component(v_x, v_y, v_XX, v_YY, v_plots[0], v_labels[0], plot_params, sim_wguide, ival, comp)
            elif comp=='uy'   and EM_AC=='AC':   plot_component(v_x, v_y, v_XX, v_YY, v_plots[1], v_labels[1], plot_params, sim_wguide, ival, comp)
            elif comp=='uz'   and EM_AC=='AC':   plot_component(v_x, v_y, v_XX, v_YY, v_plots[5], v_labels[5], plot_params, sim_wguide, ival, comp)
            elif comp=='uabs' and EM_AC=='AC':   plot_component(v_x, v_y, v_XX, v_YY, v_plots[6], v_labels[6], plot_params, sim_wguide, ival, comp)
            elif comp=='ut'   and EM_AC=='AC':   plot_component(v_x_q, v_y_q, v_XX, v_YY, vq_plots, '$(u_x, u_y)$', plot_params, sim_wguide, ival, comp)
     


#### Plot mesh #############################################
def plot_msh(x_arr, prefix_str='', suffix_str=''):
    """ Plot EM mode fields.

        Args:
            sim_wguide : A ``Struct`` instance that has had calc_modes calculated

        Keyword Args:
            n_points  (int): The number of points across unitcell to \
                interpolate the field onto.
    """

    plt.clf()
    plt.figure(figsize=(13,13))
    ax = plt.subplot(1,1,1)
    for node in range(np.shape(x_arr)[1]):
        plt.plot(x_arr[0,node], x_arr[1,node], 'og')
    ax.set_aspect('equal')
    plt.savefig('%(pre)smsh_%(add)s.pdf' %
        {'pre' : prefix_str, 'add' : suffix_str}, bbox_inches='tight')
    plt.close()



### Plot nodal arrangement on mesh triangle.
# plt.figure(figsize=(13,13))
# el = 1
# plt.clf()
# for i in range(0,6):
#     print table_nod[i][el] - 1
#     x = x_arr[0,table_nod[i][el] - 1]
#     y = x_arr[1,table_nod[i][el] - 1]
#     print 'x1 = ', x_arr[0,table_nod[i][el] - 1]
#     print 'y1 = ', x_arr[1,table_nod[i][el] - 1]
#     plt.plot(x, y, 'o')
#     plt.text(x+0.001, y+0.001, str(i))
# plt.savefig('triangle_%i.png' %el)
# plt.close()
