"""
    plotting.py is a subroutine of NumBAT that contains numerous plotting
    routines.

    Copyright (C) 2016  Bjorn Sturmberg
"""

import os
import numpy as np
from scipy import sqrt
import subprocess
from matplotlib.mlab import griddata
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from mpl_toolkits.axes_grid1 import make_axes_locatable


# font = {'family' : 'normal',
#         'weight' : 'bold',
#         'size'   : 18}
# matplotlib.rc('font', **font)
linesstrength = 2.5
title_font = 25


#### Natural constants ########################################################
ASTM15_tot_I   = 900.084            # Integral ASTM 1.5 solar irradiance W/m**2
Plancks_h      = 6.62606957*1e-34   # Planck's constant
speed_c        = 299792458          # Speed of light in vacuum
charge_e       = 1.602176565*1e-19  # Charge of an electron
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


#### Standard plotting of spectra #############################################
def plot_EM_modes(sim_wguide, n_points=500, EM_AC='EM', add_name=''):
    """ Plot EM mode fields.

        Args:
            sim_wguide : A :Struct: instance that has had calc_modes calculated

        Keyword Args:
            n_points  (int): The number of points across unitcell to \
                interpolate the field onto.
    """

    EM_mode_fields = sim_wguide.sol1
    print np.shape(EM_mode_fields)
    print repr(EM_mode_fields[-1,-1,-1,-1])

    # field mapping
    v_x=np.zeros(n_points**2)
    v_y=np.zeros(n_points**2)
    i=0
    x_min=0.0;x_max=1.0
    y_min=-1.0;y_max=0.0
    for x in np.linspace(x_min,x_max,n_points):
        for y in np.linspace(y_min,y_max,n_points):
            v_x[i] = x
            v_y[i] = y
            i+=1
    v_x = np.array(v_x)
    v_y = np.array(v_y)

    # unrolling data for the interpolators
    table_nod = sim_wguide.table_nod.T
    x_arr = sim_wguide.x_arr.T
    # print repr(table_nod)
    # print repr(x_arr)

    for ival in [0]:
    # for ival in range(len(sim_wguide.Eig_value)):
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
                v_Ex6p[i] = EM_mode_fields[0,i_node,ival,i_el]
                v_Ey6p[i] = EM_mode_fields[1,i_node,ival,i_el]
                v_Ez6p[i] = EM_mode_fields[2,i_node,ival,i_el]

                i += 1

        v_E6p = np.sqrt(np.abs(v_Ex6p)**2 +
                        np.abs(v_Ey6p)**2 +
                        np.abs(v_Ez6p)**2)
        # print np.shape(v_x6p)
        # print np.shape(v_Ex6p)
        # print triangles
        # print v_x6p[-1] 
        # print v_y6p[-1] 
        # print v_Ex6p[-1]
        # print v_Ey6p[-1]
        # print v_Ez6p[-1]

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

        # print v_x6p,v_y6p,v_triang6p
        # print x_arr[:,0],x_arr[:,1],v_triang1p

        # building interpolators: triang1p for the finder, triang6p for the values
        finder = matplotlib.tri.TrapezoidMapTriFinder(triang1p)
        ReEx = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ex6p.real,trifinder=finder)
        ImEx = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ex6p.imag,trifinder=finder)
        ReEy = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ey6p.real,trifinder=finder)
        ImEy = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ey6p.imag,trifinder=finder)
        ReEz = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ez6p.real,trifinder=finder)
        ImEz = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ez6p.imag,trifinder=finder)
        AbsE = matplotlib.tri.LinearTriInterpolator(triang6p,v_E6p,trifinder=finder)

        ### plotting
        # interpolated fields
        m_ReEx = ReEx(v_x,v_y).reshape(n_points,n_points)
        m_ReEy = ReEy(v_x,v_y).reshape(n_points,n_points)
        m_ReEz = ReEz(v_x,v_y).reshape(n_points,n_points)
        m_ImEx = ImEx(v_x,v_y).reshape(n_points,n_points)
        m_ImEy = ImEy(v_x,v_y).reshape(n_points,n_points)
        m_ImEz = ImEz(v_x,v_y).reshape(n_points,n_points)
        m_AbsE = AbsE(v_x,v_y).reshape(n_points,n_points)
        v_plots = [m_ReEx,m_ReEy,m_ReEz,m_ImEx,m_ImEy,m_ImEz,m_AbsE]
        v_labels = ["ReEx","ReEy","ReEz","ImEx","ImEy","ImEz","AbsE"]
        print m_ReEx

        # field plots
        plt.clf()
        plt.figure(figsize=(13,13))
        for i_p,plot in enumerate(v_plots):
            ax = plt.subplot(3,3,i_p+1)
            im = plt.imshow(plot.T,cmap='jet');
            # no ticks
            plt.xticks([])
            plt.yticks([])
            # titles
            plt.title(v_labels[i_p],fontsize=title_font)
            # colorbar
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.1)
            cbar = plt.colorbar(im, cax=cax)
            cbar.ax.tick_params(labelsize=title_font-10)
        # plt.tight_layout(1)

        if EM_AC=='EM':
            n_eff = sim_wguide.Eig_value[ival] * sim_wguide.wl_norm() / (2*np.pi)
            if np.imag(sim_wguide.Eig_value[ival]) < 0:
                k_str = r'k$_z = %(re_k)f6 %(im_k)f6 i$'% \
                    {'re_k' : np.real(sim_wguide.Eig_value[ival]), 
                    'im_k' : np.imag(sim_wguide.Eig_value[ival])}
                n_str = r'n$_{eff} = %(re_k)f6 %(im_k)f6 i$'% \
                    {'re_k' : np.real(n_eff), 'im_k' : np.imag(n_eff)}
            else:
                k_str = r'k$_z = %(re_k)f6 + %(im_k)f6 i$'% \
                    {'re_k' : np.real(sim_wguide.Eig_value[ival]), 
                    'im_k' : np.imag(sim_wguide.Eig_value[ival])}
                n_str = r'n$_{eff} = %(re_k)f6 + %(im_k)f6 i$'% \
                    {'re_k' : np.real(n_eff), 'im_k' : np.imag(n_eff)}
            plt.text(10, 0.3, n_str, fontsize=title_font)
        elif EM_AC=='AC':
            if np.imag(sim_wguide.Eig_value[ival]) < 0:
                k_str = r'$\Omega = %(re_k)f6 %(im_k)f6 i$'% \
                    {'re_k' : np.real(sim_wguide.Eig_value[ival]), 
                    'im_k' : np.imag(sim_wguide.Eig_value[ival])}
            else:
                k_str = r'$\Omega = %(re_k)f6 + %(im_k)f6 i$'% \
                    {'re_k' : np.real(sim_wguide.Eig_value[ival]), 
                    'im_k' : np.imag(sim_wguide.Eig_value[ival])}
        else:
            raise ValueError, "EM_AC must be either 'AC' or 'EM'."
        plt.text(10, 0.5, k_str, fontsize=title_font)
        
        plt.savefig('E_field_%(i)i%(add)s.png' % 
            {'i' : ival, 'add' : add_name}, bbox_inches='tight')


#### Plot mesh #############################################
def plot_msh(x_arr, add_name=''):
    """ Plot EM mode fields.

        Args:
            sim_wguide : A :Struct: instance that has had calc_modes calculated

        Keyword Args:
            n_points  (int): The number of points across unitcell to \
                interpolate the field onto.
    """

    plt.clf()
    plt.figure(figsize=(13,13))
    ax = plt.subplot(1,1,1)
    for node in range(np.shape(x_arr)[1]):
        plt.plot(x_arr[0,node], x_arr[1,node], 'o')
    ax.set_aspect('equal')
    plt.savefig('msh_%(add)s.pdf' % 
        {'add' : add_name}, bbox_inches='tight')



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
