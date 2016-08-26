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
def plt_mode_fields(sim_wguide, n_points=1000, xlim=None, ylim=None,
                  EM_AC='EM', pdf_png='png', add_name=''):
    """ Plot EM mode fields.
    NOTE: z component of EM field needs comes scaled by 1/(i beta), 
    which must be reintroduced!

        Args:
            sim_wguide : A :Struct: instance that has had calc_modes calculated

        Keyword Args:
            n_points  (int): The number of points across unitcell to \
                interpolate the field onto.
    """

    if EM_AC is not 'EM' and EM_AC is not 'AC':
        raise ValueError, "EM_AC must be either 'AC' or 'EM'."

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
    v_x=np.zeros(n_pts_x*n_pts_y)
    v_y=np.zeros(n_pts_x*n_pts_y)
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

    # for ival in [0]:
    for ival in range(len(sim_wguide.Eig_value)):
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
                v_Ex6p[i] = sim_wguide.sol1[0,i_node,ival,i_el]
                v_Ey6p[i] = sim_wguide.sol1[1,i_node,ival,i_el]
                if EM_AC == 'EM':
    # Note physical z-comp of EM modes is -i beta E_z, where E_z is FEM output sol
                    v_Ez6p[i] = sim_wguide.sol1[2,i_node,ival,i_el]*-1j*sim_wguide.Eig_value[ival]
                else:
                    v_Ez6p[i] = sim_wguide.sol1[2,i_node,ival,i_el]
                i += 1

        v_E6p = np.sqrt(np.abs(v_Ex6p)**2 +
                        np.abs(v_Ey6p)**2 +
                        np.abs(v_Ez6p)**2)

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
        m_ReEx = ReEx(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ReEy = ReEy(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ReEz = ReEz(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ImEx = ImEx(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ImEy = ImEy(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ImEz = ImEz(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_AbsE = AbsE(v_x,v_y).reshape(n_pts_x,n_pts_y)
        v_plots = [m_ReEx,m_ReEy,m_ReEz,m_ImEx,m_ImEy,m_ImEz,m_AbsE]
        if EM_AC=='EM':
            v_labels = ["Re(E_x)","Re(E_y)","Re(E_z)","Im(E_x)","Im(E_y)","Im(E_z)","Abs(E)"]
        else:
            v_labels = ["Re(u_x)","Re(u_y)","Re(u_z)","Im(u_x)","Im(u_y)","Im(u_z)","Abs(u)"]


        # field plots
        plt.clf()
        plt.figure(figsize=(13,13))
        for i_p,plot in enumerate(v_plots):
            ax = plt.subplot(3,3,i_p+1)
            # im = plt.imshow(plot.T,cmap='viridis');
            im = plt.imshow(plot.T,cmap='inferno');
            # ax.set_aspect('equal')
            # no ticks
            plt.xticks([])
            plt.yticks([])
            # limits
            if xlim:
                ax.set_xlim(xlim*n_points,(1-xlim)*n_points)
            if ylim:
                ax.set_ylim(ylim*n_points,(1-ylim)*n_points)
            # titles
            plt.title(v_labels[i_p],fontsize=title_font)
            # colorbar
            divider = make_axes_locatable(ax)
            cax = divider.append_axes("right", size="5%", pad=0.1)
            cbar = plt.colorbar(im, cax=cax)
            cbar.ax.tick_params(labelsize=title_font-10)

        q_step = 100
        v_x_q = v_x.reshape(n_pts_x,n_pts_y)
        v_y_q = v_y.reshape(n_pts_x,n_pts_y)
        v_x_q = v_x_q[0::q_step,0::q_step]
        v_y_q = v_y_q[0::q_step,0::q_step]
        m_ReEx_q = m_ReEx[0::q_step,0::q_step]
        m_ReEy_q = m_ReEy[0::q_step,0::q_step]
        m_ImEx_q = m_ImEx[0::q_step,0::q_step]
        m_ImEy_q = m_ImEy[0::q_step,0::q_step]
        ax = plt.subplot(3,3,i_p+2)
        plt.quiver(v_x_q, v_y_q, (m_ReEx_q+m_ImEx_q), (m_ReEy_q+m_ImEy_q),      # data
                   np.sqrt(np.real((m_ReEx_q+1j*m_ImEx_q)*(m_ReEx_q-1j*m_ImEx_q)+(m_ReEy_q+1j*m_ImEy_q)*(m_ReEy_q-1j*m_ImEy_q))),  #colour the arrows based on this array
                   cmap='inferno',     # colour map
                   pivot='mid',
                   headlength=5)        # length of the arrows
        # no ticks
        ax.set_aspect('equal')
        plt.xticks([])
        plt.yticks([])
        # limits
        if xlim:
            ax.set_xlim(xlim*n_points,(1-xlim)*n_points)
        if ylim:
            ax.set_ylim(ylim*n_points,(1-ylim)*n_points)
        # titles
        plt.title('Transverse',fontsize=title_font)
        # divider = make_axes_locatable(ax)
        # cax = divider.append_axes("right", size="5%", pad=0.1)
        # cbar = plt.colorbar(im, cax=cax)



        if EM_AC=='EM':
            n_eff = sim_wguide.Eig_value[ival] * sim_wguide.wl_m / (2*np.pi)
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
            # plt.text(10, 0.3, n_str, fontsize=title_font)
        else:
            n_str = ''
            if np.imag(sim_wguide.Eig_value[ival]) < 0:
                k_str = r'$\Omega/2\pi = %(re_k)f6 %(im_k)f6 i$ GHz'% \
                    {'re_k' : np.real(sim_wguide.Eig_value[ival]*1e-9),
                    'im_k' : np.imag(sim_wguide.Eig_value[ival]*1e-9)}
            else:
                k_str = r'$\Omega/2\pi = %(re_k)f6 + %(im_k)f6 i$ GHz'% \
                    {'re_k' : np.real(sim_wguide.Eig_value[ival]*1e-9),
                    'im_k' : np.imag(sim_wguide.Eig_value[ival]*1e-9)}
        # plt.text(10, 0.5, k_str, fontsize=title_font)
        plt.suptitle(k_str + n_str, fontsize=title_font)

        if not os.path.exists("fields"):
            os.mkdir("fields")
        if pdf_png=='png':
            plt.savefig('fields/%(s)s_field_%(i)i%(add)s.png' %
                {'s' : EM_AC, 'i' : ival, 'add' : add_name}) 
                #, bbox_inches='tight') - this caused error in Q calc... ?
        elif pdf_png=='pdf':
            plt.savefig('fields/%(s)s_field_%(i)i%(add)s.pdf' %
                {'s' : EM_AC, 'i' : ival, 'add' : add_name}, bbox_inches='tight')
        else:
            raise ValueError, "pdf_png must be either 'png' or 'pdf'."
        plt.close()


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
