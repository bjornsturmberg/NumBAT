# mode_calcs.py is a subroutine of NumBAT that contains methods to
# calculate the EM and Acoustic modes of a structure.

# Copyright (C) 2017  Bjorn Sturmberg, Kokou Dossou,

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


import time
import numpy as np
from scipy import interpolate
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# from mpl_toolkits.axes_grid1 import make_axes_locatable
import csv

import plotting
from fortran import NumBAT


def gain_and_qs(sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC,
                EM_ival_pump=0, EM_ival_Stokes=0, AC_ival=0, fixed_Q=None, typ_select_out=None):
    r""" Calculate interaction integrals and SBS gain.

        Implements Eqs. 33, 41, 45, 91 of
        Wolff et al. PRA 92, 013836 (2015) doi/10.1103/PhysRevA.92.013836
        These are for Q_photoelastic, Q_moving_boundary, the Acoustic loss "alpha",
        and the SBS gain respectively.

        Note there is a sign error in published Eq. 41. Also, in implementing Eq. 45 we use integration by parts, with a 
        boundary integral term set to zero on physical grounds, and filled in some missing subscripts. We prefer to express
        Eq. 91 with the Lorentzian explicitly visible, which makes it clear how to transform to frequency space.
        The final integrals are

        .. math:: 
            Q^{\rm PE} = -\varepsilon_0 \int_A {\rm d}^2r \sum_{ijkl} \varepsilon^2_r e^{(s)\star}_i e^{(p)}_j p_{ijkl} \partial_k u_l^{*},\\

            Q^{\rm MB} =  \int_C {\rm d \mathbf{r} (\mathbf{u}^{*} \cdot \hat n}) \big[ (\varepsilon_a - \varepsilon_b)  
            \varepsilon_0 ({\rm \hat n \times \mathbf{e}}) \cdot ({\rm \hat n \times \mathbf{e}}) - 
            (\varepsilon_a^{-1} - \varepsilon_b^{-1})  \varepsilon_0^{-1} ({\rm \hat n \cdot \mathbf{d}}) 
            \cdot ({\rm \hat n \cdot \mathbf{d}}) \big],\\

            \alpha = \frac{\Omega^2}{\mathcal{E}_{ac}} \int {\rm d}^2r \sum_{ijkl} \partial_i u_j^{*} \eta_{ijkl} \partial_k u_l,\\

            \Gamma =  \frac{2 \omega \Omega {\rm Re} (Q_1 Q_1^*)}{P_p P_s \mathcal{E}_{ac}} \frac{1}{\alpha} \frac{\alpha^2}{\alpha^2 + \kappa^2}.
  

        Args:
            sim_EM_pump  (``Simmo`` object): Contains all info on pump EM modes

            sim_EM_Stokes  (``Simmo`` object): Contains all info on Stokes EM modes

            sim_AC  (``Simmo`` object): Contains all info on AC modes

            k_AC  (float): Propagation constant of acoustic modes.

        Keyword Args:
            EM_ival_pump  (int/string): Specify mode number of EM mode 1 (pump mode)
                to calculate interactions for.
                Numbering is python index so runs from 0 to num_EM_modes-1,
                with 0 being fundamental mode (largest prop constant).
                Can also set to 'All' to include all modes.

            EM_ival_Stokes  (int/string): Specify mode number of EM mode 2 (stokes mode)
                to calculate interactions for.
                Numbering is python index so runs from 0 to num_EM_modes-1,
                with 0 being fundamental mode (largest prop constant).
                Can also set to 'All' to include all modes.

            AC_ival  (int/string): Specify mode number of AC mode
                to calculate interactions for.
                Numbering is python index so runs from 0 to num_AC_modes-1,
                with 0 being fundamental mode (largest prop constant).
                Can also set to 'All' to include all modes.

            fixed_Q  (int): Specify a fixed Q-factor for the AC modes, rather than
                calculating the acoustic loss (alpha).

        Returns:
            SBS_gain  : The SBS gain including both photoelastic and moving boundary contributions. 
                        Note this will be negative for backwards SBS because gain is expressed as 
                        gain in power as move along z-axis in positive direction, but the Stokes
                        waves experience gain as they propagate in the negative z-direction.
                        Dimensions = [num_modes_EM_Stokes,num_modes_EM_pump,num_modes_AC].

            SBS_gain_PE  : The SBS gain for only the photoelastic effect.
                           The comment about negative gain (see SBS_gain above) holds here also.
                           Dimensions = [num_modes_EM_Stokes,num_modes_EM_pump,num_modes_AC].
            
            SBS_gain_MB  : The SBS gain for only the moving boundary effect. 
                           The comment about negative gain (see SBS_gain above) holds here also.
                           Dimensions = [num_modes_EM_Stokes,num_modes_EM_pump,num_modes_AC].

            alpha  : The acoustic loss for each mode. Dimensions = [num_modes_AC].
    """

    # Notes about internals of fortran integration
    # Calc overlap of basis functions (and PE tensor and epsilon)
    # Then use this multiple times for calc of each mode field values

    # phi is values of Lagrange polynomials (1-6) at that node.
    # grad is value of gradient of Lagrange polynomials (1-6) at that node.
    # i variables refer to E field
    # j variables refer to H field
    # ww weight function
    # coeff numerical integration


    if EM_ival_pump == 'All':
        EM_ival_pump_fortran = -1
    else:
        EM_ival_pump_fortran = EM_ival_pump+1  # convert back to Fortran indexing
    if EM_ival_Stokes == 'All':
        EM_ival_Stokes_fortran = -1
    else:
        EM_ival_Stokes_fortran = EM_ival_Stokes+1  # convert back to Fortran indexing
    if AC_ival == 'All':
        AC_ival_fortran = -1
    else:
        AC_ival_fortran = AC_ival+1  # convert back to Fortran indexing

    Fortran_debug = 0
    ncomps = 3
    nnodes = 6
    num_modes_EM_pump = sim_EM_pump.num_modes
    num_modes_EM_Stokes = sim_EM_Stokes.num_modes
    num_modes_AC = sim_AC.num_modes
    n_msh_el_AC = sim_AC.n_msh_el
    trimmed_EM_pump_field = np.zeros((ncomps,nnodes,num_modes_EM_pump,n_msh_el_AC), dtype=complex)
    trimmed_EM_Stokes_field = np.zeros((ncomps,nnodes,num_modes_EM_Stokes,n_msh_el_AC), dtype=complex)
    for el in range(n_msh_el_AC):
        new_el = sim_AC.el_convert_tbl[el]
        for n in range(nnodes):
            for x in range(ncomps):
                for ival in range(num_modes_EM_pump):
                    trimmed_EM_pump_field[x,n,ival,el] = sim_EM_pump.sol1[x,n,ival,new_el]
                for ival in range(num_modes_EM_Stokes):
                    trimmed_EM_Stokes_field[x,n,ival,el] = sim_EM_Stokes.sol1[x,n,ival,new_el]

    # sim_EM_pump.sol1 = trimmed_EM_pump_field
    # sim_EM_pump.n_msh_el = sim_AC.n_msh_el
    # sim_EM_pump.n_msh_pts = sim_AC.n_msh_pts
    # sim_EM_pump.type_el = sim_AC.type_el
    # sim_EM_pump.table_nod = sim_AC.table_nod
    # sim_EM_pump.x_arr = sim_AC.x_arr
    # plotting.plt_mode_fields(sim_EM_pump, EM_AC='EM', prefix_str='int_test-', suffix_str='trim')

    relevant_eps_effs =[]
    for el_typ in range(sim_EM_pump.structure.nb_typ_el):
        if el_typ+1 in sim_AC.typ_el_AC:
            relevant_eps_effs.append(sim_EM_pump.n_list[el_typ]**2)

    print("\n-----------------------------------------------")
    if fixed_Q is None:
        # Calc alpha (loss) Eq. 45
        print("Acoustic loss calc")
        start = time.time()
        try:
            if sim_EM_pump.structure.inc_shape in sim_EM_pump.structure.linear_element_shapes:
                alpha = NumBAT.ac_alpha_int_v2(sim_AC.num_modes,
                    sim_AC.n_msh_el, sim_AC.n_msh_pts, nnodes,
                    sim_AC.table_nod, sim_AC.type_el, sim_AC.x_arr,
                    sim_AC.structure.nb_typ_el_AC, sim_AC.structure.eta_tensor,
                    k_AC, sim_AC.Omega_AC, sim_AC.sol1,
                    # sim_AC.AC_mode_power) # appropriate for alpha in [1/m]
                    sim_AC.AC_mode_energy_elastic) # appropriate for alpha in [1/s]
            else:
                if sim_EM_pump.structure.inc_shape not in sim_EM_pump.structure.curvilinear_element_shapes:
                    print("Warning: ac_alpha_int - not sure if mesh contains curvi-linear elements", 
                        "\n using slow quadrature integration by default.\n\n")
                alpha = NumBAT.ac_alpha_int(sim_AC.num_modes,
                    sim_AC.n_msh_el, sim_AC.n_msh_pts, nnodes,
                    sim_AC.table_nod, sim_AC.type_el, sim_AC.x_arr,
                    sim_AC.structure.nb_typ_el_AC, sim_AC.structure.eta_tensor,
                    k_AC, sim_AC.Omega_AC, sim_AC.sol1,
                    # sim_AC.AC_mode_power, Fortran_debug) # appropriate for alpha in [1/m]
                    sim_AC.AC_mode_energy_elastic, Fortran_debug) # appropriate for alpha in [1/s]
        except KeyboardInterrupt:
            print("\n\n Routine ac_alpha_int interrupted by keyboard.\n\n")
        alpha = np.real(alpha)
        # Q_factors = 0.5*(k_AC/alpha)*np.ones(num_modes_AC) # appropriate for alpha in [1/m]
        Q_factors = 0.5*(sim_AC.Omega_AC/alpha)*np.ones(num_modes_AC) # appropriate for alpha in [1/s]
        end = time.time()
        print("     time (sec.)", (end - start))
    else:
        # factor of a 1/2 because alpha is for power!
        # alpha [1/m] = Omega_AC/(2*vg*fixed_Q) = k_AC/fixed_Q
        # alpha [1/s] = vg * alpha [1/m]
        # alpha [1/s] = Omega_AC/(2*fixed_Q)
        # alpha = 0.5*(k_AC/fixed_Q)*np.ones(num_modes_AC) # appropriate for alpha in [1/m]
        alpha = 0.5*(sim_AC.Omega_AC/fixed_Q)*np.ones(num_modes_AC) # appropriate for alpha in [1/s]
        Q_factors = fixed_Q*np.ones(num_modes_AC)

    linewidth_Hz = alpha/np.pi # SBS linewidth of each resonance in [Hz]

    # Calc Q_photoelastic Eq. 33
    print("Photoelastic calc")
    start = time.time()
    try:
        if sim_EM_pump.structure.inc_shape in sim_EM_pump.structure.linear_element_shapes:
            Q_PE = NumBAT.photoelastic_int_v2(
                sim_EM_pump.num_modes, sim_EM_Stokes.num_modes, sim_AC.num_modes, EM_ival_pump_fortran,
                EM_ival_Stokes_fortran, AC_ival_fortran, sim_AC.n_msh_el,
                sim_AC.n_msh_pts, nnodes,
                sim_AC.table_nod, sim_AC.type_el, sim_AC.x_arr,
                sim_AC.structure.nb_typ_el_AC, sim_AC.structure.p_tensor,
                k_AC, trimmed_EM_pump_field, trimmed_EM_Stokes_field, sim_AC.sol1,
                relevant_eps_effs, Fortran_debug)
        else:
            if sim_EM_pump.structure.inc_shape not in sim_EM_pump.structure.curvilinear_element_shapes:
                print("Warning: photoelastic_int - not sure if mesh contains curvi-linear elements", 
                    "\n using slow quadrature integration by default.\n\n")
            Q_PE = NumBAT.photoelastic_int(
                sim_EM_pump.num_modes, sim_EM_Stokes.num_modes, sim_AC.num_modes, EM_ival_pump_fortran,
                EM_ival_Stokes_fortran, AC_ival_fortran, sim_AC.n_msh_el,
                sim_AC.n_msh_pts, nnodes,
                sim_AC.table_nod, sim_AC.type_el, sim_AC.x_arr,
                sim_AC.structure.nb_typ_el_AC, sim_AC.structure.p_tensor,
                k_AC, trimmed_EM_pump_field, trimmed_EM_Stokes_field, sim_AC.sol1,
                relevant_eps_effs, Fortran_debug)
    except KeyboardInterrupt:
        print("\n\n Routine photoelastic_int interrupted by keyboard.\n\n")
    end = time.time()
    print("     time (sec.)", (end - start))


    # Calc Q_moving_boundary Eq. 41
    typ_select_in = 1 # first element in relevant_eps_effs list, in fortan indexing
    if len(relevant_eps_effs) == 2: typ_select_out = 2
    elif typ_select_out is None: typ_select_out = -1
    print("Moving boundary calc")
    start = time.time()
    try:
        Q_MB = NumBAT.moving_boundary(sim_EM_pump.num_modes, sim_EM_Stokes.num_modes,
            sim_AC.num_modes, EM_ival_pump_fortran, EM_ival_Stokes_fortran,
            AC_ival_fortran, sim_AC.n_msh_el,
            sim_AC.n_msh_pts, nnodes, sim_AC.table_nod, 
            sim_AC.type_el, sim_AC.x_arr,
            sim_AC.structure.nb_typ_el_AC, typ_select_in, typ_select_out,
            trimmed_EM_pump_field, trimmed_EM_Stokes_field, sim_AC.sol1,
            relevant_eps_effs, Fortran_debug)
    except KeyboardInterrupt:
        print("\n\n Routine moving_boundary interrupted by keyboard.\n\n")
    end = time.time()
    print("     time (sec.)", (end - start))
    print("-----------------------------------------------")

    Q = Q_PE + Q_MB

    # Note: sim_EM_pump.omega_EM is the optical angular freq in units of Hz
    # Note: sim_AC.Omega_AC is the acoustic angular freq in units of Hz
    gain = 2*sim_EM_pump.omega_EM*sim_AC.Omega_AC*np.real(Q*np.conj(Q))
    gain_PE = 2*sim_EM_pump.omega_EM*sim_AC.Omega_AC*np.real(Q_PE*np.conj(Q_PE))
    gain_MB = 2*sim_EM_pump.omega_EM*sim_AC.Omega_AC*np.real(Q_MB*np.conj(Q_MB))
    normal_fact = np.zeros((num_modes_EM_Stokes, num_modes_EM_pump, num_modes_AC), dtype=complex)
    for i in range(num_modes_EM_Stokes):
        P1 = sim_EM_Stokes.EM_mode_power[i]
        for j in range(num_modes_EM_pump):
            P2 = sim_EM_pump.EM_mode_power[j]
            for k in range(num_modes_AC):
                # P3 = sim_AC.AC_mode_power[k]
                P3 = sim_AC.AC_mode_energy_elastic[k]
                normal_fact[i, j, k] = P1*P2*P3*alpha[k]
    SBS_gain = np.real(gain/normal_fact)
    SBS_gain_PE = np.real(gain_PE/normal_fact)
    SBS_gain_MB = np.real(gain_MB/normal_fact)

    return SBS_gain, SBS_gain_PE, SBS_gain_MB, linewidth_Hz, Q_factors, alpha


#### Categorise modes by their symmetries #############################################
def symmetries(sim_wguide, n_points=10, negligible_threshold=1e-5):
    """ Plot EM mode fields.

        Args:
            sim_wguide : A ``Struct`` instance that has had calc_modes calculated

        Keyword Args:
            n_points  (int): The number of points across unitcell to \
                interpolate the field onto.
    """

    mode_fields = sim_wguide.sol1

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

    sym_list = []

    for ival in range(len(sim_wguide.Eig_values)):
        # dense triangulation with multiple points
        v_x6p = np.zeros(6*sim_wguide.n_msh_el)
        v_y6p = np.zeros(6*sim_wguide.n_msh_el)
        v_Ex6p = np.zeros(6*sim_wguide.n_msh_el, dtype=np.complex128)
        v_Ey6p = np.zeros(6*sim_wguide.n_msh_el, dtype=np.complex128)
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
                v_Ex6p[i] = mode_fields[0,i_node,ival,i_el]
                v_Ey6p[i] = mode_fields[1,i_node,ival,i_el]
                i += 1

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

        ### plotting
        # interpolated fields
        m_ReEx = ReEx(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ReEy = ReEy(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ImEx = ImEx(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ImEy = ImEy(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_Ex = m_ReEx + 1j*m_ImEx
        m_Ey = m_ReEy + 1j*m_ImEy

        if np.max(np.abs(m_Ex[~np.isnan(m_Ex)])) < negligible_threshold:
            m_Ex = np.zeros(np.shape(m_Ex))
        if np.max(np.abs(m_Ey[~np.isnan(m_Ey)])) < negligible_threshold:
            m_Ey = np.zeros(np.shape(m_Ey))


        m_Ex_ymirror = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        m_Ex_xmirror = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        m_Ex_rotated = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        m_Ey_ymirror = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        m_Ey_xmirror = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        m_Ey_rotated = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        Ex_sigma_y = 0
        Ey_sigma_y = 0
        Ex_sigma_x = 0
        Ey_sigma_x = 0
        Ex_C_2 = 0
        Ey_C_2 = 0

        for ix in range(n_pts_x):
            for iy in range(n_pts_y):
                m_Ex_ymirror[ix,iy] = (m_Ex[ix,n_pts_y-iy-1])
                m_Ey_ymirror[ix,iy] = -1*(m_Ey[ix,n_pts_y-iy-1])
                m_Ex_xmirror[ix,iy] = -1*(m_Ex[n_pts_x-ix-1,iy])
                m_Ey_xmirror[ix,iy] = (m_Ey[n_pts_x-ix-1,iy])
                m_Ex_rotated[ix,iy] = -1*(m_Ex[n_pts_x-ix-1,n_pts_y-iy-1])
                m_Ey_rotated[ix,iy] = -1*(m_Ey[n_pts_x-ix-1,n_pts_y-iy-1])

        Ex_sigma_y = np.sum(np.abs(m_Ex - m_Ex_ymirror))
        Ey_sigma_y = np.sum(np.abs(m_Ey - m_Ey_ymirror))
        Ex_sigma_x = np.sum(np.abs(m_Ex - m_Ex_xmirror))
        Ey_sigma_x = np.sum(np.abs(m_Ey - m_Ey_xmirror))
        Ex_C_2 = np.sum(np.abs(m_Ex - m_Ex_rotated))
        Ey_C_2 = np.sum(np.abs(m_Ey - m_Ey_rotated))
        sigma_y = (Ex_sigma_y + Ey_sigma_y)/(n_pts_x*n_pts_y)
        sigma_x = (Ex_sigma_x + Ey_sigma_x)/(n_pts_x*n_pts_y)
        C_2 = (Ex_C_2 + Ey_C_2)/(n_pts_x*n_pts_y)

        if abs(C_2) > 0.2:
            C_2_print = -1
        else:
            C_2_print = 1
        if abs(sigma_y) > 0.1:
            sigma_y_print = -1
        else:
            sigma_y_print = 1
        if abs(sigma_x) > 0.1:
            sigma_x_print = -1
        else:
            sigma_x_print = 1
        sym_list.append([C_2_print, sigma_y_print, sigma_x_print])

        # v_plots = [np.real(m_Ex_ymirror),np.real(m_Ey_ymirror),np.real(m_Ez_ymirror),
        #     np.imag(m_Ex_ymirror),np.imag(m_Ey_ymirror),np.imag(m_Ez_ymirror)]
        # # field plots
        # plt.clf()
        # plt.figure(figsize=(13,13))
        # for i_p,plot in enumerate(v_plots):
        #     ax = plt.subplot(3,3,i_p+1)
        #     im = plt.imshow(plot.T,cmap='inferno');
        #     # colorbar
        #     divider = make_axes_locatable(ax)
        #     cax = divider.append_axes("right", size="5%", pad=0.1)
        #     cbar = plt.colorbar(im, cax=cax)
        # plt.savefig('fields/field_%(i)i-ymirror.pdf' %
        #         {'i' : ival}, bbox_inches='tight')
        # plt.close()
        # v_plots = [np.real(m_Ex_xmirror),np.real(m_Ey_xmirror),np.real(m_Ez_xmirror),
        #     np.imag(m_Ex_xmirror),np.imag(m_Ey_xmirror),np.imag(m_Ez_xmirror)]
        # # field plots
        # plt.clf()
        # plt.figure(figsize=(13,13))
        # for i_p,plot in enumerate(v_plots):
        #     ax = plt.subplot(3,3,i_p+1)
        #     im = plt.imshow(plot.T,cmap='inferno');
        #     # colorbar
        #     divider = make_axes_locatable(ax)
        #     cax = divider.append_axes("right", size="5%", pad=0.1)
        #     cbar = plt.colorbar(im, cax=cax)
        # plt.savefig('fields/field_%(i)i-xmirror.pdf' %
        #         {'i' : ival}, bbox_inches='tight')
        # plt.close()
        # v_plots = [np.real(m_Ex_rotated),np.real(m_Ey_rotated),np.real(m_Ez_rotated),
        #     np.imag(m_Ex_rotated),np.imag(m_Ey_rotated),np.imag(m_Ez_rotated)]
        # # field plots
        # plt.clf()
        # plt.figure(figsize=(13,13))
        # for i_p,plot in enumerate(v_plots):
        #     ax = plt.subplot(3,3,i_p+1)
        #     im = plt.imshow(plot.T,cmap='inferno');
        #     # colorbar
        #     divider = make_axes_locatable(ax)
        #     cax = divider.append_axes("right", size="5%", pad=0.1)
        #     cbar = plt.colorbar(im, cax=cax)
        # plt.savefig('fields/field_%(i)i-rotated.pdf' %
        #         {'i' : ival}, bbox_inches='tight')
        # plt.close()
        # v_plots = [m_ReEx,m_ReEy,m_ReEz,
        #     m_ImEx,m_ImEy,m_ImEz]
        # # field plots
        # plt.clf()
        # plt.figure(figsize=(13,13))
        # for i_p,plot in enumerate(v_plots):
        #     ax = plt.subplot(3,3,i_p+1)
        #     im = plt.imshow(plot.T,cmap='inferno');
        #     # colorbar
        #     divider = make_axes_locatable(ax)
        #     cax = divider.append_axes("right", size="5%", pad=0.1)
        #     cbar = plt.colorbar(im, cax=cax)
        # plt.savefig('fields/field_%(i)i.pdf' %
        #         {'i' : ival}, bbox_inches='tight')
        # plt.close()

    return sym_list


def grad_u(dx, dy, u_mat, k_AC):
    """ Take the gradient of field as well as of conjugate of field.
    """

    m_ux = u_mat[0]
    m_uy = u_mat[1]
    m_uz = u_mat[2]
    del_x_ux = np.gradient(m_ux, dx, axis=0)
    del_y_ux = np.gradient(m_ux, dy, axis=1)
    del_x_uy = np.gradient(m_uy, dx, axis=0)
    del_y_uy = np.gradient(m_uy, dy, axis=1)
    del_x_uz = np.gradient(m_uz, dx, axis=0)
    del_y_uz = np.gradient(m_uz, dy, axis=1)
    del_z_ux = 1j*k_AC*m_ux
    del_z_uy = 1j*k_AC*m_uy
    del_z_uz = 1j*k_AC*m_uz
    del_x_ux_star = np.gradient(np.conj(m_ux), dx, axis=0)
    del_y_ux_star = np.gradient(np.conj(m_ux), dy, axis=1)
    del_x_uy_star = np.gradient(np.conj(m_uy), dx, axis=0)
    del_y_uy_star = np.gradient(np.conj(m_uy), dy, axis=1)
    del_x_uz_star = np.gradient(np.conj(m_uz), dx, axis=0)
    del_y_uz_star = np.gradient(np.conj(m_uz), dy, axis=1)
    del_z_ux_star = -1j*k_AC*np.conj(m_ux)
    del_z_uy_star = -1j*k_AC*np.conj(m_uy)
    del_z_uz_star = -1j*k_AC*np.conj(m_uz)

    del_u_mat = np.array([[del_x_ux, del_x_uy, del_x_uz], [del_y_ux, del_y_uy, del_y_uz], [del_z_ux, del_z_uy, del_z_uz]])
    del_u_mat_star = np.array([[del_x_ux_star, del_x_uy_star, del_x_uz_star], [del_y_ux_star, del_y_uy_star, del_y_uz_star], [del_z_ux_star, del_z_uy_star, del_z_uz_star]])

    return del_u_mat, del_u_mat_star


def comsol_fields(data_file, n_points, ival=0):
    """ Load Comsol field data on (assumed) grid mesh.
    """

    with open(data_file, 'rt', encoding='ascii') as csvfile:
        spamreader = csv.reader(csvfile, delimiter=' ')#, quotechar='|')
        for header_rows in range(9):
            next(spamreader)
        x_coord = []; y_coord = []
        u_x = []; u_y = []; u_z = []
        for row in spamreader:
            row = [_f for _f in row if _f]
            row = [float(x) for x in row]
            x_coord.append(row[0])
            y_coord.append(row[1])
            u_x.append(row[(ival*6)+2] + 1j*row[(ival*6)+3])
            u_y.append(row[(ival*6)+4] + 1j*row[(ival*6)+5])
            u_z.append(row[(ival*6)+6] + 1j*row[(ival*6)+7])

    x_coord = np.array(x_coord).reshape(n_points,n_points)
    y_coord = np.array(y_coord).reshape(n_points,n_points)
    x_coord = np.swapaxes(x_coord,0,1)
    y_coord = np.swapaxes(y_coord,0,1)
    u_x = np.array(u_x).reshape(n_points,n_points)
    u_y = np.array(u_y).reshape(n_points,n_points)
    u_z = np.array(u_z).reshape(n_points,n_points)
    u_x = np.swapaxes(u_x,0,1)
    u_y = np.swapaxes(u_y,0,1)
    u_z = np.swapaxes(u_z,0,1)
    field_mat = np.array([u_x, u_y, u_z])

    return x_coord, y_coord, field_mat


def interp_py_fields(sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC, n_points,
              EM_ival_pump=0, EM_ival_Stokes=0, AC_ival=0):
    """ Interpolate fields from FEM mesh to square grid.
    """

    # Trim EM fields to non-vacuum area where AC modes are defined
    num_modes_EM = sim_EM_pump.num_modes
    num_modes_AC = sim_AC.num_modes
    n_msh_el_AC = sim_AC.n_msh_el
    ncomps = 3
    nnodes = 6
    trimmed_EM_field_p = np.zeros((ncomps,nnodes,n_msh_el_AC), dtype=complex)
    trimmed_EM_field_S = np.zeros((ncomps,nnodes,n_msh_el_AC), dtype=complex)
    trimmed_EM_n = np.zeros((1,nnodes,n_msh_el_AC), dtype=complex)
    for el in range(n_msh_el_AC):
        new_el = sim_AC.el_convert_tbl[el]
        for n in range(nnodes):
            for x in range(ncomps):
                trimmed_EM_field_p[x,n,el] = sim_EM_pump.sol1[x,n,EM_ival_pump,new_el]
                trimmed_EM_field_S[x,n,el] = sim_EM_Stokes.sol1[x,n,EM_ival_Stokes,new_el]
            trimmed_EM_n[0,n,el] = sim_EM_pump.ls_material[0,n,new_el]

    # field mapping
    x_tmp = []
    y_tmp = []
    for i in np.arange(sim_AC.n_msh_pts):
        x_tmp.append(sim_AC.x_arr[0,i])
        y_tmp.append(sim_AC.x_arr[1,i])
    x_min = np.min(x_tmp); x_max=np.max(x_tmp)
    y_min = np.min(y_tmp); y_max=np.max(y_tmp)
    area = abs((x_max-x_min)*(y_max-y_min))
    n_pts_x = n_points
    n_pts_y = n_points
    v_x=np.zeros(n_pts_x*n_pts_y)
    v_y=np.zeros(n_pts_x*n_pts_y)
    i=0
    for x in np.linspace(x_min,x_max,n_pts_x):
        for y in np.linspace(y_max,y_min,n_pts_y):
            v_x[i] = x
            v_y[i] = y
            i+=1
    v_x = np.array(v_x)
    v_y = np.array(v_y)

    # unrolling data for the interpolators
    table_nod = sim_AC.table_nod.T
    x_arr = sim_AC.x_arr.T

    # dense triangulation with multiple points
    v_x6p = np.zeros(6*sim_AC.n_msh_el)
    v_y6p = np.zeros(6*sim_AC.n_msh_el)
    v_ux6p = np.zeros(6*sim_AC.n_msh_el, dtype=np.complex128)
    v_uy6p = np.zeros(6*sim_AC.n_msh_el, dtype=np.complex128)
    v_uz6p = np.zeros(6*sim_AC.n_msh_el, dtype=np.complex128)
    v_Ex6p_E_p = np.zeros(6*sim_AC.n_msh_el, dtype=np.complex128)
    v_Ey6p_E_p = np.zeros(6*sim_AC.n_msh_el, dtype=np.complex128)
    v_Ez6p_E_p = np.zeros(6*sim_AC.n_msh_el, dtype=np.complex128)
    v_Ex6p_E_S = np.zeros(6*sim_AC.n_msh_el, dtype=np.complex128)
    v_Ey6p_E_S = np.zeros(6*sim_AC.n_msh_el, dtype=np.complex128)
    v_Ez6p_E_S = np.zeros(6*sim_AC.n_msh_el, dtype=np.complex128)
    v_n = np.zeros(6*sim_AC.n_msh_el, dtype=np.complex128)
    v_triang6p = []

    i = 0
    for i_el in np.arange(sim_AC.n_msh_el):
        for i_node in np.arange(6):
            # index for the coordinates
            i_ex = table_nod[i_el, i_node]-1
            # values
            v_x6p[i] = x_arr[i_ex, 0]
            v_y6p[i] = x_arr[i_ex, 1]
            v_ux6p[i] = sim_AC.sol1[0,i_node,AC_ival,i_el]
            v_uy6p[i] = sim_AC.sol1[1,i_node,AC_ival,i_el]
            v_uz6p[i] = sim_AC.sol1[2,i_node,AC_ival,i_el]
            v_Ex6p_E_p[i] = trimmed_EM_field_p[0,i_node,i_el]
            v_Ey6p_E_p[i] = trimmed_EM_field_p[1,i_node,i_el]
            v_Ez6p_E_p[i] = trimmed_EM_field_p[2,i_node,i_el]
            v_Ex6p_E_S[i] = trimmed_EM_field_S[0,i_node,i_el]
            v_Ey6p_E_S[i] = trimmed_EM_field_S[1,i_node,i_el]
            v_Ez6p_E_S[i] = trimmed_EM_field_S[2,i_node,i_el]
            v_n[i] = trimmed_EM_n[0,i_node,i_el]
            i += 1

    xy = list(zip(v_x6p, v_y6p))
    grid_x, grid_y = np.mgrid[x_min:x_max:n_pts_x*1j, y_min:y_max:n_pts_y*1j]
    # pump mode
    m_ReEx_E = interpolate.griddata(xy, v_Ex6p_E_p.real, (grid_x, grid_y), method='cubic')
    m_ReEy_E = interpolate.griddata(xy, v_Ey6p_E_p.real, (grid_x, grid_y), method='cubic')
    m_ReEz_E = interpolate.griddata(xy, v_Ez6p_E_p.real, (grid_x, grid_y), method='cubic')
    m_ImEx_E = interpolate.griddata(xy, v_Ex6p_E_p.imag, (grid_x, grid_y), method='cubic')
    m_ImEy_E = interpolate.griddata(xy, v_Ey6p_E_p.imag, (grid_x, grid_y), method='cubic')
    m_ImEz_E = interpolate.griddata(xy, v_Ez6p_E_p.imag, (grid_x, grid_y), method='cubic')
    m_Ex_E = m_ReEx_E + 1j*m_ImEx_E
    m_Ey_E = m_ReEy_E + 1j*m_ImEy_E
    m_Ez_E = m_ReEz_E + 1j*m_ImEz_E
    m_Ex_E = m_Ex_E.reshape(n_pts_x,n_pts_y)
    m_Ey_E = m_Ey_E.reshape(n_pts_x,n_pts_y)
    m_Ez_E = m_Ez_E.reshape(n_pts_x,n_pts_y)
    E_mat_p = np.array([m_Ex_E, m_Ey_E, m_Ez_E])
    # Stokes mode
    m_ReEx_E = interpolate.griddata(xy, v_Ex6p_E_S.real, (grid_x, grid_y), method='cubic')
    m_ReEy_E = interpolate.griddata(xy, v_Ey6p_E_S.real, (grid_x, grid_y), method='cubic')
    m_ReEz_E = interpolate.griddata(xy, v_Ez6p_E_S.real, (grid_x, grid_y), method='cubic')
    m_ImEx_E = interpolate.griddata(xy, v_Ex6p_E_S.imag, (grid_x, grid_y), method='cubic')
    m_ImEy_E = interpolate.griddata(xy, v_Ey6p_E_S.imag, (grid_x, grid_y), method='cubic')
    m_ImEz_E = interpolate.griddata(xy, v_Ez6p_E_S.imag, (grid_x, grid_y), method='cubic')
    m_Ex_E = m_ReEx_E + 1j*m_ImEx_E
    m_Ey_E = m_ReEy_E + 1j*m_ImEy_E
    m_Ez_E = m_ReEz_E + 1j*m_ImEz_E
    m_Ex_E = m_Ex_E.reshape(n_pts_x,n_pts_y)
    m_Ey_E = m_Ey_E.reshape(n_pts_x,n_pts_y)
    m_Ez_E = m_Ez_E.reshape(n_pts_x,n_pts_y)
    E_mat_S = np.array([m_Ex_E, m_Ey_E, m_Ez_E])
    # AC mode
    m_Reux = interpolate.griddata(xy, v_ux6p.real, (grid_x, grid_y), method='cubic')
    m_Reuy = interpolate.griddata(xy, v_uy6p.real, (grid_x, grid_y), method='cubic')
    m_Reuz = interpolate.griddata(xy, v_uz6p.real, (grid_x, grid_y), method='cubic')
    m_Imux = interpolate.griddata(xy, v_ux6p.imag, (grid_x, grid_y), method='cubic')
    m_Imuy = interpolate.griddata(xy, v_uy6p.imag, (grid_x, grid_y), method='cubic')
    m_Imuz = interpolate.griddata(xy, v_uz6p.imag, (grid_x, grid_y), method='cubic')
    m_ux = m_Reux + 1j*m_Imux
    m_uy = m_Reuy + 1j*m_Imuy
    m_uz = m_Reuz + 1j*m_Imuz
    m_ux = m_ux.reshape(n_pts_x,n_pts_y)
    m_uy = m_uy.reshape(n_pts_x,n_pts_y)
    m_uz = m_uz.reshape(n_pts_x,n_pts_y)
    u_mat = np.array([m_ux, m_uy, m_uz])

    dx = grid_x[-1,0] - grid_x[-2,0]
    dy = grid_y[0,-1] - grid_y[0,-2]
    del_u_mat, del_u_mat_star = grad_u(dx, dy, u_mat, k_AC)

    m_Ren = interpolate.griddata(xy, v_n.real, (grid_x, grid_y), method='cubic')
    m_Imn = interpolate.griddata(xy, v_n.imag, (grid_x, grid_y), method='cubic')
    m_n = m_Ren + 1j*m_Imn
    m_n = m_n.reshape(n_pts_x,n_pts_y)

    return n_pts_x, n_pts_y, dx, dy, E_mat_p, E_mat_S, u_mat, del_u_mat, del_u_mat_star, m_n


def grid_integral(m_n, sim_AC_structure, sim_AC_Omega_AC, n_pts_x, n_pts_y, 
                  dx, dy, E_mat_p, E_mat_S, u_mat, del_u_mat, del_u_mat_star, AC_ival):
    """ Quadrature integration of AC energy density, AC loss (alpha), and PE gain.
    """

    # AC energy density integral
    F_AC_energy = 0
    for i in range(3):
        integrand_AC = np.conj(u_mat[i])*u_mat[i]*sim_AC_structure.rho
        # do a 1-D integral over every row
        I = np.zeros( n_pts_x )
        for r in range(n_pts_x):
            I[r] = np.trapz( np.real(integrand_AC[r,:]), dx=dy )
        # then an integral over the result
        F_AC_energy += np.trapz( I, dx=dx )
        # Adding imag comp
        I = np.zeros( n_pts_x )
        for r in range(n_pts_x):
            I[r] = np.trapz( np.imag(integrand_AC[r,:]), dx=dy )
        F_AC_energy += 1j*np.trapz( I, dx=dx )
    energy_py = 2*F_AC_energy*sim_AC_Omega_AC[AC_ival]**2

    # AC loss (alpha) integral
    F_alpha = 0
    for i in range(3):
        for k in range(3):
            for l in range(3):
                for j in range(3):
                    integrand = del_u_mat[i,j]*del_u_mat_star[k,l]*sim_AC_structure.eta_tensor[i,j,k,l]
                    I = np.zeros( n_pts_x )
                    for r in range(n_pts_x):
                        I[r] = np.trapz( np.real(integrand[r,:]), dx=dy )
                    F_alpha += np.trapz( I, dx=dx )
                    I = np.zeros( n_pts_x )
                    for r in range(n_pts_x):
                        I[r] = np.trapz( np.imag(integrand[r,:]), dx=dy )
                    F_alpha += 1j*np.trapz( I, dx=dx )
    alpha_py = np.real(F_alpha*sim_AC_Omega_AC[AC_ival]**2/energy_py)

    # PE gain integral
    eps_0 = 8.854187817e-12
    F_PE = 0
    for i in range(3):
        for k in range(3):
            for l in range(3):
                for j in range(3):
                    # integrand_PE = relevant_eps_effs[0]**2 * E_mat_p[j]*np.conj(E_mat_S[i])*sim_AC_structure.p_tensor[i,j,k,l]*del_u_mat_star[k,l]
                    integrand_PE = m_n**4 * E_mat_p[j]*np.conj(E_mat_S[i])*sim_AC_structure.p_tensor[i,j,k,l]*del_u_mat_star[k,l]
                    I = np.zeros( n_pts_x )
                    for r in range(n_pts_x):
                        I[r] = np.trapz( np.real(integrand_PE[r,:]), dx=dy )
                    F_PE += np.trapz( I, dx=dx )
                    I = np.zeros( n_pts_x )
                    for r in range(n_pts_x):
                        I[r] = np.trapz( np.imag(integrand_PE[r,:]), dx=dy )
                    F_PE += 1j*np.trapz( I, dx=dx )
    Q_PE_py = F_PE*eps_0

    return energy_py, alpha_py, Q_PE_py


def gain_python(sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC, comsol_data_file, comsol_ivals=1):
    """ Calculate interaction integrals and SBS gain in python.
        Load in acoustic mode displacement and calculate gain from this also.
    """

    num_modes_EM = sim_EM_pump.num_modes
    num_modes_AC = sim_AC.num_modes
    EM_ival_pump = 0
    EM_ival_Stokes = 0

    n_points = 100
    n_points_comsol_data = 100

    # relevant_eps_effs =[]
    # for el_typ in range(sim_EM_pump.structure.nb_typ_el):
    #     if el_typ+1 in sim_AC.typ_el_AC:
    #         relevant_eps_effs.append(sim_EM_pump.n_list[el_typ]**2)

    energy_py = np.zeros(comsol_ivals, dtype=np.complex128)
    alpha_py = np.zeros(comsol_ivals)
    Q_PE_py = np.zeros((len(sim_EM_pump.Eig_values),len(sim_EM_Stokes.Eig_values),comsol_ivals), dtype=np.complex128)
    energy_comsol = np.zeros(comsol_ivals, dtype=np.complex128)
    alpha_comsol = np.zeros(comsol_ivals)
    Q_PE_comsol = np.zeros((len(sim_EM_pump.Eig_values),len(sim_EM_Stokes.Eig_values),comsol_ivals), dtype=np.complex128)

    for AC_ival in range(comsol_ivals): # Comsol data only contains some AC modes
        # Interpolate NumBAT FEM fields onto grid
        n_pts_x, n_pts_y, dx, dy, E_mat_p, E_mat_S, u_mat, del_u_mat, del_u_mat_star, m_n = interp_py_fields(
            sim_EM_pump, sim_EM_Stokes, sim_AC, k_AC,
            n_points, EM_ival_pump=EM_ival_pump, EM_ival_Stokes=EM_ival_Stokes, AC_ival=AC_ival)

        # Carry out integration
        energy_py[AC_ival], alpha_py[AC_ival], Q_PE_py[EM_ival_pump,EM_ival_Stokes,AC_ival] = grid_integral(
                m_n, sim_AC.structure, sim_AC.Omega_AC, n_pts_x, n_pts_y, dx, dy, 
                E_mat_p, E_mat_S, u_mat, del_u_mat, del_u_mat_star, AC_ival)

        # Load Comsol FEM fields onto grid - acoustic displacement fields
        x_coord, y_coord, u_mat_comsol = comsol_fields(comsol_data_file, n_points_comsol_data, ival=AC_ival)
        dx_comsol = x_coord[-1,0] - x_coord[-2,0]
        dy_comsol = y_coord[0,-1] - y_coord[0,-2]
        del_u_mat_comsol, del_u_mat_star_comsol = grad_u(dx, dy, u_mat_comsol, k_AC)

        # Carry out integration
        n_pts_x_comsol = n_points_comsol_data
        n_pts_y_comsol = n_points_comsol_data
        energy_comsol[AC_ival], alpha_comsol[AC_ival], Q_PE_comsol[EM_ival_pump,EM_ival_Stokes,AC_ival] = grid_integral(
                m_n, sim_AC.structure, sim_AC.Omega_AC, n_pts_x_comsol, n_pts_y_comsol, 
                dx_comsol, dy_comsol, E_mat_p, E_mat_S, 
                u_mat_comsol, del_u_mat_comsol, del_u_mat_star_comsol, AC_ival)

    # Note this is only the PE contribution to gain.
    gain_PE_py = 2*sim_EM_pump.omega_EM*sim_AC.Omega_AC[:comsol_ivals]*np.real(Q_PE_py*np.conj(Q_PE_py))
    normal_fact_py = np.zeros((num_modes_EM, num_modes_EM, comsol_ivals), dtype=complex)
    gain_PE_comsol = 2*sim_EM_pump.omega_EM*sim_AC.Omega_AC[:comsol_ivals]*np.real(Q_PE_comsol*np.conj(Q_PE_comsol))
    normal_fact_comsol = np.zeros((num_modes_EM, num_modes_EM, comsol_ivals), dtype=complex)
    for i in range(num_modes_EM):
        P1 = sim_EM_pump.EM_mode_power[i]
        for j in range(num_modes_EM):
            P2 = sim_EM_Stokes.EM_mode_power[j]
            for k in range(comsol_ivals):
                P3_py = energy_py[k]
                normal_fact_py[i, j, k] = P1*P2*P3_py*alpha_py[k]
                P3_comsol = energy_comsol[k]
                normal_fact_comsol[i, j, k] = P1*P2*P3_comsol*alpha_comsol[k]
    SBS_gain_PE_py = np.real(gain_PE_py/normal_fact_py)
    SBS_gain_PE_comsol = np.real(gain_PE_comsol/normal_fact_comsol)

    return SBS_gain_PE_py, alpha_py, SBS_gain_PE_comsol, alpha_comsol