"""
    mode_calcs.py is a subroutine of NumBAT that contains methods to
    calculate the EM and Acoustic modes of a structure.

    Copyright (C) 2016  Bjorn Sturmberg, Kokou Dossou
"""

import time
import numpy as np
from scipy import interpolate
import matplotlib
matplotlib.use('pdf')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
# from mpl_toolkits.axes_grid1 import make_axes_locatable

import plotting
from fortran import NumBAT


def gain_and_qs(sim_EM_wguide, sim_AC_wguide, k_AC,
                EM_ival1=0, EM_ival2=0, AC_ival=0, fixed_Q=None, typ_select_out=None):
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
            Q^{\rm PE} = \varepsilon_0 \int_A {\rm d}^2r \sum_{ijkl} \varepsilon^2_r e_i e_j p_{ijkl} \partial_k u_l^{*},\\

            Q^{\rm MB} =  \int_C {\rm d \mathbf{r} (\mathbf{u}^{*} \cdot \hat n}) \big[ (\varepsilon_a - \varepsilon_b)  
            \varepsilon_0 ({\rm \hat n \times \mathbf{e}}) \cdot ({\rm \hat n \times \mathbf{e}}) - 
            (\varepsilon_a^{-1} - \varepsilon_b^{-1})  \varepsilon_0^{-1} ({\rm \hat n \cdot \mathbf{d}}) 
            \cdot ({\rm \hat n \cdot \mathbf{d}}) \big],\\

            \alpha = \frac{\Omega^2}{P_b} \int {\rm d}^2r \sum_{ijkl} \partial_i u_j^{*} \eta_{ijkl} \partial_k u_l,\\

            \Gamma =  \frac{2 \omega \Omega {\rm Re} (Q_1 Q_1^*)}{P_e P_e P_{ac}} \frac{1}{\alpha} \frac{\alpha^2}{\alpha^2 + \kappa^2}.
  

        Args:
            sim_EM_wguide  (:Simmo: object): Contains all info on EM modes

            sim_AC_wguide  (:Simmo: object): Contains all info on AC modes

            k_AC  (float): Propagation constant of acoustic modes.

        Keyword Args:
            EM_ival1  (int/string): Specify mode number of EM mode 1 (pump mode)
                to calculate interactions for.
                Numbering is python index so runs from 0 to num_EM_modes-1,
                with 0 being fundamental mode (largest prop constant).
                Can also set to 'All' to include all modes.

            EM_ival2  (int/string): Specify mode number of EM mode 2 (stokes mode)
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
            SBS_gain  (num_modes_EM,num_modes_EM,num_modes_AC): The SBS gain including both photoelastic and moving boundary contributions.
            SBS_gain_PE  (num_modes_EM,num_modes_EM,num_modes_AC): The SBS gain for only the photoelastic effect.
            SBS_gain_MB  (num_modes_EM,num_modes_EM,num_modes_AC): The SBS gain for only the moving boundary effect.
            alpha  (num_modes_AC): The acoustic loss for each mode.
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


    if EM_ival1 == 'All':
        EM_ival1_fortran = -1
    else:
        EM_ival1_fortran = EM_ival1+1  # convert back to Fortran indexing
    if EM_ival2 == 'All':
        EM_ival2_fortran = -1
    else:
        EM_ival2_fortran = EM_ival2+1  # convert back to Fortran indexing
    if AC_ival == 'All':
        AC_ival_fortran = -1
    else:
        AC_ival_fortran = AC_ival+1  # convert back to Fortran indexing

    Fortran_debug = 0
    ncomps = 3
    nnodes = 6
    num_modes_EM = sim_EM_wguide.num_modes
    num_modes_AC = sim_AC_wguide.num_modes
    n_msh_el_AC = sim_AC_wguide.n_msh_el
    trimmed_EM_field = np.zeros((ncomps,nnodes,num_modes_EM,n_msh_el_AC), dtype=complex)
    for el in range(n_msh_el_AC):
        new_el = sim_AC_wguide.el_convert_tbl[el]
        for ival in range(num_modes_EM):
            for n in range(nnodes):
                for x in range(ncomps):
                    trimmed_EM_field[x,n,ival,el] = sim_EM_wguide.sol1[x,n,ival,new_el]
    # sim_EM_wguide.sol1 = trimmed_EM_field
    # sim_EM_wguide.n_msh_el = sim_AC_wguide.n_msh_el
    # sim_EM_wguide.n_msh_pts = sim_AC_wguide.n_msh_pts
    # sim_EM_wguide.type_el = sim_AC_wguide.type_el
    # sim_EM_wguide.table_nod = sim_AC_wguide.table_nod
    # sim_EM_wguide.x_arr = sim_AC_wguide.x_arr
    # plotting.plt_mode_fields(sim_EM_wguide, EM_AC='EM', add_name='trim')

    relevant_eps_effs =[]
    for el_typ in range(sim_EM_wguide.structure.nb_typ_el):
        if el_typ+1 in sim_AC_wguide.typ_el_AC:
            relevant_eps_effs.append(sim_EM_wguide.n_list[el_typ]**2)

    print("\n -----------------------------------------------")
    if fixed_Q is None:
        # Calc alpha (loss) Eq. 45
        print("Acoustic loss calc")
        start = time.time()
        try:
            if sim_EM_wguide.structure.inc_shape in sim_EM_wguide.structure.linear_element_shapes:
                alpha = NumBAT.ac_alpha_int_v2(sim_AC_wguide.num_modes,
                    sim_AC_wguide.n_msh_el, sim_AC_wguide.n_msh_pts, nnodes,
                    sim_AC_wguide.table_nod, sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
                    sim_AC_wguide.structure.nb_typ_el_AC, sim_AC_wguide.structure.eta_tensor,
                    k_AC, sim_AC_wguide.Omega_AC, sim_AC_wguide.sol1,
                    sim_AC_wguide.AC_mode_overlap)
            else:
                if sim_EM_wguide.structure.inc_shape != 'circular':
                    print("Warning: ac_alpha_int - not sure if mesh contains curvi-linear elements", 
                        "\n using slow quadrature integration by default.\n\n")
                alpha = NumBAT.ac_alpha_int(sim_AC_wguide.num_modes,
                    sim_AC_wguide.n_msh_el, sim_AC_wguide.n_msh_pts, nnodes,
                    sim_AC_wguide.table_nod, sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
                    sim_AC_wguide.structure.nb_typ_el_AC, sim_AC_wguide.structure.eta_tensor,
                    k_AC, sim_AC_wguide.Omega_AC, sim_AC_wguide.sol1,
                    sim_AC_wguide.AC_mode_overlap, Fortran_debug)
        except KeyboardInterrupt:
            print("\n\n Routine ac_alpha_int interrupted by keyboard.\n\n")
        alpha = np.real(alpha)
        end = time.time()
        print("     time (sec.)", (end - start))
    else:
        alpha = (np.pi*k_AC/fixed_Q)*np.ones(num_modes_AC)


    # Calc Q_photoelastic Eq. 33
    print("Photoelastic calc")
    start = time.time()
    try:
        if sim_EM_wguide.structure.inc_shape in sim_EM_wguide.structure.linear_element_shapes:
            Q_PE = NumBAT.photoelastic_int_v2(
                sim_EM_wguide.num_modes, sim_AC_wguide.num_modes, EM_ival1_fortran,
                EM_ival2_fortran, AC_ival_fortran, sim_AC_wguide.n_msh_el,
                sim_AC_wguide.n_msh_pts, nnodes,
                sim_AC_wguide.table_nod, sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
                sim_AC_wguide.structure.nb_typ_el_AC, sim_AC_wguide.structure.p_tensor,
                k_AC, trimmed_EM_field, sim_AC_wguide.sol1,
                relevant_eps_effs, sim_EM_wguide.Eig_values, Fortran_debug)
        else:
            if sim_EM_wguide.structure.inc_shape != 'circular':
                print("Warning: photoelastic_int - not sure if mesh contains curvi-linear elements", 
                    "\n using slow quadrature integration by default.\n\n")
            Q_PE = NumBAT.photoelastic_int(
                sim_EM_wguide.num_modes, sim_AC_wguide.num_modes, EM_ival1_fortran,
                EM_ival2_fortran, AC_ival_fortran, sim_AC_wguide.n_msh_el,
                sim_AC_wguide.n_msh_pts, nnodes,
                sim_AC_wguide.table_nod, sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
                sim_AC_wguide.structure.nb_typ_el_AC, sim_AC_wguide.structure.p_tensor,
                k_AC, trimmed_EM_field, sim_AC_wguide.sol1,
                relevant_eps_effs, sim_EM_wguide.Eig_values, Fortran_debug)
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
        Q_MB = NumBAT.moving_boundary(sim_EM_wguide.num_modes,
            sim_AC_wguide.num_modes, EM_ival1_fortran, EM_ival2_fortran,
            AC_ival_fortran, sim_AC_wguide.n_msh_el,
            sim_AC_wguide.n_msh_pts, nnodes, sim_AC_wguide.table_nod, 
            sim_AC_wguide.type_el, sim_AC_wguide.x_arr,
            sim_AC_wguide.structure.nb_typ_el_AC, typ_select_in, typ_select_out,
            sim_EM_wguide.Eig_values, trimmed_EM_field, sim_AC_wguide.sol1,
            relevant_eps_effs, Fortran_debug)
    except KeyboardInterrupt:
        print("\n\n Routine moving_boundary interrupted by keyboard.\n\n")
    end = time.time()
    print("     time (sec.)", (end - start))
    print("-----------------------------------------------")

    Q = Q_PE + Q_MB

    # Note: sim_EM_wguide.omega_EM is the optical angular freq in units of Hz
    # Note: sim_AC_wguide.Omega_AC is the acoustic angular freq in units of Hz
    gain = 2*sim_EM_wguide.omega_EM*sim_AC_wguide.Omega_AC*np.real(Q*np.conj(Q))
    gain_PE = 2*sim_EM_wguide.omega_EM*sim_AC_wguide.Omega_AC*np.real(Q_PE*np.conj(Q_PE))
    gain_MB = 2*sim_EM_wguide.omega_EM*sim_AC_wguide.Omega_AC*np.real(Q_MB*np.conj(Q_MB))
    normal_fact = np.zeros((num_modes_EM, num_modes_EM, num_modes_AC), dtype=complex)
    for i in range(num_modes_EM):
        P1 = sim_EM_wguide.EM_mode_overlap[i]
        for j in range(num_modes_EM):
            P2 = sim_EM_wguide.EM_mode_overlap[j]
            for k in range(num_modes_AC):
                P3 = sim_AC_wguide.AC_mode_overlap[k]
                normal_fact[i, j, k] = P1*P2*P3
    SBS_gain = np.real(gain/normal_fact)
    SBS_gain_PE = np.real(gain_PE/normal_fact)
    SBS_gain_MB = np.real(gain_MB/normal_fact)

    return SBS_gain, SBS_gain_PE, SBS_gain_MB, alpha


#### Categorise modes by their symmetries #############################################
def symmetries(sim_wguide, n_points=10):
    """ Plot EM mode fields.

        Args:
            sim_wguide : A :Struct: instance that has had calc_modes calculated

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
        # v_Ez6p = np.zeros(6*sim_wguide.n_msh_el, dtype=np.complex128)
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
                # v_Ez6p[i] = mode_fields[2,i_node,ival,i_el]
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
        # ReEz = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ez6p.real,trifinder=finder)
        # ImEz = matplotlib.tri.LinearTriInterpolator(triang6p,v_Ez6p.imag,trifinder=finder)

        ### plotting
        # interpolated fields
        m_ReEx = ReEx(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ReEy = ReEy(v_x,v_y).reshape(n_pts_x,n_pts_y)
        # m_ReEz = ReEz(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ImEx = ImEx(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_ImEy = ImEy(v_x,v_y).reshape(n_pts_x,n_pts_y)
        # m_ImEz = ImEz(v_x,v_y).reshape(n_pts_x,n_pts_y)
        m_Ex = m_ReEx + 1j*m_ImEx
        m_Ey = m_ReEy + 1j*m_ImEy
        # m_Ez = m_ReEz + 1j*m_ImEz

        m_Ex_ymirror = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        m_Ex_xmirror = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        m_Ex_rotated = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        m_Ey_ymirror = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        m_Ey_xmirror = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        m_Ey_rotated = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        # m_Ez_ymirror = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        # m_Ez_xmirror = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        # m_Ez_rotated = np.zeros((n_pts_x,n_pts_y), dtype=np.complex128)
        Ex_sigma_y = 0
        Ey_sigma_y = 0
        # Ez_sigma_y = 0
        Ex_sigma_x = 0
        Ey_sigma_x = 0
        # Ez_sigma_x = 0
        Ex_C_2 = 0
        Ey_C_2 = 0
        # Ez_C_2 = 0
        # max_E = max(np.max(np.abs(m_Ex)), np.max(np.abs(m_Ey)), np.max(np.abs(m_Ez)))

        for ix in range(n_pts_x):
            for iy in range(n_pts_y):
                m_Ex_ymirror[ix,iy] = (m_Ex[ix,n_pts_y-iy-1])
                m_Ey_ymirror[ix,iy] = -1*(m_Ey[ix,n_pts_y-iy-1])
                # m_Ez_ymirror[ix,iy] = m_Ez[ix,n_pts_y-iy-1]
                m_Ex_xmirror[ix,iy] = -1*(m_Ex[n_pts_x-ix-1,iy])
                m_Ey_xmirror[ix,iy] = (m_Ey[n_pts_x-ix-1,iy])
                # m_Ez_xmirror[ix,iy] = m_Ez[n_pts_x-ix-1,iy]
                m_Ex_rotated[ix,iy] = -1*(m_Ex[n_pts_x-ix-1,n_pts_y-iy-1])
                m_Ey_rotated[ix,iy] = -1*(m_Ey[n_pts_x-ix-1,n_pts_y-iy-1])
                # m_Ez_rotated[ix,iy] = m_Ez[n_pts_x-ix-1,iy]
                # Ex_sigma_y += np.abs(m_Ex[ix,iy] - m_Ex_ymirror[ix,iy])/max_E
                # Ey_sigma_y += np.abs(m_Ey[ix,iy] - m_Ey_ymirror[ix,iy])/max_E
                # Ez_sigma_y += np.abs(m_Ez[ix,iy] - m_Ez_ymirror[ix,iy])/max_E
                # Ex_sigma_x += np.abs(m_Ex[ix,iy] - m_Ex_xmirror[ix,iy])/max_E
                # Ey_sigma_x += np.abs(m_Ey[ix,iy] - m_Ey_xmirror[ix,iy])/max_E
                # Ez_sigma_x += np.abs(m_Ez[ix,iy] - m_Ez_xmirror[ix,iy])/max_E
                # Ex_C_2 += np.abs(m_Ex[ix,iy] - m_Ex_rotated[ix,iy])/max_E
                # Ey_C_2 += np.abs(m_Ey[ix,iy] - m_Ey_rotated[ix,iy])/max_E
                # Ez_C_2 += np.abs(m_Ez[ix,iy] - m_Ez_rotated[ix,iy])/max_E

        Ex_sigma_y = np.sum(np.abs(m_Ex - m_Ex_ymirror))
        Ey_sigma_y = np.sum(np.abs(m_Ey - m_Ey_ymirror))
        # Ez_sigma_y = np.sum(np.abs(m_Ez - m_Ez_ymirror))
        Ex_sigma_x = np.sum(np.abs(m_Ex - m_Ex_xmirror))
        Ey_sigma_x = np.sum(np.abs(m_Ey - m_Ey_xmirror))
        # Ez_sigma_x = np.sum(np.abs(m_Ez - m_Ez_xmirror))
        Ex_C_2 = np.sum(np.abs(m_Ex - m_Ex_rotated))
        Ey_C_2 = np.sum(np.abs(m_Ey - m_Ey_rotated))
        # Ez_C_2 = np.sum(np.abs(m_Ez - m_Ez_rotated))
        # sigma_y = (Ex_sigma_y + Ey_sigma_y + Ez_sigma_y)/(n_pts_x*n_pts_y)
        # sigma_x = (Ex_sigma_x + Ey_sigma_x + Ez_sigma_x)/(n_pts_x*n_pts_y)
        # C_2 = (Ex_C_2 + Ey_C_2 + Ez_C_2)/(n_pts_x*n_pts_y)
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
        # print '------'
        # print ival
        # print 'C_2', C_2_print
        # print 'C_2', C_2
        # print 'sigma_x', sigma_x_print
        # print 'sigma_x', sigma_x
        # print 'sigma_y', sigma_y_print
        # print 'sigma_y', sigma_y
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
