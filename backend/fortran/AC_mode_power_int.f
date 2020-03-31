C Calculate the overlap integral of an AC mode with itself using
C numerical quadrature. 
C
      subroutine AC_mode_power_int (nval, 
     *  nel, npt, nnodes, table_nod, type_el, x,
     *  nb_typ_el, c_tensor_z, beta_AC, Omega_AC, soln_AC,
     *  debug, overlap)
c
      implicit none
      integer*8 nval, ival
      integer*8 nel, npt, nnodes, nb_typ_el
      integer*8 type_el(nel), debug
      integer*8 table_nod(nnodes,nel)
      double precision x(2,npt)
c      complex*16 x(2,npt)
      complex*16 soln_AC(3,nnodes,nval,nel)
      complex*16 Omega_AC(nval)
      complex*16 beta_AC
      complex*16, dimension(nval) :: overlap
      complex*16 c_tensor_z(3,3,3,nb_typ_el)

c     Local variables
      integer*8 nnodes0
      parameter (nnodes0 = 6)
      integer*8 nod_el_p(nnodes0)
      double precision xel(2,nnodes0)
      complex*16 basis_overlap(3*nnodes0,3,3*nnodes0)
      complex*16 U, Ustar
      integer*8 i, j, k, l, j1, typ_e
      integer*8 iel, ind_ip, i_eq, k_eq
      integer*8 ltest, ind_lp, l_eq
      integer*8 itrial, ui
      complex*16 z_tmp1, ii
      double precision mat_B(2,2), mat_T(2,2)
c
c     NQUAD: The number of quadrature points used in each element.
      integer*8 nquad, nquad_max, iq
      parameter (nquad_max = 25)
      double precision wq(nquad_max)
      double precision xq(nquad_max), yq(nquad_max)
      double precision xx(2), xx_g(2), ww, det
      integer*8 info_curved, n_curved
      double precision r_tmp1, ZERO, ONE
      parameter ( ZERO = 0.0D0, ONE = 1.0D0)
      complex*16 coeff_1, coeff_2
      double precision phi2_list(6), grad2_mat0(2,6)
      double precision grad2_mat(2,6)
C
C
Cf2py intent(in) nval, nel, npt, nnodes, table_nod
Cf2py intent(in) type_el, x, nb_typ_el, c_tensor_z, beta_AC 
Cf2py intent(in) soln_AC, debug, Omega_AC
C
Cf2py depend(table_nod) nnodes, nel
Cf2py depend(type_el) npt
Cf2py depend(x) npt
Cf2py depend(soln_AC) nnodes, nval, nel
Cf2py depend(c_tensor_z) nb_typ_el
Cf2py depend(Omega_AC) nval
C
Cf2py intent(out) overlap
C
C
CCCCCCCCCCCCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCC
C
      ui = 6
      ii = cmplx(0.0d0, 1.0d0, 8)
C
      if ( nnodes .ne. 6 ) then
        write(ui,*) "AC_mode_power_int: problem nnodes = ", nnodes
        write(ui,*) "AC_mode_power_int: nnodes should be equal to 6 !"
        write(ui,*) "AC_mode_power_int: Aborting..."
        stop
      endif
C
      call quad_triangle (nquad, nquad_max, wq, xq, yq)
      if (debug .eq. 1) then
        write(ui,*) "AC_mode_power_int: nquad, nquad_max = ",
     *              nquad, nquad_max
      endif
C
      do i=1,nval
        overlap(i) = 0.0d0
      enddo
C
cccccccccccc
C Loop over elements - start
cccccccccccc
      do iel=1,nel
        typ_e = type_el(iel)
        do j=1,nnodes
          j1 = table_nod(j,iel)
          nod_el_p(j) = j1
          xel(1,j) = x(1,j1)
          xel(2,j) = x(2,j1)
        enddo
        call curved_elem_tri (nnodes, xel, info_curved, r_tmp1)
        if (info_curved .eq. 1) then
          n_curved = n_curved + 1
        endif
cccccccccc
        do i=1,3*nnodes
          do k=1,3
            do l=1,3*nnodes
              basis_overlap(i,k,l) = 0.0d0
            enddo
          enddo
        enddo
cccccccccc
C For each quadrature point evaluate overlap of Lagrange polynomials 
C or derivative of Lagrange polynomials 
        do iq=1,nquad
          xx(1) = xq(iq)
          xx(2) = yq(iq)
          ww = wq(iq)
c         xx   = coordinate on the reference triangle
c         xx_g = coordinate on the actual triangle
C         phi2_list = values of Lagrange polynomials (1-6) at each local node.
C         grad2_mat0 = gradient on the reference triangle (P2 element)
          call phi2_2d_mat(xx, phi2_list, grad2_mat0)
c
          if (info_curved .eq. 0) then
c           Rectilinear element
            call jacobian_p1_2d(xx, xel, nnodes,
     *               xx_g, det, mat_B, mat_T)
          else
c           Isoparametric element
            call jacobian_p2_2d(xx, xel, nnodes, phi2_list,
     *               grad2_mat0, xx_g, det, mat_B, mat_T)
          endif
           if(abs(det) .lt. 1.0d-20) then
             write(*,*)
             write(*,*) "   ???"
             write(*,*) "AC_m_en_int: det = 0 : iel, det = ", iel, det
             write(*,*) "AC_m_en_int: Aborting..."
             stop
           endif
c          grad_i  = gradient on the actual triangle
c          grad_i  = Transpose(mat_T)*grad_i0
c          Calculation of the matrix-matrix product:
          call DGEMM('Transpose','N', 2, 6, 2, ONE, mat_T, 2,
     *           grad2_mat0, 2, ZERO, grad2_mat, 2)
          coeff_1 = ww * abs(det)
C Calculate overlap of basis functions at quadrature point, 
C which is a superposition of P2 polynomials for each function (field).
          do itrial=1,nnodes0
            do i_eq=1,3
              ind_ip = i_eq + 3*(itrial-1)
C             Gradient of transverse components of basis function
              do k_eq=1,2
                do ltest=1,nnodes0
                  do l_eq=1,3
                    ind_lp = l_eq + 3*(ltest-1)
                    z_tmp1 = phi2_list(itrial) * grad2_mat(k_eq,ltest)
                    coeff_2 = c_tensor_z(i_eq,k_eq,l_eq,typ_e)
                    basis_overlap(ind_ip,k_eq,ind_lp) =
     *                basis_overlap(ind_ip,k_eq,ind_lp) 
     *                + coeff_1 * coeff_2 * z_tmp1
                  enddo
                enddo
              enddo
C             Gradient of longitudinal components of basis function,
C             which is i*beta*phi because field is assumed to be of
C             form e^{i*beta*z} phi.
              k_eq=3
              do ltest=1,nnodes0
                do l_eq=1,3
                  ind_lp = l_eq + 3*(ltest-1)
                  z_tmp1 = phi2_list(itrial)
     *                    * phi2_list(ltest) * ii * beta_AC
                  coeff_2 = c_tensor_z(i_eq,k_eq,l_eq,typ_e)
                  basis_overlap(ind_ip,k_eq,ind_lp) =
     *              basis_overlap(ind_ip,k_eq,ind_lp)
     *                + coeff_1 * coeff_2 * z_tmp1
                enddo
              enddo
            enddo
          enddo
        enddo
cccccccccc
C Having calculated overlap of basis functions on element
C now multiply by specific field values for modes of interest.
        do ival=1,nval
          do itrial=1,nnodes0
            do i_eq=1,3
              ind_ip = i_eq + 3*(itrial-1)
              Ustar = conjg(soln_AC(i_eq,itrial,ival,iel))
              do ltest=1,nnodes0
                do l_eq=1,3
                  ind_lp = l_eq + 3*(ltest-1)
                  U = soln_AC(l_eq,ltest,ival,iel)
                  do k_eq=1,3
                    z_tmp1 = basis_overlap(ind_ip,k_eq,ind_lp)
                    overlap(ival) = overlap(ival) + Ustar * U * z_tmp1
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
cccccccccccc
C Loop over elements - end
cccccccccccc
      enddo
C Multiply through prefactor
      do i=1,nval
        overlap(i) = -2.0 * ii * Omega_AC(i) * overlap(i)
      enddo
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      end subroutine AC_mode_power_int
