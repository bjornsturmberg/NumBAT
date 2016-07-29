C Calculate the overlap integral of two EM modes and an AC mode
C using numerical quadrature.
C
      subroutine photoelastic_int (nval_EM, nval_AC, ival1,
     *  ival2, ival3, nel, npt, nnodes, table_nod, type_el, x,
     *  nb_typ_el, p_tensor, beta_AC, soln_EM, soln_AC, eps_lst, 
     *  debug, overlap, basis_overlap)
c
      implicit none
      integer*8 nval_EM, nval_AC, ival1, ival2, ival3
      integer*8 nel, npt, nnodes, nb_typ_el
      integer*8 type_el(nel), debug
      integer*8 table_nod(nnodes,nel)
      double precision x(2,npt)
c      complex*16 x(2,npt)
      complex*16 soln_EM(3,nnodes,nval_EM,nel)
      complex*16 soln_AC(3,nnodes,nval_AC,nel)
      complex*16 overlap, beta_AC
      complex*16 p_tensor(3,3,3,3,nb_typ_el)

c     Local variables
      integer*8 nnodes0
      parameter (nnodes0 = 6)
      double precision xel(2,nnodes0)
      complex*16 basis_overlap(3*nnodes0,3*nnodes0,3,3*nnodes0)
      complex*16 E1star, E2, Ustar, eps
      integer*8 i, j, k, l, j1, typ_e
      integer*8 iel, ind_ip, i_eq
      integer*8 jtest, ind_jp, j_eq, k_eq
      integer*8 ltest, ind_lp, l_eq
      integer*8 itrial, ui
      complex*16 eps_lst(nb_typ_el)
      complex*16 z_tmp1, ii
      double precision mat_B(2,2), mat_T(2,2), eps_0
c
c     NQUAD: The number of quadrature points used in each element.
      integer*8 nquad, nquad_max, iq
      parameter (nquad_max = 16) ! Limit to P2 polynomials
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
Cf2py intent(in) nval_EM, nval_AC, ival1, ival2, ival3, nb_typ_el
Cf2py intent(in) nel, npt, nnodes, table_nod, p_tensor, beta_AC , debug
Cf2py intent(in) type_el, x, soln_EM, soln_AC, eps_lst
C
Cf2py depend(table_nod) nnodes, nel
Cf2py depend(type_el) npt
Cf2py depend(x) npt
Cf2py depend(soln_EM) nnodes, nval_EM, nel
Cf2py depend(soln_AC) nnodes, nval_AC, nel
Cf2py depend(p_tensor) nb_typ_el
C
Cf2py intent(out) overlap, basis_overlap
C
C
CCCCCCCCCCCCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCC
C
      ui = 6
      eps_0 = 8.854187817d-12
      ii = cmplx(0.0d0, 1.0d0)
C
      if ( nnodes .ne. 6 ) then
        write(ui,*) "photoelastic_int: problem nnodes = ", nnodes
        write(ui,*) "photoelastic_int: nnodes should be equal to 6 !"
        write(ui,*) "photoelastic_int: Aborting..."
        stop
      endif
C
      overlap = 0.0d0
      write(*,*) "call quad_triangle begin"
      call quad_triangle (nquad, nquad_max, wq, xq, yq)
      if (debug .eq. 1) then
        write(ui,*) "photoelastic_int: nquad, nquad_max = ",
     *              nquad, nquad_max
      endif
        write(*,*) "call quad_triangle end"
cccccccccccc
C Loop over elements - start
cccccccccccc
      do iel=1,nel
        typ_e = type_el(iel)
        write(*,*) "iel", iel
        do j=1,nnodes
          j1 = table_nod(j,iel)
          xel(1,j) = x(1,j1)
          xel(2,j) = x(2,j1)
        enddo
        write(*,*) "call curved_elem_tri"
        call curved_elem_tri (nnodes, xel, info_curved, r_tmp1)
        if (info_curved .eq. 1) then
          n_curved = n_curved + 1
        endif
        write(*,*) "call curved_elem_tri end"
cccccccccc
        do i=1,3*nnodes
          do j=1,3*nnodes
            do k=1,3
              do l=1,3*nnodes
                basis_overlap(i,j,k,l) = 0.0d0
              enddo
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
          write(*,*) "call phi2_2d_mat"
          call phi2_2d_mat(xx, phi2_list, grad2_mat0)
          write(*,*) "call phi2_2d_mat end"
c
          if (info_curved .eq. 0) then
c           Rectilinear element
            call jacobian_p1_2d(xx, xel, nnodes,
     *               xx_g, det, mat_B, mat_T)
            if(det .le. 0 .and. debug .eq. 2 .and. iq .eq. 1) then
              write(*,*) "   !!!"
              write(*,*) "PE_int: det <= 0: iel, det ", iel, det
            endif
          else
c           Isoparametric element
            call jacobian_p2_2d(xx, xel, nnodes, phi2_list,
     *               grad2_mat0, xx_g, det, mat_B, mat_T)
          endif
C            if(abs(det) .lt. 1.0d-10) then
           if(abs(det) .lt. 1.0d-20) then
             write(*,*)
             write(*,*) "   ???"
             write(*,*) "PE_int: det = 0 : iel, det = ", iel, det
             write(*,*) "PE_int: Aborting..."
             stop
           endif
c          grad_i  = gradient on the actual triangle
c          grad_i  = Transpose(mat_T)*grad_i0
c          Calculation of the matrix-matrix product:
          write(*,*) "call DGEMM"
          call DGEMM('Transpose','N', 2, 6, 2, ONE, mat_T, 2,
     *           grad2_mat0, 2, ZERO, grad2_mat, 2)
          coeff_1 = ww * abs(det)
C Calculate overlap of basis functions at quadrature point, 
C which is a superposition of P2 polynomials for each function (field).
          do itrial=1,nnodes0
            do i_eq=1,3
              ind_ip = i_eq + 3*(itrial-1)
              do jtest=1,nnodes0
                do j_eq=1,3
                  ind_jp = j_eq + 3*(jtest-1)
C                 Gradient of transverse components of basis function
                  do k_eq=1,2
                    do ltest=1,nnodes0
                      do l_eq=1,3
                        ind_lp = l_eq + 3*(ltest-1)
                        z_tmp1 = phi2_list(itrial) * phi2_list(jtest)
     *                          * grad2_mat(k_eq,ltest)
                        coeff_2 = p_tensor(i_eq,j_eq,k_eq,l_eq,typ_e)
                        eps = eps_lst(typ_e)
                        z_tmp1 = coeff_1 * coeff_2 * eps**2 * z_tmp1
                        basis_overlap(ind_ip,ind_jp,k_eq,ind_lp) =
     *                    basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)
     *                          + z_tmp1
                      enddo
                    enddo
                  enddo
C                 Gradient of longitudinal components of basis function,
C                 which is i*beta*phi because field is assumed to be of
C                 form e^{i*beta*z} phi.
                  k_eq=3
                  do ltest=1,nnodes0
                    do l_eq=1,3
                      ind_lp = l_eq + 3*(ltest-1)
                      z_tmp1 = phi2_list(itrial) * phi2_list(jtest)
     *                        * phi2_list(ltest) * ii * beta_AC
                      coeff_2 = p_tensor(i_eq,j_eq,k_eq,l_eq,typ_e)
                      eps = eps_lst(typ_e)
                      z_tmp1 = coeff_1 * coeff_2 * eps**2 * z_tmp1
                      basis_overlap(ind_ip,ind_jp,k_eq,ind_lp) =
     *                  basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)
     *                        + z_tmp1
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
        enddo
cccccccccc
C Having calculated overlap of basis functions on element
C now multiply by specific field values for modes of interest.
        do itrial=1,nnodes0
          do i_eq=1,3
            ind_ip = i_eq + 3*(itrial-1)
            E1star = conjg(soln_EM(i_eq,itrial,ival1,iel))
            do jtest=1,nnodes0
              do j_eq=1,3
                ind_jp = j_eq + 3*(jtest-1)
                E2 = soln_EM(j_eq,jtest,ival2,iel)
                do ltest=1,nnodes0
                  do l_eq=1,3
                    ind_lp = l_eq + 3*(ltest-1)
                    Ustar = conjg(soln_AC(l_eq,ltest,ival3,iel))
                    do k_eq=1,3
                      z_tmp1 = basis_overlap(ind_ip,ind_jp,k_eq,ind_lp)
                      z_tmp1 = E1star * E2 * Ustar * z_tmp1
                      overlap = overlap + z_tmp1
                    enddo
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
C Apply scaling that sits outside of integration.
      overlap = overlap * eps_0
      if(debug .eq. 1) then
        write(*,*) "PE_int: overlap"
        write(*,*) overlap
      endif
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      end subroutine photoelastic_int
