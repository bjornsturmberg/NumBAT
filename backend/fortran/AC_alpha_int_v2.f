C Calculate the overlap integral of an AC mode with itself using
C Direct integration
C
      subroutine AC_alpha_int_v2 (nval, 
     *  nel, npt, nnodes, table_nod, type_el, x,
     *  nb_typ_el, eta_tensor, beta_AC, Omega_AC, soln_AC,
     *  AC_mode_energy_elastic, overlap)
c
      implicit none
      integer*8 nval, ival
      integer*8 nel, npt, nnodes, nb_typ_el
      integer*8 type_el(nel)
      integer*8 table_nod(nnodes,nel)
      double precision x(2,npt)
C       complex*16 x(2,npt)
      complex*16 soln_AC(3,nnodes,nval,nel)
      complex*16 Omega_AC(nval)
      complex*16 beta_AC, AC_mode_energy_elastic(nval)
      complex*16, dimension(nval) :: overlap
      complex*16 eta_tensor(3,3,3,3,nb_typ_el)

c     Local variables
      integer*8 nnodes0
      parameter (nnodes0 = 6)
      integer*8 nod_el_p(nnodes0)
      double precision xel(2,nnodes0)
      complex*16 basis_overlap(3*nnodes0,3,3,3*nnodes0)
      complex*16 U, Ustar
      integer*8 i, j, j1, typ_e
      integer*8 iel, ind_ip, i_eq, k_eq
      integer*8 ltest, ind_lp, l_eq, j_eq
      integer*8 itrial, ui
      complex*16 z_tmp1, ii

      double precision p2_p2(6,6), p2_p2x(6,6), p2_p2y(6,6)
      double precision p2x_p2x(6,6), p2y_p2y(6,6), p2x_p2y(6,6)
      double precision mat_B(2,2), mat_T(2,2), mat_T_tr(2,2)
      double precision det_b
c
      double precision ZERO, ONE
      parameter ( ZERO = 0.0D0, ONE = 1.0D0)
      complex*16 coeff
C
C
Cf2py intent(in) nval, nel, npt, nnodes, table_nod
Cf2py intent(in) type_el, x, nb_typ_el, eta_tensor, beta_AC 
Cf2py intent(in) soln_AC, debug, Omega_AC, AC_mode_energy_elastic
C
Cf2py depend(table_nod) nnodes, nel
Cf2py depend(type_el) npt
Cf2py depend(x) npt
Cf2py depend(soln_AC) nnodes, nval, nel
Cf2py depend(eta_tensor) nb_typ_el
Cf2py depend(Omega_AC) nval
Cf2py depend(AC_mode_energy_elastic) nval
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
        write(ui,*) "AC_alpha_int_v2: problem nnodes = ", nnodes
        write(ui,*) "AC_alpha_int_v2: nnodes should be equal to 6 !"
        write(ui,*) "AC_alpha_int_v2: Aborting..."
        stop
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

      do i=1,2
        do j=1,2
          mat_B(j,i) = xel(j,i+1) - xel(j,1)
        enddo
      enddo
      det_b = mat_B(1,1) * mat_B(2,2) - mat_B(1,2) * mat_B(2,1)
      if (abs(det_b) .le. 1.0d-22) then
        write(*,*) '?? AC_alpha_int_v2: Determinant = 0 :', det_b
        write(*,*) "xel = ", xel
        write(*,*) 'Aborting...'
        stop
      endif
C
c     mat_T = Inverse of mat_B
      mat_T(1,1) = mat_B(2,2) / det_b
      mat_T(2,2) = mat_B(1,1) / det_b
      mat_T(1,2) = -mat_B(1,2) / det_b
      mat_T(2,1) = -mat_B(2,1) / det_b
c
c	mat_T_tr = Tanspose(mat_T)
      mat_T_tr(1,1) = mat_T(1,1)
      mat_T_tr(1,2) = mat_T(2,1)
      mat_T_tr(2,1) = mat_T(1,2)
      mat_T_tr(2,2) = mat_T(2,2)
C
      call mat_p2_p2(p2_p2, det_b)
      call mat_p2_p2x (p2_p2x, mat_T_tr, det_b)
      call mat_p2_p2y (p2_p2y, mat_T_tr, det_b)
      call mat_p2x_p2x (p2x_p2x, mat_T_tr, det_b)
      call mat_p2x_p2y (p2x_p2y, mat_T_tr, det_b)
      call mat_p2y_p2y (p2y_p2y, mat_T_tr, det_b)
cccccccccc
C Calculate overlap of basis functions
C which is a superposition of P2 polynomials for each function (field).
          do itrial=1,nnodes0
            do i_eq=1,3
              ind_ip = i_eq + 3*(itrial-1)
              do j_eq=1,3
                do k_eq=1,3
                  do ltest=1,nnodes0
                    do l_eq=1,3
                      ind_lp = l_eq + 3*(ltest-1)
c                     See Eq. (45) of C. Wolff et al. PRB (2015)
                      if(j_eq == 1 .and. k_eq == 1) then
                        z_tmp1 = p2x_p2x(itrial,ltest)
                      elseif(j_eq == 1 .and. k_eq == 2) then
                        z_tmp1 = p2x_p2y(itrial,ltest)
                      elseif(j_eq == 1 .and. k_eq == 3) then
                        z_tmp1 = p2_p2x(ltest,itrial)
                        z_tmp1 = z_tmp1 * (ii * beta_AC)
cccccccccccccccccccccc
                      elseif(j_eq == 2 .and. k_eq == 1) then
                        z_tmp1 = p2x_p2y(ltest,itrial)
                      elseif(j_eq == 2 .and. k_eq == 2) then
                        z_tmp1 = p2y_p2y(itrial,ltest)
                      elseif(j_eq == 2 .and. k_eq == 3) then
                        z_tmp1 = p2_p2y(ltest,itrial)
                        z_tmp1 = z_tmp1 * (ii * beta_AC)
cccccccccccccccccccccc
                      elseif(j_eq == 3 .and. k_eq == 1) then
                        z_tmp1 = p2_p2x(itrial,ltest)
                        z_tmp1 = z_tmp1 * (-ii * beta_AC)
                      elseif(j_eq == 3 .and. k_eq == 2) then
                        z_tmp1 = p2_p2y(itrial,ltest)
                        z_tmp1 = z_tmp1 * (-ii * beta_AC)
                      elseif(j_eq == 3 .and. k_eq == 3) then
                        z_tmp1 = p2_p2(itrial,ltest)
                        z_tmp1 = z_tmp1 *  beta_AC**2
                      endif
                      coeff = eta_tensor(i_eq,j_eq,k_eq,l_eq,typ_e)
                      basis_overlap(ind_ip,j_eq,k_eq,ind_lp) =
     *                            coeff * z_tmp1
                    enddo
                  enddo
                enddo
              enddo
            enddo
          enddo
cccccccccc
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
                  do j_eq=1,3
                    do k_eq=1,3
                      z_tmp1 = basis_overlap(ind_ip,j_eq,k_eq,ind_lp)
                      z_tmp1 = Ustar * U * z_tmp1
                      overlap(ival) = overlap(ival) + z_tmp1
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
C Multiply through prefactor
      do i=1,nval
C         z_tmp1 = -1.0 * Omega_AC(i)**2 / AC_mode_energy_elastic(i)
C       Flipped sign as assuming did not do integration by parts - going off CW advice.
        z_tmp1 = Omega_AC(i)**2 / AC_mode_energy_elastic(i)
        overlap(i) = z_tmp1 * overlap(i)
      enddo

C       open (unit=26,file="Output/overlap_alpha_v2.txt")
C       do i=1,nval
C         write(26,*) i, Omega_AC(i), abs(overlap(i)), overlap(i), 
C      *              AC_mode_energy_elastic(i)
C       enddo
C       close (unit=26)

C       open (unit=26,file="Output/basis_overlap_v2.txt")
C       do i=1,3*nnodes
C         do j=1,3
C           do k=1,3
C             do l=1,3*nnodes
C             write(26,*)  i,j,k,l, basis_overlap(i,j,k,l)
C             enddo
C           enddo
C         enddo
C       enddo
C       close (unit=26)

C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
      end subroutine AC_alpha_int_v2
