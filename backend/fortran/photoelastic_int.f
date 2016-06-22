C   Calculate the Overlap integral of two EM modes and an AC mode.
C
      subroutine photoelastic_int (lambda, nval_EM, nval_AC, ival1,
     *  ival2, ival3, nel, npt, nnodes, table_nod, type_el, x,
     *  nb_typ_el, p_tensor, betas, soln_EM, soln_AC, eps_lst, overlap)
c
      implicit none
      integer*8 nval_EM, nval_AC, ival1, ival2, ival3
      integer*8 nel, npt, nnodes, nb_typ_el
      integer*8 type_el(nel)
      integer*8 table_nod(nnodes,nel)
      double precision x(2,npt)
      complex*16 soln_EM(3,nnodes,nval_EM,nel)
      complex*16 soln_AC(3,nnodes,nval_AC,nel)
      complex*16 beta1, beta2, overlap
      complex*16 betas(nval_EM)
      complex*16 p_tensor(6,6,nb_typ_el)
      double precision k_0, pi, lambda

c     Local variables
      integer*8 nnodes_0
      parameter (nnodes_0 = 6)
      integer*8 nod_el_p(nnodes_0)
      double precision xel(2,nnodes_0)
      complex*16 E_field_el1(3,nnodes_0)
      complex*16 E_field_el2(3,nnodes_0)
      complex*16 U_field_el(3,nnodes_0)
      double precision p2_p2(nnodes_0,nnodes_0)
      integer*8 i, j, j1, typ_e
      integer*8 iel
      integer*8 itrial, jtest, ui
      complex*16 vec_i(3), vec_j(3), epsilon, eps_lst(nb_typ_el)
      complex*16 z_tmp1, z_tmp2, z_tmp3, ii
      double precision mat_B(2,2), mat_T(2,2), det_b
C
C
Cf2py intent(in) lambda, nval_EM, nval_AC, ival1, ival2, ival3
Cf2py intent(in) nel, npt, nnodes, table_nod, p_tensor
Cf2py intent(in) type_el, x, betas, soln_EM, soln_AC
C
Cf2py depend(table_nod) nnodes, nel
Cf2py depend(type_el) npt
Cf2py depend(x) npt
Cf2py depend(betas) nval_EM
Cf2py depend(soln_EM) nnodes, nval_EM, nel
Cf2py depend(soln_AC) nnodes, nval_AC, nel
Cf2py depend(p_tensor) nb_typ_el
C
Cf2py intent(out) overlap
C
C
CCCCCCCCCCCCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCC
C
      ui = 6
      pi = 3.141592653589793d0
      k_0 = 2.0d0*pi/lambda
      ii = cmplx(0.0d0, 1.0d0)
C
      if ( nnodes .ne. 6 ) then
        write(ui,*) "EM_mode_energy_int: problem nnodes = ", nnodes
        write(ui,*) "EM_mode_energy_int: nnodes should be equal to 6 !"
        write(ui,*) "EM_mode_energy_int: Aborting..."
        stop
      endif
C
      overlap = 0.0d0
      beta1 = betas(ival1)
      beta2 = betas(ival2)
      do iel=1,nel
        typ_e = type_el(iel)
        epsilon = eps_lst(typ_e)
        do j=1,nnodes
          j1 = table_nod(j,iel)
          nod_el_p(j) = j1
          xel(1,j) = x(1,j1)
          xel(2,j) = x(2,j1)
        enddo
cccccccc
c       The geometric transformation (x,y) -> (x_g,y_g) = mat_B*(x,y)^t + (x_0, y_0, z_0)^t
c       maps the current triangle to the reference triangle.
        do i=1,2
          do j=1,2
            mat_B(j,i) = xel(j,i+1) - xel(j,1)
          enddo
        enddo
        det_b = mat_B(1,1) * mat_B(2,2) - mat_B(1,2) * mat_B(2,1)
        if (abs(det_b) .le. 1.0d-22) then  ! TEMPORARY CHANGE
cc        if (abs(det_b) .le. 1.0d-8) then
          write(*,*) '?? get_H_field: Determinant = 0 :', det_b
          write(*,*) "xel = ", xel
          write(*,*) 'Aborting...'
          stop
        endif
c       We also need, is the matrix mat_T of the reverse transmation
c                (from reference to current triangle):
c       mat_T = inverse matrix of de mat_B
        mat_T(1,1) =  mat_B(2,2) / det_b
        mat_T(2,2) =  mat_B(1,1) / det_b
        mat_T(1,2) = -mat_B(1,2) / det_b
        mat_T(2,1) = -mat_B(2,1) / det_b
c       Note that if grad_i_0 is the gradient on the reference triangle,
c       then the gradient on the actual triangle is:
c       grad_i  = Transpose(mat_T)*grad_i0
cccccccc
        do i=1,nnodes
c         The components (E_x,E_y) of the mode ival
          do j=1,2
            z_tmp1 = soln_EM(j,i,ival1,iel)
            E_field_el1(j,i) = z_tmp1
            z_tmp2 = soln_EM(j,i,ival2,iel)
            E_field_el2(j,i) = z_tmp2
            z_tmp3 = soln_AC(j,i,ival3,iel)
            U_field_el(j,i) = z_tmp3
          enddo
c         The component E_z of the mode ival. The FEM code uses the scaling:
c         E_z = ii * beta1 * \hat{E}_z
          j=3
            z_tmp1 = soln_EM(j,i,ival1,iel)
            E_field_el1(j,i) = z_tmp1 * ii * beta1
            z_tmp2 = soln_EM(j,i,ival2,iel)
            E_field_el2(j,i) = z_tmp2 * ii * beta2
            z_tmp3 = soln_AC(j,i,ival3,iel)
            U_field_el(j,i) = z_tmp3
        enddo
c       The matrix p2_p2 contains the overlap integrals between the P2-polynomial basis functions
        call mat_p2_p2(p2_p2, det_b)
        do itrial=1,nnodes_0
          do i=1,3
            z_tmp1 = E_field_el1(i,itrial)
            vec_i(i) = conjg(z_tmp1)
          enddo
          do jtest=1,nnodes_0
            do i=1,3
              z_tmp1 = E_field_el2(i,jtest)
              vec_j(i) = z_tmp1
            enddo
c           Cross-product Z.(E^* X H) of E^*=vec_i and H=vec_j
            z_tmp1 = vec_i(1) * vec_j(2) - vec_i(2) * vec_j(1)
            overlap = overlap + z_tmp1 * p2_p2(itrial, jtest)
          enddo
        enddo
      enddo
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      end subroutine photoelastic_int
