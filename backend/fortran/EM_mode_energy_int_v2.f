C   Calculate the Overlap integral of an EM mode with itself.
C
      subroutine EM_mode_energy_int_v2 (lambda, nval, nel, npt,
     *  nnodes, table_nod,
     *  type_el, x, betas, soln_k1, overlap)
c
      implicit none
      integer*8 nval, nel, npt, nnodes      
      integer*8 type_el(nel)
      integer*8 table_nod(nnodes,nel)
      double precision x(2,npt)
      complex*16 soln_k1(3,nnodes+7,nval,nel)
      complex*16 beta1, overlap1
      complex*16 betas(nval)
      complex*16, dimension(nval) :: overlap
      double precision k_0, pi, lambda

c     Local variables

      integer*8 nnodes_0
      parameter (nnodes_0 = 6)
      integer*8 nod_el_p(nnodes_0)
      double precision xel(2,nnodes_0)
      complex*16 E_field_el(3,nnodes_0)
      complex*16 H_field_el(3,nnodes_0)
      double precision p2_p2(nnodes_0,nnodes_0)
      integer*8 i, j, j1, typ_e
      integer*8 iel, ival
      integer*8 itrial, jtest, ui
      complex*16 vec_i(3), vec_j(3)
      complex*16 z_tmp1

      double precision mat_B(2,2), mat_T(2,2), det_b
C
C
Cf2py intent(in) lambda, nval, nel, npt,
Cf2py intent(in) nnodes, table_nod
Cf2py intent(in) type_el, x, betas, soln_k1
C
Cf2py depend(table_nod) nnodes, nel
Cf2py depend(type_el) npt
Cf2py depend(x) npt
Cf2py depend(betas) nval
C
Cf2py intent(out) overlap
C
C
CCCCCCCCCCCCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCC
C
      ui = 6
      pi = 3.141592653589793d0
      k_0 = 2.0d0*pi/lambda
C
      if ( nnodes .ne. 6 ) then
        write(ui,*) "EM_mode_energy_int: problem nnodes = ", nnodes
        write(ui,*) "EM_mode_energy_int: nnodes should be equal to 6 !"
        write(ui,*) "EM_mode_energy_int: Aborting..."
        stop
      endif
C
c

      do ival=1,nval
      overlap1 = 0.0d0
      beta1 = betas(ival)
      do iel=1,nel
        typ_e = type_el(iel)
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
c         The components (E_x,E_y,E_z) the mode ival
          do j=1,3
            z_tmp1 = soln_k1(j,i,ival,iel)
            E_field_el(j,i) = z_tmp1 * beta1
          enddo
        enddo
        call get_H_field (nnodes, k_0, beta1, mat_T, 
     *    E_field_el, H_field_el)
c       The matrix p2_p2 contains the overlap integrals between the P2-polynomial basis functions
        call mat_p2_p2(p2_p2, det_b)
        do itrial=1,nnodes_0
          do i=1,3
            z_tmp1 = E_field_el(i,itrial)
            vec_i(i) = conjg(z_tmp1)
          enddo
          do jtest=1,nnodes_0
            do i=1,3
              z_tmp1 = H_field_el(i,jtest)
              vec_j(i) = z_tmp1
            enddo
c           Cross-product Z.(E^* X H) of E^*=vec_i and H=vec_j
            z_tmp1 = vec_i(1) * vec_j(2) - vec_i(2) * vec_j(1)
            overlap1 = overlap1 + z_tmp1 * p2_p2(itrial, jtest)
          enddo
        enddo
      enddo
      overlap(ival) = overlap1
      enddo
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      end subroutine EM_mode_energy_int_v2
