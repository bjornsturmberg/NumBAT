c
C Calculate the H-field soln_H1 from the E-field soln_k1 of a mode
C  The z-component of the E-field is not normalised
C
      subroutine H_mode_field_Ez (k_0, nval, nel, npt,
     *  nnodes_P2, table_nod,
     *  x, betas, soln_k1, soln_H1)
c
      implicit none
      integer*8 nval, nel, npt, nnodes_P2
      integer*8 table_nod(nnodes_P2,nel)
      double precision x(2,npt)
c      complex*16 x(2,npt)
      complex*16 soln_k1(3,nnodes_P2+7,nval,nel)
      complex*16 soln_H1(3,nnodes_P2,nval,nel)
      complex*16 beta1
      complex*16 betas(nval)
      double precision k_0

c     Local variables

      integer*8 nnodes_P2_0, nnodes_P3_0
      parameter (nnodes_P2_0 = 6, nnodes_P3_0 = 10)
      integer*8 nod_el_p(nnodes_P2_0)
      double precision xel(2,nnodes_P2_0)
      complex*16 E_field_el(3,nnodes_P2_0)
      complex*16 H_field_el(3,nnodes_P2_0)
C       !  P3 Ez-field
      complex*16 Ez_field_el_P3(nnodes_P3_0)  
      integer*8 i, j, j1
      integer*8 iel, ival, inod
      integer*8 ui
      complex*16 z_tmp1, ii

      double precision mat_B(2,2), mat_T(2,2), det_b
C
C
Cf2py intent(in) k_0, nval, nel, npt,
Cf2py intent(in) nnodes_P2, table_nod
Cf2py intent(in) x, betas, soln_k1
C
Cf2py depend(table_nod) nnodes_P2, nel
Cf2py depend(x) npt
Cf2py depend(betas) nval
Cf2py depend(soln_k1) nnodes_P2, nval, nel
C
Cf2py intent(out) soln_H1
C
C
CCCCCCCCCCCCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCC
C
      ui = 6
      ii = cmplx(0.0d0, 1.0d0, 8)
C
      if ( nnodes_P2 .ne. 6 ) then
        write(ui,*) "EM_mode_en_int_v2: problem nnodes = ", nnodes_P2
        write(ui,*) "EM_mode_en_int_v2: nnodes should be equal to 6 !"
        write(ui,*) "EM_mode_en_int_v2: Aborting..."
        stop
      endif
C
      do ival=1,nval
        beta1 = betas(ival)
        do iel=1,nel
          do j=1,nnodes_P2
            j1 = table_nod(j,iel)
            nod_el_p(j) = j1
            xel(1,j) = x(1,j1)
            xel(2,j) = x(2,j1)
          enddo
cccccccc
c         The geometric transformation (x,y) -> (x_g,y_g) = mat_B*(x,y)^t + (x_0, y_0, z_0)^t
c         maps the current triangle to the reference triangle.
          do i=1,2
            do j=1,2
              mat_B(j,i) = xel(j,i+1) - xel(j,1)
            enddo
          enddo
          det_b = mat_B(1,1) * mat_B(2,2) - mat_B(1,2) * mat_B(2,1)
          if (abs(det_b) .le. 1.0d-22) then 
cc          if (abs(det_b) .le. 1.0d-8) then
            write(*,*) '?? H_mode_field_Ez: Deter. = 0 :', det_b
            write(*,*) "xel = ", xel
            write(*,*) 'Aborting...'
            stop
          endif
c         We also need, is the matrix mat_T of the reverse transmation
c                  (from reference to current triangle):
c         mat_T = inverse matrix of de mat_B
          mat_T(1,1) =  mat_B(2,2) / det_b
          mat_T(2,2) =  mat_B(1,1) / det_b
          mat_T(1,2) = -mat_B(1,2) / det_b
          mat_T(2,1) = -mat_B(2,1) / det_b
c         Note that if grad_i_0 is the gradient on the reference triangle,
c         then the gradient on the actual triangle is:
c         grad_i  = Transpose(mat_T)*grad_i0
cccccccc
          do inod=1,nnodes_P2
c           The components (E_x,E_y) of the mode ival
            do j=1,2
              z_tmp1 = soln_k1(j,inod,ival,iel)
              E_field_el(j,inod) = z_tmp1
            enddo
c           The component E_z of the mode ival. The FEM code uses the scaling:
c           E_z = ii * beta1 * \hat{E}_z
            j=3
              z_tmp1 = soln_k1(j,inod,ival,iel)
              E_field_el(j,inod) = z_tmp1
          enddo
c         E_z-field:
          do inod=1,3
c           The longitudinal component at the vertices (P3 elements)
            j=3
            z_tmp1 = soln_k1(j,inod,ival,iel)
            Ez_field_el_P3(inod) = z_tmp1
          enddo
          do inod=4,nnodes_P3_0
c           The longitudinal component at the edge nodes and interior node (P3 elements)
            j=3
            z_tmp1 = soln_k1(j,inod+nnodes_P2-3,ival,iel)
            Ez_field_el_P3(inod) = z_tmp1
          enddo
c
          call get_H_field_p3 (nnodes_P2, k_0, beta1, mat_T,
     *      E_field_el, Ez_field_el_P3, H_field_el)
c
          do inod=1,nnodes_P2_0
            do j=1,3
              z_tmp1 = H_field_el(j,inod)
              soln_H1(j,inod,ival,iel) = z_tmp1
            enddo
          enddo
        enddo
      enddo
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      end subroutine H_mode_field_Ez
