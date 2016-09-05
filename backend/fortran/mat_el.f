c
c***********************************************************************
c
c
c***********************************************************************
c
      subroutine mat_el (xel, beta, c_tensor_el, rho_el,
     *  mat_K, mat_M)
c     *  mat_M, mat_K)
c
c***********************************************************************
c
      implicit none
      double precision xel(2,6)
      complex*16 beta
      complex*16 mat_K(18,18), mat_M(18,18)
      complex*16 c_tensor_el(6,6), rho_el

c     Local variables

      double precision p2_p2(6,6), p2_p2x(6,6), p2_p2y(6,6)
      double precision p2x_p2x(6,6), p2y_p2y(6,6), p2x_p2y(6,6)
      double precision mat_B(2,2), mat_T(2,2), mat_T_tr(2,2)
      double precision det_b
      complex*16 ii, z_tmp1
      integer*8 i, j, i_p, j_p, i_xyz,  j_xyz
      integer*8 debug

c    Compute the Affine mappings from the current triangle to the
c     reference unit triangle.
c    Integration will be performed on the reference unit triangle
c
c
ccccccccccccccccccccccccccccccccccccccc
c
      debug = 0

c  ii = sqrt(-1)
      ii = dcmplx(0.0d0, 1.0d0)

      do i=1,2
        do j=1,2
          mat_B(j,i) = xel(j,i+1) - xel(j,1)
        enddo
      enddo
      det_b = mat_B(1,1) * mat_B(2,2) - mat_B(1,2) * mat_B(2,1)
      if (abs(det_b) .le. 1.0d-22) then  ! TEMPORARY CHANGE
cc      if (abs(det_b) .le. 1.0d-8) then
        write(*,*) '?? mat_el: Determinant = 0 :', det_b
        write(*,*) "xel = ", xel
        write(*,*) 'Aborting...'
        stop
      endif
c     mat_T = Inverse of mat_B
      mat_T(1,1) = mat_B(2,2) / det_b
      mat_T(2,2) = mat_B(1,1) / det_b
      mat_T(1,2) = -mat_B(1,2) / det_b
      mat_T(2,1) = -mat_B(2,1) / det_b
c
c
c	mat_T_tr = Tanspose(mat_T)
c
      mat_T_tr(1,1) = mat_T(1,1)
      mat_T_tr(1,2) = mat_T(2,1)
      mat_T_tr(2,1) = mat_T(1,2)
      mat_T_tr(2,2) = mat_T(2,2)

      call mat_p2_p2(p2_p2, det_b)
      call mat_p2_p2x (p2_p2x, mat_T_tr, det_b)
      call mat_p2_p2y (p2_p2y, mat_T_tr, det_b)
      call mat_p2x_p2x (p2x_p2x, mat_T_tr, det_b)
      call mat_p2x_p2y (p2x_p2y, mat_T_tr, det_b)
      call mat_p2y_p2y (p2y_p2y, mat_T_tr, det_b)
c
      do i=1,18
        do j=1,18
          mat_K(j,i) = 0
          mat_M(j,i) = 0
        enddo
      enddo

c=================  Construction of the matrix mat_M =================
c     Integral [rho * P(i) * P(i)]
      do i=1,6
        do i_xyz=1,3  ! The components x, y and z
          i_p = 3*(i-1) + i_xyz
          do j=1,6
            j_xyz = i_xyz
            j_p = 3*(j-1) + j_xyz
            z_tmp1 = p2_p2(i,j) * rho_el
            mat_M(i_p,j_p) = mat_M(i_p,j_p) + z_tmp1
          enddo
        enddo
      enddo
c
c=================  Construction of the matrix mat_K =================
c     Integral [K_{ij} = Gradient_s(conjg(P_k(i))) x c_tensor x Gradient_s(P_k(j))], where k=x,y,z
c     Reference: see Eqs. (7) and (8) in:
c     A.-C. Hladky-Hennion
c     "Finite element analysis of the propagation of acoustic waves in waveguides," 
c     Journal of Sound and Vibration, vol. 194, no. 2, pp. 119-136, 1996. 
      do i=1,6
        do i_xyz=1,3
          i_p = 3*(i-1) + i_xyz
          do j=1,6
            j_xyz = i_xyz
            j_p = 3*(j-1) + j_xyz
            if (i_xyz == 1) then
c             Overlap: row 1 of [C]*[B] ###########
              z_tmp1 = p2x_p2x(i,j)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,i_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2x(j,i) * beta * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2y(i,j)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 5 of [C]*[B]  ###########
              z_tmp1 = p2_p2x(i,j) * (-beta) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,i_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2(i,j) * beta**2
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2y(i,j) * (-beta) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 6 of [C]*[B]  ###########
              z_tmp1 = p2x_p2y(j,i)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+5,i_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2y(j,i) * beta * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+5,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2y_p2y(i,j)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+5,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            elseif (i_xyz == 2) then
c             Overlap: row 2 of [C]*[B]  ###########
              z_tmp1 = p2y_p2y(i,j)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,i_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2y(j,i) * beta * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2y(j,i)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 4 of [C]*[B]  ###########
              z_tmp1 = p2_p2y(i,j) * (-beta) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,i_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2(i,j) * beta**2
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2x(i,j) * (-beta) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 6 of [C]*[B]  ###########
              z_tmp1 = p2x_p2y(i,j)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,i_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2x(j,i) * (-beta) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2x(i,j)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            elseif (i_xyz == 3) then
c             Overlap: column 3 of [C]*[B]  ###########
              z_tmp1 = p2_p2(i,j) * beta**2
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,i_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2y(i,j) * (-beta) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2x(i,j) * (-beta) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 4 of [C]*[B]  ###########
              z_tmp1 = p2_p2y(j,i) * beta * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+1,i_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2y_p2y(i,j)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+1,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2y(j,i)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+1,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 5 of [C]*[B]  ###########
              z_tmp1 = p2_p2x(j,i) * beta * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,i_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2y(i,j)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2x(i,j)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            endif
          enddo
        enddo
      enddo
      do i=1,6
        do i_xyz=1,3
          i_p = 3*(i-1) + i_xyz
          do j=1,6
            if (i_xyz == 1) then
c             Overlap: row 1 of [C]*[B]  ###########
              j_xyz = 2
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2x_p2y(i,j)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2x(j,i) * beta * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2x(i,j)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 5 of [C]*[B]  ###########
              j_xyz = 2
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2_p2y(i,j) * (-beta) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2(i,j) * beta**2
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2x(i,j) * (-beta) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 6 of [C]*[B]  ###########
              j_xyz = 2
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2y_p2y(i,j)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+5,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2y(j,i) * beta * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+5,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2y(j,i)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+5,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            elseif (i_xyz == 2) then
c             Overlap: row 2 of [C]*[B]  ###########
              j_xyz = 3
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2_p2y(j,i) * (-beta)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2y_p2y(i,j) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2y(j,i) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 4 of [C]*[B]  ###########
              j_xyz = 3
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2_p2(i,j) * beta**2 * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2y(i,j) * beta
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2x(i,j) * beta
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 6 of [C]*[B]  ###########
              j_xyz = 3
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2_p2x(j,i) * (-beta)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2y(i,j) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2x(i,j) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            elseif (i_xyz == 3) then
c             Overlap: row 3 of [C]*[B]  ###########
              j_xyz = 1
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2_p2x(i,j) * (-beta)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2(i,j) * (-beta**2) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2y(i,j) * (-beta)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 4 of [C]*[B]  ###########
              j_xyz = 1
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2x_p2y(j,i) * (-ii)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+1,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2y(j,i) * beta
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+1,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2y_p2y(i,j) * (-ii)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+1,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 5 of [C]*[B]  ###########
              j_xyz = 1
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2x_p2x(i,j) * (-ii)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2x(j,i) * beta
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2y(i,j) * (-ii)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            endif
          enddo
        enddo
      enddo
      do i=1,6
        do i_xyz=1,3
          i_p = 3*(i-1) + i_xyz
          do j=1,6
            if (i_xyz == 1) then
c             Overlap: row 1 of [C]*[B]  ###########
              j_xyz = 3
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2_p2x(j,i) * (-beta)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2y(i,j) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2x(i,j) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 5 of [C]*[B]  ###########
              j_xyz = 3
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2_p2(i,j) * beta**2 * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2y(i,j) * beta
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2x(i,j) * beta
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 6 of [C]*[B]  ###########
              j_xyz = 3
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2_p2y(j,i) * (-beta)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+5,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2y_p2y(i,j) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+5,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2y(j,i) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+5,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            elseif (i_xyz == 2) then
c             Overlap: row 2 of [C]*[B]  ###########
              j_xyz = 1
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2x_p2y(j,i)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2y(j,i) * beta * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2y_p2y(i,j)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 4 of [C]*[B]  ###########
              j_xyz = 1
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2_p2x(i,j) * (-beta) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2(i,j) * beta**2
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2y(i,j) * (-beta) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 6 of [C]*[B]  ###########
              j_xyz = 1
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2x_p2x(i,j)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2x(j,i) * beta * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,5)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2y(i,j)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+4,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            elseif (i_xyz == 3) then
c             Overlap: row 3 of [C]*[B]  ###########
              j_xyz = 2
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2_p2y(i,j) * (-beta)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2(i,j) * (-beta**2) * ii
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2x(i,j) * (-beta)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 4 of [C]*[B]  ###########
              j_xyz = 2
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2y_p2y(i,j) * (-ii)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+1,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2y(j,i) * beta
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+1,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2y(j,i) * (-ii)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+1,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
c             Overlap: row 5 of [C]*[B]  ###########
              j_xyz = 2
              j_p = 3*(j-1) + j_xyz
              z_tmp1 = p2x_p2y(i,j) * (-ii)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,j_xyz)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2_p2x(j,i) * beta
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,4)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
              z_tmp1 = p2x_p2x(i,j) * (-ii)
              z_tmp1 = z_tmp1 * c_tensor_el(i_xyz+2,6)
              mat_K(i_p,j_p) = mat_K(i_p,j_p) + z_tmp1
            endif
          enddo
        enddo
      enddo

cc
cc
      if (debug .eq. 1) then
        open(4,file="Output/mat_K.txt",status='unknown')
        do i=1,18
          do j=1,18
            write(4,"(2(I6),6(e20.10))") i,j, mat_K(j,i),
     *                mat_K(i,j), mat_K(i,j) - conjg(mat_K(j,i))
          enddo
        enddo
        close(4)
      endif
cc
      if (debug .eq. 1) then
        open(4,file="Output/mat_M.txt",status='unknown')
        do i=1,18
          do j=1,18
            write(4,"(2(I6),6(e20.10))") i,j, mat_M(j,i),
     *                mat_M(i,j), mat_M(i,j) - conjg(mat_M(j,i))
          enddo
        enddo
        close(4)
      endif
cc
      if (debug .eq. 1) then
        open(4,file="Output/rho_el.txt",status="unknown")
          write(4,"(2(e20.10))") rho_el
        close(4)
        open(4,file="Output/c_tensor_el.txt",status="unknown")
          do j=1,6
            do i=1,6
              write(4,"(2(I6),2(e20.10))") i,j, c_tensor_el(i,j)
            enddo
          enddo
        close(4)
      endif
cc
cc
      if (debug .eq. 1) then
        open(4,file="Output/xel.txt",status='unknown')
        do i=1,6
          write(4,"(I6,2(e20.10))") i, xel(1,i), xel(2,i)
        enddo
        close(4)
      endif

cc      write(*,*)
cc      write(*,*) "mat_el: Aborting..."
cc      stop

c
      return
      end

