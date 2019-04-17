
c
c***********************************************************************
c
c	Construction of the matrix of power flow
c	(integral of the z-component of the acoustic Poynting vector)
c
c***********************************************************************
c
      subroutine mat_el_energy_rho (xel, rho_el, mat_P)
c
c***********************************************************************
c
      implicit none
      double precision xel(2,6)
      complex*16 beta
      complex*16 mat_P(6,6)
      complex*16 rho_el

c     Local variables
      double precision p2_p2(6,6)
      double precision mat_B(2,2), mat_T(2,2), mat_T_tr(2,2)
      double precision det_b
      complex*16 ii, z_tmp1, z_tmp2
      integer*8 i, j
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
      if (abs(det_b) .le. 1.0d-22) then
cc      if (abs(det_b) .le. 1.0d-8) then
        write(*,*) '?? mat_el_energy_rho: Determinant = 0 :', det_b
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
c	mat_T_tr = Tanspose(mat_T)
c
      mat_T_tr(1,1) = mat_T(1,1)
      mat_T_tr(1,2) = mat_T(2,1)
      mat_T_tr(2,1) = mat_T(1,2)
      mat_T_tr(2,2) = mat_T(2,2)

      call mat_p2_p2(p2_p2, det_b)
c
      do i=1,6
        do j=1,6
          mat_P(j,i) = 0
        enddo
      enddo
c=================  Construction of the matrix of power flow   ==================
c                   (integral of the z-component of the acoustic Poynting vector)
      do i=1,6
        do j=1,6
          z_tmp1 = p2_p2(i,j) * rho_el
          mat_P(i,j) = mat_P(i,j) + z_tmp1
        enddo
      enddo
cc
cc
cc      if (debug .eq. 1) then
cc        open(4,file="Output/mat_P_rho.txt",status='unknown')
cc        write(4,*) "rho_el = ", rho_el
cc        do i=1,6
cc          do j=1,6
cc            write(4,"(2(I6),6(e20.10))") i,j, mat_P(j,i),
cc     *                mat_P(i,j), mat_P(i,j) - conjg(mat_P(j,i))
cc          enddo
cc        enddo
cc        close(4)
cc      endif
cc
cc
cc      if (debug .eq. 1) then
cc        open(4,file="Output/xel_rho.txt",status='unknown')
cc        do i=1,6
cc          write(4,"(I6,2(e20.10))") i, xel(1,i), xel(2,i)
cc        enddo
cc        close(4)
cc      endif


c
      return
      end

