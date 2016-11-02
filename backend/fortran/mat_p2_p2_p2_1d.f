c
c***********************************************************************
c
c      Compute the integral: m_ijk = Integrate[P(i) * P(j) * P(k), over an interval]
c
c***********************************************************************
c
      subroutine mat_p2_p2_p2_1d (mat, det_b)
c
c***********************************************************************
c
      implicit none
      double precision mat(3,3,3), det_b

      double precision dt, factor, factor2
      integer*8 i, j, k

      dt=dabs(det_b)

      mat(1, 1, 1) = 39
      mat(1, 1, 2) = -3
      mat(1, 1, 3) = 20
      mat(1, 2, 1) = -3
      mat(1, 2, 2) = -3
      mat(1, 2, 3) = -8
      mat(1, 3, 1) = 20
      mat(1, 3, 2) = -8
      mat(1, 3, 3) = 16
      mat(2, 1, 1) = -3
      mat(2, 1, 2) = -3
      mat(2, 1, 3) = -8
      mat(2, 2, 1) = -3
      mat(2, 2, 2) = 39
      mat(2, 2, 3) = 20
      mat(2, 3, 1) = -8
      mat(2, 3, 2) = 20
      mat(2, 3, 3) = 16
      mat(3, 1, 1) = 20
      mat(3, 1, 2) = -8
      mat(3, 1, 3) = 16
      mat(3, 2, 1) = -8
      mat(3, 2, 2) = 20
      mat(3, 2, 3) = 16
      mat(3, 3, 1) = 16
      mat(3, 3, 2) = 16
      mat(3, 3, 3) = 192

      factor = 420
      factor2 = dt / factor
      do k=1,3
        do j=1,3
          do i=1,3
            mat(i,j,k) = mat(i,j,k) * factor2
          enddo
        enddo
      enddo
c
      return
      end



