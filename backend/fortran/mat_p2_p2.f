c
c***********************************************************************
c
c      Compute the integral: m_ij = Integrate[P(i) * P(j), over a triangle]
c
c***********************************************************************
c
      subroutine mat_p2_p2 (mat, det_b)
c
c***********************************************************************
c
      implicit none
      double precision mat(6,6), det_b
      double precision dt, factor, factor2
      integer*8 i, j

      dt=dabs(det_b)

      mat(1,1) = 6
      mat(1,2) = -1
      mat(1,3) = -1
      mat(1,4) = 0
      mat(1,5) = -4
      mat(1,6) = 0
      mat(2,1) = -1
      mat(2,2) = 6
      mat(2,3) = -1
      mat(2,4) = 0
      mat(2,5) = 0
      mat(2,6) = -4
      mat(3,1) = -1
      mat(3,2) = -1
      mat(3,3) = 6
      mat(3,4) = -4
      mat(3,5) = 0
      mat(3,6) = 0
      mat(4,1) = 0
      mat(4,2) = 0
      mat(4,3) = -4
      mat(4,4) = 32
      mat(4,5) = 16
      mat(4,6) = 16
      mat(5,1) = -4
      mat(5,2) = 0
      mat(5,3) = 0
      mat(5,4) = 16
      mat(5,5) = 32
      mat(5,6) = 16
      mat(6,1) = 0
      mat(6,2) = -4
      mat(6,3) = 0
      mat(6,4) = 16
      mat(6,5) = 16
      mat(6,6) = 32

      factor = 360
      factor2 = dt / factor
      do j=1,6
        do i=1,6
          mat(i,j) = mat(i,j) * factor2
        enddo
      enddo
c
      return
      end

