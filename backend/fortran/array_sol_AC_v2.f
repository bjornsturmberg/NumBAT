C Difference from array_sol_AC.f is that the u_z field is multiplied by i
C which gives you the correct physical displacement field.

c   sol_0(*,i) : contains the imaginary and real parts of the solution for points such that ineq(i) != 0
c   sol(i) : contains solution for all points
c   The dimension of the geometric domain is : dim_32 = 2
c   The dimension of the vector field is : dim2 = 3
c

      subroutine array_sol_AC (nval, nel, npt, neq,
     *     nnodes, index, table_nod, type_el, ineq,
     *     x, v_cmplx, v_tmp, mode_pol, sol_0, sol)


      implicit none
      integer*8 nval, nel, npt, neq, nnodes
      integer*8 n_core(2), type_el(nel)
      integer*8 ineq(3,npt), index(*)
      integer*8 table_nod(nnodes,nel)
      complex*16 sol_0(neq,nval)
      double precision x(2,npt)
c     sol(3, 1..nnodes,nval, nel)          contains the values of the 3 components at P2 interpolation nodes
      complex*16 sol(3,nnodes,nval,nel)
      complex*16 v_cmplx(nval), v_tmp(nval), mode_pol(4,nval)



c     Local variables
      integer*8 nnodes_0
c     32-but integers for BLAS and LAPACK
      parameter (nnodes_0 = 6)
c
      double precision mode_comp(4)
      integer*8 nod_el_p(nnodes_0)
      double precision xel(2,nnodes_0)
      complex*16 sol_el(3,nnodes_0)

      double precision ZERO, ONE
      parameter ( ZERO = 0.0D0, ONE = 1.0D0)
c
      integer*8 j, i1, j1, inod, typ_e, debug
      integer*8 iel, ival, ival2, jtest, jp, ind_jp, j_eq
      complex*16 ii, z_tmp1, z_tmp2, z_sol_max

      double precision x_min, x_max, y_min, y_max
      double precision x_mid, y_mid
      double precision dx, dy, x_0, y_0
      double precision lx, ly, rx, ry

      integer*8 i_sol_max, i_sol_max_tmp
      integer*8 i_component, i_component_tmp

c
c  ii = sqrt(-1)
      ii = cmplx(0.0d0, 1.0d0, 8)
c
      debug = 0
c
      if ( nnodes .ne. 6 ) then
        write(*,*) "array_sol: problem nnodes = ", nnodes
        write(*,*) "array_sol: nnodes should be equal to 6 !"
        write(*,*) "array_sol: Aborting..."
        stop
      endif


      do ival=1,nval
        do iel=1,nel
          do inod=1,nnodes
            do j=1,3
              sol(j,inod,ival,iel) = 0
            enddo
          enddo
        enddo
      enddo

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      x_min = x(1,1)
      x_max = x(1,1)
      do j=1,npt
        x_0 = x(1,j)
        if(x_0 .lt. x_min) x_min = x_0
        if(x_0 .gt. x_max) x_max = x_0
      enddo
      y_min = x(2,1)
      y_max = x(2,1)
      do j=1,npt
        y_0 = x(2,j)
        if(y_0 .lt. y_min) y_min = y_0
        if(y_0 .gt. y_max) y_max = y_0
      enddo

      x_mid = (x_min + x_max) / 2.0d0
      y_mid = (y_min + y_max) / 2.0d0

C       !  Length in the x direction
      lx = x_max - x_min  
C       !  Length in the y direction
      ly = y_max - y_min  

      rx = sqrt(nel * lx / (2.0d0 * ly))
      ry = rx * ly/lx
c     rx and ry and such that : rx * ry = 2 * nel


      dx = lx / rx
      dy = ly / ry


cc      lx = x_max - x_min  !  Length in the x direction
cc      ly = y_max - y_min  !  Length in the y direction

      if (debug .eq. 1) then
        write(*,*)
        write(*,*) "array_sol: x_min = ", x_min
        write(*,*) "array_sol: x_max = ", x_max
        write(*,*) "array_sol: y_min = ", y_min
        write(*,*) "array_sol: y_max = ", y_max
        write(*,*) "array_sol: x_mid = ", x_mid
        write(*,*) "array_sol: y_mid = ", y_mid
        write(*,*) "array_sol:    dx = ", dx
        write(*,*) "array_sol:    dy = ", dy
        write(*,*) "array_sol:    rx = ", rx
        write(*,*) "array_sol:    ry = ", ry
        write(*,*) "array_sol:   nel = ", nel
        write(*,*)
      endif




c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c
      do j=1,nval
        j1=index(j)
        v_tmp(j) = v_cmplx(j1)
      enddo
      do j=1,nval
        v_cmplx(j) = v_tmp(j)
      enddo
c
      do ival=1,nval
        ival2 = index(ival)
        do j=1,4
          mode_pol(j,ival) = 0.0d0
        enddo

        z_sol_max = 0.0d0
        i_sol_max = 0
        i_component = 0

        do iel=1,nel
          typ_e = type_el(iel)
          do j=1,4
            mode_comp(j) = 0.0d0
          enddo
          do inod=1,nnodes
            j = table_nod(inod,iel)
            nod_el_p(inod) = j
            xel(1,inod) = x(1,j)
            xel(2,inod) = x(2,j)
          enddo
          do inod=1,nnodes
            jp = table_nod(inod,iel) 
            do j_eq=1,3
              ind_jp = ineq(j_eq,jp)
              if (ind_jp .gt. 0) then
                z_tmp1 = sol_0(ind_jp, ival2)
                sol_el(j_eq,inod) = z_tmp1
              else
                sol_el(j_eq,inod) = 0
              endif
            enddo
c           The z-compoenent must be multiplied by ii in order to get the un-normalised z-compoenent
            j_eq=3
                sol_el(j_eq,inod) = ii * sol_el(j_eq,inod)
            do j=1,3
              z_tmp2 = sol_el(j,inod)
              sol(j,inod,ival,iel) = z_tmp2
              if (abs(z_sol_max) .lt. abs(z_tmp2)) then
                z_sol_max = z_tmp2
c               We want to normalise such the the z-component is purely imaginary complex number
                if (j == 3) z_sol_max = - ii * z_sol_max
                i_sol_max = table_nod(inod,iel)
                i_component = j
              endif
            enddo
c           Contribution of the element iel to the mode component
            do j=1,3
              z_tmp2 = abs(sol_el(j,inod))**2
              mode_comp(j) = mode_comp(j) + z_tmp2
            enddo
          enddo
c         Avarage values
          do j=1,3
            mode_comp(j) = mode_comp(j)/dble(nnodes)
c            mode_comp(j) = abs(det)*mode_comp(j)/dble(nnodes)
          enddo
c         Add the contribution of the element iel to the mode component
          do j=1,3
            mode_pol(j,ival) = mode_pol(j,ival) + mode_comp(j)
          enddo
          if (typ_e .eq. n_core(1) .or. typ_e .eq. n_core(2)) then
            z_tmp2 = mode_comp(1) + mode_comp(2)
     *        + mode_comp(3)
            mode_pol(4,ival) = mode_pol(4,ival) + z_tmp2
          endif
        enddo
c       Total energy and normalization
        z_tmp2 = mode_pol(1,ival) + mode_pol(2,ival)
     *        + mode_pol(3,ival)
        if (abs(z_tmp2) .lt. 1.0d-10) then
          write(*,*) "array_sol: the total energy ",
     *       "is too small : ", z_tmp2
          write(*,*) "array_sol: ival ival2 = ", ival, ival2
          write(*,*) "array_sol: zero eigenvector; aborting..."
          stop
        endif
        do j=1,3
          mode_pol(j,ival) = mode_pol(j,ival) / z_tmp2
        enddo
        j=4
          mode_pol(j,ival) = mode_pol(j,ival) / z_tmp2
c       Check if the eigenvector is nonzero
        if (abs(z_sol_max) .lt. 1.0d-10) then
          z_sol_max = z_tmp2
          write(*,*) "array_sol: z_sol_max is too small"
          write(*,*) "array_sol: z_sol_max = ", z_sol_max
          write(*,*) "ival, ival2, nval = ", ival, ival2, nval
          write(*,*) "array_sol: zero eigenvector; aborting..."
          stop
        endif
        if (debug .eq. 1) then
          write(*,*) "array_sol (A): "
          write(*,*) "                           ival = ", ival
          write(*,*) "z_sol_max     = ", z_sol_max, abs(z_sol_max)
        endif
c       Normalization so that the maximum field component is 1
        do iel=1,nel
          do inod=1,nnodes
            i1 = table_nod(inod,iel)
            do j=1,3
              z_tmp1 = sol(j,inod,ival,iel)/z_sol_max
              sol(j,inod,ival,iel) = z_tmp1
            enddo
            i1 = table_nod(inod,iel)
            if (i1 .eq. i_sol_max_tmp .and. debug .eq. 1) then
              write(*,*) "array_sol (B):"
              write(*,*) "ival, i1, iel = ", ival, i1, iel
              write(*,*) "array_sol: Field normalisaion point:"
              write(*,*) "x = ", dble(x(1,i1))
              write(*,*) "y = ", dble(x(2,i1))
              write(*,*) "i_sol_max = ", i_sol_max
              write(*,*) "i_component_tmp = ", i_component_tmp
              write(*,*) ival, i1, iel,
     *                   (dble(sol(j,inod,ival,iel)),j=1,3)
              write(*,*) ival, i1, iel,
     *                   (imag(sol(j,inod,ival,iel)),j=1,3)
              write(*,*) ival, i1, iel,
     *                   (abs(sol(j,inod,ival,iel)),j=1,3)
            endif
          enddo
        enddo
        do j=1,neq
          z_tmp1 = sol_0(j,ival2)/z_sol_max
          sol_0(j,ival2) = z_tmp1
        enddo

        if (debug .eq. 1) then
          write(*,*)
        endif
      enddo

      if (debug .eq. 1) then
        write(*,*) "Exiting subroutine array_sol"
      endif


c
      return
      end
