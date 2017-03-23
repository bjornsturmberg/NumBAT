
      subroutine asmbly_AC (i_base, nel, npt, neq, nnodes, 
     *  shift, beta, nb_typ_el, rho, c_tensor, 
     *  table_nod, type_el, ineq, 
     *  x, nonz, row_ind, col_ptr, 
     *  mat1_re, mat1_im, mat2, i_work, symmetry_flag, debug)


      implicit none
      integer*8 nel, npt, neq, nnodes
      integer*8 i_base, nb_typ_el, nonz
      integer*8 row_ind(nonz), col_ptr(neq+1)
      integer*8 type_el(nel)
      integer*8 table_nod(nnodes,nel), ineq(3,npt)
      integer*8 i_work(3*npt), symmetry_flag
      complex*16 shift, beta
      double precision x(2,npt)
c     rho: density
      complex*16 rho(nb_typ_el), c_tensor(6,6,nb_typ_el)
      complex*16 mat2(nonz)
      double precision mat1_re(nonz), mat1_im(nonz)

c     Local variables
      integer*8 nnodes_0
      parameter (nnodes_0 = 6)
      integer*8 nod_el(nnodes_0)
      double precision xel(2,nnodes_0)

      integer*8 i_base2, debug, ui, typ_e
      integer*8 i, j, k, iel, j1
      integer*8 jtest, jp, ind_jp, j_eq
      integer*8 itrial, ip, ind_ip, i_eq
      integer*8 col_start, col_end

      complex*16 mat_K(18,18), mat_M(18,18)
      complex*16 ii, c_tensor_el(6,6), rho_el
      complex*16 z_tmp1, z_tmp2

c
ccccccccccccccccccccccccccccccccccccccc
c
c      write(*,*) "B: nnodes = ", nnodes
c
c      write(*,*) "i_base, nel, npt, neq, nnodes = ",
c     *  i_base, nel, npt, neq, nnodes

      ui = 6

      if (debug .eq. 1) then
        write(ui,*) "asmbly_AC: beta = ", beta
      endif

c
c     The CSC indexing, i.e., col_ptr, is 1-based 
c      But valpr.f may have changed the CSC indexing to 0-based indexing)
      if (i_base .eq. 0) then
        i_base2 = 1
      else
        i_base2 = 0
      endif
c
      if (nnodes .ne. 6 ) then
        write(ui,*) "asmbly_AC: problem nnodes = ", nnodes
        write(ui,*) "asmbly_AC: nnodes should be equal to 6 !"
        write(ui,*) "asmbly_AC: Aborting..."
        stop
      endif
c
ccccccccccccccccccccccccccccccccccccccc
c
c  ii = sqrt(-1)
      ii = dcmplx(0.0d0, 1.0d0)

      do i=1,nonz
        mat1_re(i) = 0.d0
        mat1_im(i) = 0.d0
        mat2(i) = 0.d0
      enddo

C
C#####################  Loop over Elements  ##########################
C
      do iel=1,nel
cc        write(ui,*) "asmbly_AC: iel, nel", iel, nel
        typ_e = type_el(iel)
        do j=1,nnodes
          j1 = table_nod(j,iel)
          nod_el(j) = j1
          xel(1,j) = x(1,j1)
          xel(2,j) = x(2,j1)
        enddo
        rho_el = rho(typ_e)
        do j=1,6
          do i=1,6
            c_tensor_el(i,j) = c_tensor(i,j,typ_e)
          enddo
        enddo
C       If c_tensor has regular symmetries use more efficient formulation
        if (symmetry_flag .eq. 1) then
          call mat_el_v2 (xel,beta,c_tensor_el,rho_el,mat_K,mat_M,debug)
        elseif (symmetry_flag .eq. 0) then
          call mat_el_v3 (xel,beta,c_tensor_el,rho_el,mat_K,mat_M,debug)
        else
          write(ui,*) "asmbly_AC: problem with symmetry_flag "
          write(ui,*) "asmbly_AC: c_tensor = ", symmetry_flag
          write(ui,*) "asmbly_AC: Aborting..."
          stop
        endif

        do jtest=1,nnodes
          jp = table_nod(jtest,iel)
          do j_eq=1,3
            ind_jp = ineq(j_eq,jp)
            if (ind_jp .gt. 0) then
              col_start = col_ptr(ind_jp) + i_base2
              col_end = col_ptr(ind_jp+1) - 1 + i_base2
c             unpack row into i_work
              do i=col_start,col_end
                i_work(row_ind(i) + i_base2) = i
              enddo
              do itrial=1,nnodes
                do i_eq=1,3
                  ip = table_nod(itrial,iel)
                  ind_ip = ineq(i_eq,ip)
                  if (ind_ip .gt. 0) then
                      z_tmp1 = mat_K(3*(itrial-1) + i_eq,
     *                               3*(jtest-1) + j_eq)
                      z_tmp2 = mat_M(3*(itrial-1) + i_eq,
     *                               3*(jtest-1) + j_eq)
                      z_tmp1 = z_tmp1 - shift*z_tmp2
                      k = i_work(ind_ip)
                      if (k .gt. 0 .and. k .le. nonz) then
                        mat1_re(k) = mat1_re(k) + dble(z_tmp1)
                        mat1_im(k) = mat1_im(k) + imag(z_tmp1)
                        mat2(k) = mat2(k) + z_tmp2
                      else
                        write(ui,*) "asmbly_AC: problem with row_ind !!"
                        write(ui,*) "asmbly_AC: k, nonz = ", k, nonz
                        write(ui,*) "asmbly_AC: Aborting..."
                        stop
                      endif
                  endif
                enddo
              enddo
            endif
          enddo
        enddo
      enddo
C
C#########################  End Element Loop  ########################
C

      if (debug .eq. 1) then
        write(ui,*) "asmbly_AC: shift = ", shift
      endif
c
      return
      end
