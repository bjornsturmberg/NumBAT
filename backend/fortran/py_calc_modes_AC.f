      subroutine calc_AC_modes(
c     Explicit inputs
     *    beta_in, nval,
     *    debug, mesh_file, npt, nel,
     *    nb_typ_el,  c_tensor, rho, d_in_m, shift,
     *    i_cond, itermax, tol,
     *    plot_modes, cmplx_max, real_max,
     *    int_max, supplied_geo_flag, type_nod, symmetry_flag,
c     Inputs and outputs
     *    table_nod, type_el, x_arr,
c     Outputs
     *    beta1, sol1, mode_pol)

C***********************************************************************
C
C  Program:
C     FEM solver of Acoustic waveguide problems.
C     This subroutine is compiled by f2py & called in mode_calcs.py
C
C  Authors:
C    Bjorn Sturmberg & Kokou B. Dossou
C
C***********************************************************************
C
      implicit none
C  Local parameters:
C       ! Propagation constant
      complex*16 beta_in 
      integer*8 int_max, cmplx_max, int_used, cmplx_used
      integer*8 real_max, real_used, plot_modes
      integer :: alloc_stat=0
C       !  (int_max)
      integer*8, dimension(:), allocatable :: a   
C       !  (cmplx_max)
      complex*16, dimension(:), allocatable :: b   
C       !  (real_max)
      double precision, dimension(:), allocatable :: c   
      integer*8 supplied_geo_flag, symmetry_flag

C  Declare the pointers of the integer super-vector
      integer*8 ip_eq
      integer*8 ip_visite

C  Declare the pointers of the real super-vector
      integer*8 jp_x, jp_mat2
      integer*8 jp_vect1, jp_vect2, jp_workd, jp_resid, jp_vschur
      integer*8 jp_trav, jp_vp, jp_rhs
      integer*8 jp_eigenval_tmp, jp_eigen_pol

c     Declare the pointers of the real super-vector
      integer*8 kp_mat1_re, kp_mat1_im
      integer*8 kp_rhs_re, kp_rhs_im, kp_lhs_re, kp_lhs_im

c     Declare the pointers of for sparse matrix storage
      integer*8 ip_col_ptr, ip_row
      integer*8 ip_work, ip_work_sort, ip_work_sort2
      integer*8 nonz, nonz_max, max_row_len

      integer*8 nb_typ_el
      complex*16 c_tensor(6,6,nb_typ_el)
      complex*16 rho(nb_typ_el)

      integer*8 i, j, ip
      integer*8 nnodes, ui, debug, namelength
      integer*8 nel, npt, i_cond, neq

C     ! Number of nodes per element
      parameter(nnodes=6)
      integer*8 type_nod(npt), type_el(nel), table_nod(nnodes, nel)

      double precision pi
      double precision lat_vecs(2,2)
      double precision lx, ly, d_in_m

      complex*16 shift
      integer*8  i_base
      complex*16 ii

C  Variable used by valpr
      integer*8 ltrav, n_conv
      double precision ls_data(10)
      complex*16 z_beta, z_tmp, z_tmp0
      integer*8, dimension(:), allocatable :: index
c     variable used by UMFPACK
      double precision control (20), info_umf (90)
      integer*8 numeric, status, filenum

      double precision time1, time2
      character*(8) start_date, end_date
      character*(10) start_time, end_time

C  Variable used by valpr
      integer*8 nval, nvect, itermax
      double precision tol

C  Names and Controls
      character mesh_file*1000, gmsh_file*1000, log_file*1000
      character gmsh_file_pos*1000, dir_name*1000
      character*1000 tchar

c     new breed of variables to prise out of a, b and c
      double precision x_arr(2,npt)
C       complex*16, target :: sol1(3,nnodes+7,nval,nel)
      complex*16, target :: sol1(3,nnodes,nval,nel)
      complex*16, target :: beta1(nval)
      complex*16 mode_pol(4,nval)

      integer*8 ival, iel, inod

Cf2py intent(in) beta_in, nval
Cf2py intent(in) debug, mesh_file, npt, nel
Cf2py intent(in) d_in_m, shift
Cf2py intent(in) i_cond, itermax, tol
Cf2py intent(in) plot_modes, c_tensor, rho
Cf2py intent(in) cmplx_max, real_max, int_max
Cf2py intent(in) nb_typ_el, supplied_geo_flag,
Cf2py intent(in) type_nod, table_nod, type_el, x_arr, symmetry_flag

C  Note: the dependent variables must be listed AFTER the
C  independent variables that they depend on in the function call!
Cf2py depend(c_tensor) nb_typ_el
Cf2py depend(rho) nb_typ_el
Cf2py depend(type_nod) npt
Cf2py depend(table_nod) nnodes, nel
Cf2py depend(type_el) nel
Cf2py depend(x_arr) npt

Cf2py intent(out) beta1
Cf2py intent(out) sol1, mode_pol, table_nod, type_el, x_arr

C
CCCCCCCCCCCCCCCCCCCC  Start Program - get parameters  CCCCCCCCCCCCCCCCCC
C
C     Set parameter for the super-vectors of integer and real numbers
C
C       !ui = Unite dImpression
      ui = 6     
C      nnodes = 6 ! Number of nodes per element
      pi = 3.141592653589793d0
c     ii = sqrt(-1)
      ii = cmplx(0.0d0, 1.0d0, 8)

C       nvect = 2*nval + nval/2 +3
      nvect = 3*nval + 3
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      allocate(b(cmplx_max), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "The allocation is unseccesfull"
        write(*,*) "alloc_stat = ", alloc_stat
        write(*,*) "Not enough memory for the complex array b"
        write(*,*) "cmplx_max = ", cmplx_max
        write(*,*) "Aborting..."
        stop
      endif

      allocate(c(real_max), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "The allocation is unseccesfull"
        write(*,*) "alloc_stat = ", alloc_stat
        write(*,*) "Not enough memory for the real array c"
        write(*,*) "real_max = ", real_max
        write(*,*) "Aborting..."
        stop
      endif

      allocate(a(int_max), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "The allocation is unseccesfull"
        write(*,*) "alloc_stat = ", alloc_stat
        write(*,*) "Not enough memory for the integer array a"
        write(*,*) "int_max = ", int_max
        write(*,*) "Aborting..."
        stop
      endif

      allocate(index(nval), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "The allocation is unsuccessful"
        write(*,*) "alloc_stat = ", alloc_stat
        write(*,*) "Not enough memory for index"
        write(*,*) "nval = ", nval
        write(*,*) "Aborting..."
        stop
      endif
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C     clean mesh_format
      namelength = len_trim(mesh_file)
      gmsh_file = mesh_file(1:namelength-5)//'.msh'
      gmsh_file_pos = mesh_file(1:namelength)
      log_file = mesh_file(1:namelength-5)//'-AC.log'
      if (debug .eq. 1) then
        write(*,*) "mesh_file = ", mesh_file
        write(*,*) "gmsh_file = ", gmsh_file
      endif

C       ! initial time  in unit = sec.
      call cpu_time(time1)  
      call date_and_time ( start_date, start_time )
C
C      tol = 0.0 ! ARPACK accuracy (0.0 for machine precision)
C      lx=1.0 ! Diameter of unit cell. Default, lx = 1.0.
C      ly=1.0 ! NOTE: currently requires ly=lx, ie rectangular unit cell
C     ToDo: sort out what's going on here - in EM lx=ly=1
      lx = d_in_m
      ly = d_in_m
      shift = (2*pi*shift)**2
C
C####################  Start FEM PRE-PROCESSING  #######################
C
      if (debug .eq. 1) then
        write(ui,*)
        write(ui,*) "lx,ly = ", lx, ly
        write(ui,*) "npt, nel, nnodes = ", npt, nel, nnodes
        write(ui,*) "mesh_file = ", mesh_file
        write(ui,*)
      endif
C
      if ((3*npt+nel+nnodes*nel) .gt. int_max) then
         write(ui,*) "py_calc_modes_AC: "
         write(ui,*) "(3*npt+nel+nnodes*nel) + npt > int_max : ",
     *    (3*npt+nel+nnodes*nel), int_max
         write(ui,*) "py_calc_modes_AC: increase the size of int_max"
         write(ui,*) "py_calc_modes_AC: Aborting..."
         stop
      endif
      if ((7*npt) .gt. cmplx_max) then
         write(ui,*) "py_calc_modes_AC: (7*npt) > cmplx_max : ",
     *    (7*npt), cmplx_max
         write(ui,*) "py_calc_modes_AC: increase the size of cmplx_max"
         write(ui,*) "py_calc_modes_AC: Aborting..."
         stop
      endif

C       ! pointer to FEM connectivity table
      ip_visite = 1 
      ip_eq = ip_visite + npt
      jp_x = 1
C
      if (supplied_geo_flag .eq. 0) then
        call geometry (nel, npt, nnodes, nb_typ_el,
     *     lx, ly, type_nod, type_el, table_nod,
     *     x_arr, mesh_file)
      endif

      call lattice_vec (npt, x_arr, lat_vecs, debug)

C       if (debug .eq. 1) then
C       open (unit=64, file="msh_check.txt",
C      *         status="unknown")
C         do i=1,nel
C           write(64,*) i, type_el(i)
C         enddo
C         write(64,*)
C         write(64,*)
C         write(64,*)
C         do i=1,nel
C           do j=1,nnodes
C             write(64,*) i, j, table_nod(j,i)
C           enddo
C         enddo
C         write(64,*)
C         write(64,*)
C         write(64,*)
C         do j=1,nnodes
C           write(64,*) j, type_nod(j)
C         enddo
C       close(63)
C       endif



      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: npt, nel = ", npt, nel
      endif

      call bound_cond_AC (i_cond, npt, neq, type_nod,
     *       a(ip_eq), debug)
C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c     Sparse matrix storage
      ip_col_ptr = ip_eq + 3*npt

      call csr_max_length_AC (nel, npt, neq, nnodes,
     *  table_nod, a(ip_eq), a(ip_col_ptr), nonz_max)

      ip = ip_col_ptr + neq + 1
      if (ip .gt. int_max) then
         write(ui,*) "py_calc_modes_AC: ip > int_max : ",
     *    ip, int_max
         write(ui,*) "py_calc_modes_AC: nonz_max = ", nonz_max
         write(ui,*) "py_calc_modes_AC: increase the size of int_max"
         write(ui,*) "py_calc_modes_AC: Aborting..."
         stop
      endif
c
      ip_row = ip_col_ptr + neq + 1

      call csr_length_AC (nel, npt, neq, nnodes, table_nod,
     *  a(ip_eq), a(ip_row), a(ip_col_ptr), nonz_max,
     *  nonz, max_row_len, ip, int_max, debug)

      ip_work = ip_row + nonz
      ip_work_sort = ip_work + 3*npt
      ip_work_sort2 = ip_work_sort + max_row_len

c     sorting csr ...
      call sort_csr (neq, nonz, max_row_len, a(ip_row),
     *  a(ip_col_ptr), a(ip_work_sort), a(ip_work),
     *  a(ip_work_sort2))

      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: nonz_max = ", nonz_max
        write(ui,*) "py_calc_modes_AC: nonz = ", nonz
        write(ui,*) "py_calc_modes_AC: cmplx_max/nonz = ",
     *    dble(cmplx_max)/dble(nonz)
      endif

      int_used = ip_work_sort2 + max_row_len

      if (int_max .lt. int_used) then
        write(ui,*)
        write(ui,*) "The size of the integer supervector is too small"
        write(ui,*) "integer super-vec: int_max  = ", int_max
        write(ui,*) "integer super-vec: int_used = ", int_used
        write(ui,*) "Aborting..."
        stop
      endif

      jp_rhs = jp_x + 2*npt
c     jp_rhs will also be used (in gmsh_post_process) to store a solution
      jp_mat2 = jp_rhs + max(neq, 3*npt)
      jp_vect1 = jp_mat2 + nonz
      jp_vect2 = jp_vect1 + neq
      jp_workd = jp_vect2 + neq
      jp_resid = jp_workd + 3*neq
      jp_eigenval_tmp = jp_resid + 3*nnodes*nval*nel
C       ! Eigenvectors
      jp_vschur = jp_eigenval_tmp + nval + 1     
      jp_eigen_pol = jp_vschur + neq*nvect
      jp_trav = jp_eigen_pol + nval*4

      ltrav = 3*nvect*(nvect+2)
      jp_vp = jp_trav + ltrav
      cmplx_used = jp_vp + neq*nval

      if (cmplx_max .lt. cmplx_used)  then
         write(ui,*) "The size of the real supervector is too small"
         write(ui,*) "real super-vec: cmplx_max  = ", cmplx_max
         write(ui,*) "real super-vec: cmplx_used = ", cmplx_used
         write(ui,*) "Aborting..."
         stop
      endif

c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      kp_rhs_re = 1
      kp_rhs_im = kp_rhs_re + neq
      kp_lhs_re = kp_rhs_im + neq
      kp_lhs_im = kp_lhs_re + neq
      kp_mat1_re = kp_lhs_im + neq
      kp_mat1_im = kp_mat1_re + nonz
      real_used = kp_mat1_im + nonz

      if (real_max .lt. real_used) then
        write(ui,*)
        write(ui,*) "The size of the real supervector is too small"
        write(ui,*) "2*nonz  = ", 2*nonz
        write(ui,*) "real super-vec: real_max  = ", real_max
        write(ui,*) "real super-vec: real_used = ", real_used
        write(ui,*) "Aborting..."
        stop
      endif

c
c###############################################
c
c       ----------------------------------------------------------------
c       convert from 1-based to 0-based
c       ----------------------------------------------------------------
c
        do 60 j = 1, neq+1
            a(j+ip_col_ptr-1) = a(j+ip_col_ptr-1) - 1
60      continue
        do 70 j = 1, nonz
            a(j+ip_row-1) = a(j+ip_row-1) - 1
70      continue
c
c
c     The CSC indexing, i.e., ip_col_ptr, is 1-based
c       (but valpr.f will change the CSC indexing to 0-based indexing)
      i_base = 0

C#####################  End FEM PRE-PROCESSING  #########################
C
      write(ui,*)
      write(ui,*) "-----------------------------------------------"
C       write(ui,*) " AC FEM, k_AC : ", real(beta_in), " 1/m"
C       write(ui,*) "-----------------------------------------------"
C       write(ui,*)

C       if (debug .eq. 1) then
C         write(ui,*) "py_calc_modes_AC: call to asmbly"
C       endif
      write(ui,*) "AC FEM, assembling linear system"
      call cpu_time(time1)
C     Assemble the coefficient matrix K and M of the finite element equations
      call asmbly_AC (i_base, nel, npt, neq, nnodes,
     *  shift, beta_in, nb_typ_el, rho, c_tensor,
     *  table_nod, type_el, a(ip_eq),
     *  x_arr, nonz, a(ip_row), a(ip_col_ptr),
     *  c(kp_mat1_re), c(kp_mat1_im), b(jp_mat2), a(ip_work), 
     *  symmetry_flag, debug)
      call cpu_time(time2)
      write(ui,*) '    assembly time (sec.)  = ', (time2-time1)
C
C       if (debug .eq. 1) then
C         write(ui,*) "py_calc_modes_AC: call to valpr"
C       endif
      write(ui,*) "AC FEM, solving linear system"
      call cpu_time(time1)
      call valpr_64_AC (i_base, nvect, nval, neq, itermax, ltrav,
     *  tol, nonz, a(ip_row), a(ip_col_ptr), c(kp_mat1_re),
     *  c(kp_mat1_im), b(jp_mat2), b(jp_vect1), b(jp_vect2),
     *  b(jp_workd), b(jp_resid), b(jp_vschur), beta1,
     *  b(jp_trav), b(jp_vp), c(kp_rhs_re), c(kp_rhs_im),
     *  c(kp_lhs_re), c(kp_lhs_im), n_conv, ls_data,
     *  numeric, filenum, status, control, info_umf, debug)
      call cpu_time(time2)
      write(ui,*) '    eigensolver time (sec.)  = ', (time2-time1)

      if (n_conv .ne. nval) then
         write(ui,*) "py_calc_modes_AC: convergence problem in valpr_64"
         write(ui,*) "py_calc_modes_AC: n_conv != nval : ",
     *    n_conv, nval
         write(ui,*) "py_calc_modes_AC: Aborting..."
c         stop
      endif
C
      do i=1,nval
        z_tmp0 = beta1(i)
        z_tmp = 1.0d0/z_tmp0+shift
        z_beta = sqrt(z_tmp) / (2.0d0 * pi)
C       Frequency (z_beta) should always be positive.
        if (dble(z_beta) .lt. 0) z_beta = -z_beta
        beta1(i) = z_beta
      enddo
c
      call z_indexx_AC (nval, beta1, index)
C
C       The eigenvectors will be stored in the array sol1
C       The eigenvalues and eigenvectors will be renumbered
C                 using the permutation vector index
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: call to array_sol"
      endif
        call array_sol_AC (nval, nel, npt, neq, nnodes,
     *   index, table_nod, type_el, a(ip_eq), x_arr, beta1,
     *   b(jp_eigenval_tmp), mode_pol, b(jp_vp), sol1)

      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: array_sol returns call"
      endif
C
      if(debug .eq. 1) then
        write(ui,*) 'index = ', (index(i), i=1,nval)
      endif
      if(debug .eq. 1) then
        write(ui,*)
C         write(ui,*) "lambda, 1/lambda = ", lambda, 1.0d0/lambda
C         write(ui,*) "sqrt(shift)/(2*pi) = ", sqrt(shift) / (2.0d0 * pi)
        do i=1,nval
          write(ui,"(i4,2(g22.14),2(g18.10))") i, beta1(i)
        enddo
      endif
      
C    Save Original solution
      if (plot_modes .eq. 1) then
        dir_name = "AC_fields"
C        call write_sol_AC (nval, nel, nnodes, lambda,
C      *       beta1, sol1, mesh_file, dir_name)
C        call write_param (lambda, npt, nel, i_cond,
C    *       nval, nvect, itermax, tol, shift, lx, ly,
C    *       mesh_file, n_conv, dir_name)
        tchar = "AC_fields/All_plots_png_abs2_eE.geo"
        open (unit=34,file=tchar)
          do i=1,nval
            call gmsh_post_process_AC (i, nval, nel, npt,
     *         nnodes, table_nod, type_el,
     *         x_arr, beta1, sol1, b(jp_rhs), a(ip_visite),
     *         gmsh_file_pos, dir_name, d_in_m, debug)
          enddo
        close (unit=34)
      endif
C
C
C#########################  End Calculations  ###########################
C
      call date_and_time ( end_date, end_time )
      call cpu_time(time2)
C
      if (debug .eq. 1) then
        write(ui,*)
        write(ui,*) 'Total CPU time (sec.)  = ', (time2-time1)
C
        open (unit=26,file=log_file)
        write(26,*)
        write(26,*) "Date and time formats = ccyymmdd ; hhmmss.sss"
        write(26,*) "Start date and time   = ", start_date,
     *    " ; ", start_time
        write(26,*) "End date and time     = ", end_date,
     *    " ; ", end_time
        write(26,*) "Total CPU time (sec.) = ",  (time2-time1)
        write(26,*)
        write(26,*) "beta_in = ", beta_in
        write(26,*) "shift = ", shift
        write(26,*)
        write(26,*) "npt, nel, nnodes  = ", npt, nel, nnodes
        write(26,*) "neq, i_cond = ", neq, i_cond
        write(26,*) " lat_vecs:  = "
        write(26,"(2(f18.10))") lat_vecs
        write(26,*) "mesh_file = ", mesh_file
        write(26,*) "gmsh_file = ", gmsh_file
        write(26,*) "log_file  = ", log_file
        close(26)
C
        write(ui,*) "   .      .      ."
        write(ui,*) "   .      .      ."
        write(ui,*) "   .      . (d=",d_in_m,")"
        write(ui,*) "  and   we're   done!"
      endif

      write(ui,*) "-----------------------------------------------"
      write(ui,*)
C
      deallocate(a,b,c,index)

      end subroutine calc_AC_modes
