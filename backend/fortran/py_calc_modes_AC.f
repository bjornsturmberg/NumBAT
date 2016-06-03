      program main

C************************************************************************
C
C  Program:
C    Finite Element Method - Acoustic waveguide problems
C
C  Authors:
C    Bjorn Sturmberg & Kokou B. Dossou
C
C************************************************************************
C
      implicit none
C  Local parameters:
      complex*16 beta ! Propagation constant
      integer*8 int_max, cmplx_max, int_used, cmplx_used
      integer*8 real_max, real_used, n_64
      integer alloc_stat

c     32-bit integer supervector for UMFPACK and ARPACK
      integer*4 int_max_32, int_used_32
      integer*4, dimension(:), allocatable :: a_32  !  (int_max_32)

      integer*8, dimension(:), allocatable :: a   !  (int_max)
      complex*16, dimension(:), allocatable :: b   !  (cmplx_max)
      double precision, dimension(:), allocatable :: c   !  (real_max)


C  Declare the pointers of the integer super-vector
      integer*8 ip_type_nod, ip_type_el, ip_table_nod
      integer*8 ip_eq
      integer*8 ip_visite

C  Declare the pointers of the real super-vector
      integer*8 jp_x, jp_mat2
      integer*8 jp_vect1, jp_vect2, jp_workd, jp_resid, jp_vschur
      integer*8 jp_eigenval, jp_sol, jp_trav, jp_vp, jp_rhs
      integer*8 jp_eigenval_tmp, jp_eigen_pol

c     Declare the pointers of the real super-vector
      integer*8 kp_mat1_re, kp_mat1_im
      integer*8 kp_rhs_re, kp_rhs_im, kp_lhs_re, kp_lhs_im

c     Declare the pointers of for sparse matrix storage
      integer*8 ip_col_ptr, ip_row
      integer*8 ip_work, ip_work_sort, ip_work_sort2
      integer*8 nonz, nonz_max, max_row_len

      integer*8 max_typ_el, nb_typ_el
      parameter (max_typ_el=10)
      complex*16 pp(max_typ_el),  qq(max_typ_el)
      complex*16 c_tensor(6,6,max_typ_el)
c     rho: density
      complex*16 rho(max_typ_el)

      integer*8 i, j, k, ip, Lambda_count
      integer*8 nnodes, ui, debug, PrintSolution
      integer*8 i_lambda, n_lambda
      integer*8 nel, npt, i_cond, neq

      double precision pi, theta, phi
      double precision lambda, lambda_1, lambda_2, d_lambda
      double precision d_freq, freq, lat_vecs(2,2)
      double precision lx, ly, d_in_nm

      complex*16 shift_0, shift

      integer*8 n_lambda_1, Loss
      parameter (n_lambda_1 = 101)
      complex*16 Complex_refract1(2,n_lambda_1)
      complex*16 Complex_refract_sub(2,n_lambda_1)


      integer*8  i_base

C  Variable used by valpr
      integer*8 ltrav, n_conv
      double precision ls_data(10)
      complex*16 z_beta, z_tmp, z_tmp0
      integer*8 index(1000)
c     variable used by UMFPACK
      double precision control (20), info_umf (90)
      integer*8 numeric, symbolic, status, sys, filenum

      integer*4 ip_col_ptr_32, ip_row_32

c     32-but integer variables used by UMFPACK
      integer*4 numeric_32, symbolic_32, status_32
      integer*4 sys_32, filenum_32
c     32-but integer variables used by UMFPACK and ARPACK
      integer*4 neq_32, nonz_32, debug_32
      integer*4 nval_32, nvect_32, itermax_32, ltrav_32
      integer*4 n_conv_32, i_base_32

      double precision time1, time2
      double precision time_Lambda_end, time_Lambda_start
      character*(8) start_date, end_date
      character*(10) start_time, end_time

C  Variable used by valpr
      integer*8 nval, nvect, itermax
      double precision tol

C  Names and Controls
      character mesh_file*100, gmsh_file*100, log_file*100
      character gmsh_file_pos*100, dir_name*100
      character*100 tchar

C
CCCCCCCCCCCCCCCCCCCC  Start Program - get parameters  CCCCCCCCCCCCCCCCCCC
C 
C     Set parameter for the super-vectors of integer and real numbers
C
      ui = 6     !ui = Unite dImpression
      nnodes = 6 ! Number of nodes per element
      pi = 3.141592653589793d0

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

      alloc_stat = 0

      n_64 = 2
      cmplx_max=n_64**24 * 2
      real_max=n_64**20 * 2
      int_max=n_64**20 * 2
      int_max_32=n_64**20

      write(*,*) "cmplx_max = ", cmplx_max
      write(*,*) "real_max = ", real_max
      write(*,*) "int_max = ", int_max
      write(*,*) "int_max_32 = ", int_max_32
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

      allocate(a_32(int_max_32), STAT=alloc_stat)
      if (alloc_stat /= 0) then
        write(*,*) "The allocation is unseccesfull"
        write(*,*) "alloc_stat = ", alloc_stat
        write(*,*) "Not enough memory for the integer array a_32"
        write(*,*) "int_max_32 = ", int_max_32
        write(*,*) "Aborting..."
        stop
      endif

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
C
      open(3,file="Parameters/controls.txt",status="old")
        read(3,*) debug
        read(3,*) PrintSolution
      close(3)
C
      open(4,file="Parameters/angles.txt",status="old")
        read(4,*) theta
        read(4,*) phi
      close(4)
c
      call cpu_time(time1)  ! initial time  in unit = sec.
      call date_and_time ( start_date, start_time )
C      
      if (debug .eq. 1) then
        write(ui,*)
        write(ui,*) "start_date = ", start_date
        write(ui,*) "start_time = ", start_time
        write(ui,*)
      endif
C
      lx = d_in_nm
      ly = d_in_nm
      call get_param(n_lambda, lambda_1, lambda_2, 
     *      npt, nel, i_cond, nval, nvect, itermax, tol,
     *      shift_0, lx, ly, mesh_file, gmsh_file,
     *      gmsh_file_pos, log_file, d_in_nm, debug)
C
C####################  Start FEM PRE-PROCESSING  ########################
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
         write(ui,*) "py_calc_modes_AC: (3*npt+nel+nnodes*nel) + npt > int_max : ",
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
      ip_type_nod = 1
      ip_type_el = ip_type_nod + npt
      ip_table_nod = ip_type_el + nel ! pointer to FEM connectivity table

      ip_visite = ip_table_nod + nnodes*nel
      ip_eq = ip_visite + npt

      jp_x = 1
C
      lx = d_in_nm
      ly = d_in_nm
      call geometry (nel, npt, nnodes, nb_typ_el,
     *     lx, ly, a(ip_type_nod), a(ip_type_el), a(ip_table_nod), 
     *     b(jp_x), mesh_file)

      call lattice_vec (npt, b(jp_x), lat_vecs)

      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: npt, nel = ", npt, nel
      endif

      call get_tensor (nb_typ_el, c_tensor, rho)

      if (debug .eq. 1) then
        open(4,file="Output/c_tensor.txt",status="unknown")
        do k=1,nb_typ_el
          do j=1,6
            do i=1,6
              write(4,"(3(I6),2(f20.10))") i,j,k, c_tensor(i,j,k)
            enddo
          enddo
        enddo
        close(4)
      endif
C
      if (PrintSolution .ge. 1) then
C  Export the mesh to gmsh format
        call mail_to_gmsh (nel, npt, nnodes, a(ip_type_el), 
     *    a(ip_type_nod), a(ip_table_nod), 
     *    nb_typ_el, b(jp_x), gmsh_file)
      endif

        call bound_cond_AC (i_cond, npt, neq, a(ip_type_nod), 
     *       a(ip_eq), debug)

C
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
c
c     Sparse matrix storage
c
      ip_col_ptr = ip_eq + 3*npt
c       ip_col_ptr = ip_index_pw_inv + neq_PW 

      call csr_max_length (nel, npt, neq, nnodes, 
     *  a(ip_table_nod), a(ip_eq), a(ip_col_ptr), nonz_max)

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

      call csr_length (nel, npt, neq, nnodes, a(ip_table_nod), 
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
      jp_sol = jp_resid + neq
      jp_eigenval = jp_sol + 3*nnodes*nval*nel
      jp_eigenval_tmp = jp_eigenval + nval + 1
      jp_vschur = jp_eigenval_tmp + nval + 1     ! Eigenvectors
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

      nvect_32 = nvect
      nval_32 = nval
      neq_32 = neq
      itermax_32 = itermax
      ltrav_32 = ltrav
      nonz_32 = nonz
      debug_32 = debug

c
c###############################################
c
c     SOME 32 bit integers for UMFPACK AND ARPACK
      ip_col_ptr_32 = 1
      ip_row_32 = ip_col_ptr_32 + neq + 1
      int_used_32 = ip_row_32 + nonz

      if (int_max_32 .lt. int_used_32) then
        write(ui,*)
        write(ui,*) "The size of the integer_32 supervector ",
     *              "is too small"
        write(ui,*) "integer super-vec: int_max_32  = ", int_max_32
        write(ui,*) "integer super-vec: int_used_32 = ", int_used_32
        write(ui,*) "int_used_32/int_max_32 = ", 
     *                dble(int_used_32)/dble(int_max_32)
        write(ui,*) "Aborting..."
        stop
      endif
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
c
c     SOME 32 bit integers for UMFPACK AND ARPACK
      i_base_32 = i_base

      do i=1,neq+1
        a_32(ip_col_ptr_32+i-1) = a(ip_col_ptr+i-1)
      enddo

      do i=1,nonz
          a_32(ip_row_32+i-1) = a(ip_row+i-1)
      enddo

cc      do i=1,nb_typ_el
cc        rho(i) = rho(i) * d_in_nm**3
cc      enddo

C
C#####################  End FEM PRE-PROCESSING  #########################
C

      if(n_lambda .gt. 1) then
        d_lambda = (lambda_2 - lambda_1)/dble(n_lambda-1)
        d_freq   = (1.0d0/lambda_2 - 1.0d0/lambda_1) /
     *             dble(n_lambda-1)
      else
        d_lambda = 0.0d0
        d_freq   = 0.0d0
      endif

      Lambda_count = 1    !  index wavelength loops for A_and_W_Lambda

C
C#####################  Loop over Wavelengths  ##########################
C
      do i_lambda=1,n_lambda
      call cpu_time(time_Lambda_start)

        write(ui,*) 
        write(ui,*) "--------------------------------------------",
     *     "-------"
        write(ui,*) "py_calc_modes_AC: Slice Progress  ", 
     *     i_lambda,"/", n_lambda

C
C  For uniformly spaced wavelengths
        lambda = lambda_1 + (i_lambda-1)*d_lambda
C
C  For uniformly spaced frequencies
C        freq = 1.0d0/lambda_1 + (i_lambda-1)*d_freq
C        lambda = 1.0d0/freq
cc        freq = 1.0d0/lambda

        write(ui,*) "py_calc_modes_AC: lambda = ", lambda


c       Propagation constant: beta = k_z
cc        beta = 28.9d0 * d_in_nm
cc        beta = 28.9d0 / d_in_nm

        beta = 1.49175d7

        shift = shift_0

cc        shift = 12.0d9
cc        shift = (2.0d0 * pi * shift)**2


cc        shift = 11.0d9
cc        shift = (2.0d0 * pi * shift)**2
cc
cc        shift = 11.0d9
cc        shift = (2.0d0 * pi * shift)**2



C
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: Asmbly: call to asmbly"
      endif


C     Assemble the coefficient matrix K and M of the finite element equations
      call asmbly_AC (i_base, nel, npt, neq, nnodes, 
     *  shift, beta, nb_typ_el, rho, c_tensor, 
     *  a(ip_table_nod), a(ip_type_el), a(ip_eq),
     *  b(jp_x), nonz, a(ip_row), a(ip_col_ptr), 
     *  c(kp_mat1_re), c(kp_mat1_im), b(jp_mat2), a(ip_work))
C
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: call to valpr"
      endif


      call valpr_32 ( i_base_32, nvect_32, nval_32, neq_32, 
     *  itermax_32, ltrav_32, tol, nonz_32, 
     *  a_32(ip_row_32), a_32(ip_col_ptr_32), c(kp_mat1_re),
     *  c(kp_mat1_im), b(jp_mat2), b(jp_vect1), b(jp_vect2),
     *  b(jp_workd), b(jp_resid), b(jp_vschur), b(jp_eigenval),
     *  b(jp_trav), b(jp_vp), c(kp_rhs_re), c(kp_rhs_im),
     *  c(kp_lhs_re), c(kp_lhs_im), n_conv_32, ls_data,
     *  numeric_32, filenum_32, status_32, control, 
     *  info_umf, debug_32)
      n_conv = n_conv_32

      if (n_conv .ne. nval) then
         write(ui,*) "py_calc_modes_AC: convergence problem with valpr_64"
         write(ui,*) "py_calc_modes_AC: n_conv != nval : ",
     *    n_conv, nval
         write(ui,*) "py_calc_modes_AC: Aborting..."
cc         stop
      endif

C
      do i=1,nval
        z_tmp0 = b(jp_eigenval+i-1)
        z_tmp = 1.0d0/z_tmp0+shift
        z_beta = sqrt(z_tmp) / (2.0d0 * pi)
C       Mode classification - we want the forward propagating mode
        if (abs(imag(z_beta)/z_beta) .lt. 1.0d-8) then
C         re(z_beta) > 0 for forward propagating mode
          if (dble(z_beta) .lt. 0) z_beta = -z_beta
        else
C         im(z_beta) > 0 for forward decaying evanescent mode
          if (imag(z_beta) .lt. 0) z_beta = -z_beta
        endif
        b(jp_eigenval+i-1) = z_beta
      enddo
c
      call z_indexx (nval, b(jp_eigenval), index)
C
C       The eigenvectors will be stored in the array b(jp_sol)
C       The eigenvalues and eigenvectors will be renumbered  
C                 using the permutation vector index
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: call to array_sol"
      endif
        call array_sol (i_cond, nval, nel, npt, neq, nnodes, 
     *   index, a(ip_table_nod), 
     *   a(ip_type_el), a(ip_eq), 
     *   b(jp_x), b(jp_eigenval), 
     *   b(jp_eigenval_tmp), b(jp_eigen_pol), b(jp_vp), b(jp_sol))


      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: array_sol returns call"
      endif
C
      if(debug .eq. 1) then
        write(ui,*) 'index = ', (index(i), i=1,nval)
      endif
      if(debug .eq. 1) then
        write(ui,*)
        write(ui,*) "lambda, 1/lambda = ", lambda, 1.0d0/lambda
        write(ui,*) "sqrt(shift)/(2*pi) = ", sqrt(shift) / (2.0d0 * pi)
        do i=1,nval
          write(ui,"(i4,2(g22.14),2(g18.10))") i, 
     *       b(jp_eigenval+i-1)
        enddo
      endif

      if (debug .eq. 1) then
          open (unit=343, file='Matrices/evals.txt',
     *         status='unknown')
          write(343,*) "sqrt(shift) = ", sqrt(shift)
          do i=1,nval
            write(343,"(i4,2(g22.14))") i,
     *       b(jp_eigenval+i-1)
          enddo
          close(343)
      endif

C    Save Original solution
      if (PrintSolution .eq. 1) then
        dir_name = "Output/Fields"
        call write_sol (nval, nel, nnodes, lambda,
     *       b(jp_eigenval), b(jp_sol), mesh_file, dir_name)
        call write_param (lambda, npt, nel, i_cond,
     *       nval, nvect, itermax, tol, shift, lx, ly, 
     *       mesh_file, n_conv, 
     *       dir_name)
cc  nb_typ_el, eps_eff, mesh_format, 
      tchar = "Output/FieldsPNG/All_plots_png_abs2_eE.geo"
      open (unit=34,file=tchar)
        do i=1,nval
          call gmsh_post_process (i, nval, nel, npt, 
     *       nnodes, a(ip_table_nod), a(ip_type_el), 
     *       b(jp_x), b(jp_eigenval), b(jp_sol),
     *       b(jp_rhs), a(ip_visite), gmsh_file_pos, dir_name)
        enddo 
      close (unit=34)
      endif

      Lambda_count = Lambda_count + 1
C
      call cpu_time(time_Lambda_end)
      write(ui,*) "py_calc_modes_AC: Time for Lambda loop:    ", 
     *                time_Lambda_end-time_Lambda_start
      write(ui,*)"-----------------------------------Time-Remaining==",
     *  (n_lambda+1-Lambda_count)
     *  *(time_Lambda_end-time_Lambda_start)

      enddo
C
C#########################  End Wavelength Loop  ########################
C
C
C
C#########################  End Calculations  ###########################
C
      call date_and_time ( end_date, end_time )
      call cpu_time(time2)
C
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
        write(26,*) "beta = ", beta
        write(26,*) "shift = ", shift
        write(26,*)
        write(26,*) "npt, nel, nnodes  = ", npt, nel, nnodes
        write(26,*) "neq, i_cond = ", neq, i_cond

        write(26,*) " lat_vecs:  = "
        write(26,"(2(f18.10))") lat_vecs
        write(26,*) "shift             = ", shift
        write(26,*) "mesh_file = ", mesh_file
        write(26,*) "gmsh_file = ", gmsh_file
        write(26,*) "log_file  = ", log_file
      close(26)

C
      write(ui,*) "   .      .      ."
      write(ui,*) "   .      .      ."
      write(ui,*) "   .      . (d=",d_in_nm,")"
      write(ui,*) "  and   we're   done!"
C
      stop
      end 
