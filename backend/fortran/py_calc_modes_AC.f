      subroutine calc_AC_modes(
c     Explicit inputs
     *    lambda, beta_in, nval,
     *    debug, mesh_file, npt, nel,
     *    nb_typ_el,  c_tensor, rho, d_in_n, shift,
     *    i_cond, itermax,
     *    plot_modes,
     *    cmplx_max, real_max, int_max,
c     Outputs
     *    beta1, sol1)

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
      complex*16 beta_in ! Propagation constant
      integer*8 int_max, cmplx_max, int_used, cmplx_used
      integer*8 real_max, real_used, plot_modes !, n_64
      integer :: alloc_stat=0
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

      integer*8 nb_typ_el
      complex*16 pp(nb_typ_el),  qq(nb_typ_el)
      complex*16 c_tensor(6,6,nb_typ_el)
c     rho: density
      complex*16 rho(nb_typ_el)

      integer*8 i, j, k, ip!, Lambda_count
      integer*8 nnodes, ui, debug!, PrintSolution
      integer*8 nel, npt, i_cond, neq

C     ! Number of nodes per element
      parameter(nnodes=6)
      integer*8 type_nod(npt), type_el(nel), table_nod(nnodes, nel)

      double precision pi!, theta, phi
      double precision lambda!, lambda_1, lambda_2, d_lambda
      double precision d_freq, freq, lat_vecs(2,2)
      double precision lx, ly, d_in_n

      complex*16 shift
      integer*8  i_base

C  Variable used by valpr
      integer*8 ltrav, n_conv
      double precision ls_data(10)
      complex*16 z_beta, z_tmp, z_tmp0
C      integer*8 index(1000)
      integer*8, dimension(:), allocatable :: index
c     variable used by UMFPACK
      double precision control (20), info_umf (90)
      integer*8 numeric, symbolic, status, sys, filenum

      double precision time1, time2
      character*(8) start_date, end_date
      character*(10) start_time, end_time

C  Variable used by valpr
      integer*8 nval, nvect, itermax
      double precision tol

C  Names and Controls
      character mesh_file*100, gmsh_file*100, log_file*100
      character gmsh_file_pos*100, dir_name*100
      character*100 tchar

c     new breed of variables to prise out of a, b and c
      double precision x_arr(2,npt)
      complex*16, target :: sol1(3,nnodes+7,nval,nel)
      complex*16, target :: beta1(nval)
C      complex*16 mode_pol(4,nval)

Cf2py intent(in) lambda, beta_in, nval
Cf2py intent(in) debug, mesh_file, npt, nel
Cf2py intent(in) d_in_n, shift
Cf2py intent(in) i_cond, itermax
Cf2py intent(in) plot_modes, c_tensor, rho
Cf2py intent(in) cmplx_max, real_max, int_max
Cf2py intent(in) nb_typ_el

C  Note: the dependent variables must be listed AFTER the 
C  independent variables that they depend on in the function call!
Cf2py depend(c_tensor) nb_typ_el
Cf2py depend(rho) nb_typ_el

Cf2py intent(out) beta1
Cf2py intent(out) sol1

C
CCCCCCCCCCCCCCCCCCCC  Start Program - get parameters  CCCCCCCCCCCCCCCCCC
C
C     Set parameter for the super-vectors of integer and real numbers
C
      ui = 6     !ui = Unite dImpression
C      nnodes = 6 ! Number of nodes per element
      pi = 3.141592653589793d0

      nvect = 2*nval + nval/2 +3
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
      call cpu_time(time1)  ! initial time  in unit = sec.
      call date_and_time ( start_date, start_time )
C
      tol = 0.0 ! ARPACK accuracy (0.0 for machine precision)
C      lx=1.0 ! Diameter of unit cell. Default, lx = 1.0.
C      ly=1.0 ! NOTE: currently requires ly=lx, ie rectangular unit cell
C     ToDo: sort out what's going on here - in EM lx=ly=1
      lx = d_in_n
      ly = d_in_n
C      call get_param(n_lambda, lambda_1, lambda_2,
C     *      npt, nel, i_cond, nval, nvect, itermax, tol,
C     *      shift_0, lx, ly, mesh_file, gmsh_file,
C     *      gmsh_file_pos, log_file, d_in_n, debug)
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
      ip_type_nod = 1
      ip_type_el = ip_type_nod + npt
      ip_table_nod = ip_type_el + nel ! pointer to FEM connectivity table

      ip_visite = ip_table_nod + nnodes*nel
      ip_eq = ip_visite + npt

      jp_x = 1
C
C      lx = d_in_n
C      ly = d_in_n
C      call geometry (nel, npt, nnodes, nb_typ_el,
C     *     lx, ly, a(ip_type_nod), a(ip_type_el), a(ip_table_nod),
C     *     b(jp_x), mesh_file)
      call geometry (nel, npt, nnodes, nb_typ_el,
     *     lx, ly, type_nod, type_el, table_nod,
     *     x_arr, mesh_file)

C      call lattice_vec (npt, b(jp_x), lat_vecs)
      call lattice_vec (npt, x_arr, lat_vecs, debug)

      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: npt, nel = ", npt, nel
      endif
C
C      if (PrintSolution .ge. 1) then
C  Export the mesh to gmsh format
C        call mail_to_gmsh (nel, npt, nnodes, a(ip_type_el),
C     *    a(ip_type_nod), a(ip_table_nod),
C     *    nb_typ_el, b(jp_x), gmsh_file)
C      endif

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

C      ToDo: is this the same as in EM prob?
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
c       Propagation constant: beta = k_z
cc        beta = 28.9d0 * d_in_n
cc        beta = 28.9d0 / d_in_n

C        beta = 1.49175d7
C
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: Asmbly: call to asmbly"
      endif
C     Assemble the coefficient matrix K and M of the finite element equations
      call asmbly_AC (i_base, nel, npt, neq, nnodes,
     *  shift, beta_in, nb_typ_el, rho, c_tensor,
     *  a(ip_table_nod), type_el, a(ip_eq),
     *  b(jp_x), nonz, a(ip_row), a(ip_col_ptr),
     *  c(kp_mat1_re), c(kp_mat1_im), b(jp_mat2), a(ip_work))
C
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: call to valpr"
      endif
      call valpr_64_AC (i_base, nvect, nval, neq, itermax, ltrav,
     *  tol, nonz, a(ip_row), a(ip_col_ptr), c(kp_mat1_re),
     *  c(kp_mat1_im), b(jp_mat2), b(jp_vect1), b(jp_vect2),
     *  b(jp_workd), b(jp_resid), b(jp_vschur), beta1,
     *  b(jp_trav), b(jp_vp), c(kp_rhs_re), c(kp_rhs_im),
     *  c(kp_lhs_re), c(kp_lhs_im), n_conv, ls_data,
     *  numeric, filenum, status, control, info_umf, debug)

      if (n_conv .ne. nval) then
         write(ui,*) "py_calc_modes_AC: convergence problem in valpr_64"
         write(ui,*) "py_calc_modes_AC: n_conv != nval : ",
     *    n_conv, nval
         write(ui,*) "py_calc_modes_AC: Aborting..."
cc         stop
      endif

C
      do i=1,nval
        z_tmp0 = beta1(i)
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
        beta1(i) = z_beta
      enddo
c
      call z_indexx_AC (nval, beta1, index)
C
C       The eigenvectors will be stored in the array b(jp_sol)
C       The eigenvalues and eigenvectors will be renumbered
C                 using the permutation vector index
      if (debug .eq. 1) then
        write(ui,*) "py_calc_modes_AC: call to array_sol"
      endif
        call array_sol_AC (i_cond, nval, nel, npt, neq, nnodes,
     *   index, table_nod,
     *   type_el, a(ip_eq),
     *   x_arr, b(jp_eigenval),
     *   b(jp_eigenval_tmp), beta1, b(jp_vp), sol1)


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
     *       beta1(i)
        enddo
      endif

C      if (debug .eq. 1) then
C          open (unit=343, file='Matrices/evals.txt',
C     *         status='unknown')
C          write(343,*) "sqrt(shift) = ", sqrt(shift)
C          do i=1,nval
C            write(343,"(i4,2(g22.14))") i,
C     *       beta1(i)
C          enddo
C          close(343)
C      endif

C    Save Original solution
      if (plot_modes .eq. 1) then
        dir_name = "AC_fields"
C        call write_sol_AC (nval, nel, nnodes, lambda,
C    *       b(jp_eigenval), b(jp_sol), mesh_file, dir_name)
C        call write_param (lambda, npt, nel, i_cond,
C    *       nval, nvect, itermax, tol, shift, lx, ly,
C    *       mesh_file, n_conv, dir_name)
        tchar = "AC_fields/All_plots_png_abs2_eE.geo"
        open (unit=34,file=tchar)
          do i=1,nval
            call gmsh_post_process (i, nval, nel, npt,
     *         nnodes, a(ip_table_nod), a(ip_type_el),
     *         b(jp_x), beta1, b(jp_sol),
     *         b(jp_rhs), a(ip_visite), gmsh_file_pos, dir_name)
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
      write(ui,*) "   .      . (d=",d_in_n,")"
      write(ui,*) "  and   we're   done!"
C
      deallocate(a,b,c,index)

      end subroutine calc_AC_modes
