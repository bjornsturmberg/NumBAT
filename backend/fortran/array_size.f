
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine array_size (nel, nval,
     *     int_size, cmplx_size, real_size)

      implicit none
      integer*8 nel
      integer*8 int_size, cmplx_size, real_size
c      integer*8 


c     Local variables
      integer*8 nnodes, npt, npt_p3
      integer*8 int_max, cmplx_max
      integer*8 real_max, int_max_32
      integer*8 nval, nvect, ordre_ls
      integer*8 nonz, nonz_max, max_row_len
      integer*8 neq, neq_PW
      integer*8 n_edge, n_face, n_ddl, n_ddl_max
      integer*8 ltrav
c
C  Declare the pointers of the integer super-vector
      integer*8 ip_type_nod, ip_type_el, ip_table_nod
      integer*8 ip_table_E, ip_table_N_E_F, ip_visite
      integer*8 ip_type_N_E_F, ip_eq
      integer*8 ip_period_N, ip_nperiod_N
      integer*8 ip_period_N_E_F, ip_nperiod_N_E_F
      integer*8 ip_index_pw_inv
C  Declare the pointers of the real super-vector
      integer*8 jp_x, jp_x_N_E_F, jp_rhs
c      integer*8 jp_matD, jp_matL, jp_matU
c      integer*8 jp_matD2, jp_matL2, jp_matU2
      integer*8 jp_vect1, jp_vect2, jp_workd, jp_resid, jp_vschur
      integer*8 jp_eigenval_tmp, jp_trav, jp_vp, jp_eigen_pol
      integer*8 jp_overlap_L, jp_overlap_J, jp_overlap_J_dagger
      integer*8 jp_flux
      integer*8 jp_overlap_K, jp_X_mat
      integer*8 jp_sol, jp_sol1, jp_sol2, jp_sol1b
      integer*8 jp_sol1_H, jp_sol1b_H
      integer*8 jp_eigenval, jp_eigenval1, jp_eigenval2
      integer*8 jp_T, jp_R, jp_T12, jp_R12, jp_T21, jp_R21
      integer*8 jp_T_Lambda, jp_R_Lambda
      integer*8 jp_X_mat_b
c     Declare the pointers of the real super-vector
      integer*8 kp_rhs_re, kp_rhs_im, kp_lhs_re, kp_lhs_im
      integer*8 kp_mat1_re, kp_mat1_im

c     Declare the pointers of for sparse matrix storage
      integer*8 ip_col_ptr, ip_row
      integer*8 jp_mat2
      integer*8 ip_work, ip_work_sort, ip_work_sort2

C      integer*4 ip_col_ptr_32, ip_row_32

Cf2py intent(in) nel, nval,

Cf2py intent(out) int_size, cmplx_size, real_size

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      nnodes = 6
      ordre_ls = 0
c
c     For most of the FEM meshes I have used, I have observed that:
c     npt is approximately equal to nel * 2.1
c
      npt = nel * 3
      n_edge = (npt + nel) / 2
      npt_p3 = npt + n_edge + nel
      nvect = 2*nval + nval/2 +3

c     Euler's polyhedron formula (no holes):
c     V - E + F = 2
c     V, E, and F are respectively the numbers of vertices (corners), edges and faces (triangles) 
c
c     Since V - E = 2 - nel and V + E = npt, we have E = n_edge = (npt+nel-2)/2
c

      n_ddl = n_edge + nel + npt_p3
      neq = 3 * (n_edge + nel) + npt_p3
      neq_PW = (2*ordre_ls+1)**2
c     For most of the FEM meshes I have used, I have observed that:
c     nonz = 34.25 * neq
      nonz = 40 * neq
      nonz_max = nonz

c     I have observed that: max_row_len < 200
      max_row_len = 200
cccccc
      ip_type_nod = 1
      ip_type_el = ip_type_nod + npt
C       ! pointer to FEM connectivity table
      ip_table_nod = ip_type_el + nel 
      ip_table_N_E_F = ip_table_nod + nnodes*nel

      n_ddl_max = npt + nel
      ip_visite =  ip_table_N_E_F  + 14*nel 
      ip_table_E = ip_visite + n_ddl_max

      ip_type_N_E_F = ip_table_E + 4*n_edge


        ip_period_N = ip_type_N_E_F + 2*n_ddl
        ip_nperiod_N = ip_period_N + npt
        ip_period_N_E_F = ip_nperiod_N + npt
        ip_nperiod_N_E_F = ip_period_N_E_F + n_ddl
        ip_eq = ip_nperiod_N_E_F + n_ddl

      ip_index_pw_inv = ip_eq + 3*n_ddl
      ip_col_ptr = ip_index_pw_inv + neq_PW 
      ip_row = ip_col_ptr + neq + 1

      ip_work = ip_row + nonz
      ip_work_sort = ip_work + 3*n_ddl
      ip_work_sort2 = ip_work_sort + max_row_len

      int_size = ip_work_sort2 + max_row_len

cccccc
      jp_x = 1
      jp_x_N_E_F = jp_x + 2*npt

      jp_rhs = jp_x_N_E_F + 3*n_ddl
c     jp_rhs will also be used (in gmsh_post_process) to store a solution
      jp_mat2 = jp_rhs + max(neq, 3*npt)

      jp_vect1 = jp_mat2 + nonz
      jp_vect2 = jp_vect1 + neq
      jp_workd = jp_vect2 + neq
      jp_resid = jp_workd + 3*neq

      jp_sol1 = jp_resid + neq
      jp_sol1b = jp_sol1 + 3*(nnodes+7)*nval*nel
      jp_sol2 = jp_sol1b + 3*(nnodes+7)*nval*nel
      jp_sol1_H = jp_sol2 + 3*(nnodes+7)*nval*nel
      jp_sol1b_H = jp_sol1_H + 3*nnodes*nval*nel
      jp_eigenval1 = jp_sol1b_H + 3*nnodes*nval*nel
      jp_eigenval2 = jp_eigenval1 + nval + 1
C       ! Eigenvectors
      jp_vschur = jp_eigenval2 + nval + 1     
      jp_eigenval = jp_vschur + neq*nvect
      jp_eigen_pol = jp_eigenval + nval + 1
      jp_eigenval_tmp = jp_eigen_pol + nval*4
      jp_trav = jp_eigenval_tmp + nval + 1

      ltrav = 3*nvect*(nvect+2)
      jp_vp = jp_trav + ltrav
      jp_overlap_L = jp_vp + neq*nval
      jp_flux = jp_overlap_L + nval*nval
      jp_overlap_J = jp_flux + nval
      jp_overlap_J_dagger = jp_overlap_J + 2*neq_PW*nval
      jp_overlap_K = jp_overlap_J_dagger + 2*neq_PW*nval
      jp_X_mat = jp_overlap_K + 2*neq_PW*nval
      jp_T12 = jp_X_mat + 2*neq_PW*2*neq_PW  
      jp_T21 = jp_T12 + 2*neq_PW*nval
      jp_R12 = jp_T21 + 2*neq_PW*nval
      jp_R21 = jp_R12 + 4*neq_PW*neq_PW
      jp_T = jp_R21 + nval*nval
      jp_R = jp_T + 2*neq_PW*2*neq_PW
      jp_T_Lambda = jp_R + 2*neq_PW*2*neq_PW
      jp_R_Lambda = jp_T_Lambda + 2*neq_PW*2*neq_PW
 
cc      if (substrate .eq. 1) then
      jp_X_mat_b = jp_R_Lambda + 2*neq_PW*2*neq_PW
      cmplx_size = jp_X_mat_b + 2*neq_PW*2*neq_PW  
cc      else
cc        cmplx_size = jp_R_Lambda + 2*neq_PW*2*neq_PW
cc      endif

cccccc

      kp_rhs_re = 1
      kp_rhs_im = kp_rhs_re + neq
      kp_lhs_re = kp_rhs_im + neq
      kp_lhs_im = kp_lhs_re + neq
      kp_mat1_re = kp_lhs_im + neq
      kp_mat1_im = kp_mat1_re + nonz
      real_size = kp_mat1_im + nonz

C cccccc
C c     SOME 32 bit integers for UMFPACK AND ARPACK
C       ip_col_ptr_32 = 1
C       ip_row_32 = ip_col_ptr_32 + neq + 1
C       int_size_32 = ip_row_32 + nonz

C cccccc
C       write(*,*) "array_size:"
C       write(*,*) "int_size = ", int_size
C       write(*,*) "cmplx_size = ", cmplx_size
C       write(*,*) "real_size = ", real_size
C       write(*,*) "int_size_32 = ", int_size_32

C       open (unit=26, file="Output/array_size.txt", status='unknown')
C         write(26,*) "int_size = ", int_size
C         write(26,*) "cmplx_size = ", cmplx_size
C         write(26,*) "real_size = ", real_size
C         write(26,*) "int_size_32 = ", int_size_32
C         write(26,*)
C         write(26,*) "nel = ", nel
C         write(26,*) "nval = ", nval
C         write(26,*) "nvect = ", nvect
C         write(26,*) "ordre_ls = ", ordre_ls
C         write(26,*)
C         write(26,*) "npt = ", npt
C         write(26,*) "n_edge = ", n_edge
C         write(26,*) "npt_p3 = ", npt_p3
C         write(26,*) "n_ddl = ", n_ddl
C         write(26,*) "neq = ", neq
C         write(26,*) "neq_PW = ", neq_PW
C         write(26,*) "nonz_max = ", nonz_max
C         write(26,*) "nonz = ", nonz
C         write(26,*) "max_row_len = ", max_row_len
C       close(26)
c
      end subroutine array_size

