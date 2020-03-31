C Calculate the energy (not power) overlap integral of an EM mode with itself using
C numerical quadrature.  
C
      subroutine EM_mode_E_energy_int (nval, nel, npt,
     *  nnodes, table_nod, type_el, nb_typ_el, n_lst,
     *  x, soln_EM, overlap)
C
      implicit none
      integer*8 nval, nel, npt, nnodes
      integer*8 table_nod(nnodes,nel), nb_typ_el
      integer*8 type_el(nel)
      double precision x(2,npt)
      complex*16 soln_EM(3,nnodes+7,nval,nel)
      complex*16 n_lst(nb_typ_el), eps_lst(nb_typ_el)
      complex*16, dimension(nval) :: overlap

c     Local variables
      integer*8 nnodes_0
      parameter (nnodes_0 = 6)
      integer*8 nod_el_p(nnodes_0)
      complex*16 basis_overlap(nnodes_0)
      integer*8 i, j, j1, typ_e
      integer*8 iel, ival
      integer*8 itrial, i_eq
      integer*8 info_curved, n_curved, debug, ui
      double precision xel(2,nnodes_0)
      double precision phi2_list(6), grad2_mat0(2,6)
      double precision grad2_mat(2,6)
      double precision ZERO, ONE, r_tmp1
      parameter ( ZERO = 0.0D0, ONE = 1.0D0)
      complex*16 coeff_1
      complex*16 E, Estar
c
c     NQUAD: The number of quadrature points used in each element.
      integer*8 nquad, nquad_max, iq
      parameter (nquad_max = 25)
      double precision wq(nquad_max)
      double precision xq(nquad_max), yq(nquad_max)
      double precision xx(2), xx_g(2), ww, det
      double precision mat_B(2,2), mat_T(2,2), eps_0
C
C
Cf2py intent(in) nval, nel, npt,
Cf2py intent(in) nnodes, table_nod, type_el, nb_typ_el
Cf2py intent(in) x, soln_EM, n_lst
C
Cf2py depend(table_nod) nnodes, nel
Cf2py depend(x) npt
Cf2py depend(soln_EM) nnodes, nval, nel
Cf2py depend(n_lst) nb_typ_el
Cf2py depend(type_el) nel
C
Cf2py intent(out) overlap
C
C
CCCCCCCCCCCCCCCCCCCCC Start Program CCCCCCCCCCCCCCCCCCCCCCCC
C
      ui = 6
      debug = 0
      eps_0 = 8.854187817d-12
C
      if ( nnodes .ne. 6 ) then
        write(ui,*) "EM_mode_E_energy_int: problem nnodes = ", nnodes
        write(ui,*)"EM_mode_E_energy_int: nnodes should be equal to 14!"
        write(ui,*) "EM_mode_E_energy_int: Aborting..."
        stop
      endif
c
      call quad_triangle (nquad, nquad_max, wq, xq, yq)
      if (debug .eq. 1) then
        write(ui,*) "EM_mode_E_energy_int: nquad, nquad_max = ",
     *              nquad, nquad_max
      endif

c     Calculate permittivity
      do i = 1, int(nb_typ_el)
        eps_lst(i) = n_lst(i)**2
      enddo
C
      do i=1,nval
        overlap(i) = 0.0d0
      enddo
c
c      n_curved = 0
      do iel=1,nel
        typ_e = type_el(iel)
        do j=1,nnodes
          j1 = table_nod(j,iel)
          nod_el_p(j) = j1
          xel(1,j) = x(1,j1)
          xel(2,j) = x(2,j1)
        enddo
        call curved_elem_tri (nnodes, xel, info_curved, r_tmp1)
        if (info_curved .eq. 1) then
          n_curved = n_curved + 1
        endif
cccccccccc
        do i=1,nnodes
          basis_overlap(i) = 0.0d0
        enddo
cccccccccc
        do iq=1,nquad
          xx(1) = xq(iq)
          xx(2) = yq(iq)
          ww = wq(iq)
c         xx   = coordinate on the reference triangle
c         xx_g = coordinate on the actual triangle
c         We will also need the gradients of the P1 element
c          grad2_mat0 = gradient on the reference triangle (P2 element)
           call phi2_2d_mat(xx, phi2_list, grad2_mat0)
c
          if (info_curved .eq. 0) then
c           Rectilinear element
            call jacobian_p1_2d(xx, xel, nnodes,
     *               xx_g, det, mat_B, mat_T)
c            if (det .le. 0) then
            if (det .le. 0 .and. debug .eq. 2 .and. iq .eq. 1) then
              write(*,*) "   !!!"
              write(*,*) "EM_m_e_en_int: det <= 0: iel, det ", iel, det
            endif
          else
c           Isoparametric element
            call jacobian_p2_2d(xx, xel, nnodes, phi2_list,
     *               grad2_mat0, xx_g, det, mat_B, mat_T)
          endif
           if(abs(det) .lt. 1.0d-20) then
             write(*,*)
             write(*,*) "   ???"
             write(*,*) "EM_m_e_en_int: det = 0 : iel, det = ", iel, det
             write(*,*) "EM_m_e_en_int: Aborting..."
             stop
           endif
c          grad_i  = gradient on the actual triangle
c          grad_i  = Transpose(mat_T)*grad_i0
c          Calculation of the matrix-matrix product:
          call DGEMM('Transpose','N', 2, 6, 2, ONE, mat_T, 2,
     *      grad2_mat0, 2, ZERO, grad2_mat, 2)
          coeff_1 = ww * abs(det)
C Calculate overlap of basis functions at quadrature point, 
C which is a superposition of P2 polynomials for each function (field).
          do itrial=1,nnodes_0
            basis_overlap(itrial) = basis_overlap(itrial) + 
     *        coeff_1 * phi2_list(itrial) * phi2_list(itrial)
          enddo
        enddo
cccccccccc
C Having calculated overlap of basis functions on element
C now multiply by specific field values for modes of interest.
        do ival=1,nval
          do itrial=1,nnodes_0
            do i_eq=1,3
              Estar = conjg(soln_EM(i_eq,itrial,ival,iel))
              E = soln_EM(i_eq,itrial,ival,iel)
              overlap(ival) = overlap(ival) + 
     *          eps_lst(typ_e) * Estar * E * basis_overlap(itrial)
            enddo
          enddo
        enddo
cccccccccccc
C Loop over elements - end
cccccccccccc
      enddo
C Multiply through prefactor
      do i=1,nval
        overlap(i) = 2.0 * overlap(i) * eps_0
      enddo
C
CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
C
C
      end subroutine EM_mode_E_energy_int
