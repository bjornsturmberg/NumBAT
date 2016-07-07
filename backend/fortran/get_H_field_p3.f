
c     For a waveguide mode: compute the H-field from the E-field
C     This program used Maxwell's equation to compute the magnetic field H from the electric field E
C
      subroutine get_H_field_p3 (nnodes_P2, k_0, beta1, mat_T, 
     *    E_field_el, Ez_field_el_P3, H_field_el)
C
      implicit none
C
      integer*8 nnodes_P2
      double precision k_0, mat_T(2,2)
      complex*16 beta1
      complex*16 E_field_el(3,nnodes_P2)
      complex*16 Ez_field_el_P3(10)  !  P3 Ez-field
      complex*16 H_field_el(3,nnodes_P2)

c     Local variables

      integer*8 nnodes_P2_0, nnodes_P3_0
      parameter (nnodes_P2_0 = 6, nnodes_P3_0 = 10)
      double precision vec_grad_P2(2,nnodes_P2_0)
      double precision vec_grad_P3(2,nnodes_P3_0)
      integer*8 j, inod, jnod
      complex*16 ii, z_tmp1, z_tmp2
      complex*16 Maxwell_coeff
C
CCCCCCCCCCCCCCCCCCCCCCCCC
C
c  ii = sqrt(-1)
      ii = dcmplx(0.0d0, 1.0d0)
C
c      By applying the Maxwell's equations to the E-field of a waveguide mode, we get:
c      H_x = [-beta*E_y + D(E_z,y)] * Coefficient
c      H_y = [ beta*E_x - D(E_z,x)] * Coefficient
c      H_z = [ D(E_y,x) - D(E_x,y)] * Coefficient

      do inod=1,nnodes_P2
c       The components (H_x,H_y,H_z) the mode ival
        do j=1,3
          H_field_el(j,inod) = 0
        enddo
      enddo
      do inod=1,nnodes_P2
        z_tmp1 = -beta1 * E_field_el(2,inod) * ii
        z_tmp2 =  beta1 * E_field_el(1,inod) * ii
        H_field_el(1,inod) = H_field_el(1,inod) + z_tmp1
        H_field_el(2,inod) = H_field_el(2,inod) + z_tmp2
      enddo
      do inod=1,nnodes_P2
c       vec_grad_p2: contains the gradients of all 6 basis polynomials at the node inod
        call phi2_grad(inod, nnodes_P2, mat_T, vec_grad_P2)
        call phi3_grad_p2(inod, nnodes_P3_0, mat_T, vec_grad_P3)
        do jnod=1,nnodes_P3_0
          z_tmp1 =  vec_grad_P3(2,jnod) * Ez_field_el_P3(jnod)
          z_tmp2 = -vec_grad_P3(1,jnod) * Ez_field_el_P3(jnod)
          H_field_el(1,inod) = H_field_el(1,inod) + z_tmp1
          H_field_el(2,inod) = H_field_el(2,inod) + z_tmp2
        enddo
        do jnod=1,nnodes_P2
          z_tmp1 = -vec_grad_P2(2,jnod) * E_field_el(1,jnod)
          z_tmp2 =  vec_grad_P2(1,jnod) * E_field_el(2,jnod)
          H_field_el(3,inod) = H_field_el(3,inod) + z_tmp1+z_tmp2
        enddo
      enddo
c     The curl of the E-field must be multiplied by a coefficient in order to get the H-field
c     For example: Maxwell_coeff = 1/ (i * k0 * mu) 
      Maxwell_coeff = 1.0d0 / (ii * k_0)
      do inod=1,nnodes_P2
        do j=1,3
          H_field_el(j,inod) = H_field_el(j,inod) * Maxwell_coeff
        enddo
      enddo


C
CCCCCCCCCCCCCCCCCCCCCCCCC
C
      return
      end 

