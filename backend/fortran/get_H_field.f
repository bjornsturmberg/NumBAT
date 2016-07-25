
c     For a waveguide mode: compute the H-field from the E-field
C     This program used Maxwell's equation to compute the magnetic field H from the electric field E
C
      subroutine get_H_field (nnodes, k_0, beta1, mat_T, 
     *    E_field_el, H_field_el)
C
      implicit none
C
      integer*8 nnodes
      double precision k_0, mat_T(2,2)
      complex*16 beta1
      complex*16 E_field_el(3,nnodes)
      complex*16 H_field_el(3,nnodes)

c     Local variables

      integer*8 nnodes_0
      parameter (nnodes_0 = 6)
      double precision vec_grad(2,nnodes_0), speed_c, omega
      integer*8 j, inod, jnod
      complex*16 ii, z_tmp1, z_tmp2
      complex*16 Maxwell_coeff
C
CCCCCCCCCCCCCCCCCCCCCCCCC
C
c  ii = sqrt(-1)
      ii = dcmplx(0.0d0, 1.0d0)
      speed_c = 299792458
C
c      By applying the Maxwell's equations to the E-field of a waveguide mode, we get:
c      H_x = [-beta*E_y + D(E_z,y)] * Coefficient
c      H_y = [ beta*E_x - D(E_z,x)] * Coefficient
c      H_z = [ D(E_y,x) - D(E_x,y)] * Coefficient

      do inod=1,nnodes
c       The components (H_x,H_y,H_z) the mode ival
        do j=1,3
          H_field_el(j,inod) = 0
        enddo
      enddo
      do inod=1,nnodes
        z_tmp1 = -beta1 * E_field_el(2,inod)
        z_tmp2 =  beta1 * E_field_el(1,inod)
        H_field_el(1,inod) = H_field_el(1,inod) + z_tmp1
        H_field_el(2,inod) = H_field_el(2,inod) + z_tmp2
      enddo
      do inod=1,nnodes
c       vec_grad: contains the gradients of all 6 basis polynomials at the node inod
        call phi2_grad(inod, nnodes, mat_T, vec_grad)
        do jnod=1,nnodes
          z_tmp1 =  vec_grad(2,jnod) * E_field_el(3,jnod)
          z_tmp2 = -vec_grad(1,jnod) * E_field_el(3,jnod)
          H_field_el(1,inod) = H_field_el(1,inod) + z_tmp1
          H_field_el(2,inod) = H_field_el(2,inod) + z_tmp2
        enddo
        do jnod=1,nnodes
          z_tmp1 = -vec_grad(2,jnod) * E_field_el(1,jnod)
          z_tmp2 =  vec_grad(1,jnod) * E_field_el(2,jnod)
          H_field_el(3,inod) = H_field_el(3,inod) + z_tmp1+z_tmp2
        enddo
      enddo
c     The curl of the E-field must be multiplied by a coefficient in order to get the H-field
c     For example: Maxwell_coeff = 1/ (i * k0 * mu) 
      omega = k_0 * speed_c
      Maxwell_coeff = 1.0d0 / (ii * omega)
C       Maxwell_coeff = 1.0d0 / (ii * k_0)
      do inod=1,nnodes
        do j=1,3
          H_field_el(j,inod) = H_field_el(j,inod) * Maxwell_coeff
        enddo
      enddo


C
CCCCCCCCCCCCCCCCCCCCCCCCC
C
      return
      end 

