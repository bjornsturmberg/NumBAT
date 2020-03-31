      subroutine array_material_AC (npt, nel, nb_typ_el, type_el,
     *  rho, c_tensor, p_tensor, eta_tensor, ls_material)

      implicit none
      integer*8 nnodes_0
      parameter (nnodes_0 = 6)
      integer*8 nel, npt, nb_typ_el
      integer*8 type_el(nel)
      complex*16 rho(nb_typ_el), c_tensor(6,6,nb_typ_el)
      complex*16 p_tensor(3,3,3,3,nb_typ_el)
      complex*16 eta_tensor(3,3,3,3,nb_typ_el)
      complex*16 ls_material(10,nnodes_0,nel)

c     Local variables
      integer*8 debug, k_typ
      integer*8 i, iel, inod
      complex*16 rho_el, c_11, c_12, c_44
      complex*16 p_11, p_12, p_44, eta_11, eta_12, eta_44

Cf2py intent(in) npt, nel, nb_typ_el, type_el
Cf2py intent(in) rho, c_tensor, p_tensor
Cf2py intent(in) eta_tensor

Cf2py depend(type_el) nel
Cf2py depend(rho) nb_typ_el
Cf2py depend(c_tensor) nb_typ_el
Cf2py depend(p_tensor) nb_typ_el
Cf2py depend(eta_tensor) nb_typ_el

Cf2py intent(out) ls_material

c
      debug = 0

C
C#####################  Loop over Elements  ##########################
C
      do iel=1,nel
        k_typ = type_el(iel)
        rho_el = rho(k_typ)
        c_11 = c_tensor(1,1,k_typ)
        c_12 = c_tensor(1,2,k_typ)
        c_44 = c_tensor(4,4,k_typ)
        p_11 = p_tensor(1,1,1,1,k_typ)
        p_12 = p_tensor(1,1,2,2,k_typ)
        p_44 = p_tensor(2,3,2,3,k_typ)
        eta_11 = eta_tensor(1,1,1,1,k_typ)
        eta_12 = eta_tensor(1,1,2,2,k_typ)
        eta_44 = eta_tensor(2,3,2,3,k_typ)
        do inod=1,nnodes_0
          ls_material(1,inod,iel) = rho_el
          ls_material(2,inod,iel) = c_11
          ls_material(3,inod,iel) = c_12
          ls_material(4,inod,iel) = c_44
          ls_material(5,inod,iel) = p_11
          ls_material(6,inod,iel) = p_12
          ls_material(7,inod,iel) = p_44
          ls_material(8,inod,iel) = eta_11
          ls_material(9,inod,iel) = eta_12
          ls_material(10,inod,iel) = eta_44
        enddo
      enddo
C
C#########################  End Element Loop  ########################
C

      if (debug .eq. 1) then
      open (unit=63, file="Output/ls_material.txt",
     *         status="unknown")
        do iel=1,nel
          do inod=1,nnodes_0
            write(63,*) iel, inod,
     *         (ls_material(i,inod,iel),i=1,4)
          enddo
          write(63,*)
        enddo
      close(63)
      endif
c
      end subroutine array_material_AC