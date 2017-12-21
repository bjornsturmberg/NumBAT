

      subroutine array_material_AC (nel, npt, 
     *  nb_typ_el, rho, c_tensor, 
     *  type_el, ls_material)

      implicit none
      integer*8 nnodes_0
      parameter (nnodes_0 = 6)
      integer*8 nel, npt, nb_typ_el
      integer*8 type_el(nel)
      complex*16 rho(nb_typ_el), c_tensor(6,6,nb_typ_el)
      complex*16 ls_material(4,nnodes_0,nel)

c     Local variables
      integer*8 debug, k_typ
      integer*8 i, j, iel, inod
      complex*16 rho_0, c_11, c_12, c_44
      complex*16 c_tensor_el(6,6), rho_el

c
      debug = 1

C
C#####################  Loop over Elements  ##########################
C
      do iel=1,nel
        k_typ = type_el(iel)
        do j=1,6
          do i=1,6
            c_tensor_el(i,j) = c_tensor(i,j,k_typ)
          enddo
        enddo
        c_11 = c_tensor(1,1,k_typ)
        c_44 = c_tensor(4,4,k_typ)
        c_12 = c_tensor(1,2,k_typ)
        rho_el = rho(k_typ)
        do inod=1,nnodes_0
          ls_material(1,inod,iel) = rho_el
          ls_material(2,inod,iel) = c_11
          ls_material(3,inod,iel) = c_44
          ls_material(4,inod,iel) = c_12
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
      return
      end
