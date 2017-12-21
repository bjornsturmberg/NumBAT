

      subroutine array_material_EM (nel, npt, 
     *  nb_typ_el, n_index, 
     *  type_el, ls_material)

      implicit none
      integer*8 nnodes_0
      parameter (nnodes_0 = 6)
      integer*8 nel, npt, nb_typ_el
      integer*8 type_el(nel)
      complex*16 n_index(nb_typ_el)
      complex*16 ls_material(1,nnodes_0+7,nel)

c     Local variables
      integer*8 debug, k_typ
      integer*8 i, j, iel, inod
      complex*16 n_index_el

c
      debug = 0

C
C#####################  Loop over Elements  ##########################
C
      do iel=1,nel
        k_typ = type_el(iel)
        n_index_el = n_index(k_typ)
        do inod=1,nnodes_0+7
          ls_material(1,inod,iel) = n_index_el
        enddo
      enddo
C
C#########################  End Element Loop  ########################
C

      if (debug .eq. 1) then
      open (unit=63, file="Output/ls_material.txt",
     *         status="unknown")
        do iel=1,nel
          do inod=1,nnodes_0+7
            write(63,*) iel, inod,
     *         ls_material(1,inod,iel)
          enddo
          write(63,*)
        enddo
      close(63)
      endif
c
      return
      end
