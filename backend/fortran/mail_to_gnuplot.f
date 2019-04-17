c*******************************************************
c
c     mail_to_gnuplot: covert the FEM mesh format to the gnuplot mesh format
c
c*******************************************************
c
c     nnodes : Number of nodes per element (10-node second order tetrahedron)
c
c*******************************************************
c
      subroutine mail_to_gnuplot (nel, npt, nnodes, type_el, 
     *  type_nod, table_nod, nb_typ_el, n_eff, 
     *  x, gmsh_file)

c
      implicit none
      integer*8 nel, npt, nnodes, nb_typ_el
      integer*8 type_el(nel), type_nod(npt)
      integer*8 table_nod(nnodes,nel)
      integer*8 max_typ_el
      parameter (max_typ_el=10)
      complex*16 x(2,npt), n_eff(max_typ_el)
      character gmsh_file*100
c
c     Local variables
      integer*8 nvertex
      parameter (nvertex = 3)
      double precision xel(2,nvertex+1)
      integer*8 i, j, k, j1, iel
      integer*8 namelen, namelen_print, gnuplot_info

      character gnuplot_file*100
      character gnuplot_file_print*100
      character gnuplot_file_png*100

c
ccccccccccccccccccccccccccccccccccccc
c

C       ! if gnuplot_info = 1: Write the  gnuplot command line file
      gnuplot_info = 1  


      namelen = len_trim(gmsh_file)
      gnuplot_file = gmsh_file(1:namelen-4)// "_gnuplot.txt"
C       ! 1+7: Remove "Normed/" from the file name
      gnuplot_file_png = gmsh_file(1+7:namelen-4)// ".png" 

      open (unit=26,file=gnuplot_file)
      do iel=1,nel
        do j=1,nvertex
          j1 = table_nod(j,iel)
          xel(1,j) = x(1,j1)
          xel(2,j) = x(2,j1)
        enddo
C         ! The nodes xel(:,j) form a closed hexagon
        k=nvertex+1  
        j=1
          j1 = table_nod(j,iel)
          xel(1,k) = x(1,j1)
          xel(2,k) = x(2,j1)
        do j=1,nvertex+1
          write(26,*) xel(1,j), xel(2,j)
        enddo
          write(26,*)
      enddo
      close(26)

      if(gnuplot_info == 1) then
        open (unit=26,file="Normed/gnuplot_commands.txt")
          write(26,*)
          write(26,*) " # You can type the following command line ",
     *    "to start gnuplot: gnuplot < gnuplot_info.txt"
          write(26,*)
          write(26,*) " # To plot the data in a file:"
          write(26,*)
          write(26,*) " set style data lines"
          write(26,*) " set size ratio -1"
          namelen = len_trim(gnuplot_file)
          gnuplot_file_print = " plot """//
C      ! Remove "Normed/" from the file name
     *    gnuplot_file(1+7:namelen)//"""" 
          namelen_print = len_trim(gnuplot_file)
          write(26,*) gnuplot_file_print(1:namelen_print)
          write(26,*)
          write(26,*) " # To create a PNG format of the gnuplot image:"
          write(26,*) 
          write(26,*) " set term png"
          namelen = len_trim(gnuplot_file_png)
          gnuplot_file_print = " set output """//
     *    gnuplot_file_png(1:namelen)//""""
          namelen_print = len_trim(gnuplot_file)
          write(26,*) gnuplot_file_print (1:namelen_print)
          write(26,*)
          namelen = len_trim(gnuplot_file)
          gnuplot_file_print = " plot """//
C      ! Remove "Normed/" from the file name
     *    gnuplot_file(1+7:namelen)//"""" 
          namelen_print = len_trim(gnuplot_file)
          write(26,*) gnuplot_file_print(1:namelen_print)
          write(26,*)
          write(26,*) " quit"
          write(26,*)
          write(26,*)
          write(26,*) " # Meaning of set size ratio -1:"
          write(26,*) " #    gnuplot tries to set the scales so that "
          write(26,*) " #    the unit has the same length on both "
          write(26,*) " #    the x and y axes"
        close(26)
      endif
c
ccccccccccccccccccccccccc
c


      return
      end
