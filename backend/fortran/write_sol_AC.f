c

      subroutine write_sol_AC (nval, nel, nnodes,
     *     lambda, beta, sol, mesh_file, dir_name)
c
      implicit none
      integer*8 nval, nel, nnodes
      double precision lambda
      complex*16 sol(3,nnodes,nval,nel)
      complex*16 beta(nval)
      character mesh_file*1000
c     Local variables
      integer*8 nnodes_0
      parameter (nnodes_0 = 6)

      integer*8 i, iel, ival
      integer*8 namelen, namelength, namelen2

      character*11 ivalue, jvalue
      character dir_name*1000
      character*1000 tchar1, tchar2
c
      namelen = len_trim(mesh_file)
      namelength = len_trim(dir_name)

      open (unit=63, file="Matrices/beta.txt",
     *         status="unknown")
        do ival=1,nval
          write(63,12) ival, beta(ival)
        enddo
      close(63)

      do ival=1,nval
        jvalue = ivalue(ival)
        namelen2 = len_trim(jvalue)
        tchar1 = dir_name(1:namelength)//"/r_mode_"
     *     //jvalue(1:namelen2)//".txt"
        tchar2 = dir_name(1:namelength)//"/i_mode_"
     *     //jvalue(1:namelen2)//".txt"
        open(63, file=tchar1,
     *     form="formatted")
        open(64, file=tchar2,
     *     form="formatted")
        write(63,*) mesh_file(1:namelen)
        write(64,*) mesh_file(1:namelen)
        write(63,"((g25.17),2(I7))") lambda, ival
        write(64,"((g25.17),2(I7))") lambda, ival
        write(63,"(2(g25.17))") beta(ival)
        write(64,"(2(g25.17))") beta(ival)
        do iel=1,nel
          write(63,12) iel, (dble(sol(1,i,ival,iel)),i=1,nnodes)
          write(63,12) iel, (dble(sol(2,i,ival,iel)),i=1,nnodes)
          write(63,12) iel, (dble(sol(3,i,ival,iel)),i=1,nnodes)
c
          write(64,12) iel, (imag(sol(1,i,ival,iel)),i=1,nnodes)
          write(64,12) iel, (imag(sol(2,i,ival,iel)),i=1,nnodes)
          write(64,12) iel, (imag(sol(3,i,ival,iel)),i=1,nnodes)
        enddo
        close(63)
        close(64)
      enddo
12    format(I7,13(g25.12) )
c
      return
      end
