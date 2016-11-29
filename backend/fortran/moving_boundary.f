c
cccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine moving_boundary (nval_EM, nval_AC, ival1,
     *    ival2, ival3, nel, npt, nnodes, table_nod, type_el, x,
     *    nb_typ_el, typ_select_in, typ_select_out, soln_EM, soln_AC,
     *    eps_lst, debug, overlap)
c
      implicit none
      integer*8 nel, npt, nnodes, nb_typ_el
      integer*8 type_el(nel)
      integer*8 table_nod(6,nel)
      double precision x(2,npt)
cc      complex*16 x(2,npt)
      integer*8 nval_EM, nval_AC, ival1, ival2, ival3
      integer*8 ival3s, ival2s, ival1s
      integer*8 typ_select_in, typ_select_out
      complex*16 soln_EM(3,nnodes,nval_EM,nel)
      complex*16 soln_AC(3,nnodes,nval_AC,nel)
      complex*16 eps_lst(nb_typ_el)
      complex*16 overlap(nval_EM, nval_EM, nval_AC)

c     Local variables
      integer alloc_stat, debug
C       integer*8, allocatable :: nb_visite(:) ! (npt)
C       integer*8, allocatable :: ls_edge_endpoint(:,:) ! (npt)
C       integer*8, allocatable :: edge_direction(:) ! (npt)
      integer*8 nb_visite(npt)
      integer*8 ls_edge_endpoint(2,npt)
      integer*8 edge_direction(npt)

      integer*8 iel, inod, typ_e
      integer*8 inod_1, inod_2, inod_3, ls_inod(3)
      integer*8 j, j_1, j_2, j_3, i, k
      integer*8 nb_edges, nb_interface_edges
      integer*8 edge_endpoints(2,3), opposite_node(3)
      double precision xy_1(2), xy_2(2), xy_3(2), ls_xy(2,3)
      double precision edge_vec(2), edge_perp(2), vec_0(2)
      double precision edge_length, r_tmp, zz
      double precision eps_0
      complex*16 ls_n_dot(3), ls_n_cross(3,3)
      complex*16 vec(3,3)!, z_integral
      complex*16 n_dot_d(2)
      complex*16 eps_a, eps_b, tmp1, tmp2
      double precision p2_p2_p2_1d(3,3,3)
      double precision version_number
      integer file_type, data_size
      integer physical_tag, elementary_tag
      integer element_type, number_of_tags
      integer number_of_string_tags
      integer number_of_real_tags
      integer number_of_integer_tags
C
C
Cf2py intent(in) nval_EM, nval_AC, ival1, ival2, ival3, nb_typ_el
Cf2py intent(in) nel, npt, nnodes, table_nod, debug
Cf2py intent(in) type_el, x, soln_EM, soln_AC
Cf2py intent(in) typ_select_in, typ_select_out, eps_lst, debug
C
Cf2py depend(table_nod) nnodes, nel
Cf2py depend(type_el) npt
Cf2py depend(x) npt
Cf2py depend(soln_EM) nnodes, nval_EM, nel
Cf2py depend(soln_AC) nnodes, nval_AC, nel
Cf2py depend(eps_lst) nb_typ_el
C
Cf2py intent(out) overlap
C
ccccccccccccccccccccccccccccccccccccc
c
c     typ_select_in: Only the elements iel with type_el(iel)=typ_select_in will be analysed
c     When nb_visite(j) is not zero: nb_visite(j) indicates the number of element the edge j belongs
c
ccccccccccccccccccccccccccccccccccccc
c
C     alloc_stat = 0

C     allocate(nb_visite(npt), ls_edge_endpoint(2,npt),
C    *     edge_direction(npt), STAT=alloc_stat)
C     if (alloc_stat /= 0) then
C       write(*,*)
C       write(*,*) "edge_orientation: ",
C    *     "The allocation is unsuccessful"
C       write(*,*) "alloc_stat = ", alloc_stat
C       write(*,*) "Not enough memory for the array nb_visite"
C       write(*,*) "npt = ", npt
C       write(*,*) "edge_orientation: Aborting..."
C       stop
C     endif

c     Initialisation
      do inod=1,npt
        nb_visite(inod) = 0
        ls_edge_endpoint(1,inod) = 0
        ls_edge_endpoint(2,inod) = 0
        edge_direction(inod) = 0
      enddo
c
ccccccccccccccccccccccccccccccccccccc
c
      edge_endpoints(1,1) = 1
      edge_endpoints(2,1) = 2
      edge_endpoints(1,2) = 2
      edge_endpoints(2,2) = 3
      edge_endpoints(1,3) = 3
      edge_endpoints(2,3) = 1
c
c     opposite_node(i): Node which is opposite to the edge i
C     i = 1 is inod = 4 etc
      opposite_node(1) = 3
      opposite_node(2) = 1
      opposite_node(3) = 2
c
ccccccccccccccccccccccccccccccccccccc
      do i=1,nval_EM
        do j=1,nval_EM
          do k=1,nval_AC
            overlap(i,j,k) = 0.0d0
          enddo
        enddo
      enddo
cccccccccccc
C
      eps_0 = 8.854187817d-12
C
      do iel=1,nel
        typ_e = type_el(iel)
        if(typ_e == typ_select_in) then
          do inod=4,6  ! Scan the edges
            j = table_nod(inod,iel)
            nb_visite(j) = nb_visite(j) + 1  ! Will indicate the number of
          enddo
        endif
      enddo
      nb_edges = 0
      nb_interface_edges = 0
      do inod=1,npt
        if (nb_visite(inod) >= 1) then
          nb_edges = nb_edges + 1
        endif
        if (nb_visite(inod) == 1) then
          nb_interface_edges = nb_interface_edges + 1
        endif
      enddo
      if (debug .eq. 1) then
        write(*,*)
        write(*,*) "edge_orientation: npt, nel = ", npt, nel
        write(*,*) "edge_orientation: nb_edges = ", nb_edges
        write(*,*) "nb_interface_edges = ", nb_interface_edges
      endif
c     Outward pointing normal vector to the interface edges
      do iel=1,nel
        typ_e = type_el(iel)
        if(typ_e == typ_select_in) then
          do inod=4,6  ! Scan the edges
            j = table_nod(inod,iel)
            if (nb_visite(j) == 1) then
              inod_1 = edge_endpoints(1,inod-3)
              inod_2 = edge_endpoints(2,inod-3)
              ls_edge_endpoint(1,j) = table_nod(inod_1,iel)
              ls_edge_endpoint(2,j) = table_nod(inod_2,iel)
              xy_1(1) = x(1,table_nod(inod_1,iel))
              xy_1(2) = x(2,table_nod(inod_1,iel))
              xy_2(1) = x(1,table_nod(inod_2,iel))
              xy_2(2) = x(2,table_nod(inod_2,iel))
c             edge_vec: vector parallel to the edge
              edge_vec(1) = xy_2(1) - xy_1(1)
              edge_vec(2) = xy_2(2) - xy_1(2)
c             Normalisation of edge_vec
              r_tmp = sqrt(edge_vec(1)**2+edge_vec(2)**2)
              edge_vec(1) = edge_vec(1) / r_tmp
              edge_vec(2) = edge_vec(2) / r_tmp
c             edge_vec: vector perpendicular to the edge (rotation of edge_vec by -pi/2)
              edge_perp(1) = edge_vec(2)
              edge_perp(2) = -edge_vec(1)
c             Node opposite to the edge inod
              inod_3 = opposite_node(inod-3)
              xy_3(1) = x(1,table_nod(inod_3,iel))
              xy_3(2) = x(2,table_nod(inod_3,iel))
              vec_0(1) = xy_3(1) - xy_1(1)
              vec_0(2) = xy_3(2) - xy_1(2)
c             Scalar product of edge_perp and vec_0:
              r_tmp = edge_perp(1)*vec_0(1)+edge_perp(2)*vec_0(2)
c             if r_tmp < 0: then edge_perp is oriented in the outward direction
              if( r_tmp < 0) then
                edge_direction(j) = 1
              elseif( r_tmp > 0) then
                edge_direction(j) = -1
              else
                write(*,*) "edge_orientation: illegal:"
                write(*,*) "edge_perp is perpendicular to vec_0"
                write(*,*) "edge_orientation: Aborting..."
                stop
              endif
            endif
          enddo
        endif
      enddo
c
ccccccccccccccccccccccccccccccccccccc
c
c     Numerical integration
      do iel=1,nel
        typ_e = type_el(iel)
        if(typ_e == typ_select_in) then
          eps_a = eps_lst(1)!typ_e)
          eps_b = 1.0d0!eps_lst(2)!typ_select_out)
C           write(*,*) eps_a, eps_b
          do inod=4,6  ! Scan the edges
            j = table_nod(inod,iel)
            xy_3(1) = x(1,j)
            xy_3(2) = x(2,j)
            if (ls_edge_endpoint(1,j) .ne. 0) then
C               write(*,*) "an edge"
              inod_1 = ls_edge_endpoint(1,j)
              inod_2 = ls_edge_endpoint(2,j)
              xy_1(1) = x(1,inod_1)
              xy_1(2) = x(2,inod_1)
              xy_2(1) = x(1,inod_2)
              xy_2(2) = x(2,inod_2)
c             List of the nodes coordinates
              ls_xy(1,1) = xy_1(1) ! x-coord. of node 1
              ls_xy(2,1) = xy_1(2) ! y-coord. of node 1
              ls_xy(1,2) = xy_2(1) ! x-coord. of node 2
              ls_xy(2,2) = xy_2(2) ! y-coord. of node 2
              ls_xy(1,3) = xy_3(1) ! x-coord. of mid-edge node
              ls_xy(2,3) = xy_3(2) ! y-coord. of mid-edge node
c
              edge_vec(1) = ls_xy(1,2) - ls_xy(1,1)
              edge_vec(2) = ls_xy(2,2) - ls_xy(2,1)
c             Normalisation of edge_vec
              r_tmp = sqrt(edge_vec(1)**2+edge_vec(2)**2)
              edge_vec(1) = -1*edge_direction(j)*edge_vec(1) / r_tmp
              edge_vec(2) = -1*edge_direction(j)*edge_vec(2) / r_tmp
c             edge_vec: vector perpendicular to the edge (rotation of edge_vec by -pi/2)
              edge_perp(1) = -1*edge_vec(2)
              edge_perp(2) = edge_vec(1)
c
              r_tmp = (ls_xy(1,2) - ls_xy(1,1))**2
     *              + (ls_xy(2,2) - ls_xy(2,1))**2
              edge_length = sqrt(r_tmp)
              call mat_p2_p2_p2_1d (p2_p2_p2_1d, edge_length)
c             Identification number of the two end-points and mid-edge point
              ls_inod(1) = edge_endpoints(1,inod-3)
              ls_inod(2) = edge_endpoints(2,inod-3)
              ls_inod(3) = inod
c
C If only want overlap of one given combination of EM modes and AC mode.
              if (ival1 .ge. 0 .and. ival2 .ge. 0 .and. 
     *            ival3 .ge. 0) then
c             Nodes of the edge
              do j_1=1,3
c               (x,y,z)-components of the electric field
                vec(1,1) = soln_EM(1,ls_inod(j_1),ival1,iel)
                vec(2,1) = soln_EM(2,ls_inod(j_1),ival1,iel)
                vec(3,1) = soln_EM(3,ls_inod(j_1),ival1,iel)
c               ls_n_dot(1): Normal component of vec(:,1)
                ls_n_dot(1) = vec(1,1) * edge_perp(1)
     *              + vec(2,1) * edge_perp(2)
                ls_n_cross(1,1) = vec(3,1) * edge_perp(2)
                ls_n_cross(2,1) = -1*vec(3,1) * edge_perp(1)
                ls_n_cross(3,1) = vec(2,1) * edge_perp(1)
     *              - vec(1,1) * edge_perp(2)
                do j_2=1,3
c                 (x,y,z)-components of the electric field
                  vec(1,2)=soln_EM(1,ls_inod(j_2),ival2,iel)
                  vec(2,2)=soln_EM(2,ls_inod(j_2),ival2,iel)
                  vec(3,2)=soln_EM(3,ls_inod(j_2),ival2,iel)
c                 ls_n_dot(2): Normal component of vec(:,2)
                  ls_n_dot(2) = vec(1,2) * edge_perp(1)
     *                + vec(2,2) * edge_perp(2)
                  ls_n_cross(1,2) = vec(3,2) * edge_perp(2)
                  ls_n_cross(2,2) = -1*vec(3,2) * edge_perp(1)
                  ls_n_cross(3,2) = vec(2,2) * edge_perp(1)
     *                - vec(1,2) * edge_perp(2)
                  do j_3=1,3
c                   (x,y,z)-components of the acoustic field
                    vec(1,3) = soln_AC(1,ls_inod(j_3),ival3,iel)
                    vec(2,3) = soln_AC(2,ls_inod(j_3),ival3,iel)
                    vec(3,3) = soln_AC(3,ls_inod(j_3),ival3,iel)
c                   ls_n_dot(3): scalar product of vec(:,3) and normal vector edge_perp
                    ls_n_dot(3) = vec(1,3) * edge_perp(1)
     *                  + vec(2,3) * edge_perp(2)
                    tmp1 = (eps_a - eps_b)*eps_0
                    tmp1 = tmp1*((ls_n_cross(1,1))*ls_n_cross(1,2)
     *                    + (ls_n_cross(2,1))*ls_n_cross(2,2)
     *                    + (ls_n_cross(3,1))*ls_n_cross(3,2))
                    n_dot_d(1) = eps_0*eps_a * ls_n_dot(1)
                    n_dot_d(2) = eps_0*eps_a * ls_n_dot(2)
                    tmp2 = (1.0d0/eps_b - 1.0d0/eps_a)*(1.0d0/eps_0)
                    tmp2 = tmp2*(n_dot_d(1))*n_dot_d(2)
                    r_tmp = p2_p2_p2_1d(j_1, j_2, j_3)
              overlap(ival1,ival2,ival3) = overlap(ival1,ival2,ival3) +
     *           r_tmp*conjg(ls_n_dot(3))*(tmp1 - tmp2)
                  enddo
                enddo
              enddo
C
C If want overlap of given EM mode 1 and 2 and all AC modes.
              else if (ival1 .ge. 0 .and. ival2 .ge. 0 .and.
     *                 ival3 .eq. -1) then
c             Nodes of the edge
              do j_1=1,3
c               (x,y,z)-components of the electric field
                vec(1,1) = soln_EM(1,ls_inod(j_1),ival1,iel)
                vec(2,1) = soln_EM(2,ls_inod(j_1),ival1,iel)
                vec(3,1) = soln_EM(3,ls_inod(j_1),ival1,iel)
c               ls_n_dot(1): Normal component of vec(:,1)
                ls_n_dot(1) = vec(1,1) * edge_perp(1)
     *              + vec(2,1) * edge_perp(2)
c               ls_n_cross(1): Tangential component of vec(:,1)
C                   ls_n_cross(1) = vec(1,1) * edge_perp(2)
C      *                - vec(2,1) * edge_perp(1)
                ls_n_cross(1,1) = vec(3,1) * edge_perp(2)
                ls_n_cross(2,1) = -1*vec(3,1) * edge_perp(1)
                ls_n_cross(3,1) = vec(2,1) * edge_perp(1)
     *              - vec(1,1) * edge_perp(2)
                do j_2=1,3
c                 (x,y,z)-components of the electric field
                  vec(1,2)=soln_EM(1,ls_inod(j_2),ival2,iel)
                  vec(2,2)=soln_EM(2,ls_inod(j_2),ival2,iel)
                  vec(3,2)=soln_EM(3,ls_inod(j_2),ival2,iel)
c                 ls_n_dot(2): Normal component of vec(:,2)
                  ls_n_dot(2) = vec(1,2) * edge_perp(1)
     *                + vec(2,2) * edge_perp(2)
c                 ls_n_cross(2): Tangential component of vec(:,2)
C                   ls_n_cross(2) = vec(1,2) * edge_perp(2)
C      *                - vec(2,2) * edge_perp(1)
                  ls_n_cross(1,2) = vec(3,2) * edge_perp(2)
                  ls_n_cross(2,2) = -1*vec(3,2) * edge_perp(1)
                  ls_n_cross(3,2) = vec(2,2) * edge_perp(1)
     *                - vec(1,2) * edge_perp(2)

                  do ival3s = 1,nval_AC
                  do j_3=1,3
c                   (x,y,z)-components of the acoustic field
                    vec(1,3) = soln_AC(1,ls_inod(j_3),ival3s,iel)
                    vec(2,3) = soln_AC(2,ls_inod(j_3),ival3s,iel)
                    vec(3,3) = soln_AC(3,ls_inod(j_3),ival3s,iel)
c                   ls_n_dot(3): scalar product of vec(:,3) and normal vector edge_perp
                    ls_n_dot(3) = vec(1,3) * edge_perp(1)
     *                  + vec(2,3) * edge_perp(2)
                    tmp1 = (eps_a - eps_b)*eps_0
                    tmp1 = tmp1*((ls_n_cross(1,1))*ls_n_cross(1,2)
     *                    + (ls_n_cross(2,1))*ls_n_cross(2,2)
     *                    + (ls_n_cross(3,1))*ls_n_cross(3,2))
                    n_dot_d(1) = eps_0*eps_a * ls_n_dot(1)
                    n_dot_d(2) = eps_0*eps_a * ls_n_dot(2)
                    tmp2 = (1.0d0/eps_b - 1.0d0/eps_a)*(1.0d0/eps_0)
                    tmp2 = tmp2*(n_dot_d(1))*n_dot_d(2)
                    r_tmp = p2_p2_p2_1d(j_1, j_2, j_3)
             overlap(ival1,ival2,ival3s) = overlap(ival1,ival2,ival3s) +
     *           r_tmp*conjg(ls_n_dot(3))*(tmp1 - tmp2)
C                     write(*,*) "tmp1", tmp1
C                     write(*,*) "tmp2", tmp2
C                     write(*,*) "r_tmp", r_tmp
C                     write(*,*) "ls_n_dot", ls_n_dot(3)
C                     write(*,*) "fie", r_tmp*ls_n_dot(3)*(tmp1 - tmp2)
                  enddo
                  enddo
                enddo
              enddo
              endif
C cc
            endif
          enddo
        endif
      enddo
C       if (debug .eq. 1) then
C         write(*,*)
C         write(*,*) "nb_interface_edges: z_integral = ", z_integral
C         write(*,*)
C       endif
c
ccccccccccccccccccccccccccccccccccccc
c
C       open (unit=26,file="Output/edge_data.txt")
C       write(26,*)
C       write(26,*) "typ_select_in = ", typ_select_in
C       write(26,*) "npt, nel = ", npt, nel
C       write(26,*) "nb_edges = ", nb_edges
C       write(26,*) "nb_interface_edges = ", nb_interface_edges
C       write(26,*) "nb_interface_edges: z_integral = ", z_integral
C       j = 0
C       do inod=1,npt
C         if (ls_edge_endpoint(1,inod) .ne. 0) then
C           j = j + 1
C           write(26,*) j, inod, ls_edge_endpoint(1,inod),
C      *              ls_edge_endpoint(2,inod),
C      *              edge_direction(inod)
C         endif
C       enddo
C       close(26)
c
ccccccccccccccccccccccccccccccccccccc
c
      debug = 1
      if (debug .eq. 1) then
        version_number = 2.2
        file_type = 0  ! An integer equal to 0 in the ASCII file format
        data_size = 8 ! An integer equal to the size of the floating point numbers used in the file
        open (unit=27,file="../Output/edge_data.msh")
        write(27,'(a11)') "$MeshFormat"
        write(27,'((f4.1,1x,I1,1x,I1,1x))') version_number,
     *            file_type, data_size
        write(27,'(a14)') "$EndMeshFormat"
        write(27,'(a6)') "$Nodes"
        write(27,'(I0.1)') nb_interface_edges
        zz = 0.0d0
        j = 0
        do inod=1,npt
          if (ls_edge_endpoint(1,inod) .ne. 0) then
              xy_1(1) = 100*x(1,inod)
              xy_1(2) = 100*x(2,inod)
            j = j + 1
            write(27,*) j, xy_1(1), xy_1(2), zz
          endif
        enddo
        write(27,'(a9)') "$EndNodes"
        write(27,'(a9)') "$Elements"
        write(27,'(I0.1)') nb_interface_edges
        element_type = 15  ! 1-node point
        number_of_tags = 2
        j = 0
        do inod=1,npt
          if (ls_edge_endpoint(1,inod) .ne. 0) then
            j = j + 1
          physical_tag = j
          elementary_tag = j
          write(27,'(100(I0.1,2x))') j, element_type,
     *      number_of_tags, physical_tag, elementary_tag,
     *      j
          endif
        enddo
        write(27,'(a12)') "$EndElements"
        number_of_string_tags = 1
        number_of_real_tags = 1
        number_of_integer_tags = 3
        write(27,'(a9)') "$NodeData"
        write(27,*) number_of_string_tags
        write(27,*) " ""View of tangential vector"" "
        write(27,*) number_of_real_tags
        write(27,*) 0.0
        write(27,*) number_of_integer_tags
        write(27,*) 0 ! the time step (0; time steps always start at 0)
        write(27,*) 3 ! 3-component (vector) field
        write(27,*) nb_interface_edges ! Number of associated nodal values
c        node-number value
        zz = 0.0d0
        j = 0
        do inod=1,npt
          if (ls_edge_endpoint(1,inod) .ne. 0) then
            inod_1 = ls_edge_endpoint(1,inod)
            inod_2 = ls_edge_endpoint(2,inod)
            xy_1(1) = x(1,inod_1)
            xy_1(2) = x(2,inod_1)
            xy_2(1) = x(1,inod_2)
            xy_2(2) = x(2,inod_2)
            edge_vec(1) = xy_2(1) - xy_1(1)
            edge_vec(2) = xy_2(2) - xy_1(2)
c            Normalisation of edge_vec
            r_tmp = sqrt(edge_vec(1)**2+edge_vec(2)**2)
            edge_vec(1) = -1*edge_direction(inod)*edge_vec(1) / r_tmp
            edge_vec(2) = -1*edge_direction(inod)*edge_vec(2) / r_tmp
            j = j + 1
            write(27,*) j, edge_vec(1), edge_vec(2), zz
          endif
        enddo
        write(27,'(a12)') "$EndNodeData"
c
ccccccccccccccccccccccccccccccccccccc
c
        write(27,'(a9)') "$NodeData"
        write(27,*) number_of_string_tags
        write(27,*) " ""View of the normal vector"" "
        write(27,*) number_of_real_tags
        write(27,*) 0.0
        write(27,*) number_of_integer_tags
        write(27,*) 0 ! the time step (0; time steps always start at 0)
        write(27,*) 3 ! 3-component (vector) field
        write(27,*) nb_interface_edges ! Number of associated nodal values
c        node-number value
        zz = 0.0d0
        j = 0
        do inod=1,npt
          if (ls_edge_endpoint(1,inod) .ne. 0) then
            inod_1 = ls_edge_endpoint(1,inod)
            inod_2 = ls_edge_endpoint(2,inod)
            xy_1(1) = x(1,inod_1)
            xy_1(2) = x(2,inod_1)
            xy_2(1) = x(1,inod_2)
            xy_2(2) = x(2,inod_2)
            edge_vec(1) = xy_2(1) - xy_1(1)
            edge_vec(2) = xy_2(2) - xy_1(2)
c            Normalisation of edge_vec
            r_tmp = sqrt(edge_vec(1)**2+edge_vec(2)**2)
            edge_vec(1) = -1*edge_direction(inod)*edge_vec(1) / r_tmp
            edge_vec(2) = -1*edge_direction(inod)*edge_vec(2) / r_tmp
c            edge_vec: vector perpendicular to the edge (rotation of edge_vec by -pi/2)
            edge_perp(1) = -edge_vec(2)
            edge_perp(2) = edge_vec(1)
C             edge_perp(1) = edge_perp(1) * edge_direction(inod)
C             edge_perp(2) = edge_perp(2) * edge_direction(inod)
            j = j + 1
            write(27,*) j, edge_perp(1), edge_perp(2), zz
          endif
        enddo
        write(27,'(a12)') "$EndNodeData"
        close(27)
      endif
c
ccccccccccccccccccccccccccccccccccccc
c
C       return
C       end
      end subroutine moving_boundary


