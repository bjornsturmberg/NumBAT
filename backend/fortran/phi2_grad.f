c
ccccccccccccccccccccccccccccccc
c
c       Note that if grad_i_0 is the gradient on the reference triangle, 
c       then the gradient on the actual triangle is:
c       grad_i  = Transpose(mat_T)*grad_i0
c
ccccccccccccccccccccccccccccccc
c
      subroutine phi2_grad(inode, nnodes, c, vec_grad)

cc  dt, mat_jac
c
      implicit none
      integer*8 inode, nnodes
      double precision vec_grad(2,nnodes)
      double precision c(2,2)
cc      double precision mat_jac(2,2), dt 

c     Local variables
      integer*8 nnodes_0
      parameter (nnodes_0 = 6)
      double precision xel_0(2,nnodes_0)
      double precision phi_xi, phi_yi
      double precision phi0_xi, phi0_yi
      double precision x, y
      integer*8 i

c      integer*8 ui, debug
c
c     Coordinates (x,y)= xel_0(1..2,inode) of the P2 Lagrange interpolaion nodes
      xel_0(1,1) = 0
      xel_0(2,1) = 0
      xel_0(1,2) = 1
      xel_0(2,2) = 0
      xel_0(1,3) = 0
      xel_0(2,3) = 1
      xel_0(1,4) = 0.5d0
      xel_0(2,4) = 0
      xel_0(1,5) = 0.5d0
      xel_0(2,5) = 0.5d0
      xel_0(1,6) = 0
      xel_0(2,6) = 0.5d0

cc      if(abs(dt) .gt. 1.0d-12) then
cc        c(1,1) = mat_jac(2,2)/dt
cc        c(1,2) = -mat_jac(1,2)/dt
cc        c(2,1) = -mat_jac(2,1)/dt
cc        c(2,2) = mat_jac(1,1)/dt
cc      else
cc        write(*,*) 'phi2_grad: ATTENTION: determinant is zero !'
cc        write(*,*) 'dt = ', dt
cc        write(*,*) 'Aborting...'
cc        stop
cc      endif
      x = xel_0(1,inode)
      y = xel_0(2,inode)
c
      i = 1
        phi0_xi = 4.0d0*(x+y) - 3.0d0  ! x-derivative over the reference triangle
        phi0_yi = 4.0d0*(x+y) - 3.0d0 ! y-derivative over the reference triangle
        phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi ! x-derivative over the current triangle
        phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi ! y-derivative over the current triangle
        vec_grad(1,i) = phi_xi
        vec_grad(2,i) = phi_yi
      i = 2
        phi0_xi = 4.0d0*x - 1.0d0
        phi0_yi = 0.0d0
        phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
        phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
        vec_grad(1,i) = phi_xi
        vec_grad(2,i) = phi_yi
      i = 3
        phi0_xi = 0.0d0
        phi0_yi = 4.0d0*y - 1.0d0
        phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
        phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
        vec_grad(1,i) = phi_xi
        vec_grad(2,i) = phi_yi
      i = 4
        phi0_xi = 4.0d0*(1.0d0 - 2.0d0*x - y)
        phi0_yi = -4.0d0*x
        phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
        phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
        vec_grad(1,i) = phi_xi
        vec_grad(2,i) = phi_yi
      i = 5
        phi0_xi = 4.0d0*y
        phi0_yi = 4.0d0*x
        phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
        phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
        vec_grad(1,i) = phi_xi
        vec_grad(2,i) = phi_yi
      i = 6
        phi0_xi = -4.0d0*y
        phi0_yi = 4.0d0*(1.0d0 - x - 2.0d0*y)
        phi_xi =c(1,1)*phi0_xi+c(1,2)*phi0_yi
        phi_yi =c(2,1)*phi0_xi+c(2,2)*phi0_yi
        vec_grad(1,i) = phi_xi
        vec_grad(2,i) = phi_yi
ccccccc
c
      return
      end
