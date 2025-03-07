program main
   use mpmodule
   use bspline_mod
   implicit none

   integer :: d, n, iter, nprime, maxit
   type(mp_real) :: amin, amax, Z, C, kappa, eps, w
   type(mp_real), dimension(:, :), allocatable :: A, B
   type(mp_real), dimension(:), allocatable :: work, wa

   !----------------------------------------------------------------
   ! Define Important Variables
   d = 4 ! Order of Mathemathica + 1
   n = 6 ! Number of B-Splines
   Z = mpreal(0.d0)
   C = 137
   kappa = mpreal(-1.d0)
   amin = mpreal(0.1d0)
   amax = mpreal(10.d0)

   nprime = n - 2
   allocate (A(2*nprime, 2*nprime))
   allocate (B(2*nprime, 2*nprime))

   !----------------------------------------------------------------
   ! Generate A and B matrices
   call matrixAB(d, n, Z, kappa, C, amin, amax, A, B, .true.)

   !----------------------------------------------------------------
   ! Get Eigenvalues

   print *, 'Calculating Eigenvalues'

   open(1, file='eigenvalues.dat', status='replace')

   allocate (wa(3*2*nprime + 1))
   allocate (work(2*nprime))
   work(1) = mpreal(1.d0)
   work = mpreal(0.d0)
   w = mpreal(10.d4)
   eps = mpreal(0.d-28)
   maxit = 10

   print *, 'Working'
   do iter = 1, 3
      call invsg(A, B, 2*nprime, w, work, eps, 0, maxit, wa)
      call mpwrite(1, 20, 10, w)
   end do
   close (1)

   deallocate (A, B, work, wa)

end program main

function epsilonn(alpha)
   !> @brief Calculate the machine epsilon
   !> @param alpha The value to calculate the machine epsilon
   USE mpmodule
   implicit type(mp_real) (a - h, o - z)

   ten = '10.d0'
   epsilonn = ten**(-mpipl)

   return
end function epsilonn
