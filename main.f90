program main
   use mpmodule
   use bspline_mod
   implicit none

   integer :: d, n, iter, nprime, maxit, rn, i_tmp, j_tmp
   type(mp_real) :: amin, amax, Z, C, kappa, eps, w, ga, ext, kappa_abs, zqvdn
   type(mp_real), dimension(:, :), allocatable :: A, B
   type(mp_real), dimension(:), allocatable :: work, wa, avec, bvec

   !----------------------------------------------------------------
   ! Define Important Variables
   d = 7 ! Order of Mathemathica + 1
   n = 40 ! Number of B-Splines
   Z = mpreal(1.d0)
   C = mpreal(137.d0)
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
   ! Generate Exact Solution

   print *, 'Calculating Exact Solution'

   rn = mpreal(1.d0)
   kappa_abs = abs(kappa)
   ga = (kappa_abs**mpreal(2.d0) - (Z/C)**mpreal(2.d0))**mpreal(0.5d0)
   ext = (c**mpreal(2.d0))*(1 + (Z/(C*(rn - kappa_abs + ga)))**mpreal(2.d0))**mpreal(-0.5d0) - (c**mpreal(2.d0))

   zqvdn = mpreal(0.99d0)
   w = ext*zqvdn

   !----------------------------------------------------------------
   ! Get Eigenvalues

   print *, 'Calculating Eigenvalues'

   open (2, file='eigenvalues.dat', status='replace')

   allocate (avec(2*nprime*(2*nprime + 1)/2))
   allocate (bvec(2*nprime*(2*nprime + 1)/2))

   avec = mpreal(0.d0)
   bvec = mpreal(0.d0)

   do i_tmp = 1, 2*nprime
      do j_tmp = 1, i_tmp
         avec((i_tmp*(i_tmp - 1))/2 + j_tmp) = A(i_tmp, j_tmp)
         bvec((i_tmp*(i_tmp - 1))/2 + j_tmp) = B(i_tmp, j_tmp)
      end do
   end do

   allocate (wa(3*2*nprime + 1))
   allocate (work(2*nprime))
   work = mpreal(0.d0)
   work(1) = mpreal(1.d0)
   eps = mpreal(0.d-28)
   maxit = 10

   print *, 'Working'
   do iter = 1, 2*nprime
      call invsg(avec, bvec, 2*nprime, w, work, eps, 0, maxit, wa)
      call mpwrite(2, 30, 10, w)
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
