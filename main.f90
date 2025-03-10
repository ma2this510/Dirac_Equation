program main
   use mpmodule
   use bspline_mod
   implicit none

   integer :: d, n, nprime, i_tmp, ierr
   type(mp_real) :: amin, amax, Z, C, kappa
   type(mp_real), dimension(:, :), allocatable :: A, B, vect
   type(mp_real), dimension(:), allocatable :: w, fv1, fv2

   !----------------------------------------------------------------
   ! Define Important Variables
   d = 9 ! Order of Mathemathica + 1
   n = 80 ! Number of B-Splines
   Z = mpreal(2.d0)
   C = mpreal(137.036d0) ! check CODATA 1986
   kappa = mpreal(-1.d0)
   amin = mpreal(0.002d0)
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

   allocate (w(2*nprime))
   allocate (fv1(2*nprime))
   allocate (fv2(2*nprime))
   allocate (vect(2*nprime, 2*nprime))

   call rsg(2*nprime, 2*nprime, A, B, w, 0, vect, fv1, fv2, ierr)

   open (2, file='eigenvalues.dat', status='replace')

   do i_tmp = 1, 2*nprime
      call mpwrite(2, 30, 10, w(i_tmp))
   end do
   
   close (2)

   deallocate (A, B, fv1, fv2, w)

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
