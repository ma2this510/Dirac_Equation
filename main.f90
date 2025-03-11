program main
   use mpmodule
   use bspline_mod
   implicit none

   integer :: d, n, nprime, i_tmp, ierr, method, n_remove
   type(mp_real) :: amin, amax, Z, C, kappa, clt
   type(mp_real), dimension(:, :), allocatable :: A, B, vect
   type(mp_real), dimension(:), allocatable :: w, fv1, fv2

   type(mp_real) :: zero, one

   zero = '0.d0'
   one = '1.d0'

   !-----------------------------------------------------------------!
   ! Define Important Variables
   d = 7 ! Order of Mathemathica + 1
   n = 200 ! Number of B-Splines
   n_remove = 1 ! Number of Bspline to remove at the start and the end
   Z = '2.d0'
   C = '137.0359895d0' ! check CODATA 1986
   kappa = '-1.d0'
   amin = '2.d-3'
   amax = '4.d1'
   method = 0
   clt = '1.025d-1' ! Clustering factor useful only if method = 1 and become h if method == 2
   !-----------------------------------------------------------------!

   nprime = n - 2*n_remove
   allocate (A(2*nprime, 2*nprime))
   allocate (B(2*nprime, 2*nprime))

   !----------------------------------------------------------------
   ! Generate A and B matrices
   call matrixAB(d, n, n_remove, Z, kappa, C, amin, amax, method, clt, A, B, .true., 40, 20)

   !----------------------------------------------------------------
   ! Get Eigenvalues

   print *, 'Calculating Eigenvalues'

   allocate (w(2*nprime))
   allocate (fv1(2*nprime))
   allocate (fv2(2*nprime))
   allocate (vect(2*nprime, 2*nprime))

   call rsg(2*nprime, 2*nprime, A, B, w, 0, vect, fv1, fv2, ierr)

   print *, "Error Status : ", ierr

   open (2, file='eigenvalues.dat', status='replace')

   do i_tmp = 1, 2*nprime
      call mpwrite(2, 45, 25, w(2*nprime - i_tmp + 1))
   end do
   
   close (2)

   deallocate (A, B, fv1, fv2, w, vect)

   print *, 'Done'

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
