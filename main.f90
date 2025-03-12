program main
   use mpmodule
   use bspline_mod
   implicit none

   integer :: d, nmin, nmax, n_remove, n, step, i_tmp
   type(mp_real) :: amin, amax, Z, C, kappa

   type(mp_real) :: zero, one

   zero = '0.d0'
   one = '1.d0'

   !-----------------------------------------------------------------!
   ! Define Important Variables
   d = 4 ! Order of Mathemathica + 1
   nmin = 20 ! Number minimal of B-Spline
   nmax = 1000 ! Number minimal of B-Spline
   step = 10 ! Number of Intermedier n
   n_remove = 2 ! Number of Bspline to remove at the start and the end
   Z = '2.d0'
   C = '137.0359895d0' ! check CODATA 1986
   kappa = '-1.d0'
   amin = '2.d-3'
   amax = '4.d1'
   !-----------------------------------------------------------------!

   do i_tmp = 1, step
      n = nint(nmin*(nmax/nmin)**(real(i_tmp - 1)/real(step - 1)))
      print *, 'Starting Process for n =', n, '& d =', d
      
      call get_eigen(d, n, n_remove, Z, kappa, C, amin, amax, .false., 45, 25)

   end do
   print *, 'All Process Done and writen to ./result'

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
