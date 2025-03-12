program main
   use mpmodule
   use bspline_mod
   implicit none

   integer :: d, n_remove, n
   type(mp_real) :: amin, amax, Z, C, kappa

   type(mp_real) :: zero, one

   zero = '0.d0'
   one = '1.d0'

   !-----------------------------------------------------------------!
   ! Define Important Variables
   d = 7 ! Order of Mathemathica + 1
   n = 60 ! Number of Bspline
   n_remove = 2 ! Number of Bspline to remove at the start and the end
   Z = '2.d0'
   C = '137.0359895d0' ! check CODATA 1986
   kappa = '1.d0'
   amin = '2.d-3'
   amax = '4.d1'
   !-----------------------------------------------------------------!

   call get_eigen(d, n, n_remove, Z, kappa, C, amin, amax, .true., 45, 25)

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
