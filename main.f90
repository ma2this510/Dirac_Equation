program main
   use mpmodule
   use bspline_mod
   implicit none

   integer :: d, n_remove, n
   type(mp_real) :: amin, amax, Z, C, kappa, zero, one

   type(mp_real), dimension(:, :), allocatable :: Gauche, Droite

   zero = '0.d0'
   one = '1.d0'

   !-----------------------------------------------------------------!
   ! Define Important Variables
   d = 4 ! Order of Mathemathica + 1
   n = 8 ! Number of Bspline
   n_remove = 2 ! Number of Bspline to remove at the start and the end
   Z = '2.d0'
   C = '137.0359895d0' ! check CODATA 1986
   kappa = '-1.d0'
   amin = '1.d-1'
   amax = '1.d1'
   !-----------------------------------------------------------------!

   allocate(Gauche(n - 2*n_remove, n - 2*n_remove))
   allocate(Droite(n - 2*n_remove, n - 2*n_remove))

   call matrixAB(d, n, n_remove, Z, kappa, C, amin, amax, Gauche, Droite, .true., 45, 25)

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
