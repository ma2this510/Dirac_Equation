program main
   use bspline_mod
   implicit none

   integer :: d, n, i_tmp, j_tmp
   integer, dimension(:), allocatable :: i_range
   real, dimension(:), allocatable :: knot
   real, dimension(:, :, :), allocatable :: bspline

   real, dimension(:, :), allocatable :: ovrlp_mat, deriv_mat, potential_mat, aprime_mat
   real :: amin, amax, Z, C, kappa, result

   !----------------------------------------------------------------
   ! Format for printing and logging
   open(1, file="log.txt", status="replace")
10 format(250f15.8)
   !----------------------------------------------------------------
   ! Define Important Variables
   d = 10 ! Order of Mathemathica + 1
   n = 200 ! Number of B-Splines
   Z = 1.0
   C = 137
   kappa = 1.0
   amin = 0.001
   amax = 40.0
   !----------------------------------------------------------------
   ! Generate Knots

   allocate (knot(n + d ))

   knot = exp_knot(size(knot), amax, amin)
   !----------------------------------------------------------------
   ! Generate Bsplines

   allocate (i_range(n))
   i_range = [(i_tmp, i_tmp=1, n)]

   allocate (bspline(i_range(n), size(knot), d))

   do i_tmp = 1, n
      call init_bspine(d, i_range(i_tmp), knot, bspline(i_tmp, :, :), .false.)
   end do

   print *, "Number of BSplines: ", n
   write (1,'(a,i4)') "Number of BSplines: ", n
   print *, "Order of BSplines: ", d, " and number of knots: ", size(knot)
   write (1,'(a,i4)') "Order of BSplines: ", d, " and number of knots: ", size(knot)
   print *, "Knots: "
   write (1,'(a)') "Knots: "
   write (*, 10) knot
   write (1, 10) knot

   !----------------------------------------------------------------
   ! Generate Overlap Matrix (C)

   print *, "----------------------------------------------------------------"
   write (1,'(a)') "----------------------------------------------------------------"
   print *, "BSplines Overlaps : "
   write (1,'(a)') "BSplines Overlaps : "

   allocate (ovrlp_mat(n, n))

   do i_tmp = 1, n
      do j_tmp = 1, n
         result = 0.0
         call int_overlp(bspline(i_tmp, :, :), bspline(j_tmp, :, :), knot, result)
         ovrlp_mat(i_tmp, j_tmp) = result
      end do
      write (*, 10) ovrlp_mat(i_tmp, :)
      write (1, 10) ovrlp_mat(i_tmp, :)
   end do

   !----------------------------------------------------------------
   ! Generate Derivative Matrix (D)

   print *, "----------------------------------------------------------------"
   write (1,'(a)') "----------------------------------------------------------------"
   print *, "BSplines Derivatives : "
   write (1,'(a)') "BSplines Derivatives : "

   allocate (deriv_mat(n, n))

   do i_tmp = 1, n
      do j_tmp = 1, n
         result = 0.0
         call int_deriv(bspline(i_tmp, :, :), bspline(j_tmp, :, :), knot, result)
         deriv_mat(i_tmp, j_tmp) = result
      end do
      write (*, 10) deriv_mat(i_tmp, :)
      write (1, 10) deriv_mat(i_tmp, :)
   end do

   !----------------------------------------------------------------
   ! Generate Potential Matrix (V)

   print *, "----------------------------------------------------------------"
   write (1,'(a)') "----------------------------------------------------------------"
   print *, "BSplines Potential : "
   write (1,'(a)') "BSplines Potential : "

   allocate (potential_mat(n, n))

   do i_tmp = 1, n
      do j_tmp = 1, n
         result = 0.0
         call int_potential(bspline(i_tmp, :, :), bspline(j_tmp, :, :), Z, knot, result)
         potential_mat(i_tmp, j_tmp) = result
      end do
      write (*, 10) potential_mat(i_tmp, :)
      write (1, 10) potential_mat(i_tmp, :)
   end do

   !----------------------------------------------------------------
   ! Generate Boundary Matrix (A-Prime)

   print *, "----------------------------------------------------------------"
   write (1,'(a)') "----------------------------------------------------------------"
   print *, "Boundary Matrix : "
   write (1,'(a)') "Boundary Matrix : "

   allocate (aprime_mat(2*n, 2*n))

   call boundary_cond(n, C, kappa, aprime_mat)
   do i_tmp = 1, 2*n
      write (*, 10) aprime_mat(i_tmp, :)
      write (1, 10) aprime_mat(i_tmp, :)
   end do

   !----------------------------------------------------------------
   ! Test

!    write (1, '(a)') "----------------------------------------------------------------"
!    write (1, '(a)') "Bspline 1 :"
!    call print_table(d, knot, bspline(1, :, :), 1)
!    write (1, '(a)') "----------------------------------------------------------------"
!    write (1, '(a)') "Bspline 2 :"
!    call print_table(d, knot, bspline(2, :, :), 1)
!    write (1, '(a)') "----------------------------------------------------------------"
!    write (1, '(a)') "Integral Overlap of Bspline 1 and Bspline 2 :"
!    call int_overlp(bspline(1, :, :), bspline(2, :, :), knot, result)
!    write (1, '(f10.4)') result
!    write (1, '(a)') "----------------------------------------------------------------"
!    write (1, '(a)') "Integral Potential of Bspline 1 and Bspline 2 :"
!    call int_potential(bspline(1, :, :), bspline(2, :, :), Z, knot, result)
!    write (1, '(f10.4)') result
!    write (1, '(a)') "----------------------------------------------------------------"

   close(1)
end program main
