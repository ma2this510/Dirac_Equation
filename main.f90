! program main
!     use bspline_mod
!     implicit none

!     integer :: d, n
!     integer, allocatable :: i(:)
!     real, dimension(16) :: knot
!     real, dimension(:, :, :), allocatable :: result
!     real, dimension(:, :), allocatable :: primitive

!     integer :: i_tmp

!     d = 4 ! Order of Mathemathica + 1
!     n = size(knot) - d - 3

!     knot = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 6.0, 9.0, 12.0, 15.0, 15.0, 15.0, 15.0]

!     allocate(i(n))
!     i = [(i_tmp, i_tmp = 1, n)]

!     allocate(result(n, size(knot), d))

!     do i_tmp = 1, n
!         print *, i_tmp
!         call init_bspine(d, i(i_tmp), knot, result(i_tmp, :, :), .false.)
!     end do

!     call print_table(d, knot, result(1, :, :))
!     call print_table(d, knot, result(2, :, :))

!     allocate(primitive(size(result, 2), 2*size(result, 3)))

!     call int_overlp(result(1, :, :), result(2, :, :), knot, primitive)

!     call print_table(2*d, knot, primitive)
    
! end program main

program main
    use bspline_mod
    implicit none

    integer :: d
    real, dimension(16) :: knot
    real, dimension(:, :), allocatable :: b1, b2

    real, dimension(:, :), allocatable :: result

    d = 4 ! Order of Mathemathica + 1

    knot = [0.0, 0.0, 0.0, 0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 6.0, 9.0, 12.0, 15.0, 15.0, 15.0, 15.0]

    allocate(b1(size(knot), d))
    allocate(b2(size(knot), d))

    allocate(result(size(b1, 1), 2*size(b1, 2)))

    call init_bspine(d, 2, knot, b1, .false.)

    call init_bspine(d, 3, knot, b2, .false.)

    print *, "B-Spline number 1"
    print *, "-----------------"
    call print_table(d, knot, b1)

    print *, "B-Spline number 2"
    print *, "-----------------"
    call print_table(d, knot, b2)

    call int_overlp(b1, b2, knot, result)

    print *, "Primitive"
    print *, "---------"
    call print_table(2*d, knot, result)

end program main