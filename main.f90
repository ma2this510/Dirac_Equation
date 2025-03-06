program main
    use bspline_mod
    implicit none

    integer :: d, n, i_tmp, j_tmp
    integer, dimension(:), allocatable :: i_range
    real, dimension(:), allocatable :: knot
    real, dimension(:, :, :), allocatable :: bspline

    real, dimension(:,:), allocatable :: ovrlp_mat
    real :: result

    d = 7 ! Order of Mathemathica + 1
    n = 40
    allocate(knot(n + d + 3))

    knot = exp_knot(size(knot), 10.0, 0.1)

    allocate(i_range(n))
    i_range = [(i_tmp, i_tmp = 1, n)]

    allocate(bspline(i_range(n), size(knot), d))

    do i_tmp = 1, n
        call init_bspine(d, i_range(i_tmp), knot, bspline(i_tmp, :, :), .false.)
    end do

    print *, "Number of BSplines: ", n
    print *, "Order of BSplines: ", d, " and number of knots: ", size(knot)
    print *, "Knots: "
    write(*, 10) knot
    print *, "----------------------------------------------------------------"
    print *, "BSplines Overlaps : "

    allocate(ovrlp_mat(n,n))

    do i_tmp = 1, n
        do j_tmp = 1, n
            result = 0.0
            call int_overlp(bspline(i_tmp, :, :), bspline(j_tmp, :, :), knot, result)
            ovrlp_mat(i_tmp, j_tmp) = result
        end do
        write(*, 10) ovrlp_mat(i_tmp, :)
        10 format(50f7.1)
    end do

end program main