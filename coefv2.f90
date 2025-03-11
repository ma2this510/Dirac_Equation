module bspline_mod
   use mpmodule
   implicit none

   type(mp_real), save :: zero, one

   private

   public :: init_bspine
   public :: fusion_coef
   public :: print_table
   public :: int_overlp
   public :: exp_knot
   public :: int_deriv
   public :: int_potential
   public :: boundary_cond
   public :: matrixAB

contains

   function fusion_coef(coef1, coef2)
      !> This function calculates the product of two polynoms
      !>
      !> @param coef1 : real(:) : the coef of the first polynom by increasing order
      !> @param coef2 : real(:) : the coef of the second polynom by increasing order
      !> @return fusion_coef : real(:) : the coef of the product of the two polynoms by increasing order
      type(mp_real), intent(in) :: coef1(:), coef2(:)
      type(mp_real) :: fusion_coef(size(coef1) + size(coef2) - 1)

      integer :: s1, s2, i, j

      fusion_coef = zero

      s1 = size(coef1)
      s2 = size(coef2)

      do i = 1, s1
         do j = 1, s2
            fusion_coef(i + j - 1) = fusion_coef(i + j - 1) + coef1(i)*coef2(j)
         end do
      end do

   end function fusion_coef

   recursive subroutine rec_coef(d, i, knot, table, sol_int, tot, index)
      !> @brief Recursive function to calculate the coefficients of the B-spline
      !> @param d : integer : the degree of the B-spline
      !> @param i : integer : the index of the B-spline
      !> @param knot : real(:) : the knot vector
      !> @param table : real(:,:) : the coef of the B-spline
      !> @param sol_int : real(:) : Variable to keep track of the filling of tot.
      !> @param tot : real(:) : Variable TOTAL fill when reach end of recursion
      !> @param index : integer : the index of the current B-spline

      integer, intent(in) :: d, i
      type(mp_real), intent(in), dimension(:) :: knot
      type(mp_real), intent(in), dimension(:) :: table
      type(mp_real), intent(out), dimension(:, :) :: tot
      integer, intent(out) :: sol_int
      integer, intent(out), dimension(:) :: index

      type(mp_real) :: coef1(2), coef2(2)
      type(mp_real) :: denum1, denum2

      type(mp_real), dimension(:), allocatable :: res1, res2

      denum1 = knot(i + d - 1) - knot(i)
      denum2 = knot(i + d) - knot(i + 1)

      if (denum1 == zero) then
         coef1(1) = zero
         coef1(2) = zero
      else
         coef1(1) = 1/denum1
         coef1(2) = -knot(i)/denum1
      end if

      if (denum2 == zero) then
         coef2(1) = zero
         coef2(2) = zero
      else
         coef2(1) = -1/denum2
         coef2(2) = knot(i + d)/denum2
      end if

      allocate (res1(size(table) + 1), res2(size(table) + 1))

      res1 = fusion_coef(table, coef1)
      res2 = fusion_coef(table, coef2)

      if (d > 2) then
         call rec_coef(d - 1, i, knot, res1, sol_int, tot, index)
         call rec_coef(d - 1, i + 1, knot, res2, sol_int, tot, index)
      else
         sol_int = sol_int + 1
         index(sol_int) = i
         tot(sol_int, :) = res1

         sol_int = sol_int + 1
         index(sol_int) = i + 1
         tot(sol_int, :) = res2
      end if

   end subroutine rec_coef

   subroutine print_table(d, knot, table, file)
      !> @brief Print the polynome between each node of a given spline
      !> @param d : integer : the degree of the B-spline
      !> @param knot : real(:) : the knot vector
      !> @param table : real(:,:) : the coef of the B-spline
      !> @param file : integer : the file to write the result
      implicit none
      integer, intent(in) :: d
      type(mp_real), intent(in) :: knot(:)
      type(mp_real), intent(in), dimension(size(knot), d) :: table
      integer, intent(in), optional :: file

      integer :: i_tmp, j_tmp

      if (present(file)) then
         do i_tmp = 1, size(knot) - 1
            write (file, *) "------------------------------------------------------------------"
            write (file, *) "for node", i_tmp, knot(i_tmp), "<= x < ", knot(i_tmp + 1)
            do j_tmp = 1, d
               write (file, *) "x^", d - j_tmp, " : ", table(i_tmp, j_tmp)
            end do
         end do
      else
         do i_tmp = 1, size(knot) - 1
            print *, "------------------------------------------------------------------"
            print *, "for node", i_tmp, knot(i_tmp), "<= x < ", knot(i_tmp + 1)
            do j_tmp = 1, d
               print *, "x^", d - j_tmp, " : ", table(i_tmp, j_tmp)
            end do
         end do
      end if
   end subroutine print_table

   function calcul_double(table, index, s) result(total)
      !> @brief Sum the different branch linked to the same internal node and return the final coefs per node
      !> @param table : real(:,:) : the un-proccessed coef of the B-spline
      !> @param index : integer(:) : the index of the current contribution
      !> @param s : integer : the number of nodes
      !> @return total : real(:,:) : the proccesed coef of the B-spline
      implicit none
      type(mp_real), intent(in) :: table(:, :)
      integer, intent(in) :: index(:)
      integer, intent(in) :: s
      type(mp_real), allocatable :: total(:, :)

      integer :: i, j

      allocate (total(s, size(table, 2)))

      total = zero

      do i = 1, size(table, 1)
         do j = 1, size(table, 2)
            total(index(i), j) = total(index(i), j) + table(i, j)
         end do
      end do

   end function calcul_double

   subroutine init_bspine(d, i, knot, result, display)
      !> @brief Main function to calculate the B-spline coefficients.
      !> @warning The degree and the index are the Mathematica values + 1
      !> @param d : integer : the degree of the B-spline
      !> @param i : integer : the index of the B-spline
      !> @param knot : mp_real(:) : the knot vector
      !> @param result : mp_real(:,:) : the final coef of the B-spline
      !> @param display : logical : display the result
      implicit none
      integer, intent(in) :: d, i
      type(mp_real), intent(in) :: knot(:)
      type(mp_real), intent(inout), dimension(size(knot), d) :: result
      logical, intent(in) :: display

      type(mp_real), dimension(1) :: table
      integer :: sol_int
      type(mp_real), dimension(2**(d - 1), d) :: tot
      integer, dimension(2**(d - 1)) :: index

      zero = '0.d0'
      one = '1.d0'

      if (d == 1) then
         print *, "The degree of the B-spline must be greater than 1"
         stop
      else if (i + d > size(knot)) then
         print *, "For the given degree and index, the knot vector needs to be at least of size ", i + d
         stop
      end if

      table = zero
      table(1) = one
      sol_int = 0
      tot = zero
      index = 0

      call rec_coef(d, i, knot, table, sol_int, tot, index)

      result = calcul_double(tot, index, size(knot))

      if (display) then
         call print_table(d, knot, result)
      end if

   end subroutine init_bspine

   subroutine int_overlp(b1, b2, knot, result)
      !> @brief Calculate the integral of the product of two B-splines
      !> @param b1 : real(:,:) : the coef of the first B-spline
      !> @param b2 : real(:,:) : the coef of the second B-spline
      !> @param knot : real(:) : the knot vector
      !> @param result : real : the result of the integral
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      type(mp_real), intent(in) :: knot(:)
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(size(b1, 1)) :: int_per_knot
      type(mp_real) :: int_init, int_final

      integer :: i_tmp, j_tmp
      type(mp_real), dimension(:, :), allocatable :: produit, primitive

      zero = '0.d0'
      one = '1.d0'

      allocate (produit(size(b1, 1), 2*size(b1, 2) - 1))
      allocate (primitive(size(b1, 1), 2*size(b1, 2)))

      produit = zero
      primitive = zero
      int_per_knot = zero

      do i_tmp = 1, size(b1, 1) ! Loop over the number of B-splines
         produit(i_tmp, :) = fusion_coef(b1(i_tmp, :), b2(i_tmp, :))
         primitive(i_tmp, 1:2*size(b1, 2) - 1) = produit(i_tmp, :)

      end do

      do j_tmp = 1, size(primitive, 2) - 1 ! Loop over the different orders
         do i_tmp = 1, size(primitive, 1)
            primitive(i_tmp, j_tmp) = primitive(i_tmp, j_tmp)/(size(primitive, 2) - j_tmp)
         end do
      end do

      do i_tmp = 1, size(b1, 1) - 1 ! Loop over the number of B-splines
         int_init = zero
         int_final = zero
         do j_tmp = 1, size(primitive, 2) - 1
            int_init = int_init + primitive(i_tmp, j_tmp)*knot(i_tmp)**(size(primitive, 2) - j_tmp)
            int_final = int_final + primitive(i_tmp, j_tmp)*knot(i_tmp + 1)**(size(primitive, 2) - j_tmp)
         end do
         int_per_knot(i_tmp) = int_final - int_init
      end do

      result = zero
      do i_tmp = 1, size(b1, 1)
         result = result + int_per_knot(i_tmp)
      end do

   end subroutine int_overlp

   function exp_knot(n, a_max, a_min, d, method, clt) result(result)
      !> @brief Generate a knot vector with exponential distribution
      !> @param n : integer : the number of knots
      !> @param a_max : mp_real : the maximum value of the knot vector
      !> @param a_min : mp_real : the minimum value of the knot vector
      !> @param d : integer : the degree of the B-spline
      !> @param method : integer : which methd will be used, 0 for lin, 1 for exp, 2 for fisher
      !> @param clt : mp_real : the clustering factor
      !> @return result : mp_real(n) : the knot vector
      integer, intent(in) :: n, d, method
      type(mp_real), intent(in) :: a_max, a_min, clt
      type(mp_real), dimension(n) :: result

      integer :: i_tmp
      zero = '0.d0'
      one = '1.d0'

      do i_tmp = 1, d
         result(i_tmp) = zero
      end do
      do i_tmp = d + 1, n - d + 1
         if (method == 0) then
            result(i_tmp) = a_min*(a_max/a_min)**(((i_tmp - (d + 1))*one/(n - 2*d)*one))
         else if (method == 1) then
            result(i_tmp) = a_min + (a_max - a_min)*(exp(clt*(i_tmp - d)/(n - 2*d + 1)) - 1)/(exp(clt) - 1)
         end if
      end do
      if (method == 2) then
         result(d + 1 : n - d + 1) = fisher_knot(n-2*d+1, a_max, a_min, clt)
      end if
      do i_tmp = n - d + 2, n
         result(i_tmp) = a_max
      end do

   end function exp_knot

   function fisher_knot(n, a_max, a_min, h) result(result)
      integer, intent(in) :: n
      type(mp_real), intent(in) :: a_max, a_min, h
      type(mp_real), dimension(n) :: result

      integer :: i_tmp

      result(1) = h
      do i_tmp = 2, n - 1
         if (result(i_tmp - 1)/5 < h) then
            result(i_tmp) = 1.2d0*result(i_tmp - 1)
         else
            result(i_tmp) = result(i_tmp - 1) + h
         end if
      end do
      result(n) = a_max

   end function fisher_knot

   function deriv(b1) result(result)
      !> @brief Calculate the derivative of a B-spline
      !> @param b1 : real(size(knot),d) : the coef of the B-spline
      !> @return result : real(size(knot),d) : the coef of the derivative of the B-spline
      type(mp_real), intent(in) :: b1(:, :)
      type(mp_real), dimension(size(b1, 1), size(b1, 2)) :: result

      integer :: i_tmp, j_tmp

      zero = '0.d0'
      one = '1.d0'

      result = zero
      result(:, 2:size(b1, 2)) = b1(:, 1:size(b1, 2) - 1) ! need to check

      do i_tmp = 1, size(b1, 1)
         do j_tmp = 1, size(b1, 2)
            result(i_tmp, j_tmp) = result(i_tmp, j_tmp)*(size(b1, 2) - j_tmp + 1)
         end do
      end do
   end function deriv

   subroutine int_deriv(b1, b2, knot, result)
      !> @brief Calculate the integral of the product of a B-spline and the derivative of another B-spline
      !> @param b1 : real(:,:) : the coef of the first B-spline
      !> @param b2 : real(:,:) : the coef of the second B-spline
      !> @param knot : real(:) : the knot vector
      !> @param result : real : the result of the integral
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      type(mp_real), intent(in) :: knot(:)
      type(mp_real), intent(out) :: result

      call int_overlp(b1, deriv(b2), knot, result)

   end subroutine int_deriv

   subroutine int_potential(b1, b2, Z, knot, result)
      !> @brief Calculate the integral of the product of a B-spline and the potential of another B-spline
      !> @param b1 : real(:,:) : the coef of the first B-spline
      !> @param b2 : real(:,:) : the coef of the second B-spline
      !> @param Z : real : the potential
      !> @param knot : real(:) : the knot vector
      !> @param result : real : the result of the integral
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      type(mp_real), intent(in) :: Z
      type(mp_real), intent(in) :: knot(:)
      type(mp_real), intent(out) :: result
      type(mp_real), dimension(size(b1, 1)) :: int_per_knot
      type(mp_real) :: int_init, int_final

      integer :: i_tmp, j_tmp
      type(mp_real), dimension(:, :), allocatable :: produit, primitive

      allocate (produit(size(b1, 1), 2*size(b1, 2) - 1))
      allocate (primitive(size(b1, 1), 2*size(b1, 2) - 1))

      zero = '0.d0'
      one = '1.d0'

      produit = zero
      primitive = zero
      int_per_knot = zero

      do i_tmp = 1, size(b1, 1) ! Loop over the number of B-splines
         produit(i_tmp, :) = fusion_coef(b1(i_tmp, :), b2(i_tmp, :))
      end do

      primitive = produit

      do j_tmp = 1, size(primitive, 2) - 1 ! Loop over the different orders
         ! Remove the last order (x^-1) as there is no factor in front of it
         do i_tmp = 1, size(primitive, 1)
            primitive(i_tmp, j_tmp) = primitive(i_tmp, j_tmp)/(size(primitive, 2) - j_tmp)
         end do
      end do

      do i_tmp = 1, size(b1, 1) - 1 ! Loop over the number of B-splines
         int_init = zero
         int_final = zero
         do j_tmp = 1, size(primitive, 2) - 1 ! Loop over the different orders
            int_init = int_init + primitive(i_tmp, j_tmp)*knot(i_tmp)**(size(primitive, 2) - j_tmp)
            int_final = int_final + primitive(i_tmp, j_tmp)*knot(i_tmp + 1)**(size(primitive, 2) - j_tmp)
         end do
         ! Integrate the last order (x^-1) to log(x)
         if (primitive(i_tmp, size(primitive, 2)) /= zero) then
            int_init = int_init + primitive(i_tmp, size(primitive, 2))*log(knot(i_tmp))
            int_final = int_final + primitive(i_tmp, size(primitive, 2))*log(knot(i_tmp + 1))
         end if
         int_per_knot(i_tmp) = int_final - int_init
      end do

      result = zero
      do i_tmp = 1, size(b1, 1) - 1
         result = result + int_per_knot(i_tmp)
      end do

      result = -Z*result

   end subroutine int_potential

   subroutine boundary_cond(n, C, kappa, aprime)
      !> @brief Calculate the boundary conditions
      !> @param n : integer : the number of B-splines
      !> @param C : real : speed of light
      !> @param kappa : real : Relativistic quantum number
      !> @param aprime : real(n,n) : the boundary condition
      implicit none
      integer, intent(in) :: n
      type(mp_real), intent(in) :: C, kappa
      type(mp_real), intent(out) :: aprime(2*n, 2*n)

      zero = '0.d0'
      one = '1.d0'

      aprime = zero
      ! if (kappa > zero) then
      !    aprime(1, 1) = 2*C**2
      !    aprime(1, n + 1) = -C/2
      !    aprime(n + 1, 1) = -C/2
      !    aprime(n, n) = C/2
      !    aprime(2*n, 2*n) = -C/2
      ! else if (kappa < zero) then
      !    aprime(1, 1) = C
      !    aprime(1, n + 1) = -C/2
      !    aprime(n + 1, 1) = -C/2
      !    aprime(n, n) = C/2
      !    aprime(2*n, 2*n) = -C/2
      ! end if
   end subroutine boundary_cond

   subroutine write_lists(r1s, i1, i2, i3)
      !> @brief Write a list of mp_real
      !> @param r1s : mp_real(:) : the list of mp_real
      !> @param i1 : integer : the unit number
      !> @param i2 : integer : the field width
      !> @param i3 : integer : the number of decimal
      implicit none
      integer, intent(in) :: i1, i2, i3
      type(mp_real), intent(in), dimension(:) :: r1s

      character(i2) :: str_tmp(1)
      integer :: i_tmp

      do i_tmp = 1, size(r1s)
         str_tmp = ""
         call mpeform(r1s(i_tmp), i2, i3, str_tmp)
         write (i1, '(a)', advance='no') str_tmp
      end do
      write (i1, *)  ! End the line after the list is printed

   end subroutine write_lists

   subroutine matrixAB(d, n, n_remove, Z, kappa, C, amin, amax, method, clt, A, B, log_bool, i2, i3)
      !> @brief Generate the matrix A and B
      !> @param d : integer : the degree of the B-spline
      !> @param n : integer : the number of B-splines
      !> @param n : integer : the number of B-splines to be remove at the start and the end
      !> @param Z : mp_real : the potential constant
      !> @param kappa : mp_real : the relativistic quantum number
      !> @param C : mp_real : the speed of light
      !> @param amin : mp_real : the minimum value of the knot vector
      !> @param amax : mp_real : the maximum value of the knot vector
      !> @param A : mp_real(2*nprime,2*nprime) : the matrix A
      !> @param B : mp_real(2*nprime,2*nprime) : the matrix B
      !> @param log_bool : logical : display the result
      !> @param i2 : integer : the field width
      !> @param i3 : integer : the number of decimal
      implicit none
      integer, intent(in) :: d, n, method, n_remove
      type(mp_real), intent(in) :: Z, kappa, C, amin, amax, clt
      type(mp_real), intent(out), allocatable, dimension(:, :) :: A, B
      logical, intent(in), optional :: log_bool
      integer, intent(in), optional :: i2, i3

      integer :: nprime
      type(mp_real) :: zero, one

      type(mp_real), dimension(n + d) :: knot
      type(mp_real), dimension(n, size(knot), d) :: bspline

      type(mp_real), dimension(:, :), allocatable :: ovrlp_mat, deriv_mat, potential_mat, k_mat
      type(mp_real), allocatable, dimension(:, :) :: aprime_mat

      integer :: i_tmp, j_tmp

      zero = '0.d0'
      one = '1.d0'

      ! Number of B-splines without the two first and last one
      nprime = n - 2*n_remove

      ! Generate the knot vector
      print *, "Generate the knot vector"
      knot = exp_knot(n + d, amax, amin, d, method, clt)

      ! Generate the B-splines
      print *, "Generate the B-splines"
      !OMP PARALLEL DO PRIVATE(i_tmp) SHARED(d, knot, bspline, n)
      do i_tmp = 1, n
         call init_bspine(d, i_tmp, knot, bspline(i_tmp, :, :), .false.)
      end do
      !OMP END PARALLEL DO

      ! Generate the overlap matrix
      print *, "Generate the overlap matrix"
      allocate (ovrlp_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, knot, ovrlp_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
            call int_overlp(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), knot, ovrlp_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO

      ! Generate the derivative matrix
      print *, "Generate the derivative matrix"
      allocate (deriv_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, knot, deriv_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
            call int_deriv(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), knot, deriv_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO

      ! Generate the potential matrix
      print *, "Generate the potential matrix"
      allocate (potential_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, Z, knot, potential_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
            call int_potential(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), Z, knot, potential_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO

      ! Generate K matrix
      print *, "Generate the K matrix"
      allocate (k_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, kappa, knot, k_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
            call int_potential(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), -kappa, knot, k_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO

      ! Generate the boundary condition matrix
      print *, "Generate the boundary condition matrix"
      allocate (aprime_mat(2*nprime, 2*nprime))
      call boundary_cond(nprime, C, kappa, aprime_mat)

      ! Generate the matrix A
      print *, "Generate the matrix A"
      allocate (A(2*nprime, 2*nprime))
      A = zero
      A(1:nprime, 1:nprime) = potential_mat

      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
            A(i_tmp, nprime + j_tmp) = C*(deriv_mat(i_tmp, j_tmp) - k_mat(i_tmp, j_tmp))
            A(nprime + i_tmp, j_tmp) = -C*(deriv_mat(i_tmp, j_tmp) + k_mat(i_tmp, j_tmp))
            A(nprime + i_tmp, nprime + j_tmp) = potential_mat(i_tmp, j_tmp) - 2*(C**2)*ovrlp_mat(i_tmp, j_tmp)
         end do
      end do

      do i_tmp = 1, 2*nprime
         do j_tmp = 1, 2*nprime
            A(i_tmp, j_tmp) = A(i_tmp, j_tmp) + aprime_mat(i_tmp, j_tmp)
         end do
      end do

      ! Generate the matrix B
      print *, "Generate the matrix B"
      allocate (B(2*nprime, 2*nprime))
      B = zero
      B(1:nprime, 1:nprime) = ovrlp_mat
      B(nprime + 1:2*nprime, nprime + 1:2*nprime) = ovrlp_mat

      if (present(log_bool) .and. log_bool) then
         print *, "Writing Logs"
         open (1, file="log_2.txt", status="replace")

         write (1, '(a,i4)') "Number of BSplines: ", n
         write (1, '(a,i4)') "Order of BSplines: ", d, " and number of knots: ", size(knot)
         write (1, '(a)') "Knots: "
         call write_lists(knot, 1, i2, i3)
         write (1, '(a)') "----------------------------------------------------------------"
         write (1, '(a)') "BSplines Overlaps Matrix : "
         do i_tmp = 1, nprime
            call write_lists(ovrlp_mat(i_tmp, :), 1, i2, i3) ! This is where the Segmentation fault occurs from valgrind
         end do
         write (1, '(a)') "----------------------------------------------------------------"
         write (1, '(a)') "BSplines Derivatives Matrix : "
         do i_tmp = 1, nprime
            call write_lists(deriv_mat(i_tmp, :), 1, i2, i3)
         end do
         write (1, '(a)') "----------------------------------------------------------------"
         write (1, '(a)') "BSplines Potentials Matrix : "
         do i_tmp = 1, nprime
            call write_lists(potential_mat(i_tmp, :), 1, i2, i3)
         end do
         write (1, '(a)') "----------------------------------------------------------------"
         write (1, '(a)') "K/r Matrix : "
         do i_tmp = 1, nprime
            call write_lists(k_mat(i_tmp, :), 1, i2, i3)
         end do
         write (1, '(a)') "----------------------------------------------------------------"
         write (1, '(a)') "Boundary Condition Matrix : "
         do i_tmp = 1, 2*nprime
            call write_lists(aprime_mat(i_tmp, :), 1, i2, i3)
         end do
         write (1, '(a)') "----------------------------------------------------------------"
         write (1, '(a)') "Matrix A : "
         do i_tmp = 1, 2*nprime
            call write_lists(A(i_tmp, :), 1, i2, i3)
         end do
         write (1, '(a)') "----------------------------------------------------------------"
         write (1, '(a)') "Matrix B : "
         do i_tmp = 1, 2*nprime
            call write_lists(B(i_tmp, :), 1, i2, i3)
         end do
         close (1)
         print *, "Logs written"
      end if
   end subroutine matrixAB
end module bspline_mod
