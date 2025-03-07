module bspline_mod
   implicit none

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
      real, intent(in) :: coef1(:), coef2(:)
      real :: fusion_coef(size(coef1) + size(coef2) - 1)

      integer :: s1, s2, i, j

      fusion_coef = 0.0

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
      real, intent(in), dimension(:) :: knot
      real, intent(in), dimension(:) :: table
      real, intent(out), dimension(:, :) :: tot
      integer, intent(out) :: sol_int
      integer, intent(out), dimension(:) :: index

      real :: coef1(2), coef2(2)
      real :: denum1, denum2

      real, dimension(:), allocatable :: res1, res2

      denum1 = knot(i + d - 1) - knot(i)
      denum2 = knot(i + d) - knot(i + 1)

      if (denum1 == 0.0) then
         coef1(1) = 0.0
         coef1(2) = 0.0
      else
         coef1(1) = 1/denum1
         coef1(2) = -knot(i)/denum1
      end if

      if (denum2 == 0.0) then
         coef2(1) = 0.0
         coef2(2) = 0.0
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
      real, intent(in) :: knot(:)
      real, intent(in), dimension(size(knot), d) :: table
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
      real, intent(in) :: table(:, :)
      integer, intent(in) :: index(:)
      integer, intent(in) :: s
      real, allocatable :: total(:, :)

      integer :: i

      allocate (total(s, size(table, 2)))

      total = 0.0

      do i = 1, size(table, 1)
         total(index(i), :) = total(index(i), :) + table(i, :)
      end do

   end function calcul_double

   subroutine init_bspine(d, i, knot, result, display)
      !> @brief Main function to calculate the B-spline coefficients.
      !> @warning The degree and the index are the Mathematica values + 1
      !> @param d : integer : the degree of the B-spline
      !> @param i : integer : the index of the B-spline
      !> @param knot : real(:) : the knot vector
      !> @param result : real(:,:) : the final coef of the B-spline
      !> @param display : logical : display the result
      implicit none
      integer, intent(in) :: d, i
      real, intent(in) :: knot(:)
      real, intent(inout), dimension(size(knot), d) :: result
      logical, intent(in) :: display

      real, dimension(1) :: table
      integer :: sol_int
      real, dimension(2**(d - 1), d) :: tot
      integer, dimension(2**(d - 1)) :: index

      if (d == 1) then
         print *, "The degree of the B-spline must be greater than 1"
         stop
      else if (i + d > size(knot)) then
         print *, "For the given degree and index, the knot vector needs to be at least of size ", i + d
         stop
      end if

      table = 0.0
      table(1) = 1.0
      sol_int = 0
      tot = 0.0
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
      real, intent(in), dimension(:, :) :: b1, b2
      real, intent(in) :: knot(:)
      real, intent(out) :: result
      real, dimension(size(b1, 1)) :: int_per_knot
      real :: int_init, int_final

      integer :: i_tmp, j_tmp
      real, dimension(:, :), allocatable :: produit, primitive

      allocate (produit(size(b1, 1), 2*size(b1, 2) - 1))
      allocate (primitive(size(b1, 1), 2*size(b1, 2)))

      produit = 0.0
      primitive = 0.0

      do i_tmp = 1, size(b1, 1) ! Loop over the number of B-splines
         produit(i_tmp, :) = fusion_coef(b1(i_tmp, :), b2(i_tmp, :))
         primitive(i_tmp, 1:2*size(b1, 2) - 1) = produit(i_tmp, :)

      end do

      do j_tmp = 1, size(primitive, 2) - 1 ! Loop over the different orders
         primitive(:, j_tmp) = primitive(:, j_tmp)/(size(primitive, 2) - j_tmp)
      end do

      do i_tmp = 1, size(b1, 1) - 1 ! Loop over the number of B-splines
         int_init = 0.0
         int_final = 0.0
         do j_tmp = 1, size(primitive, 2) - 1
            int_init = int_init + primitive(i_tmp, j_tmp)*knot(i_tmp)**(size(primitive, 2) - j_tmp)
            int_final = int_final + primitive(i_tmp, j_tmp)*knot(i_tmp + 1)**(size(primitive, 2) - j_tmp)
         end do
         int_per_knot(i_tmp) = int_final - int_init
      end do

      result = 0.0
      do i_tmp = 1, size(b1, 1)
         result = result + int_per_knot(i_tmp)
      end do

   end subroutine int_overlp

   function exp_knot(n, a_max, a_min) result(result)
      !> @brief Generate a knot vector with exponential distribution
      !> @param n : integer : the number of knots
      !> @param a_max : real : the maximum value of the knot vector
      !> @param a_min : real : the minimum value of the knot vector
      !> @return result : real(n) : the knot vector
      integer, intent(in) :: n
      real, intent(in) :: a_max, a_min
      real, dimension(n) :: result

      integer :: i_tmp

      do i_tmp = 1, n
         result(i_tmp) = a_min*(a_max/a_min)**((real(i_tmp - 1)/real(n - 1)))
      end do

   end function exp_knot

   function deriv(b1) result(result)
      !> @brief Calculate the derivative of a B-spline
      !> @param b1 : real(size(knot),d) : the coef of the B-spline
      !> @return result : real(size(knot),d-1) : the coef of the derivative of the B-spline
      real, intent(in) :: b1(:, :)
      real, dimension(size(b1, 1), size(b1, 2)) :: result

      integer :: i_tmp

      result = 0.0
      result(:, 2:size(b1, 2) - 1) = b1(:, 1:size(b1, 2) - 1)

      do i_tmp = 1, size(b1, 1)
         result(i_tmp, :) = result(i_tmp, :)*(size(b1, 2) - [(i_tmp - 1, i_tmp=1, size(b1, 2))])
      end do
   end function deriv

   subroutine int_deriv(b1, b2, knot, result)
      !> @brief Calculate the integral of the product of a B-spline and the derivative of another B-spline
      !> @param b1 : real(:,:) : the coef of the first B-spline
      !> @param b2 : real(:,:) : the coef of the second B-spline
      !> @param knot : real(:) : the knot vector
      !> @param result : real : the result of the integral
      implicit none
      real, intent(in), dimension(:, :) :: b1, b2
      real, intent(in) :: knot(:)
      real, intent(out) :: result

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
      real, intent(in), dimension(:, :) :: b1, b2
      real, intent(in) :: Z
      real, intent(in) :: knot(:)
      real, intent(out) :: result
      real, dimension(size(b1, 1)) :: int_per_knot
      real :: int_init, int_final

      integer :: i_tmp, j_tmp
      real, dimension(:, :), allocatable :: produit, primitive

      allocate (produit(size(b1, 1), 2*size(b1, 2) - 1))
      allocate (primitive(size(b1, 1), 2*size(b1, 2) - 1))

      produit = 0.0
      primitive = 0.0

      do i_tmp = 1, size(b1, 1) ! Loop over the number of B-splines
         produit(i_tmp, :) = fusion_coef(b1(i_tmp, :), b2(i_tmp, :))
      end do

      primitive = produit

      do j_tmp = 1, size(primitive, 2) - 1 ! Loop over the different orders
         ! Remove the last order (x^-1) as there is no factor in front of it
         primitive(:, j_tmp) = primitive(:, j_tmp)/(size(primitive, 2) - j_tmp)
      end do

      do i_tmp = 1, size(b1, 1) - 1 ! Loop over the number of B-splines
         int_init = 0.0
         int_final = 0.0
         do j_tmp = 1, size(primitive, 2) - 1 ! Loop over the different orders
            int_init = int_init + primitive(i_tmp, j_tmp)*knot(i_tmp)**(size(primitive, 2) - j_tmp)
            int_final = int_final + primitive(i_tmp, j_tmp)*knot(i_tmp + 1)**(size(primitive, 2) - j_tmp)
         end do
         ! Integrate the last order (x^-1) to log(x)
         int_init = int_init + primitive(i_tmp, size(primitive, 2) )*log(knot(i_tmp))
         int_final = int_final + primitive(i_tmp, size(primitive, 2))*log(knot(i_tmp + 1))
         int_per_knot(i_tmp) = int_final - int_init
      end do

      result = 0.0
      do i_tmp = 1, size(b1, 1)
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
      real, intent(in) :: C, kappa
      real, intent(out) :: aprime(2*n, 2*n)

      aprime = 0.0
      if (kappa > 0.0) then
         aprime(1, 1) = 2*C**2
         aprime(1, n+1) = -C/2
         aprime(n+1, 1) = -C/2
         aprime(n, n) = C/2
         aprime(2*n, 2*n) = -C/2
      else if (kappa < 0.0) then
         aprime(1, 1) = C
         aprime(1, n+1) = -C/2
         aprime(n+1, 1) = -C/2
         aprime(n, n) = C/2
         aprime(2*n, 2*n) = -C/2
      end if
   end subroutine boundary_cond

   subroutine matrixAB(d, n, Z, kappa, C, amin, amax, A, B, log)
      !> @brief Generate the matrix A and B
      !> @param d : integer : the degree of the B-spline
      !> @param n : integer : the number of B-splines
      !> @param Z : real : the potential constant
      !> @param kappa : real : the relativistic quantum number
      !> @param C : real : the speed of light
      !> @param amin : real : the minimum value of the knot vector
      !> @param amax : real : the maximum value of the knot vector
      !> @param A : real(2n,2n) : the matrix A
      !> @param B : real(2n,2n) : the matrix B
      !> @param log : logical : display the result
      implicit none
      integer, intent(in) :: d, n
      real, intent(in) :: Z, kappa, C, amin, amax
      real, intent(out), allocatable, dimension(:, :) :: A, B
      logical, intent(in), optional :: log

      integer :: nprime

      real, dimension(n+d) :: knot
      real, dimension(n, size(knot), d) :: bspline

      real, dimension(n, n) :: ovrlp_mat, deriv_mat, potential_mat, k_mat
      real, allocatable, dimension(:, :) :: aprime_mat

      integer :: i_tmp, j_tmp

      ! Number of B-splines without the first and last one
      nprime = n-2

      ! Generate the knot vector
      knot = exp_knot(n + d, amax, amin)

      ! Generate the B-splines
      do i_tmp = 1, n
         call init_bspine(d, i_tmp, knot, bspline(i_tmp, :, :), .false.)
      end do

      ! Generate the overlap matrix
      do i_tmp = 1, n
         do j_tmp = 1, n
            call int_overlp(bspline(i_tmp, :, :), bspline(j_tmp, :, :), knot, ovrlp_mat(i_tmp, j_tmp))
         end do
      end do

      ! Generate the derivative matrix
      do i_tmp = 1, n
         do j_tmp = 1, n
            call int_deriv(bspline(i_tmp, :, :), bspline(j_tmp, :, :), knot, deriv_mat(i_tmp, j_tmp))
         end do
      end do

      ! Generate the potential matrix
      do i_tmp = 1, n
         do j_tmp = 1, n
            call int_potential(bspline(i_tmp, :, :), bspline(j_tmp, :, :), Z, knot, potential_mat(i_tmp, j_tmp))
         end do
      end do

      ! Generate K matrix
      do i_tmp = 1, n
         do j_tmp = 1, n
            call int_potential(bspline(i_tmp, :, :), bspline(j_tmp, :, :), -kappa, knot, k_mat(i_tmp, j_tmp))
         end do
      end do

      ! Generate the boundary condition matrix
      allocate (aprime_mat(2*nprime, 2*nprime))
      call boundary_cond(nprime, C, kappa, aprime_mat)

      ! Generate the matrix A
      allocate (A(2*nprime, 2*nprime))
      A = 0.0
      A(1:nprime, 1:nprime) = potential_mat(2:n-1, 2:n-1)
      A(1:nprime, nprime+1:2*nprime) = C * (deriv_mat(2:n-1, 2:n-1) + k_mat(2:n-1, 2:n-1))
      A(nprime+1:2*nprime, 1:nprime) = -C * (deriv_mat(2:n-1, 2:n-1) + k_mat(2:n-1, 2:n-1))
      A(nprime+1:2*nprime, nprime+1:2*nprime) = potential_mat(2:n-1, 2:n-1) - 2*(C**2) * ovrlp_mat(2:n-1, 2:n-1)
      A = A + aprime_mat

      ! Generate the matrix B
      allocate (B(2*nprime, 2*nprime))
      B = 0.0
      B(1:nprime, 1:nprime) = ovrlp_mat(2:n-1, 2:n-1)
      B(nprime+1:2*nprime, nprime+1:2*nprime) = ovrlp_mat(2:n-1, 2:n-1)

      if (present(log)) then
         open(1, file="log_2.txt", status="replace")
10       format(500f15.8)
         write (1,'(a,i4)') "Number of BSplines: ", n
         write (1,'(a,i4)') "Order of BSplines: ", d, " and number of knots: ", size(knot)
         write (1,'(a)') "Knots: "
         write (1, 10) knot
         write (1,'(a)') "----------------------------------------------------------------"
         write (1,'(a)') "BSplines Overlaps Matrix : "
         do i_tmp = 1, n
            write (1, 10) ovrlp_mat(i_tmp, :)
         end do
         write (1,'(a)') "----------------------------------------------------------------"
         write (1,'(a)') "BSplines Derivatives Matrix : "
         do i_tmp = 1, n
            write (1, 10) deriv_mat(i_tmp, :)
         end do
         write (1,'(a)') "----------------------------------------------------------------"
         write (1,'(a)') "BSplines Potentials Matrix : "
         do i_tmp = 1, n
            write (1, 10) potential_mat(i_tmp, :)
         end do
         write (1,'(a)') "----------------------------------------------------------------"
         write (1,'(a)') "K/r Matrix : "
         do i_tmp = 1, n
            write (1, 10) k_mat(i_tmp, :)
         end do
         write (1,'(a)') "----------------------------------------------------------------"
         write (1,'(a)') "Boundary Condition Matrix : "
         do i_tmp = 1, 2*nprime
            write (1, 10) aprime_mat(i_tmp, :)
         end do
         write (1,'(a)') "----------------------------------------------------------------"
         write (1,'(a)') "Matrix A : "
         do i_tmp = 1, 2*nprime
            write (1, 10) A(i_tmp, :)
         end do
         write (1,'(a)') "----------------------------------------------------------------"
         write (1,'(a)') "Matrix B : "
         do i_tmp = 1, 2*nprime
            write (1, 10) B(i_tmp, :)
         end do
         close(1)
      end if
   end subroutine matrixAB
end module bspline_mod
