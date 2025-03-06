module bspline_mod
   implicit none

   private

   public :: init_bspine
   public :: fusion_coef
   public :: print_table
   public :: int_overlp
   public :: exp_knot

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

   subroutine print_table(d, knot, table)
      !> @brief Print the polynome between each node of a given spline
      !> @param d : integer : the degree of the B-spline
      !> @param knot : real(:) : the knot vector
      !> @param table : real(:,:) : the coef of the B-spline
      implicit none
      integer, intent(in) :: d
      real, intent(in) :: knot(:)
      real, intent(in), dimension(size(knot), d) :: table

      integer :: i_tmp, j_tmp

      do i_tmp = 1, size(knot) - 1
         print *, "------------------------------------------------------------------"
         print *, "for node", i_tmp, knot(i_tmp), "<= x < ", knot(i_tmp + 1)
         do j_tmp = 1, d
            print *, "x^", d - j_tmp, " : ", table(i_tmp, j_tmp)
         end do
      end do

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

      do i_tmp = 1, size(b1, 1)
         produit(i_tmp, :) = fusion_coef(b1(i_tmp, :), b2(i_tmp, :))
         primitive(i_tmp, 1:2*size(b1, 2) - 1) = produit(i_tmp, :)

      end do

      do j_tmp = 1, size(primitive, 2) - 1
         primitive(:, j_tmp) = primitive(:, j_tmp)/(size(primitive, 2) - j_tmp)
      end do

      do i_tmp = 1, size(b1, 1) - 1
         int_init = 0.0
         int_final = 0.0
         do j_tmp = 1, size(primitive, 2) - 1
            int_init = int_init + primitive(i_tmp, j_tmp)*knot(i_tmp)**(size(primitive, 2) - j_tmp)
            int_final = int_final + primitive(i_tmp, j_tmp)*knot(i_tmp + 1)**(size(primitive, 2) - j_tmp)
         end do
         int_per_knot(i_tmp) = int_final - int_init
      end do

      do i_tmp = 1, size(b1, 1)
         result = result + int_per_knot(i_tmp)
      end do

   end subroutine int_overlp

   function exp_knot(n, a_max, a_min) result(result)
      integer, intent(in) :: n
      real, intent(in) :: a_max, a_min
      real, dimension(n) :: result

      integer :: i_tmp

      do i_tmp = 1, n
         result(i_tmp) = a_min*(a_max/a_min)**((real(i_tmp - 1)/real(n - 1)))
      end do

   end function exp_knot
end module bspline_mod
