module bspline_mod
   use mpmodule
   use bspline_gen
   use tools_mp
   implicit none

   type(mp_real), save :: zero, one

   private

   public :: int_overlp
   public :: exp_knot
   public :: int_deriv
   public :: int_potential

contains

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
      !> @param method : integer : which methd will be used, 0 for lin, 1 for exp
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

      do i_tmp = n - d + 2, n
         result(i_tmp) = a_max
      end do

   end function exp_knot

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

   

end module bspline_mod
