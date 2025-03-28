module bspline_mod
   use mpmodule
   use bspline_gen
   use tools_mp
   implicit none

   type(mp_real), save :: zero, one

   private

   public :: theoric_val
   public :: exp_knot
   public :: matrixAB
   public :: get_eigen

contains

   function theoric_val(n, Z, kappa, C)
      implicit none
      integer, intent(in) :: n
      type(mp_real), intent(in) :: Z, kappa, C
      type(mp_real) :: theoric_val

      theoric_val = (C**2)/sqrt(1 + ((Z/C)**2)/(n - abs(kappa) + sqrt(kappa**2 - (Z/C)**2))**2)
   end function theoric_val

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

   subroutine integral(b1, order1, b2, order2, knot, result)
      !> @brief Calculate the integral of the product of two Piecewise polynomial of any order
      !> @param b1 : real(:,:) : the coef of the first Piecewise polynomial
      !> @param order1 : integer : the max order of the first Piecewise polynomial
      !> @param b2 : real(:,:) : the coef of the second Piecewise polynomial
      !> @param order2 : integer : the max order of the second Piecewise polynomial
      !> @param knot : real(:) : the knot vector
      !> @param result : real : the result of the integral
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      integer, intent(in) :: order1, order2
      type(mp_real), intent(in) :: knot(:)
      type(mp_real), intent(out) :: result

      type(mp_real), dimension(size(b1, 1)) :: int_per_knot

      integer :: order_prod, i_tmp, j_tmp
      type(mp_real), dimension(:, :), allocatable :: produit, int_init, int_final

      zero = '0.d0'
      one = '1.d0'

      allocate (produit(size(b1, 1), 2*size(b1, 2) - 1))
      allocate (int_final(size(b1, 1), 2*size(b1, 2)))
      allocate (int_init(size(b1, 1), 2*size(b1, 2)))

      produit = zero
      order_prod = order1 + order2
      int_final = zero
      int_init = zero
      int_per_knot = zero

      do i_tmp = 1, size(b1, 1) ! Loop over the number of Piecewise polynomial
         produit(i_tmp, :) = fusion_coef(b1(i_tmp, :), b2(i_tmp, :))
      end do

      do i_tmp = 1, size(produit, 1) - 1 ! Loop over the number of Piecewise polynomial
         do j_tmp = 1, size(produit, 2) ! Loop over the different orders
            if (order_prod - j_tmp + 2 /= 0) then
               if (knot(i_tmp + 1) /= zero) then
                int_final(i_tmp, j_tmp) = produit(i_tmp, j_tmp)*(knot(i_tmp + 1)**(order_prod - j_tmp + 2))/(order_prod - j_tmp + 2)
                  ! +2 because of the integral that add one order
               end if
               if (knot(i_tmp) /= zero) then
                  int_init(i_tmp, j_tmp) = produit(i_tmp, j_tmp)*(knot(i_tmp)**(order_prod - j_tmp + 2))/(order_prod - j_tmp + 2)
               end if
            else
               if (knot(i_tmp + 1) /= zero) then
                  int_final(i_tmp, j_tmp) = produit(i_tmp, j_tmp)*log(knot(i_tmp + 1))
               end if
               if (knot(i_tmp) /= zero) then
                  int_init(i_tmp, j_tmp) = produit(i_tmp, j_tmp)*log(knot(i_tmp))
               end if
            end if
         end do
      end do

      do i_tmp = 1, size(b1, 1) - 1 ! Loop over the number of Piecewise polynomial
         int_per_knot(i_tmp) = zero
         do j_tmp = 1, size(produit, 2) ! Loop over the different orders
            int_per_knot(i_tmp) = int_per_knot(i_tmp) + int_final(i_tmp, j_tmp) - int_init(i_tmp, j_tmp)
         end do
      end do

      result = zero
      do i_tmp = 1, size(b1, 1) - 1
         result = result + int_per_knot(i_tmp)
      end do

   end subroutine integral

   function deriv_single(b, order) result(result)
      !> @brief Calculate the first derivative of a B-spline
      !> @Warning : work when divided by r or r^2
      !> @param b : real(:, :) : the coef of the B-spline
      !> @param order : integer : the max order of the coef
      !> @return result : real(:) : the coef of the first derivative of the B-spline
      implicit none
      type(mp_real), intent(in) :: b(:, :)
      integer, intent(in) :: order
      type(mp_real) :: result(size(b, 1), size(b, 2))

      integer :: i_tmp, j_tmp

      result = zero
      do i_tmp = 1, size(b, 1)
         do j_tmp = 1, size(b, 2)
            result(i_tmp, j_tmp) = b(i_tmp, j_tmp)*dble(order - j_tmp + 1)
         end do
      end do
      ! Implicite : order of the max coef is decreased by 1
   end function deriv_single

   subroutine evaluate_poly(b, order, x, result)
      !> @brief Evaluate a polynomial at a point
      !> @param b : real(:) : the coef of the polynomial
      !> @param order : integer : the max order of the coef
      !> @param x : real : the point to evaluate the polynomial
      !> @param result : real : the result of the evaluation
      implicit none
      type(mp_real), intent(in) :: b(:)
      integer, intent(in) :: order
      type(mp_real), intent(in) :: x
      type(mp_real), intent(out) :: result

      integer :: i_tmp

      result = zero
      do i_tmp = 1, size(b)
         if (order - i_tmp + 1 /= 0) then
            if (x /= zero) then
               result = result + b(i_tmp)*x**(order - i_tmp + 1)
            end if
         else 
            result = result + b(i_tmp)
         end if
      end do
   end subroutine evaluate_poly

   subroutine matrix_S(b1, b2, knot, Z, kappa, C, result)
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      type(mp_real), intent(in) :: Z, kappa, C
      type(mp_real), intent(in) :: knot(:)
      type(mp_real) :: result

      call integral(b1, size(b1, 2) - 1, b2, size(b2, 2) - 1, knot, result)
   end subroutine matrix_S

   subroutine matrix_phi(b1, b2, knot, Z, kappa, C, result)
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      type(mp_real), intent(in) :: Z, kappa, C
      type(mp_real), intent(in) :: knot(:)
      type(mp_real) :: result

      type(mp_real), dimension(size(b1, 1), size(b1, 2)) :: term1

      term1 = multiply_elem(Z, b2)

      call integral(b1, size(b1, 2) - 1, term1, size(b2, 2) - 2, knot, result)
   end subroutine matrix_phi

   subroutine matrix_T_plus(b1, b2, knot, Z, kappa, C, result)
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      type(mp_real), intent(in) :: Z, kappa, C
      type(mp_real), intent(in) :: knot(:)
      type(mp_real) :: result

      type(mp_real), dimension(size(b1, 1), size(b1, 2)) :: term1, term2
      type(mp_real) :: result1, result2
      integer :: order1, order2

      term1 = deriv_single(deriv_single(b2, size(b2, 2) - 1), size(b2, 2) - 2)
      order1 = size(b2, 2) - 3

      term2 = multiply_elem(kappa*(1 + kappa), b2)
      order2 = size(b2, 2) - 3

      call integral(b1, size(b1, 2) - 1, term1, order1, knot, result1)
      call integral(b1, size(b1, 2) - 1, term2, order2, knot, result2)

      result = (-1)*(result1 - result2)/2
   end subroutine matrix_T_plus

   subroutine matrix_T_minus(b1, b2, knot, Z, kappa, C, result)
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      type(mp_real), intent(in) :: Z, kappa, C
      type(mp_real), intent(in) :: knot(:)
      type(mp_real) :: result

      type(mp_real), dimension(size(b1, 1), size(b1, 2)) :: term1, term2
      type(mp_real) :: result1, result2
      integer :: order1, order2

      term1 = deriv_single(deriv_single(b2, size(b2, 2) - 1), size(b2, 2) - 2)
      order1 = size(b2, 2) - 3

      term2 = multiply_elem(kappa*(1 - kappa), b2)
      order2 = size(b2, 2) - 3

      call integral(b1, size(b1, 2) - 1, term1, order1, knot, result1)
      call integral(b1, size(b1, 2) - 1, term2, order2, knot, result2)

      result = (-1)*(result1 + result2)/2
   end subroutine matrix_T_minus

   subroutine matrix_W_plus(b1, b2, knot, Z, kappa, C, result)
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      type(mp_real), intent(in) :: Z, kappa, C
      type(mp_real), intent(in) :: knot(:)
      type(mp_real) :: result

      type(mp_real), dimension(size(b1, 1), size(b1, 2)) :: term1, term2, term3, term4
      type(mp_real) :: result1, result2, result3, result4
      integer :: order1, order2, order3, order4

      term1 = multiply_elem(Z, deriv_single(b2, size(b2, 2) - 1))
      order1 = size(b1, 2) - 3

      term2 = multiply_elem(Z*kappa, b2)
      order2 = size(b1, 2) - 3

      term3 = multiply_elem(Z*kappa, deriv_single(b2, size(b2, 2) - 1))
      order3 = size(b1, 2) - 4

      term4 = multiply_elem(Z*kappa*kappa, b2)
      order4 = size(b1, 2) - 4

      call integral(deriv_single(b1, size(b1, 2) - 1), size(b1, 2) - 2, term1, order1, knot, result1)
      call integral(deriv_single(b1, size(b1, 2) - 1), size(b1, 2) - 2, term2, order2, knot, result2)
      call integral(b1, size(b1, 2) - 1, term3, order3, knot, result3)
      call integral(b1, size(b1, 2) - 1, term4, order4, knot, result4)

      result = result1 + result2 + result3 + result4
   end subroutine matrix_W_plus

   subroutine matrix_W_minus(b1, b2, knot, Z, kappa, C, result)
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      type(mp_real), intent(in) :: Z, kappa, C
      type(mp_real), intent(in) :: knot(:)
      type(mp_real) :: result

      type(mp_real), dimension(size(b1, 1), size(b1, 2)) :: term1, term2, term3, term4
      type(mp_real) :: result1, result2, result3, result4
      integer :: order1, order2, order3, order4

      term1 = multiply_elem(Z, deriv_single(b2, size(b2, 2) - 1))
      order1 = size(b1, 2) - 3

      term2 = multiply_elem(-Z*kappa, b2)
      order2 = size(b1, 2) - 3

      term3 = multiply_elem(-Z*kappa, deriv_single(b2, size(b2, 2) - 1))
      order3 = size(b1, 2) - 4

      term4 = multiply_elem(Z*kappa*kappa, b2)
      order4 = size(b1, 2) - 4

      call integral(deriv_single(b1, size(b1, 2) - 1), size(b1, 2) - 2, term1, order1, knot, result1)
      call integral(deriv_single(b1, size(b1, 2) - 1), size(b1, 2) - 2, term2, order2, knot, result2)
      call integral(b1, size(b1, 2) - 1, term3, order3, knot, result3)
      call integral(b1, size(b1, 2) - 1, term4, order4, knot, result4)

      result = result1 + result2 + result3 + result4
   end subroutine matrix_W_minus

   subroutine matrix_A_plus(b1, b2, knot, Z, kappa, C, result)
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      type(mp_real), intent(in) :: Z, kappa, C
      type(mp_real), intent(in) :: knot(:)
      type(mp_real) :: result

      type(mp_real), dimension(size(b1, 1), size(b1, 2)) :: term1, term2, term3, term4
      type(mp_real) :: result1, result2, result3, result4
      integer :: order1, order2, order3, order4

      term1 = multiply_elem(Z, deriv_single(b2, size(b2, 2) - 1))
      order1 = size(b1, 2) - 3

      term2 = multiply_elem(-Z*kappa, b2)
      order2 = size(b1, 2) - 3

      term3 = multiply_elem(Z, b2)
      order3 = size(b1, 2) - 2

      term4 = multiply_elem(Z*kappa, b2)
      order4 = size(b1, 2) - 3

      call integral(b1, size(b1, 2) - 1, term1, order1, knot, result1)
      call integral(b1, size(b1, 2) - 1, term2, order2, knot, result2)
      call integral(deriv_single(b1, size(b1, 2) - 1), size(b1, 2) - 2, term3, order3, knot, result3)
      call integral(b1, size(b1, 2) - 1, term4, order4, knot, result4)

      result = result1 + result2 + result3 + result4
   end subroutine matrix_A_plus

   subroutine matrix_A_minus(b1, b2, knot, Z, kappa, C, result)
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      type(mp_real), intent(in) :: Z, kappa, C
      type(mp_real), intent(in) :: knot(:)
      type(mp_real) :: result

      type(mp_real), dimension(size(b1, 1), size(b1, 2)) :: term1, term2, term3, term4
      type(mp_real) :: result1, result2, result3, result4
      integer :: order1, order2, order3, order4

      term1 = multiply_elem(Z, deriv_single(b2, size(b2, 2) - 1))
      order1 = size(b1, 2) - 3

      term2 = multiply_elem(Z*kappa, b2)
      order2 = size(b1, 2) - 3

      term3 = multiply_elem(Z, b2)
      order3 = size(b1, 2) - 2

      term4 = multiply_elem(-Z*kappa, b2)
      order4 = size(b1, 2) - 3

      call integral(b1, size(b1, 2) - 1, term1, order1, knot, result1)
      call integral(b1, size(b1, 2) - 1, term2, order2, knot, result2)
      call integral(deriv_single(b1, size(b1, 2) - 1), size(b1, 2) - 2, term3, order3, knot, result3)
      call integral(b1, size(b1, 2) - 1, term4, order4, knot, result4)

      result = result1 + result2 + result3 + result4
   end subroutine matrix_A_minus

   subroutine matrix_B_plus(b1, b2, knot, Z, kappa, C, result)
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      type(mp_real), intent(in) :: Z, kappa, C
      type(mp_real), intent(in) :: knot(:)
      type(mp_real) :: result

      type(mp_real), dimension(size(b1, 1), size(b1, 2)) :: term1, term2, term3, term4
      type(mp_real) :: result1, result2, result3, result4
      integer :: order1, order2, order3, order4

      term1 = deriv_single(deriv_single(b2, size(b2, 2) - 1), size(b2, 2) - 2)
      order1 = size(b1, 2) - 3

      term2 = multiply_elem(kappa*(1-kappa), b2)
      order2 = size(b1, 2) - 3

      term3 = multiply_elem(kappa, deriv_single(deriv_single(b2, size(b2, 2) - 1), size(b2, 2) - 2))
      order3 = size(b1, 2) - 4

      term4 = multiply_elem(kappa*kappa*(1-kappa), b2)
      order4 = size(b1, 2) - 4

      call integral(deriv_single(b1, size(b1, 2) - 1), size(b1, 2) - 2, term1, order1, knot, result1)
      call integral(deriv_single(b1, size(b1, 2) - 1), size(b1, 2) - 2, term2, order2, knot, result2)
      call integral(b1, size(b1, 2) - 1, term3, order3, knot, result3)
      call integral(b1, size(b1, 2) - 1, term4, order4, knot, result4)

      result = (result1 + result2 + result3 + result4)/2
   end subroutine matrix_B_plus

   subroutine matrix_B_minus(b1, b2, knot, Z, kappa, C, result)
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      type(mp_real), intent(in) :: Z, kappa, C
      type(mp_real), intent(in) :: knot(:)
      type(mp_real) :: result

      type(mp_real), dimension(size(b1, 1), size(b1, 2)) :: term1, term2, term3, term4
      type(mp_real) :: result1, result2, result3, result4
      integer :: order1, order2, order3, order4

      term1 = deriv_single(deriv_single(b2, size(b2, 2) - 1), size(b2, 2) - 2)
      order1 = size(b1, 2) - 3

      term2 = multiply_elem(-kappa*(1+kappa), b2)
      order2 = size(b1, 2) - 3

      term3 = multiply_elem(-kappa, deriv_single(deriv_single(b2, size(b2, 2) - 1), size(b2, 2) - 2))
      order3 = size(b1, 2) - 4

      term4 = multiply_elem(kappa*kappa*(1+kappa), b2)
      order4 = size(b1, 2) - 4

      call integral(deriv_single(b1, size(b1, 2) - 1), size(b1, 2) - 2, term1, order1, knot, result1)
      call integral(deriv_single(b1, size(b1, 2) - 1), size(b1, 2) - 2, term2, order2, knot, result2)
      call integral(b1, size(b1, 2) - 1, term3, order3, knot, result3)
      call integral(b1, size(b1, 2) - 1, term4, order4, knot, result4)

      result = (result1 + result2 + result3 + result4)/2
   end subroutine matrix_B_minus

   subroutine matrixAB(d, n, n_remove, Z, kappa, C, amin, amax, H_mat, S_mat, log_bool, i2, i3)
      !> @brief Generate the matrix A and B
      !> @param d : integer : the degree of the B-spline
      !> @param n : integer : the number of B-splines
      !> @param n : integer : the number of B-splines to be remove at the start and the end
      !> @param Z : mp_real : the potential constant
      !> @param kappa : mp_real : the relativistic quantum number
      !> @param C : mp_real : the speed of light
      !> @param amin : mp_real : the minimum value of the knot vector
      !> @param amax : mp_real : the maximum value of the knot vector
      !> @param H_mat : mp_real(2*nprime,2*nprime) : the matrix H
      !> @param H_mat : mp_real(2*nprime,2*nprime) : the matrix S
      !> @param log_bool : logical : display the result
      !> @param i2 : integer : the field width
      !> @param i3 : integer : the number of decimal
      implicit none
      integer, intent(in) :: d, n, n_remove
      type(mp_real), intent(in) :: Z, kappa, C, amin, amax
      type(mp_real), intent(out), allocatable, dimension(:, :) :: H_mat, S_mat
      logical, intent(in), optional :: log_bool
      integer, intent(in), optional :: i2, i3

      integer :: nprime

      type(mp_real), dimension(n + d) :: knot
      type(mp_real), dimension(n, size(knot), d) :: bspline
      type(mp_real), dimension(:), allocatable :: tmp
      type(mp_real), dimension(:, :), allocatable :: Ss_mat, Phi_mat, T_plus_mat, T_minus_mat, W_plus_mat, W_minus_mat, A_plus_mat, A_minus_mat, B_plus_mat, B_minus_mat

      integer :: i_tmp, j_tmp

      zero = '0.d0'
      one = '1.d0'

      ! Number of B-splines without the two first and last one
      nprime = n - 2*n_remove

      ! Generate the knot vector
      print *, "Generate the knot vector"
      knot = exp_knot(n + d, amax, amin, d, 0, zero)

      ! Generate the B-splines
      print *, "Generate the B-splines"
      !OMP PARALLEL DO PRIVATE(i_tmp) SHARED(d, knot, bspline, n)
      do i_tmp = 1, n
         call init_bspine(d, i_tmp, knot, bspline(i_tmp, :, :), .false.)
      end do
      !OMP END PARALLEL DO


      print *, "Generate Small S matrix"
      allocate (Ss_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, knot, Z, kappa, C, Ss_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
          call matrix_S(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), knot, Z, kappa, C, Ss_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO


      print *, "Generate Phi matrix"
      allocate (Phi_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, knot, Z, kappa, C, Phi_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
          call matrix_phi(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), knot, Z, kappa, C, Phi_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO


      ! Warning for now cause singularities thus use transpose of H21 instead
      print *, "Generate T+ matrix"
      allocate (T_plus_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, knot, Z, kappa, C, T_plus_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
          call matrix_T_plus(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), knot, Z, kappa, C, T_plus_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO


      ! Warning for now cause singularities thus use transpose of H21 instead
      print *, "Generate T- matrix"
      allocate (T_minus_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, knot, Z, kappa, C, T_minus_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
          call matrix_T_minus(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), knot, Z, kappa, C, T_minus_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO


      print *, "Generate W+ matrix"
      allocate (W_plus_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, knot, Z, kappa, C, W_plus_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
          call matrix_W_plus(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), knot, Z, kappa, C, W_plus_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO


      print *, "Generate W- matrix"
      allocate (W_minus_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, knot, Z, kappa, C, W_minus_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
          call matrix_W_minus(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), knot, Z, kappa, C, W_minus_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO


      print *, "Generate A+ matrix"
      allocate (A_plus_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, knot, Z, kappa, C, A_plus_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
          call matrix_A_plus(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), knot, Z, kappa, C, A_plus_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO


      print *, "Generate A- matrix"
      allocate (A_minus_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, knot, Z, kappa, C, A_minus_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
          call matrix_A_minus(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), knot, Z, kappa, C, A_minus_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO


      print *, "Generate B+ matrix"
      allocate (B_plus_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, knot, Z, kappa, C, B_plus_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
          call matrix_B_plus(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), knot, Z, kappa, C, B_plus_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO


      print *, "Generate B- matrix"
      allocate (B_minus_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, knot, Z, kappa, C, B_minus_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
          call matrix_B_minus(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), knot, Z, kappa, C, B_minus_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO


      ! Generate the matrix H
      print *, "Generate H matrix"
      allocate (H_mat(2*nprime, 2*nprime))
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
            H_mat(i_tmp, j_tmp) = c*c*Ss_mat(i_tmp, j_tmp) + 3*T_plus_mat(i_tmp, j_tmp)/2 - Phi_mat(i_tmp, j_tmp) - W_plus_mat(i_tmp, j_tmp)/(4*c*c)
            H_mat(i_tmp, j_tmp + nprime) = (-A_plus_mat(i_tmp, j_tmp) + B_plus_mat(i_tmp, j_tmp))/(2*c)
            H_mat(i_tmp + nprime, j_tmp) = (-A_minus_mat(i_tmp, j_tmp) - B_minus_mat(i_tmp, j_tmp))/(2*c)
            H_mat(i_tmp + nprime, j_tmp + nprime) = -c*c*Ss_mat(i_tmp, j_tmp) - 3*T_minus_mat(i_tmp, j_tmp)/2 - Phi_mat(i_tmp, j_tmp) - W_minus_mat(i_tmp, j_tmp)/(4*c*c)
         end do
      end do

      ! Generate the matrix S
      print *, "Generate S matrix"
      allocate (S_mat(2*nprime, 2*nprime))
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
            S_mat(i_tmp, j_tmp) = Ss_mat(i_tmp, j_tmp) + T_plus_mat(i_tmp, j_tmp)/(2*c*c)
            S_mat(i_tmp, j_tmp + nprime) = zero
            S_mat(i_tmp + nprime, j_tmp) = zero
            S_mat(i_tmp + nprime, j_tmp + nprime) = Ss_mat(i_tmp, j_tmp) + T_minus_mat(i_tmp, j_tmp)/(2*c*c)
         end do
      end do

      if (present(log_bool) .and. log_bool) then
         print *, "Writing Logs"
         open (1, file="log_2.log", status="replace")

         allocate (tmp(2*nprime))

         write (1, '(a,i4)') "Number of BSplines: ", n
         write (1, '(a,i4)') "Order of BSplines: ", d, " and number of knots: ", size(knot)
         write (1, '(a)') "Knots: "
         call write_lists(knot, 1, i2, i3)
         write (1, '(a)') "----------------------------------------------------------------"
         write (1, '(a)') "Matrix H : "
         do i_tmp = 1, 2*nprime
            call write_lists(H_mat(i_tmp, :), 1, i2, i3)
         end do
         write (1, '(a)') "----------------------------------------------------------------"
         write (1, '(a)') "Check H symmetry : "
         do i_tmp = 1, 2*nprime
            do j_tmp = 1, 2*nprime
               tmp(j_tmp) = H_mat(j_tmp, i_tmp) - H_mat(i_tmp, j_tmp)
            end do
            call write_lists(tmp, 1, i2, i3)
         end do
         write (1, '(a)') "----------------------------------------------------------------"
         write (1, '(a)') "Matrix S : "
         do i_tmp = 1, 2*nprime
            call write_lists(S_mat(i_tmp, :), 1, i2, i3)
         end do
         close (1)
         print *, "Logs written"
      end if
   end subroutine matrixAB

   subroutine get_eigen(d, n, n_remove, Z, kappa, C, amin, amax, log_bool, i2, i3)
      integer, intent(in) :: d, n, n_remove, i2, i3
      type(mp_real), intent(in) :: Z, kappa, C, amin, amax
      logical, intent(in), optional :: log_bool

      integer :: nprime, i_tmp, ierr, method, solnum
      type(mp_real) :: clt
      type(mp_real), dimension(:, :), allocatable :: A, B, vect
      type(mp_real), dimension(:), allocatable :: w, fv1, fv2
      character(len=100) :: log_file

      zero = '0.d0'
      one = '1.d0'

      nprime = n - 2*n_remove
      method = 0
      clt = '0.5d0'

      allocate (A(2*nprime, 2*nprime))
      allocate (B(2*nprime, 2*nprime))

      call matrixAB(d, n, n_remove, Z, kappa, C, amin, amax, A, B, log_bool, i2, i3)

      print *, "Calculating the eigenvalues"

      allocate (w(2*nprime))
      allocate (vect(2*nprime, 2*nprime))
      allocate (fv1(2*nprime))
      allocate (fv2(2*nprime))

      call rsg(2*nprime, 2*nprime, A, B, w, 0, vect, fv1, fv2, ierr)

      print *, "Error code: ", ierr

      write (log_file, '(a,I4,a,I2,a)') "./result_DKB/error_", n, "_", d, "_1.txt"

      open (2, file=log_file)

      write (2, '(a, i4, a, i4, a, i4)') "Number of BSplines: ", n, ", Order of BSplines: ", d, ", Removed Bsplines", n_remove
      write (2, '(a)') "Speed of light: "
      call mpwrite(2, i2, i3, C)
      write (2, '(a)') "Relativistic quantum number: "
      call mpwrite(2, i2, i3, kappa)
      write (2, '(a)') "Potential constant: "
      call mpwrite(2, i2, i3, Z)
      write (2, '(a)') "Knot min: "
      call mpwrite(2, i2, i3, amin)
      write (2, '(a)') "Knot max: "
      call mpwrite(2, i2, i3, amax)
      write (2, '(a)') "--------------------------------------------------------------"

      solnum = 1
      do i_tmp = 1, 2*nprime
         if (w(i_tmp) - C**2 < zero .AND. abs(w(i_tmp) - c**2) < 1.d3*one) then
            call mpwrite(2, i2, i3, abs(w(i_tmp) - theoric_val(solnum, Z, kappa, C)))
            solnum = solnum + 1
         end if
      end do

      close (2)

      print *, "Errors written to ", log_file

      write (log_file, '(a,I4,a,I2,a)') "./result_DKB/eigenvalues_", n, "_", d, ".txt"

      open (2, file=log_file)

      write (2, '(a, i4, a, i4, a, i4)') "Number of BSplines: ", n, ", Order of BSplines: ", d, ", Removed Bsplines", n_remove
      write (2, '(a)') "Speed of light: "
      call mpwrite(2, i2, i3, C)
      write (2, '(a)') "Relativistic quantum number: "
      call mpwrite(2, i2, i3, kappa)
      write (2, '(a)') "Potential constant: "
      call mpwrite(2, i2, i3, Z)
      write (2, '(a)') "Knot min: "
      call mpwrite(2, i2, i3, amin)
      write (2, '(a)') "Knot max: "
      call mpwrite(2, i2, i3, amax)
      write (2, '(a)') "--------------------------------------------------------------"

      solnum = 1
      do i_tmp = 1, 2*nprime
         if (w(i_tmp) - C**2 < zero .AND. abs(w(i_tmp) - c**2) < 1.d3*one) then
            call mpwrite(2, i2, i3, w(i_tmp) - C**2)
            solnum = solnum + 1
         end if
      end do

      close (2)

      deallocate (A, B, w, vect, fv1, fv2)

      print *, "Eigenvalues written to ", log_file

   end subroutine get_eigen
end module bspline_mod
