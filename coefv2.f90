module bspline_mod
   use mpmodule
   use bspline_gen
   use tools_mp
   implicit none

   type(mp_real), save :: zero, one

   private

   public :: int_S
   public :: exp_knot
   public :: int_V
   public :: int_T
   public :: matrixAB

contains

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
                  int_final(i_tmp, j_tmp) = produit(i_tmp, j_tmp)*(knot(i_tmp + 1)**(order_prod - j_tmp + 2))/(order_prod - j_tmp + 2) ! +2 because of the integral that add one order
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

   subroutine int_S(b1, b2, knot, result)
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

   end subroutine int_S

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

   subroutine int_V(b1, b2, Z, knot, result)
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

   end subroutine int_V

   function deriv_double(b) result(result)
      !> @brief Calculate the second derivative of a B-spline
      !> @Warning : doesn't work when divided by r
      !> @param b : real(:) : the coef of the B-spline
      !> @return result : real(:) : the coef of the second derivative of the B-spline
      implicit none
      type(mp_real), intent(in) :: b(:, :)
      type(mp_real) :: result(size(b, 1), size(b, 2))

      integer :: i_tmp, j_tmp

      result = zero
      result(:, 3:size(b, 2)) = b(:, 1:size(b, 2) - 2)

      do i_tmp = 1, size(b, 1)
         do j_tmp = 1, size(b, 2)
            result(i_tmp, j_tmp) = result(i_tmp, j_tmp)*(size(b, 2) - j_tmp + 2)*(size(b, 2) - j_tmp + 1)
            ! print *, size(b, 2) - j_tmp + 2, size(b, 2) - j_tmp + 1
         end do
      end do
   end function deriv_double

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

   subroutine int_T(b1, b2, knot, result)
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      type(mp_real), intent(in) :: knot(:)
      type(mp_real), intent(out) :: result

      zero = '0.d0'
      one = '1.d0'

      call int_S(b1, deriv_double(b2), knot, result)

      result = -result/2

   end subroutine int_T

   subroutine int_Wdiag(b1, b2, Z, kappa, knot, result)
      implicit none
      type(mp_real), intent(in), dimension(:, :) :: b1, b2
      type(mp_real), intent(in) :: Z, kappa
      type(mp_real), intent(in) :: knot(:)
      type(mp_real), intent(out) :: result

      type(mp_real), dimension(size(b1, 1), size(b1, 2)) :: term1, term2, term3, term4
      type(mp_real) :: result1, result2, result3, result4
      integer :: order1, order2, order3, order4

      zero = '0.d0'
      one = '1.d0'

      term1 = deriv_single(deriv_single(b2, size(b2, 2) - 1), size(b2, 2) - 3) ! dr (dr(b2)/r)
      order1 = size(b2, 2) - 4 ! order init - 3 because divided by r and two derivative

      term2 = deriv_single(multiply_elem((1 + kappa), b2), size(b2, 2) - 3) ! dr ((1+kappa)*b2/r²)
      order2 = size(b2, 2) - 4 ! order init - 3 because divided by r² and one derivative

      term3 = multiply_elem((1 + kappa), deriv_single(b2, size(b2, 2) - 1)) ! (1+kappa)/r² * dr(b2)
      order3 = size(b2, 2) - 4 ! order init - 3 because divided by r² and one derivative

      term4 = multiply_elem((1 + kappa)**2, b2) ! (1+kappa)² * b2/r³
      order4 = size(b2, 2) - 4 ! order init - 3 because divided by r³

      call integral(b1, size(b1, 2) - 1, term1, order1, knot, result1)
      call integral(b1, size(b1, 2) - 1, term2, order2, knot, result2)
      call integral(b1, size(b1, 2) - 1, term3, order3, knot, result3)
      call integral(b1, size(b1, 2) - 1, term4, order4, knot, result4)

      result = Z*(result1 + result2 + result3 + result4)

   end subroutine int_Wdiag

   ! For WLS
   ! term1 = deriv_single(b2, size(b2, 2) - 2) ! dr (b2/r)
   ! order1 = size(b2, 2) - 3 ! order init - 2 because divided by r and derivative

   ! term2 = deriv_single(b2, size(b2, 2) - 1) ! dr(b2)/r
   ! order2 = size(b2, 2) - 3 ! order init - 2 because derivative and divided by r

   ! term3 = deriv_single(deriv_single(deriv_single(b2, size(b2, 2)-1), size(b2, 2)-2), size(b2, 2)-3) ! d^3r(b2)/r^3
   ! order3 = size(b2, 2) - 4 ! order init - 3 because three derivative

   ! term4 = deriv_single(deriv_single(b2, size(b2, 2)-2), size(b2, 2)-3) ! d^2r((1+k)*b2/r)
   ! order4 = size(b2, 2) - 4 ! order init - 3 because two derivative and divided by r

   subroutine matrixAB(d, n, n_remove, Z, kappa, C, amin, amax, A, B, log_bool, i2, i3)
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
      integer, intent(in) :: d, n, n_remove
      type(mp_real), intent(in) :: Z, kappa, C, amin, amax
      type(mp_real), intent(out), allocatable, dimension(:, :) :: A, B
      logical, intent(in), optional :: log_bool
      integer, intent(in), optional :: i2, i3

      integer :: nprime

      type(mp_real), dimension(n + d) :: knot
      type(mp_real), dimension(n, size(knot), d) :: bspline

      type(mp_real), dimension(:, :), allocatable :: ovrlp_mat, deriv_mat, potential_mat, k_mat, wdiag_mat
      type(mp_real), allocatable, dimension(:, :) :: aprime_mat

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

      ! Generate the overlap matrix
      print *, "Generate S matrix"
      allocate (ovrlp_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, knot, ovrlp_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
            call int_S(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), knot, ovrlp_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO

      ! Generate the derivative matrix
      print *, "Generate T matrix"
      allocate (deriv_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, knot, deriv_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
            call int_T(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), knot, deriv_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO

      ! Generate the potential matrix
      print *, "Generate V matrix"
      allocate (potential_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, Z, knot, potential_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
            call int_V(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), Z, knot, potential_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO

      ! Generate the matrix Wdiag
      print *, "Generate Wdiag matrix"
      allocate (wdiag_mat(nprime, nprime))
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(i_tmp, j_tmp) SHARED(bspline, Z, kappa, knot, Wdiag_mat, nprime)
      do i_tmp = 1, nprime
         do j_tmp = 1, nprime
           call int_Wdiag(bspline(i_tmp + n_remove, :, :), bspline(j_tmp + n_remove, :, :), Z, kappa, knot, wdiag_mat(i_tmp, j_tmp))
         end do
      end do
      !$OMP END PARALLEL DO

      if (present(log_bool) .and. log_bool) then
         print *, "Writing Logs"
         open (1, file="log_2.txt", status="replace")

         write (1, '(a,i4)') "Number of BSplines: ", n
         write (1, '(a,i4)') "Order of BSplines: ", d, " and number of knots: ", size(knot)
         write (1, '(a)') "Knots: "
         call write_lists(knot, 1, i2, i3)
         write (1, '(a)') "----------------------------------------------------------------"
         write (1, '(a)') "BSplines S Matrix : "
         do i_tmp = 1, nprime
            call write_lists(ovrlp_mat(i_tmp, :), 1, i2, i3)
         end do
         write (1, '(a)') "----------------------------------------------------------------"
         write (1, '(a)') "BSplines T Matrix : "
         do i_tmp = 1, nprime
            call write_lists(deriv_mat(i_tmp, :), 1, i2, i3)
         end do
         write (1, '(a)') "----------------------------------------------------------------"
         write (1, '(a)') "BSplines V Matrix : "
         do i_tmp = 1, nprime
            call write_lists(potential_mat(i_tmp, :), 1, i2, i3)
         end do
         write (1, '(a)') "----------------------------------------------------------------"
         write (1, '(a)') "BSplines Wdiag Matrix : "
         do i_tmp = 1, nprime
            call write_lists(wdiag_mat(i_tmp, :), 1, i2, i3)
         end do
         close (1)
         print *, "Logs written"
      end if
   end subroutine matrixAB
end module bspline_mod
