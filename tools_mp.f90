module tools_mp
   use mpmodule
   implicit none

   private 

   public :: write_lists
   public :: multiply_elem

   contains

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

   function multiply_elem(scal, mat) result(result)
      !> @brief Multiply a scalar with a matrix
      !> @param scal : mp_real : the scalar
      !> @param mat : mp_real(:,:) : the matrix
      !> @return result : mp_real(:,:) : the result of the multiplication
      implicit none
      type(mp_real), intent(in) :: scal
      type(mp_real), intent(in), dimension(:,:) :: mat
      type(mp_real), dimension(size(mat, 1), size(mat, 2)) :: result
      
      integer :: i_tmp, j_tmp

      do i_tmp = 1, size(mat, 1)
         do j_tmp = 1, size(mat, 2)
            result(i_tmp, j_tmp) = scal * mat(i_tmp, j_tmp)
         end do
      end do

   end function multiply_elem


end module tools_mp