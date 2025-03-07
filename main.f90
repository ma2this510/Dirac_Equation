program main
   use bspline_mod
   implicit none

   integer :: d, n, i_tmp, nprime
   real :: amin, amax, Z, C, kappa
   real, dimension(:, :), allocatable :: A, B
   real, dimension(:), allocatable :: w, work

   integer :: itype, lda, ldb, lwork, info
   character(1) :: jobz, uplo


   !----------------------------------------------------------------
   ! Define Important Variables
   d = 7 ! Order of Mathemathica + 1
   n = 100 ! Number of B-Splines
   Z = 2.0
   C = 137
   kappa = -1.0
   amin = 0.002
   amax = 40.0

   nprime = n - 2
   allocate(A(2*nprime, 2*nprime), B(2*nprime, 2*nprime), w(2*nprime))

   !----------------------------------------------------------------
   ! Generate A and B matrices
   call matrixAB(d, n, Z, kappa, C, amin, amax, A, B, .true.)

!    do i_tmp = 1, 2*n
!       print *, A(i_tmp, 1:2*n)
!    end do

!    do i_tmp = 1, 2*n
!       print *, B(i_tmp, 1:2*n)
!    end do

   !----------------------------------------------------------------
   ! Get Eigenvalues

   ! LAPACK parameters
   itype = 1      ! Solves A*x = lambda*B*x
   jobz = 'N'     ! Compute eigenvalues and eigenvectors
   uplo = 'U'     ! Upper triangular part is stored
   lda = 2*nprime      ! Leading dimension of A
   ldb = 2*nprime      ! Leading dimension of B

   ! Get optimal workspace size
   lwork = -1
   allocate(work(1))
   call dsygv(itype, jobz, uplo, 2*nprime, A, lda, B, ldb, W, work, lwork, info)

   ! Compute Eigenvalues
   lwork = int(work(1))
   deallocate(work)
   allocate(work(lwork))
   call dsygv(itype, jobz, uplo, 2*nprime, A, lda, B, ldb, W, work, lwork, info)

   ! Save Eigenvalues
   open(1, file='eigenvalues.dat', status='replace')
   do i_tmp = 1, 2*nprime
      write(1, '(f25.8)') W(i_tmp)
   end do
   close(1)

   deallocate(A, B, w, work)

end program main
