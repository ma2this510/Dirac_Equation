      subroutine reduc(nm,n,a,b,dl,ierr)
USE mpmodule
USE omp_lib
      implicit none
      integer           nm,n,ierr
      type (mp_real)    a(nm,n),b(nm,n),dl(n)
      type (mp_real), allocatable   ::  vaux(:),bl(:,:)
     

!
!     this subroutine is a translation of the algol procedure reduc1,
!     Num. Math. 11, 99-110(1968) by Martin and Wilkinson.
!     Handbook for Auto. Comp., vol.II-Linear Algebra, 303-314(1971).
!
!     this subroutine reduces the generalized symmetric eigenproblem
!     ax=(lambda)bx, where b is positive definite, to the standard
!     symmetric eigenproblem using the cholesky factorization of b.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrices a and b.  if the cholesky
!          factor l of b is already available, n should be prefixed
!          with a minus sign.
!
!        a and b contain the real symmetric input matrices.  only the
!          full upper triangles of the matrices need be supplied.  if
!          n is negative, the strict lower triangle of b contains,
!          instead, the strict lower triangle of its cholesky factor l.
!
!        dl contains, if n is negative, the diagonal elements of l.
!
!     on output
!
!        a contains in its full lower triangle the full lower triangle
!          of the symmetric matrix derived from the reduction to the
!          standard form.  the strict upper triangle of a is unaltered.
!
!        b contains in its strict lower triangle the strict lower
!          triangle of its cholesky factor l.  the full upper
!          triangle of b is unaltered.
!
!        dl contains the diagonal elements of l.
!
!        ierr is set to
!          zero       for normal return,
!          7*n+1      if b is not positive definite.
!
!     questions and comments should be directed to Burton S. Garbow,
!     Mathematics and Computer Science Div, Argonne National Laboratory
!
!     this version dated November 2000.
!
      integer           i,j,k,nn,nnum,mnum,it
      type (mp_real)    x,y,z,zero,one,red
      real              xrr
!
      zero = '0.d0'
      one = '1.d0'
!
      ierr = 0
      nn = abs(n)

      nnum = 1
      !$omp parallel
      !$ mnum = omp_get_thread_num()
      !$ if (mnum == 1) then
      !$    nnum = omp_get_num_threads()
      !$ end if 
      !$omp end parallel
      allocate(vaux(nnum),bl(nm,n))

      do i=1,n
         do j=1,n
            bl(i,j) = b(i,j)
         end do
      end do



!!     .......... form l in the arrays b and dl ..........


      do j=1,n
         
         x = zero
         do k=1,j-1
            x = x + bl(j,k)*bl(j,k)
         end do

         !bl(j,j) = sqrt(b(j,j) - x)
         bl(j,j) = sqrt(abs(b(j,j) - x))
         
         !$OMP parallel private(x)
         !$OMP do schedule(dynamic,1)
         do i=j+1,n
            x = zero
            do k=1,j-1
               x = x + bl(i,k) * bl(j,k)
            end do
            bl(i,j) = (one/bl(j,j))*(b(i,j) - x) 
         end do
         !$OMP end do 
         !$OMP end parallel 

      end do

      do i=1,n
         dl(i) = bl(i,i)
         do j=1,i-1
            b(i,j) = bl(i,j)
         end do
      end do



!     .......... form the transpose of the upper triangle of inv(l)*a
!                in the lower triangle of the array a ..........
      do i=1,nn
         y = dl(i)
         do j=i,nn
            x = a(i,j)

            red = zero
            vaux(1:nnum) = red
            mnum = 0
            !$OMP parallel private(mnum)
            !$ mnum = omp_get_thread_num()
            mnum = mnum + 1
            !$OMP do 
            do k=1,i-1
               vaux(mnum) = vaux(mnum) + b(i,k) * a(j,k)
            end do
            !$OMP end do 
            !$OMP end parallel 
            do it=1,nnum
               red=red+vaux(it)
            end do
            x = x - red

            a(j,i) = x / y
         end do
      end do
!     .......... pre-multiply by inv(l) and overwrite ..........
      do j=1,nn
         do i=j,nn
            x = a(i,j)

            red = zero
            vaux(1:nnum) = red
            mnum = 0
            !$OMP parallel private(mnum)
            !$ mnum = omp_get_thread_num()
            mnum = mnum + 1
            !$OMP do 
            do k=j,i-1
               vaux(mnum) = vaux(mnum) + a(k,j) * b(i,k)
            end do
            !$OMP end do 
            !$OMP end parallel 
            do it=1,nnum
               red=red+vaux(it)
            end do
            x = x - red
            
            red = zero
            vaux(1:nnum) = red
            mnum = 0
            !$OMP parallel private(mnum)
            !$ mnum = omp_get_thread_num()
            mnum = mnum + 1
            !$OMP do 
            do k=1,j-1
               vaux(mnum) = vaux(mnum) + a(j,k) * b(i,k)
            end do
            !$OMP end do 
            !$OMP end parallel 
            do it=1,nnum
               red=red+vaux(it)
            end do
            x = x - red
            
            a(i,j) = x / dl(i)
         end do
      end do
!
      return
      end
