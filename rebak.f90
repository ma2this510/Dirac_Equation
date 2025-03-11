      subroutine rebak(nm,n,b,dl,m,z)
USE mpmodule
USE omp_lib
      implicit none
      integer           m,n,nm
      type (mp_real)    b(nm,n),dl(n),z(nm,m)
      type (mp_real), allocatable   ::  vaux(:)
!
!     this subroutine is a translation of the algol procedure rebaka,
!     num. math. 11, 99-110(1968) by martin and wilkinson.
!     handbook for auto. comp., vol.ii-linear algebra, 303-314(1971).
!
!     this subroutine forms the eigenvectors of a generalized
!     symmetric eigensystem by back transforming those of the
!     derived symmetric matrix determined by  reduc.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix system.
!
!        b contains information about the similarity transformation
!          (cholesky decomposition) used in the reduction by  reduc
!          in its strict lower triangle.
!
!        dl contains further information about the transformation.
!
!        m is the number of eigenvectors to be back transformed.
!
!        z contains the eigenvectors to be back transformed
!          in its first m columns.
!
!     on output
!
!        z contains the transformed eigenvectors
!          in its first m columns.
!
!     questions and comments should be directed to burton s. garbow,
!     mathematics and computer science div, argonne national laboratory
!
!     this version dated November 2000.
!
      integer           i,j,k,nnum,mnum,it
      type (mp_real)    x,zero,red
!
      zero = '0.d0'
!
      nnum = 1
      !$omp parallel
      !$ mnum = omp_get_thread_num()
      !$ if (mnum == 1) then
      !$    nnum = omp_get_num_threads()
      !$ end if 
      !$omp end parallel
      allocate(vaux(nnum))
!
      do j=1,m
         do i=n,1,-1
            x = z(i,j)

            red = zero
            vaux(1:nnum) = red
            mnum = 0
            !$OMP parallel private(mnum)
            !$ mnum = omp_get_thread_num()
            mnum = mnum + 1
            !$omp do
            do k=i+1,n
               vaux(mnum) = vaux(mnum) + b(k,i) * z(k,j)
            end do
            !$OMP end do 
            !$OMP end parallel 
            do it=1,nnum
               red=red+vaux(it)
            end do
            x = x - red

            z(i,j) = x / dl(i)
         end do
      end do
!
      return
      end
