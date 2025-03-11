      subroutine tred2(nm,n,a,d,e1,z)
USE mpmodule
USE omp_lib
      implicit none
      integer           n,nm
      type (mp_real)    a(nm,n),d(n),e1(n),z(nm,n)
      type (mp_real), allocatable   ::  vaux(:)
!
!     this subroutine is a translation of the algol procedure TRED2,
!     Num. Math. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
!     Handbook for Auto. Comp., vol.II-Linear Algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix to a
!     symmetric tridiagonal matrix using and accumulating
!     orthogonal similarity transformations.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        a contains the real symmetric input matrix.  only the
!          lower triangle of the matrix need be supplied.
!
!     on output
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e1(1) is set to zero.
!
!        z contains the orthogonal transformation matrix
!          produced in the reduction.
!
!        a and z may coincide.  if distinct, a is unaltered.
!
!     questions and comments should be directed to Burton S. Garbow,
!     Mathematics and Computer Science Div, Argonne National Laboratory
!
!     this version dated November 2000.
!
!     ------------------------------------------------------------------
!
      integer               i,j,k,nnum,mnum,it
      type (mp_real)        f,g,h,hh,scl
      type (mp_real)        zero,one,red
!

      nnum = 1
      !$omp parallel
      !$ mnum = omp_get_thread_num()
      !$ if (mnum == 1) then
      !$    nnum = omp_get_num_threads()
      !$ end if 
      !$omp end parallel
      allocate(vaux(nnum))

      zero = '0.d0'
      one = '1.d0'
!
      do i=1,n
         do j=i,n
            z(j,i) = a(j,i)
         end do
         d(i) = a(n,i)
      end do
!
      do i=n,2,-1
         h = zero
         scl = zero
         if (i == 2) go to 130
!     .......... scale row (algol tol then not needed) ..........
         do k=1,i-1
            scl = scl + abs(d(k))
         end do
!
         if (scl /= zero) goto 140
  130    e1(i) = d(i-1)
!
         do j=1,i-1
            d(j) = z(i-1,j)
            z(i,j) = zero
            z(j,i) = zero
         end do
!
         go to 290
!
  140    do k=1,i-1
            d(k) = d(k) / scl
            h = h + d(k) * d(k)
         end do
!
         f = d(i-1)
         g = -sign(sqrt(h),f)
         e1(i) = scl * g
         h = h - f * g
         d(i-1) = f - g
!     .......... form a*u ..........
         do j=1,i-1
            e1(j) = zero
         end do
!
         do j=1,i-1
            f = d(j)
            z(j,i) = f
            g = e1(j) + z(j,j) * f

            red = zero
            vaux(1:nnum) = red
            mnum = 0
            !$OMP parallel private(mnum)
            !$ mnum = omp_get_thread_num()
            mnum = mnum + 1
            !$OMP do 
            do k=j+1,i-1
               vaux(mnum) = vaux(mnum) + z(k,j) * d(k)
               e1(k) = e1(k) + z(k,j) * f
            end do
            !$OMP end do 
            !$OMP end parallel 
            do it=1,nnum
               red=red+vaux(it)
            end do
            g = g + red
            
            e1(j) = g
         end do
!     .......... form p ..........
         f = zero
!
         do j=1,i-1
            e1(j) = e1(j) / h
            f = f + e1(j) * d(j)
         end do
!
         hh = f / (h + h)
!     .......... form q ..........
         do j=1,i-1
            e1(j) = e1(j) - hh * d(j)
         end do
!     .......... form reduced a ..........
         do j=1,i-1
            f = d(j)
            g = e1(j)
!
            !$omp parallel do 
            do k=j,i-1
               z(k,j) = z(k,j) - f * e1(k) - g * d(k)
            end do
            !$omp end parallel do
!
            d(j) = z(i-1,j)
            z(i,j) = zero
         end do
!
  290    d(i) = h
      end do
!     .......... accumulation of transformation matrices ..........
      do i=2,n
         z(n,i-1) = z(i-1,i-1)
         z(i-1,i-1) = one
         h = d(i)
         if (h /= zero) then
!
            do k=1,i-1
               d(k) = z(k,i) / h
            end do
!
            do j=1,i-1
               g = zero

               vaux(1:nnum) = zero
               mnum = 0
               !$OMP parallel private(mnum)
               !$ mnum = omp_get_thread_num()
               mnum = mnum + 1
               !$OMP do 
               do k=1,i-1
                  vaux(mnum) = vaux(mnum) + z(k,i) * z(k,j)
               end do
               !$OMP end do 
               !$OMP end parallel 
               do it=1,nnum
                  g=g+vaux(it)
               end do
            
               !$omp parallel do 
               do k=1,i-1
                  z(k,j) = z(k,j) - g * d(k)
               end do
               !$omp end parallel do
            end do
!
         end if
!
         do k=1,i-1
            z(k,i) = zero
         end do
!
      end do
!
      do i=1,n
         d(i) = z(n,i)
         z(n,i) = zero
      end do
!
      z(n,n) = one
      e1(1) = zero
!
      return
      end
