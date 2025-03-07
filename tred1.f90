      subroutine tred1(nm,n,a,d,e1,e2)
USE mpmodule
      implicit none
      integer           n,nm
      type (mp_real)    a(nm,n),d(n),e1(n),e2(n)
!
!     this subroutine is a translation of the algol procedure TRED1,
!     Num. Math. 11, 181-195(1968) by Martin, Reinsch, and Wilkinson.
!     Handbook for Auto. Comp., vol.II-Linear Algebra, 212-226(1971).
!
!     this subroutine reduces a real symmetric matrix
!     to a symmetric tridiagonal matrix using
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
!        a contains information about the orthogonal trans-
!          formations used in the reduction in its strict lower
!          triangle.  the full upper triangle of a is unaltered.
!
!        d contains the diagonal elements of the tridiagonal matrix.
!
!        e contains the subdiagonal elements of the tridiagonal
!          matrix in its last n-1 positions.  e1(1) is set to zero.
!
!        e2 contains the squares of the corresponding elements of e.
!          e2 may coincide with e if the squares are not needed.
!
!     questions and comments should be directed to Burton S. Garbow,
!     Mathematics and Computer Science Div, Argonne National Laboratory
!
!     this version dated November 2000.
!
      integer           i,j,k
      type (mp_real)    f,g,h, scl
      type (mp_real)    zero
!
      zero = '0.d0'
!
      do i=1,n
         d(i) = a(n,i)
         a(n,i) = a(i,i)
      end do
!
      do i=n,1,-1
         h = zero
         scl = zero
!     .......... scale row (algol tol then not needed) ..........
         do k=1,i-1
            scl = scl + abs(d(k))
         end do
!
         if (scl == zero) then
!
            do j=1,i-1
               d(j) = a(i-1,j)
               a(i-1,j) = a(i,j)
               a(i,j) = zero
            end do
!
            e1(i) = zero
            e2(i) = zero
            cycle
         end if
!
         do k=1,i-1
            d(k) = d(k) / scl
            h = h + d(k) * d(k)
         end do
!
         e2(i) = scl * scl * h
         f = d(i-1)
         g = -sign(sqrt(h),f)
         e1(i) = scl * g
         h = h - f * g
         d(i-1) = f - g
         if (i == 2) goto 200
!     .......... form a*u ..........
         do j=1,i-1
            e1(j) = zero
         end do
!
         do j=1,i-1
            f = d(j)
            g = e1(j) + a(j,j) * f
            do k=j+1,i-1
               g = g + a(k,j) * d(k)
               e1(k) = e1(k) + a(k,j) * f
            end do
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
         h = f / (h + h)
!     .......... form q ..........
         do j=1,i-1
            e1(j) = e1(j) - h * d(j)
         end do
!     .......... form reduced a ..........
         do j=1,i-1
            f = d(j)
            g = e1(j)
!
            do k=j,i-1
               a(k,j) = a(k,j) - f * e1(k) - g * d(k)
            end do
!
         end do
!
  200    do j=1,i-1
            f = d(j)
            d(j) = a(i-1,j)
            a(i-1,j) = a(i,j)
            a(i,j) = f * scl
         end do
!
      end do
!
      return
      end
