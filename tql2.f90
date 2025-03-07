      subroutine tql2(nm,n,d,e1,z,ierr)
USE mpmodule
USE omp_lib
      implicit none
      integer           n,nm,ierr
      type (mp_real)    d(n),e1(n),z(nm,n)
!
!     this subroutine is a translation of the algol procedure TQL2,
!     Num. Math. 11, 293-306(1968) by Bowdler, Martin, Reinsch, and
!     Wilkinson.
!     Handbook for Auto. Comp., vol.II-Linear Algebra, 227-240(1971).
!
!     this subroutine finds the eigenvalues and eigenvectors
!     of a symmetric tridiagonal matrix by the ql method.
!     the eigenvectors of a full symmetric matrix can also
!     be found if  tred2  has been used to reduce this
!     full matrix to tridiagonal form.
!
!     on input
!
!        nm must be set to the row dimension of two-dimensional
!          array parameters as declared in the calling program
!          dimension statement.
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e contains the subdiagonal elements of the input matrix
!          in its last n-1 positions.  e1(1) is arbitrary.
!
!        z contains the transformation matrix produced in the
!          reduction by  tred2, if performed.  if the eigenvectors
!          of the tridiagonal matrix are desired, z must contain
!          the identity matrix.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct but
!          unordered for indices 1,2,...,ierr-1.
!
!        e has been destroyed.
!
!        z contains orthonormal eigenvectors of the symmetric
!          tridiagonal (or full) matrix.  if an error exit is made,
!          z contains the eigenvectors associated with the stored
!          eigenvalues.
!
!        ierr is set to
!          zero       for normal return,
!          j          if the j-th eigenvalue has not been
!                     determined after 30 iterations.
!
!     calls pythag for  sqrt(a*a+b*b) .
!
!     questions and comments should be directed to Burton S. Garbow,
!     Mathematics and Computer Science Div, Argonne National Laboratory
!
!     this version dated November 2000.
!
      integer           i,j,k,l,m
      integer        :: jmax = 30
      type (mp_real)    c,c2,c3,dl1,el1,f,g,h,p,r,s,s2,tst1,tst2,pythag,eps
      type (mp_real)    zero, one, two, epsilonn
!
      zero = '0.d0'
      one = '1.d0'
      two = '2.d0'
      !eps = '1.d-46'!epsilonn(one)
      eps = epsilonn(one)
!
      ierr = 0
      if (n == 1) return
!
      do i=2,n
         e1(i-1) = e1(i)
      end do
!
      f = zero
      tst1 = zero
      e1(n) = zero
!
      do l=1,n
         j = 0
         h = abs(d(l)) + abs(e1(l))
         if (tst1 < h) tst1 = h
!     .......... look for small sub-diagonal element ..........
         do m=l,n
            tst2 = tst1 + abs(e1(m))
            if (tst2-tst1 < tst1*eps) exit
!     .......... e1(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
         end do
!
         if (m /= l) then
            do  ! cycle of QL transformations for the interval [l,m]
               if (j == jmax) goto 1000
               j = j + 1
!     .......... form shift ..........
               g = d(l)
               p = (d(l+1) - g) / (two * e1(l))
               r = pythag(p,one)
               d(l) = e1(l) / (p + sign(r,p))
               d(l+1) = e1(l) * (p + sign(r,p))
               dl1 = d(l+1)
               h = g - d(l)
               do i=l+2,n
                  d(i) = d(i) - h
               end do
               f = f + h
!     .......... ql transformation ..........
               p = d(m)
               c = one
               c2 = c
               el1 = e1(l+1)
               s = zero
               do i=m-1,l,-1 
                  c3 = c2
                  c2 = c
                  s2 = s
                  g = c * e1(i)
                  h = c * p
                  r = pythag(p,e1(i))
                  e1(i+1) = s * r
                  s = e1(i) / r
                  c = p / r
                  p = c * d(i) - s * g
                  d(i+1) = h + s * (c * g + s * d(i))
!     .......... form vector ..........
                  !$omp parallel do private(h) schedule(dynamic,1)
                  do k=1,n
                     h = z(k,i+1)
                     z(k,i+1) = s * z(k,i) + c * h
                     z(k,i) = c * z(k,i) - s * h
                  end do
                  !$omp end parallel do
               end do
!
               p = -s * s2 * c3 * el1 * e1(l) / dl1
               e1(l) = s * p
               d(l) = c * p
               tst2 = tst1 + abs(e1(l))
               if (tst2-tst1 < tst1*eps) exit
            end do
         end if
         d(l) = d(l) + f
      end do
!     .......... order eigenvalues and eigenvectors ..........
      do i=1,n-1
         k = i
         p = d(i)
!
         do j=i+1,n
            if (d(j) >= p) cycle
            k = j
            p = d(j)
         end do
!
         if (k == i) cycle
         d(k) = d(i)
         d(i) = p
!
         do j=1,n
            p = z(j,i)
            z(j,i) = z(j,k)
            z(j,k) = p
         end do
!
      end do
!
      return
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
      return
      end
