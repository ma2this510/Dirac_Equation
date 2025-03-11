      subroutine tqlrat(n,d,e2,ierr)
USE mpmodule
      implicit none
      integer           n,ierr
      type (mp_real)    d(n),e2(n)
!
!     this subroutine is a translation of the algol procedure TQLRAT,
!     algorithm 464, Comm. ACM 16, 689(1973) by Reinsch.
!
!     this subroutine finds the eigenvalues of a symmetric
!     tridiagonal matrix by the rational ql method.
!
!     on input
!
!        n is the order of the matrix.
!
!        d contains the diagonal elements of the input matrix.
!
!        e2 contains the squares of the subdiagonal elements of the
!          input matrix in its last n-1 positions.  e2(1) is arbitrary.
!
!      on output
!
!        d contains the eigenvalues in ascending order.  if an
!          error exit is made, the eigenvalues are correct and
!          ordered for indices 1,2,...ierr-1, but may not be
!          the smallest eigenvalues.
!
!        e2 has been destroyed.
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
      integer           i,j,l,m
      integer        :: jmax = 30
      type (mp_real)    b,c,f,g,h,p,r,s,t, pythag
      type (mp_real)    zero, one, two, eps, epsilonn
!
      zero = '0.d0'
      one = '1.d0'
      two = '2.d0'
      eps = epsilonn(one)
!
      ierr = 0
      if (n == 1) return
!
      do i=2,n
         e2(i-1) = e2(i)
      end do
!
      f = zero
      t = zero
      e2(n) = zero
!
      do l = 1, n
         j = 0
         h = abs(d(l)) + sqrt(e2(l))
         if (t < h) then
            t = h
            b = eps*t
            c = b * b
         end if
!     .......... look for small squared sub-diagonal element ..........
         do m = l, n
            if (e2(m) <= c) exit
!     .......... e2(n) is always zero, so there is no exit
!                through the bottom of the loop ..........
         end do
!
         if (m /= l) then
            do  ! cycle of QL transformations for the interval [l,m]
               if (j == jmax) goto 1000
               j = j + 1
!     .......... form shift ..........
               s = sqrt(e2(l))
               g = d(l)
               p = (d(l+1) - g) / (two * s)
               r = pythag(p,one)
               d(l) = s / (p + sign(r,p))
               h = g - d(l)
!
               do i=l+1,n
                  d(i) = d(i) - h
               end do
!
               f = f + h
!     .......... rational ql transformation ..........
               g = d(m)
               if (g == zero) g = b
               h = g
               s = zero
!     .......... for i=m-1 step -1 until l do -- ..........
               do i=m-1,l,-1
                  p = g * h
                  r = p + e2(i)
                  e2(i+1) = s * r
                  s = e2(i) / r
                  d(i+1) = h + s * (h + d(i))
                  g = d(i) - e2(i) / g
                  if (g == zero) g = b
                  h = g * p / r
               end do
!
               e2(l) = s * g
               d(l) = h
!     .......... guard against underflow in convergence test ..........
               if (h == zero) exit
               if (abs(e2(l)) <= abs(c/h)) exit
               e2(l) = h * e2(l)
               if (e2(l) == zero) exit
            end do
         end if
         p = d(l) + f
!     .......... order eigenvalues ..........
         do i=l,2,-1
            if (p >= d(i-1)) goto 100
            d(i) = d(i-1)
         end do
!
         i = 1
  100    d(i) = p
      end do
!
      return
!     .......... set error -- no convergence to an
!                eigenvalue after 30 iterations ..........
 1000 ierr = l
      return
      end
