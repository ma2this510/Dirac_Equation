      subroutine rsg(nm,n,a,b,w,matz,z,fv1,fv2,ierr)
USE mpmodule
      implicit none
      integer           n,nm,ierr,matz
      type (mp_real)    a(nm,n),b(nm,n),w(n),z(nm,n),fv1(n),fv2(n)
!
!     this subroutine calls the recommended sequence of
!     subroutines from the eigensystem subroutine package (eispack)
!     to find the eigenvalues and eigenvectors (if desired)
!     for the real symmetric generalized eigenproblem  ax = (lambda)bx.
!
!     on input
!
!        nm  must be set to the row dimension of the two-dimensional
!        array parameters as declared in the calling program
!        dimension statement.
!
!        n  is the order of the matrices  a  and  b.
!
!        a  contains a real symmetric matrix.
!
!        b  contains a positive definite real symmetric matrix.
!
!        matz  is an integer variable set equal to zero if
!        only eigenvalues are desired.  otherwise it is set to
!        any non-zero integer for both eigenvalues and eigenvectors.
!
!     on output
!
!        w  contains the eigenvalues in ascending order.
!
!        z  contains the eigenvectors if matz is not zero.
!
!        ierr  is an integer output variable set equal to an error
!           completion code described in the documentation for tqlrat
!           and tql2.  the normal completion code is zero.
!
!        fv1  and  fv2  are temporary storage arrays.
!
!     questions and comments should be directed to Burton S. Garbow,
!     Mathematics and Computer Science Div, Argonne National Laboratory
!
!     this version dated November 2000.
!
!     ------------------------------------------------------------------
!
      
      if (n > nm) then
         ierr = 10 * n
         return
      end if
!
      call reduc(nm,n,a,b,fv2,ierr)
      if (ierr /= 0) return
      if (matz == 0) then
!     .......... find eigenvalues only ..........
         call  tred1(nm,n,a,w,fv1,fv2)
         call  tqlrat(n,w,fv2,ierr)
      else
!     .......... find both eigenvalues and eigenvectors ..........
         call  tred2(nm,n,a,w,fv1,z)
         call  tql2(nm,n,w,fv1,z,ierr)
         if (ierr /= 0) return
         call  rebak(nm,n,b,fv2,n,z)
      end if
!
      return
      end
