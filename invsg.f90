      subroutine invsg(a,b,n,eig,v,eps,ijob,maxit,wa)
!-------------------------------------------------------------
!
!   latest revision     - april 23, 1992
!
!   purpose             - inverse iteration for the
!                           generalized eigenvalue problem
!                           a*x = x*b*x - symmetric storage
!                           mode.
!
!   usage               - call invsg(a,b,n,eig,v,eps,ijob,maxit,wa)
!
!   arguments    a      - linear array, on input a contains
!                           the elements of matrix a in
!                           symmetric storage mode, the length
!                           of a is n*(n+1)/2, on output a is
!                           destroyed.
!                b      - linear array, on input b contains the
!                           elements of matrix b in symmetric
!                           storage mode.
!                n      - dimension of matrices a and b.
!                eig    - on input, eig contains the initial
!                           approximation to the eigenvalue,
!                         on output, eig contains renewed eigenvalue.
!                v      - on input,  v  contains the initial
!                           approximation to the eigenvector,
!                         on output,  v  contains renewed eigenvector
!                           with unit b-norm, (x,b*x)=1.
!                eps    - parameter of regularization,
!                           if eig is one of the lowest
!                             eigenvalues, then eps must be
!                             greater than or equal to zero,
!                           if eig is one of the uppest
!                             eigenvalues, then eps must be
!                             less than or equal to zero,
!                           if the problem is well defined,
!                             then it is better to set the
!                             value of eps to zero.
!                ijob   - job option parameter,
!                           if ijob = 0, then form matrix (a-xb)
!                             and prepare it for iterations,
!                           if ijob = 1, the matrix (a-xb)
!                             was previously factored by routine
!                             leq1s.
!                maxit  - maximal number of iterations (input).
!                wa     - work array of length 3*n+1
!                           on output, wa(1) contains
!                           the number of eigenvalues below
!                           the initial value of eig.
!
!   reqd. routines      - leq1s, avms
!
!-------------------------------------------------------------
USE mpmodule
      implicit type (mp_real) (a-h,o-z)
      dimension           a(*),b(*),v(*),wa(*)
!
      zero = '0.d0'
      one = '1.d0'
      sixtn = '16.d0'
      cscale = '256.d0'
      aln256 = log(cscale)
      reps = epsilonn(reps)
!
      ndet = 2*n+1
      errest = max(sixtn*n*reps,abs(eps))
!
!                         replace matrix a with matrix
!                           (a-eig*b)
!
      eold = eig
      nsc = n
      if (ijob.ne.1) then
         ii = 0
         do i=1,n
            do j=1,i
               ii = ii+1
               a(ii) = a(ii)-eig*b(ii)
            end do
            temp = max(abs(a(ii)),abs(eig*b(ii)))
            isc = anint(log(temp)/aln256)
            wa(nsc+i) = sixtn**(-isc)
            do j=0,i-1
               a(ii-j) = a(ii-j)*wa(nsc+i)*wa(nsc+i-j)
               b(ii-j) = b(ii-j)*wa(nsc+i)*wa(nsc+i-j)
            end do
            a(ii) = a(ii)+eps*abs(a(ii))
         end do
!
         call leq1s(a,n,v,1,n,1,wa(ndet),ier)
!
      end if
!
!                         begin iterations
!
      do it=1,maxit
!
!                         multiply on matrix b and normalize
!                           vector v to unity
!
         call avms(b,n,v,wa)
         sm = dot_productt(v(1:n),v(1:n),n)
         sm = one/sqrt(sm)
         do i=1,n
            t = wa(i)*sm
            wa(i) = v(i)*sm
            v(i) = t
         end do
!
!                                                     -1
!                         multiply on matrix (a-eig*b)
!
         call leq1s(a,n,v,1,n,2,wa(ndet),ier)
!
!                         calculate refined eigenvalue
!
         sm = dot_productt(v(1:n),wa(1:n),n)
         eprev = eig
         eig = eold+one/sm
         if (abs(eprev-eig).lt.errest*abs(eig)) goto 100
!
      end do
!
!                         normalize eigenvector
!
  100 call avms(b,n,v,wa)
      sm = dot_productt(wa(1:n),v(1:n),n)
      sm = one/sqrt(sm)
      !v(1:n) = v(1:n)*wa(nsc+1:nsc+n)*sm
      do i=1,n !this loop replaces the line above
         v(i) = v(i)*wa(nsc+i)*sm
      end do

      wa(1) = wa(ndet+n)
!
!                        exit
!
      return
      end
      subroutine avms(a,n,v,vm)
!------------------------------------------------------------
!
!   purpose             - multiply a vector v on matrix a
!                           stored in symmetric mode
!                           - nucleus for invsg.
!
!   reqd. routines      - none required.
!
!-------------------------------------------------------------
USE mpmodule
      implicit type (mp_real) (a-h,o-z)
      dimension           v(*),vm(*),a(*)
!
      zero = '0.d0'

      vm(1) = v(1)*a(1)
      ii = 1
      do i=2,n
         im1 = i-1
         vm(i) = zero
         do j=1,im1
            vm(i) = vm(i)+a(ii+j)*v(j)
            vm(j) = vm(j)+a(ii+j)*v(i)
         end do
         ii = ii+i
         vm(i) = vm(i)+a(ii)*v(i)
      end do
!
      return
      end



      function dot_productt(v1,v2,n)
!------------------------------------------------------------
!
!   purpose             - dot product
!
!-------------------------------------------------------------
USE mpmodule
      implicit type (mp_real) (a-h,o-z)
      dimension           v1(*),v2(*)
     
      zero = '0.d0'
      dot_productt = zero  
      do i=1,n
         dot_productt = dot_productt + v1(i)*v2(i)
      end do
!
      return
      end
