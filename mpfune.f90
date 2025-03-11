!*****************************************************************************

!  MPFUN20-Fort: A thread-safe arbitrary precision computation package
!  Special functions module (module MPFUNE)

!  Revision date:  19 May 2022

!  AUTHOR:
!    David H. Bailey
!    Lawrence Berkeley National Lab (retired) and University of California, Davis
!    Email: dhbailey@lbl.gov

!  COPYRIGHT AND DISCLAIMER:
!    All software in this package (c) 2022 David H. Bailey.
!    By downloading or using this software you agree to the copyright, disclaimer
!    and license agreement in the accompanying file DISCLAIMER.txt.

!  PURPOSE OF PACKAGE:
!    This package permits one to perform floating-point computations (real and
!    complex) to arbitrarily high numeric precision, by making only relatively
!    minor changes to existing Fortran-90 programs.  All basic arithmetic
!    operations and transcendental functions are supported, together with several
!    special functions.

!    In addition to fast execution times, one key feature of this package is a
!    100% THREAD-SAFE design, which means that user-level applications can be
!    easily converted for parallel execution, say using a threaded parallel
!    environment such as OpenMP.

!  DOCUMENTATION:
!    A detailed description of this package, and instructions for compiling
!    and testing this program on various specific systems are included in the
!    README file accompanying this package, and, in more detail, in the
!    following technical paper:

!    David H. Bailey, "MPFUN2020: A new thread-safe arbitrary precision package,"
!    available at http://www.davidhbailey.com/dhbpapers/mpfun2020.pdf.

!  DESCRIPTION OF THIS MODULE (MPFUNE):
!    This module contains subroutines to perform special functions. Additional
!    functions will be added as they are completed.

!  NOTE ON PROGRAMMING CONVENTION FOR THIS MODULE:
!    This module is designed to facilitate easy translation (using the program
!    convfmp.f90) to MPFR calls, for use in the MPFUN-MPFR package.

module mpfune
use mpfuna
use mpfunb
use mpfunc
use mpfund

contains

!  These routines perform simple operations on the MP data structure. Those
!  listed between !> and !>> are for MPFUN20-Fort; those between !>> and !>>>
!  are for MPFUN20-MPFR. The translator selects the proper set.
!>

subroutine mpinitwds (ra, mpnw)
implicit none
integer (mpiknd) ra(0:)
integer mpnw
ra(0) = mpnw + 6
ra(1) = mpnw
ra(2) = 0
ra(3) = 0
ra(4) = 0
return
end subroutine mpinitwds

function mpwprecr (ra)
implicit none
integer (mpiknd) ra(0:)
integer mpwprecr
mpwprecr = ra(1)
return
end function mpwprecr

function mpspacer (ra)
implicit none
integer (mpiknd) ra(0:)
integer mpspacer
mpspacer = ra(0)
return
end function mpspacer

!>>

! subroutine mpabrt (ier)
! implicit none
! integer, intent(in):: ier
! write (mpldb, 1) ier
! 1 format ('*** MPABRT: Execution terminated, error code =',i4)
! stop
! end subroutine mpabrt

! subroutine mpinitwds (ra, mpnw)
! implicit none
! integer (mpiknd) ra(0:)
! integer mpnw
! ra(0) = mpnw + 6
! ra(1) = mpnw * mpnbt
! ra(2) = 1
! ra(3) = mpnan
! ra(4) = loc (ra(4)) + 8
! ra(mpnw+5) = 0
! return
! end subroutine mpinitwds

! subroutine mpfixlocr (ia)
! implicit none
! integer (mpiknd) ia(0:)
! ia(4) = loc (ia(4)) + 8
! return
! end subroutine

! function mpwprecr (ra)
! implicit none
! integer (mpiknd) ra(0:)
! integer mpwprecr
! mpwprecr = ra(1) / mpnbt
! return
! end function mpwprecr

! function mpspacer (ra)
! implicit none
! integer (mpiknd) ra(0:)
! integer mpspacer
! mpspacer = ra(0)
! return
! end function mpspacer

!>>>

subroutine mpberner (nb1, nb2, berne, mpnw)

!   This returns the array berne, containing Bernoulli numbers indexed 2*k for
!   k = 1 to n, to mpnw words precision. This is done by first computing
!   zeta(2*k), based on the following known formulas:

!   coth (pi*x) = cosh (pi*x) / sinh (pi*x)

!            1      1 + (pi*x)^2/2! + (pi*x)^4/4! + ...
!        =  ---- * -------------------------------------
!           pi*x    1 + (pi*x)^2/3! + (pi*x)^4/5! + ...

!        = 1/(pi*x) * (1 + (pi*x)^2/3 - (pi*x)^4/45 + 2*(pi*x)^6/945 - ...)

!        = 2/(pi*x) * Sum_{k >= 1} (-1)^(k+1) * zeta(2*k) * x^{2*k}

!   The strategy is to calculate the coefficients of the series by polynomial
!   operations. Polynomial division is performed by computing the reciprocal
!   of the denominator polynomial, by a polynomial Newton iteration, as follows.
!   Let N(x) be the polynomial approximation to the numerator series; let D(x) be
!   a polynomial approximation to the numerator numerator series; and let Q_k(x)
!   be polynomial approximations to R(x) = 1/D(x).  Then iterate:

!   Q_{k+1} = Q_k(x) + [1 - D(x)*Q_k(x)]*Q_k(x)

!   In these iterations, both the degree of the polynomial Q_k(x) and the
!   precision level in words are initially set to 4. When convergence is
!   achieved at this precision level, the degree is doubled, and iterations are
!   continued, etc., until the final desired degree is achieved. Then the
!   precision level is doubled and iterations are performed in a similar way,
!   until the final desired precision level is achieved. The reciprocal polynomial
!   R(x) produced by this process is then multiplied by the numerator polynomial
!   N(x) to yield an approximation to the quotient series. The even zeta values
!   are then the coefficients of this series, scaled according to the formula above.

!   Once the even integer zeta values have been computed in this way, the even
!   Bernoulli numbers are computed via the formula (for n > 0):

!   B(2*n) = (-1)^(n-1) * 2 * (2*n)! * zeta(2*n) / (2*pi)^(2*n)

!   Note: The notation in the description above is not the same as in the code below.

implicit none
integer, intent(in):: nb1, nb2, mpnw
integer i, idb, i1, ic1, j, kn, mpnw1, nwds, n, n1, nn1
integer (mpiknd), intent(out):: berne(0:nb1+5,nb2)
integer (mpiknd) c1(0:mpnw+6,0:nb2), cp2(0:mpnw+6), p1(0:mpnw+6,0:nb2), &
  p2(0:mpnw+6,0:nb2), q(0:mpnw+6,0:nb2), q1(0:mpnw+6), &
  r(0:mpnw+6,0:nb2), s(0:mpnw+6,0:nb2), t1(0:mpnw+6), t2(0:mpnw+6), &
  t3(0:mpnw+6), t4(0:mpnw+6), eps(0:mpnw+6)
real (mprknd) alog102, d1, dd1, dd2, dd3, pi
parameter (idb = 0, alog102 = 0.30102999566398119d0, pi = 3.1415926535897932385d0)
integer (mpiknd) ix8

!  End of declaration.

i1 = 1000000000

do i = 1, nb2
  i1 = min (i1, mpspacer (berne(0:nb1+5,i)))
enddo

if (mpnw < 4 .or. i1 < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBERNER: Uninitialized or inadequately sized arrays')
  call mpabrt ( 501)
endif

n = nb2
mpnw1 = mpnw + 1
nwds = mpnw1

!   Check if Pi has been precomputed.

call mpmdc (mppicon, d1, n1, mpnw1)
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw1) then
  write (mpldb, 2) mpnw1
2 format ('*** MPBERNER: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 502)
endif

if (idb > 0) write (mpldb, 3) n, mpnw
3 format ('Even Bernoulli number calculation; n, mpnw =',2i6)

call mpinitwds (cp2, nwds)
call mpinitwds (q1, nwds)
call mpinitwds (t1, nwds)
call mpinitwds (t2, nwds)
call mpinitwds (t3, nwds)
call mpinitwds (t4, nwds)
call mpinitwds (eps, nwds)

do i = 0, nb2
  call mpinitwds (c1(0:mpnw1+5,i), nwds)
  call mpinitwds (p1(0:mpnw1+5,i), nwds)
  call mpinitwds (p2(0:mpnw1+5,i), nwds)
  call mpinitwds (q(0:mpnw1+5,i), nwds)
  call mpinitwds (r(0:mpnw1+5,i), nwds)
  call mpinitwds (s(0:mpnw1+5,i), nwds)
enddo

call mpmul (mppicon, mppicon, cp2, nwds)
call mpdmc (1.d0, 0, c1(0:mpnw1+5,0), nwds)
call mpdmc (1.d0, 0, p1(0:mpnw1+5,0), nwds)
call mpdmc (1.d0, 0, p2(0:mpnw1+5,0), nwds)
call mpdmc (1.d0, 0, q(0:mpnw1+5,0), nwds)

!   Construct numerator and denominator polynomials.

do i = 1, n
  call mpdmc (0.d0, 0, c1(0:mpnw1+5,i), nwds)
  dd1 = 2.d0 * (i + 1) - 3.d0
  dd2 = dd1 + 1.d0
  dd3 = dd2 + 1.d0
  call mpmul (cp2, p1(0:mpnw1+5,i-1), t1, nwds)
  call mpdivd (t1, dd1 * dd2, p1(0:mpnw1+5,i), nwds)
  call mpmul (cp2, p2(0:mpnw1+5,i-1), t1, nwds)
  call mpdivd (t1, dd2 * dd3, p2(0:mpnw1+5,i), nwds)
  call mpdmc (0.d0, 0, q(0:mpnw1+5,i), nwds)
enddo

kn = 4
nwds = 4

call mpdmc (2.d0, 0, t1, nwds)
call mpnpwr (t1, 70 - nwds * mpnbt, eps, nwds)
if (idb > 0) then
  call mpmdc (eps, dd1, nn1, nwds)
  write (mpldb, 4) nwds, nint (alog102*nn1)
4 format ('nwds, log10eps =',2i6)
endif
call mpdmc (0.d0, 0, q1, nwds)

!   Perform Newton iterations with dynamic precision levels, using an
!   iteration formula similar to that used to evaluate reciprocals.

do j = 1, 10000
  if (idb > 0) write (mpldb, 5) j, kn, nwds
5 format ('j, kn, nwds =',3i6)

  call mppolymul (mpnw1, kn, p2, q, r, nwds)
  call mppolysub (mpnw1, kn, c1, r, s, nwds)
  call mppolymul (mpnw1, kn, s, q, r, nwds)
  call mppolyadd (mpnw1, kn, q, r, q, nwds)
  call mpsub (q(0:mpnw1+5,kn), q1, t1, nwds)

  if (idb > 0) then
    call mpmdc (t1, dd1, nn1, nwds)
    if (dd1 .eq. 0.d0) then
      write (mpldb, 6)
6     format ('Newton error = 0')
    else
      write (mpldb, 7) nint (alog102*nn1)
7     format ('Newton error = 10^',i6)
    endif
  endif

  call mpabs (t1, t2, nwds)
  call mpcpr (t2, eps, ic1, nwds)
  if (ic1 < 0) then
    if (kn == n .and. nwds == mpnw1) goto 100
    if (kn < n) then
      kn = min (2 * kn, n)
      call mpdmc (0.d0, 0, q1, nwds)
    elseif (nwds < mpnw1) then
      nwds = min (2 * nwds, mpnw1)
      call mpdmc (2.d0, 0, t1, nwds)
      call mpnpwr (t1, 70 - nwds * mpnbt, eps, nwds)
      call mpdmc (0.d0, 0, q1, nwds)
      if (idb > 0) then
        call mpmdc (eps, dd1, nn1, nwds)
        write (mpldb, 4) nwds, nint (alog102*nn1)
      endif
    endif
  else
    call mpeq (q(0:mpnw1+5,kn), q1, nwds)
  endif
enddo

write (mpldb, 8)
8 format ('*** MPBERNER: Loop end error')
call mpabrt ( 503)

100 continue

if (idb > 0) write (mpldb, 9)
9 format ('Even zeta computation complete')

!   Multiply numerator polynomial by reciprocal of denominator polynomial.

call mppolymul (mpnw1, n, p1, q, r, nwds)

!   Apply formula to produce Bernoulli numbers.

call mpdmc (-2.d0, 0, t1, nwds)
call mpdmc (1.d0, 0, t2, nwds)

do i = 1, n
  d1 = - dble (2*i-1) * dble (2*i)
  call mpmuld (t1, d1, t3, nwds)
  call mpeq (t3, t1, nwds)
  call mpmuld (cp2, 4.d0, t3, nwds)
  call mpmul (t3, t2, t4, nwds)
  call mpeq (t4, t2, nwds)
  call mpmuld (t1, 0.5d0, t3, nwds)
  call mpdiv (t3, t2, t4, nwds)
  call mpabs (r(0:mpnw1+5,i), t3, nwds)
  call mpmul (t4, t3, berne(0:nb1+5,i), mpnw)
enddo

if (idb > 0) write (mpldb, 10)
10 format ('Bernoulli number computation complete')
return
end subroutine mpberner

subroutine mppolyadd (nd1, n, a, b, c, nwds)

!   This adds two polynomials, as is required by mpberne.
!   The output array C may be the same as A or B.

implicit none
integer, intent(in):: n, nd1, nwds
integer k
integer (mpiknd), intent(in):: a(0:nd1+5,0:n), b(0:nd1+5,0:n)
integer (mpiknd), intent(out):: c(0:nd1+5,0:n)
integer (mpiknd) t1(0:nwds+5), t2(0:nwds+5)

call mpinitwds (t1, nwds)
call mpinitwds (t2, nwds)

do k = 0, n
  call mpeq (a(0:nd1+5,k), t1, nwds)
  call mpeq (b(0:nd1+5,k), t2, nwds)
  call mpadd (t1, t2, c(0:nd1+5,k), nwds)
enddo

return
end subroutine mppolyadd

subroutine mppolysub (nd1, n, a, b, c, nwds)

!   This adds two polynomials, as is required by mpberne.
!   The output array C may be the same as A or B.

implicit none
integer, intent(in):: n, nd1, nwds
integer k
integer (mpiknd), intent(in):: a(0:nd1+5,0:n), b(0:nd1+5,0:n)
integer (mpiknd), intent(out):: c(0:nd1+5,0:n)
integer (mpiknd) t1(0:nwds+5), t2(0:nwds+5)

call mpinitwds (t1, nwds)
call mpinitwds (t2, nwds)

do k = 0, n
  call mpeq (a(0:nd1+5,k), t1, nwds)
  call mpeq (b(0:nd1+5,k), t2, nwds)
  call mpsub (t1, t2, c(0:nd1+5,k), nwds)
enddo

return
end subroutine mppolysub

subroutine mppolymul (nd1, n, a, b, c, nwds)

!   This adds two polynomials (ignoring high-order terms), as is required
!   by mpberne. The output array C may not be the same as A or B.

implicit none
integer, intent(in):: n, nd1, nwds
integer j, k
integer (mpiknd), intent(in):: a(0:nd1+5,0:n), b(0:nd1+5,0:n)
integer (mpiknd), intent(out):: c(0:nd1+5,0:n)
integer (mpiknd) t0(0:nwds+5), t1(0:nwds+5), t2(0:nwds+5), t3(0:nwds+5)

call mpinitwds (t0, nwds)
call mpinitwds (t1, nwds)
call mpinitwds (t2, nwds)
call mpinitwds (t3, nwds)

do k = 0, n
  call mpdmc (0.d0, 0, t0, nwds)

  do j = 0, k
    call mpeq (a(0:nd1+5,j), t1, nwds)
    call mpeq (b(0:nd1+5,k-j), t2, nwds)
    call mpmul (t1, t2, t3, nwds)
    call mpadd (t0, t3, t2, nwds)
    call mpeq (t2, t0, nwds)
  enddo

  call mpeq (t0, c(0:nd1+5,k), nwds)
enddo

return
end subroutine mppolymul

subroutine mpbesselinr (nu, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselI (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.25.2 for modest RR,
!   and DLMF 10.40.1 for large RR, relative to precision.

implicit none
integer, intent(in):: nu, mpnw
integer ic1, itrmax, k, mpnw1, nua, n1
real (mprknd) dfrac, d1, pi
parameter (itrmax = 1000000, dfrac = 1.5d0, pi = 3.1415926535897932385d0)
integer (mpiknd), intent(in):: rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer (mpiknd) f1(0:mpnw+6), f2(0:mpnw+6), sum(0:mpnw+6), td(0:mpnw+6), &
  tn(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  rra(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (rr) < mpnw + 4 .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELINR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 504)
endif

!   Check for RR = 0.

if (mpsgn (rr) == 0) then
  write (mpldb, 2)
2 format ('*** MPBESSELINR: Second argument is zero')
  call mpabrt ( 505)
endif

!   Check if PI has been precomputed.

mpnw1 = mpnw + 1
call mpmdc (mppicon, d1, n1, mpnw1)
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw1) then
  write (mpldb, 3) mpnw1
3 format ('*** MPBESSELINR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 506)
endif

nua = abs (nu)
call mpinitwds (f1, mpnw1)
call mpinitwds (f2, mpnw1)
call mpinitwds (sum, mpnw1)
call mpinitwds (td, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (rra, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)
call mpabs (rr, rra, mpnw1)
call mpmdc (rra, d1, n1, 4)
d1 = d1 * 2.d0 ** n1

if (d1 < dfrac * mpnw1 * mpdpw) then
  call mpdmc (1.d0, 0, tn, mpnw1)
  call mpdmc (1.d0, 0, f1, mpnw1)
  call mpdmc (1.d0, 0, f2, mpnw1)
  call mpmul (rra, rra, t2, mpnw1)
  call mpmuld (t2, 0.25d0, t1, mpnw1)

  do k = 1, nua
    call mpmuld (f2, dble (k), t2, mpnw1)
    call mpeq (t2, f2, mpnw1)
  enddo

  call mpmul (f1, f2, td, mpnw1)
  call mpdiv (tn, td, t2, mpnw1)
  call mpeq (t2, sum, mpnw1)

  do k = 1, itrmax
    call mpmuld (f1, dble (k), t2, mpnw1)
    call mpeq (t2, f1, mpnw1)
    call mpmuld (f2, dble (k + nua), t2, mpnw1)
    call mpeq (t2, f2, mpnw1)
    call mpmul (t1, tn, t2, mpnw1)
    call mpeq (t2, tn, mpnw1)
    call mpmul (f1, f2, td, mpnw1)
    call mpdiv (tn, td, t2, mpnw1)
    call mpadd (sum, t2, t3, mpnw1)
    call mpeq (t3, sum, mpnw1)

    call mpabs (t2, tc1, 4)
    call mpmul (eps, sum, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 100
  enddo

  write (mpldb, 4)
  4 format ('*** MPBESSELINR: Loop end error 1')
  call mpabrt ( 507)

100 continue

  call mpmuld (rra, 0.5d0, t1, mpnw1)
  call mpnpwr (t1, nua, t2, mpnw1)
  call mpmul (sum, t2, t3, mpnw1)
else
!  sum1 = mpreal (1.d0, mpnw)
!  t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
!  t2 = mpreal (1.d0, mpnw)
!  t3 = mpreal (1.d0, mpnw)

  call mpdmc (1.d0, 0, sum, mpnw1)
  d1 = 4.d0 * dble (nua)**2
  call mpdmc (d1, 0, t1, mpnw1)
  call mpdmc (1.d0, 0, tn, mpnw1)
  call mpdmc (1.d0, 0, td, mpnw1)

  do k = 1, itrmax
!  t2 = -t2 * (t1 - (2.d0*k - 1.d0)**2)
!  t3 = t3 * dble (k) * 8.d0 * xa
!  t4 = t2 / t3
!  sum1 = sum1 + t4

    d1 = 2.d0 * k - 1.d0
    call mpdmc (d1, 0, t2, mpnw1)
    call mpmul (t2, t2, t3, mpnw1)
    call mpsub (t1, t3, t2, mpnw1)
    call mpmul (tn, t2, t3, mpnw1)
    call mpneg (t3, tn, mpnw1)
    call mpmuld (rra, 8.d0 * k, t2, mpnw1)
    call mpmul (td, t2, t3,  mpnw1)
    call mpeq (t3, td, mpnw1)
    call mpdiv (tn, td, t4, mpnw1)
    call mpadd (sum, t4, t3, mpnw1)
    call mpeq (t3, sum, mpnw1)

!   if (abs (t4) / abs (sum1) < eps) goto 110

    call mpabs (t4, tc1, 4)
    call mpmul (eps, sum, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 110
  enddo

write (mpldb, 5)
5 format ('*** MPBESSELINR: Loop end error 2')
call mpabrt ( 508)

110 continue

! t1 = exp (xa) / sqrt (2.d0 * mppi (mpnw) * xa)
! besseli = t1 * sum1

  call mpexp (rra, t1, mpnw1)
  call mpmuld (mppicon, 2.d0, t2, mpnw1)
  call mpmul (t2, rra, t3, mpnw1)
  call mpsqrt (t3, t4, mpnw1)
  call mpdiv (t1, t4, t2, mpnw1)
  call mpmul (t2, sum, t3, mpnw1)
endif

! if (x < 0.d0 .and. mod (nu, 2) /= 0) besseli = - besseli

if (mpsgn (rr) < 0 .and. mod (nu, 2) /= 0) then
  call mpneg (t3, t4, mpnw1)
  call mpeq (t4, t3, mpnw1)
endif

call mproun (t3, mpnw)
call mpeq (t3, ss, mpnw)

return
end subroutine mpbesselinr

subroutine mpbesselir (qq, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselI (QQ,RR) for QQ and RR
!   both MPR. The algorithm is DLMF formula 10.25.2 for modest RR, and
!   DLMF 10.40.1 for large RR, relative to precision.

implicit none
integer, intent(in):: mpnw
integer ic1, i1, i2, itrmax, k, mpnw1, n1
real (mprknd) dfrac, d1, pi
parameter (itrmax = 1000000, dfrac = 1.5d0, pi = 3.1415926535897932385d0)
integer (mpiknd), intent(in):: qq(0:), rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer (mpiknd) f1(0:mpnw+6), f2(0:mpnw+6), sum(0:mpnw+6), td(0:mpnw+6), &
  tn(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  rra(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (rr) < mpnw + 4 .or. mpspacer (qq) < mpnw + 4 &
  .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELIR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 509)
endif

!   Check if PI has been precomputed.

mpnw1 = mpnw + 1
call mpmdc (mppicon, d1, n1, mpnw1)
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw1) then
  write (mpldb, 2) mpnw1
2 format ('*** MPBESSELIR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 510)
endif

call mpinitwds (f1, mpnw1)
call mpinitwds (f2, mpnw1)
call mpinitwds (sum, mpnw1)
call mpinitwds (td, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (rra, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

!   If QQ is integer, call mpbesselinr; if qq < 0 and rr <= 0, then error.

call mpinfr (qq, t1, t2, mpnw)
i1 = mpsgn (qq)
i2 = mpsgn (rr)
if (mpsgn (t2) == 0) then
  call mpmdc (qq, d1, n1, mpnw)
  d1 = d1 * 2.d0**n1
  n1 = nint (d1)
  call mpbesselinr (n1, rr, t3, mpnw)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (mpldb, 3)
3 format ('*** MPBESSELIR: First argument < 0 and second argument <= 0')
  call mpabrt ( 511)
endif

call mpabs (rr, rra, mpnw1)
call mpmdc (rra, d1, n1, 4)
d1 = d1 * 2.d0 ** n1

if (d1 < dfrac * mpnw1 * mpdpw) then
  call mpdmc (1.d0, 0, tn, mpnw1)
  call mpdmc (1.d0, 0, f1, mpnw1)
  call mpadd (qq, f1, t1, mpnw1)
  call mpgammar (t1, f2, mpnw1)
  call mpmul (rra, rra, t2, mpnw1)
  call mpmuld (t2, 0.25d0, t1, mpnw1)

  call mpmul (f1, f2, td, mpnw1)
  call mpdiv (tn, td, t2, mpnw1)
  call mpeq (t2, sum, mpnw1)

  do k = 1, itrmax
    call mpmuld (f1, dble (k), t2, mpnw1)
    call mpeq (t2, f1, mpnw1)
    call mpdmc (dble (k), 0, t3, mpnw1)
    call mpadd (qq, t3, t4, mpnw1)
    call mpmul (f2, t4, t3, mpnw1)
    call mpeq (t3, f2, mpnw1)
    call mpmul (t1, tn, t2, mpnw1)
    call mpeq (t2, tn, mpnw1)
    call mpmul (f1, f2, td, mpnw1)
    call mpdiv (tn, td, t2, mpnw1)
    call mpadd (sum, t2, t3, mpnw1)
    call mpeq (t3, sum, mpnw1)

    call mpabs (t2, tc1, 4)
    call mpmul (eps, sum, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 100
  enddo

  write (mpldb, 4)
  4 format ('*** MPBESSELIR: Loop end error 1')
  call mpabrt ( 512)

100 continue

  call mpmuld (rra, 0.5d0, t1, mpnw1)
  call mppower (t1, qq, t2, mpnw1)
  call mpmul (sum, t2, t3, mpnw1)
else
!  sum1 = mpreal (1.d0, mpnw)
!  t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
!  t2 = mpreal (1.d0, mpnw)
!  t3 = mpreal (1.d0, mpnw)

  call mpdmc (1.d0, 0, sum, mpnw1)
  call mpmul (qq, qq, t2, mpnw1)
  call mpmuld (t2, 4.d0, t1, mpnw1)
  call mpdmc (1.d0, 0, tn, mpnw1)
  call mpdmc (1.d0, 0, td, mpnw1)

  do k = 1, itrmax
!  t2 = -t2 * (t1 - (2.d0*k - 1.d0)**2)
!  t3 = t3 * dble (k) * 8.d0 * xa
!  t4 = t2 / t3
!  sum1 = sum1 + t4

    d1 = 2.d0 * k - 1.d0
    call mpdmc (d1, 0, t2, mpnw1)
    call mpmul (t2, t2, t3, mpnw1)
    call mpsub (t1, t3, t2, mpnw1)
    call mpmul (tn, t2, t3, mpnw1)
    call mpneg (t3, tn, mpnw1)
    call mpmuld (rra, 8.d0 * k, t2, mpnw1)
    call mpmul (td, t2, t3,  mpnw1)
    call mpeq (t3, td, mpnw1)
    call mpdiv (tn, td, t4, mpnw1)
    call mpadd (sum, t4, t3, mpnw1)
    call mpeq (t3, sum, mpnw1)

!   if (abs (t4) / abs (sum1) < eps) goto 110

    call mpabs (t4, tc1, 4)
    call mpmul (eps, sum, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 110
  enddo

write (mpldb, 5)
5 format ('*** MPBESSELIR: Loop end error 2')
call mpabrt ( 513)

110 continue

! t1 = exp (xa) / sqrt (2.d0 * mppi (mpnw) * xa)
! besseli = t1 * sum1

  call mpexp (rra, t1, mpnw1)
  call mpmuld (mppicon, 2.d0, t2, mpnw1)
  call mpmul (t2, rra, t3, mpnw1)
  call mpsqrt (t3, t4, mpnw1)
  call mpdiv (t1, t4, t2, mpnw1)
  call mpmul (t2, sum, t3, mpnw1)
endif

120 continue

call mproun (t3, mpnw)
call mpeq (t3, ss, mpnw)

return
end subroutine mpbesselir

subroutine mpbesseljnr (nu, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselJ (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.2.2 for modest RR,
!   and DLMF 10.17.3 for large RR, relative to precision.

implicit none
integer, intent(in):: nu, mpnw
integer ic1, ic2, itrmax, k, mpnw1, nua, n1
real (mprknd) dfrac, d1, d2, pi
parameter (itrmax = 1000000, dfrac = 1.5d0, pi = 3.1415926535897932385d0)
integer (mpiknd), intent(in):: rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer (mpiknd) f1(0:2*mpnw+6), f2(0:2*mpnw+6), sum1(0:2*mpnw+6), &
  sum2(0:2*mpnw+6), td1(0:2*mpnw+6), td2(0:2*mpnw+6), tn1(0:2*mpnw+6), &
  tn2(0:2*mpnw+6), t1(0:2*mpnw+6), t2(0:2*mpnw+6), t3(0:2*mpnw+6), &
  t41(0:2*mpnw+6), t42(0:2*mpnw+6), t5(0:2*mpnw+6), rra(0:2*mpnw+6), &
  rr2(0:2*mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (rr) < mpnw + 4 .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELJNR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 514)
endif

!   Check for RR = 0.

if (mpsgn (rr) == 0) then
  write (mpldb, 2)
2 format ('*** MPBESSELJNR: Second argument is zero')
  call mpabrt ( 515)
endif

!   Check if PI has been precomputed.

mpnw1 = 2 * mpnw + 1
call mpmdc (mppicon, d1, n1, mpnw1)
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw1) then
  write (mpldb, 3) mpnw1
3 format ('*** MPBESSELJNR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 516)
endif

nua = abs (nu)
call mpinitwds (f1, mpnw1)
call mpinitwds (f2, mpnw1)
call mpinitwds (sum1, mpnw1)
call mpinitwds (sum2, mpnw1)
call mpinitwds (td1, mpnw1)
call mpinitwds (tn1, mpnw1)
call mpinitwds (td2, mpnw1)
call mpinitwds (tn2, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t41, mpnw1)
call mpinitwds (t42, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (rra, mpnw1)
call mpinitwds (rr2, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw*mpnbt, eps, 4)
call mpmdc (rr, d1, n1, 4)
d1 = abs (d1) * 2.d0 ** n1

if (d1 < dfrac * mpnw1 * mpdpw) then
  mpnw1 = min (mpnw + nint (d1 / (dfrac * mpdpw)), 2 * mpnw + 1)
  call mpabs (rr, rra, mpnw1)
  call mpdmc (1.d0, 0, tn1, mpnw1)
  call mpdmc (1.d0, 0, f1, mpnw1)
  call mpdmc (1.d0, 0, f2, mpnw1)
  call mpmul (rra, rra, t2, mpnw1)
  call mpmuld (t2, 0.25d0, t1, mpnw1)

  do k = 1, nua
    call mpmuld (f2, dble (k), t2, mpnw1)
    call mpeq (t2, f2, mpnw1)
  enddo

  call mpmul (f1, f2, td1, mpnw1)
  call mpdiv (tn1, td1, t2, mpnw1)
  call mpeq (t2, sum1, mpnw1)

  do k = 1, itrmax
    call mpmuld (f1, dble (k), t2, mpnw1)
    call mpeq (t2, f1, mpnw1)
    call mpmuld (f2, dble (k + nua), t2, mpnw1)
    call mpeq (t2, f2, mpnw1)
    call mpmul (t1, tn1, t2, mpnw1)
    call mpneg (t2, tn1, mpnw1)
    call mpmul (f1, f2, td1, mpnw1)
    call mpdiv (tn1, td1, t2, mpnw1)
    call mpadd (sum1, t2, t3, mpnw1)
    call mpeq (t3, sum1, mpnw1)

    call mpabs (t2, tc1, 4)
    call mpmul (eps, sum1, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 100
  enddo

  write (mpldb, 4)
4 format ('*** MPBESSELJNR: Loop end error 1')
  call mpabrt ( 517)

100 continue

  call mpmuld (rra, 0.5d0, t1, mpnw1)
  call mpnpwr (t1, nua, t2, mpnw1)
  call mpmul (sum1, t2, t3, mpnw1)
else
! xa2 = xa**2
! t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
! tn1 = mpreal (1.d0, mpnw)
! tn2 = (t1 - 1.d0) / 8.d0
! td1 = mpreal (1.d0, mpnw)
! td2 = xa
! sum1 = tn1 / td1
! sum2 = tn2 / t32

  mpnw1 = mpnw + 1
  call mpabs (rr, rra, mpnw1)
  call mpmul (rra, rra, rr2, mpnw1)
  d1 = 4.d0 * dble (nua)**2
  call mpdmc (d1, 0, t1, mpnw1)
  call mpdmc (1.d0, 0, tn1, mpnw1)
  call mpsub (t1, tn1, t2, mpnw1)
  call mpdivd (t2, 8.d0, tn2, mpnw1)
  call mpdmc (1.d0, 0, td1, mpnw1)
  call mpeq (rra, td2, mpnw1)
  call mpdiv (tn1, td1, sum1, mpnw1)
  call mpdiv (tn2, td2, sum2, mpnw1)

  do k = 1, itrmax
!   tn1 = -tn1 * (t1 - (2.d0*(2.d0*k-1.d0) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k-1.d0) - 1.d0)**2
    d2 = (2.d0*(2.d0*k) - 1.d0)**2
    call mpdmc (d1, 0, t2, mpnw1)
    call mpsub (t1, t2, t3, mpnw1)
    call mpdmc (d2, 0, t2, mpnw1)
    call mpsub (t1, t2, t5, mpnw1)
    call mpmul (t3, t5, t2, mpnw1)
    call mpmul (tn1, t2, t3, mpnw1)
    call mpneg (t3, tn1, mpnw1)

!   td1 = td1 * dble (2*k-1) * dble (2*k) * 64.d0 * xa2

    d1 = dble (2*k-1) * dble (2*k) * 64.d0
    call mpmuld (td1, d1, t2, mpnw1)
    call mpmul (t2, rr2, td1, mpnw1)

!   t41 = tn1 / td1
!   sum1 = sum1 + t41

    call mpdiv (tn1, td1, t41, mpnw1)
    call mpadd (sum1, t41, t2, mpnw1)
    call mpeq (t2, sum1, mpnw1)

!   tn2 = -tn2 * (t1 - (2.d0*(2.d0*k) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k+1.d0) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k) - 1.d0)**2
    d2 = (2.d0*(2.d0*k+1.d0) - 1.d0)**2
    call mpdmc (d1, 0, t2, mpnw1)
    call mpsub (t1, t2, t3, mpnw1)
    call mpdmc (d2, 0, t2, mpnw1)
    call mpsub (t1, t2, t5, mpnw1)
    call mpmul (t3, t5, t2, mpnw1)
    call mpmul (tn2, t2, t3, mpnw1)
    call mpneg (t3, tn2, mpnw1)

!   td2 = td2 * dble (2*k) * dble (2*k+1) * 64.d0 * xa2

    d1 = dble (2*k) * dble (2*k+1) * 64.d0
    call mpmuld (td2, d1, t2, mpnw1)
    call mpmul (t2, rr2, td2, mpnw1)

!  t42 = tn2 / td2
!  sum2 = sum2 + t42

    call mpdiv (tn2, td2, t42, mpnw1)
    call mpadd (sum2, t42, t2, mpnw1)
    call mpeq (t2, sum2, mpnw1)

!  if (abs (t41) / abs (sum1) < eps .and. abs (t42) / abs (sum2) < eps ) goto 110

    call mpabs (t41, tc1, 4)
    call mpmul (eps, sum1, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    call mpabs (t42, tc1, 4)
    call mpmul (eps, sum2, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic2, 4)
    if (ic1 <= 0 .and. ic2 <= 0) goto 110
  enddo

  write (mpldb, 5)
5 format ('*** MPBESSELJNR: Loop end error 2')
  call mpabrt ( 518)

110 continue

! t1 = xa - 0.5d0 * nua * pi - 0.25d0 * pi
! besselj = sqrt (2.d0 / (pi * xa)) * (cos (t1) * sum1 - sin (t1) * sum2)

  call mpmuld (mppicon, 0.5d0 * nua, t1, mpnw1)
  call mpsub (rra, t1, t2, mpnw1)
  call mpmuld (mppicon, 0.25d0, t1, mpnw1)
  call mpsub (t2, t1, t3, mpnw1)
  call mpcssnr (t3, t41, t42, mpnw1)
  call mpmul (t41, sum1, t1, mpnw1)
  call mpmul (t42, sum2, t2, mpnw1)
  call mpsub (t1, t2, t5, mpnw1)
  call mpmul (mppicon, rra, t1, mpnw1)
  call mpdmc (2.d0, 0, t2, mpnw1)
  call mpdiv (t2, t1, t3, mpnw1)
  call mpsqrt (t3, t1, mpnw1)
  call mpmul (t1, t5, t3, mpnw1)
endif

if (mod (nu, 2) /= 0) then
!  if (nu < 0 .and. x > 0.d0 .or. nu > 0 .and. x < 0.d0) besselj = - besselj

  ic1 = mpsgn (rr)
  if (nu < 0 .and. ic1 > 0 .or. nu > 0 .and. ic1 < 0) then
    call mpneg (t3, t2, mpnw1)
    call mpeq (t2, t3, mpnw1)
  endif
endif

call mproun (t3, mpnw)
call mpeq (t3, ss, mpnw)

return
end subroutine mpbesseljnr

subroutine mpbesseljr (qq, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselJ (QQ,RR) for QQ and RR
!   both MPR. The algorithm is DLMF formula 10.2.2 for modest RR,
!   and DLMF 10.17.3 for large RR, relative to precision.

implicit none
integer, intent(in):: mpnw
integer ic1, ic2, i1, i2, itrmax, k, mpnw1, n1
real (mprknd) dfrac, d1, d2, pi
parameter (itrmax = 1000000, dfrac = 1.5d0, pi = 3.1415926535897932385d0)
integer (mpiknd), intent(in):: qq(0:), rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer (mpiknd) f1(0:2*mpnw+6), f2(0:2*mpnw+6), sum1(0:2*mpnw+6), &
  sum2(0:2*mpnw+6), td1(0:2*mpnw+6), td2(0:2*mpnw+6), tn1(0:2*mpnw+6), &
  tn2(0:2*mpnw+6), t1(0:2*mpnw+6), t2(0:2*mpnw+6), t3(0:2*mpnw+6), &
  t4(0:2*mpnw+6), t41(0:2*mpnw+6), t42(0:2*mpnw+6), t5(0:2*mpnw+6), &
  rra(0:2*mpnw+6), rr2(0:2*mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (rr) < mpnw + 4 .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELJR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 519)
endif

!   Check if PI has been precomputed.

mpnw1 = 2 * mpnw + 1
call mpmdc (mppicon, d1, n1, mpnw1)
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw1) then
  write (mpldb, 2) mpnw1
2 format ('*** MPBESSELJR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 520)
endif

call mpinitwds (f1, mpnw1)
call mpinitwds (f2, mpnw1)
call mpinitwds (sum1, mpnw1)
call mpinitwds (sum2, mpnw1)
call mpinitwds (td1, mpnw1)
call mpinitwds (tn1, mpnw1)
call mpinitwds (td2, mpnw1)
call mpinitwds (tn2, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t41, mpnw1)
call mpinitwds (t42, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (rra, mpnw1)
call mpinitwds (rr2, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw*mpnbt, eps, 4)

!   If QQ is integer, call mpbesseljnr; if qq < 0 and rr <= 0, then error.

call mpinfr (qq, t1, t2, mpnw)
i1 = mpsgn (qq)
i2 = mpsgn (rr)
if (mpsgn (t2) == 0) then
  call mpmdc (qq, d1, n1, mpnw)
  d1 = d1 * 2.d0**n1
  n1 = nint (d1)
  call mpbesseljnr (n1, rr, t3, mpnw)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (mpldb, 3)
3 format ('*** MPBESSELJR: First argument < 0 and second argument <= 0')
  call mpabrt ( 521)
endif

call mpmdc (rr, d1, n1, 4)
d1 = abs (d1) * 2.d0 ** n1

if (d1 < dfrac * mpnw1 * mpdpw) then
  mpnw1 = min (mpnw + nint (d1 / (dfrac * mpdpw)), 2 * mpnw + 1)
  call mpabs (rr, rra, mpnw1)
  call mpdmc (1.d0, 0, tn1, mpnw1)
  call mpdmc (1.d0, 0, f1, mpnw1)
  call mpadd (qq, f1, t1, mpnw1)
  call mpgammar (t1, f2, mpnw1)
  call mpmul (rra, rra, t2, mpnw1)
  call mpmuld (t2, 0.25d0, t1, mpnw1)

  call mpmul (f1, f2, td1, mpnw1)
  call mpdiv (tn1, td1, t2, mpnw1)
  call mpeq (t2, sum1, mpnw1)

  do k = 1, itrmax
    call mpmuld (f1, dble (k), t2, mpnw1)
    call mpeq (t2, f1, mpnw1)
    call mpdmc (dble (k), 0, t3, mpnw1)
    call mpadd (qq, t3, t4, mpnw1)
    call mpmul (f2, t4, t3, mpnw1)
    call mpeq (t3, f2, mpnw1)
    call mpmul (t1, tn1, t2, mpnw1)
    call mpneg (t2, tn1, mpnw1)
    call mpmul (f1, f2, td1, mpnw1)
    call mpdiv (tn1, td1, t2, mpnw1)
    call mpadd (sum1, t2, t3, mpnw1)
    call mpeq (t3, sum1, mpnw1)

    call mpabs (t2, tc1, 4)
    call mpmul (eps, sum1, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 100
  enddo

  write (mpldb, 4)
4 format ('*** MPBESSELJR: Loop end error 1')
  call mpabrt ( 522)

100 continue

  call mpmuld (rra, 0.5d0, t1, mpnw1)
  call mppower (t1, qq, t2, mpnw1)
  call mpmul (sum1, t2, t3, mpnw1)
else
! xa2 = xa**2
! t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
! tn1 = mpreal (1.d0, mpnw)
! tn2 = (t1 - 1.d0) / 8.d0
! td1 = mpreal (1.d0, mpnw)
! td2 = xa
! sum1 = tn1 / td1
! sum2 = tn2 / t32

  mpnw1 = mpnw + 1
  call mpabs (rr, rra, mpnw1)
  call mpmul (rra, rra, rr2, mpnw1)
  call mpmul (qq, qq, t2, mpnw1)
  call mpmuld (t2, 4.d0, t1, mpnw1)
  call mpdmc (1.d0, 0, tn1, mpnw1)
  call mpsub (t1, tn1, t2, mpnw1)
  call mpdivd (t2, 8.d0, tn2, mpnw1)
  call mpdmc (1.d0, 0, td1, mpnw1)
  call mpeq (rra, td2, mpnw1)
  call mpdiv (tn1, td1, sum1, mpnw1)
  call mpdiv (tn2, td2, sum2, mpnw1)

  do k = 1, itrmax
!   tn1 = -tn1 * (t1 - (2.d0*(2.d0*k-1.d0) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k-1.d0) - 1.d0)**2
    d2 = (2.d0*(2.d0*k) - 1.d0)**2
    call mpdmc (d1, 0, t2, mpnw1)
    call mpsub (t1, t2, t3, mpnw1)
    call mpdmc (d2, 0, t2, mpnw1)
    call mpsub (t1, t2, t5, mpnw1)
    call mpmul (t3, t5, t2, mpnw1)
    call mpmul (tn1, t2, t3, mpnw1)
    call mpneg (t3, tn1, mpnw1)

!   td1 = td1 * dble (2*k-1) * dble (2*k) * 64.d0 * xa2

    d1 = dble (2*k-1) * dble (2*k) * 64.d0
    call mpmuld (td1, d1, t2, mpnw1)
    call mpmul (t2, rr2, td1, mpnw1)

!   t41 = tn1 / td1
!   sum1 = sum1 + t41

    call mpdiv (tn1, td1, t41, mpnw1)
    call mpadd (sum1, t41, t2, mpnw1)
    call mpeq (t2, sum1, mpnw1)

!   tn2 = -tn2 * (t1 - (2.d0*(2.d0*k) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k+1.d0) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k) - 1.d0)**2
    d2 = (2.d0*(2.d0*k+1.d0) - 1.d0)**2
    call mpdmc (d1, 0, t2, mpnw1)
    call mpsub (t1, t2, t3, mpnw1)
    call mpdmc (d2, 0, t2, mpnw1)
    call mpsub (t1, t2, t5, mpnw1)
    call mpmul (t3, t5, t2, mpnw1)
    call mpmul (tn2, t2, t3, mpnw1)
    call mpneg (t3, tn2, mpnw1)

!   td2 = td2 * dble (2*k) * dble (2*k+1) * 64.d0 * xa2

    d1 = dble (2*k) * dble (2*k+1) * 64.d0
    call mpmuld (td2, d1, t2, mpnw1)
    call mpmul (t2, rr2, td2, mpnw1)

!  t42 = tn2 / td2
!  sum2 = sum2 + t42

    call mpdiv (tn2, td2, t42, mpnw1)
    call mpadd (sum2, t42, t2, mpnw1)
    call mpeq (t2, sum2, mpnw1)

!  if (abs (t41) / abs (sum1) < eps .and. abs (t42) / abs (sum2) < eps ) goto 110

    call mpabs (t41, tc1, 4)
    call mpmul (eps, sum1, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    call mpabs (t42, tc1, 4)
    call mpmul (eps, sum2, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic2, 4)
    if (ic1 <= 0 .and. ic2 <= 0) goto 110
  enddo

  write (mpldb, 5)
5 format ('*** MPBESSELJR: Loop end error 2')
  call mpabrt ( 523)

110 continue

! t1 = xa - 0.5d0 * nua * pi - 0.25d0 * pi
! besselj = sqrt (2.d0 / (pi * xa)) * (cos (t1) * sum1 - sin (t1) * sum2)

!   call mpmuld (mppicon, 0.5d0 * nua, t1, mpnw1)

  call mpmul (mppicon, qq, t2, mpnw1)
  call mpmuld (t2, 0.5d0, t1, mpnw1)
  call mpsub (rra, t1, t2, mpnw1)
  call mpmuld (mppicon, 0.25d0, t1, mpnw1)
  call mpsub (t2, t1, t3, mpnw1)
  call mpcssnr (t3, t41, t42, mpnw1)
  call mpmul (t41, sum1, t1, mpnw1)
  call mpmul (t42, sum2, t2, mpnw1)
  call mpsub (t1, t2, t5, mpnw1)
  call mpmul (mppicon, rra, t1, mpnw1)
  call mpdmc (2.d0, 0, t2, mpnw1)
  call mpdiv (t2, t1, t3, mpnw1)
  call mpsqrt (t3, t1, mpnw1)
  call mpmul (t1, t5, t3, mpnw1)
endif

! if (mod (nu, 2) /= 0) then
!  if (nu < 0 .and. x > 0.d0 .or. nu > 0 .and. x < 0.d0) besselj = - besselj

!   ic1 = mpsgn (rr)
!   if (nu < 0 .and. ic1 > 0 .or. nu > 0 .and. ic1 < 0) then
!     call mpneg (t3, t2, mpnw1)
!     call mpeq (t2, t3, mpnw1)
!   endif
! endif

120 continue

call mproun (t3, mpnw)
call mpeq (t3, ss, mpnw)

return
end subroutine mpbesseljr

subroutine mpbesselknr (nu, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselK (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.31.1 for modest RR,
!   and DLMF 10.40.2 for large RR, relative to precision.

implicit none
integer, intent(in):: nu, mpnw
integer ic1, itrmax, k, mpnw1, nua, n1
real (mprknd) dfrac, d1, egam, pi
parameter (itrmax = 1000000, dfrac = 1.5d0, egam = 0.5772156649015328606d0, &
  pi = 3.1415926535897932385d0)
integer (mpiknd), intent(in):: rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer (mpiknd) f1(0:mpnw+6), f2(0:mpnw+6), f3(0:mpnw+6), f4(0:mpnw+6), &
  f5(0:mpnw+6), sum1(0:mpnw+6), sum2(0:mpnw+6), sum3(0:mpnw+6), td(0:mpnw+6), &
  tn(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  rra(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (rr) < mpnw + 4 .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELKNR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 524)
endif

!   Check for RR = 0.

if (mpsgn (rr) == 0) then
  write (mpldb, 2)
2 format ('*** MPBESSELKNR: Second argument is zero')
  call mpabrt ( 525)
endif

!   Check if EGAMMA and PI have been precomputed.

mpnw1 = mpnw + 1
call mpmdc (mpegammacon, d1, n1, mpnw1)
if (n1 /= -1 .or. abs (d1 * 2.d0**n1 - egam) > mprdfz &
  .or. mpwprecr (mpegammacon) < mpnw1) then
  write (mpldb, 3) mpnw
3 format ('*** MPBESSELKNR: EGAMMA must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 526)
endif
call mpmdc (mppicon, d1, n1, mpnw1)
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw1) then
  write (mpldb, 4) mpnw1
4 format ('*** MPBESSELKNR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 527)
endif

nua = abs (nu)
mpnw1 = mpnw + 1
call mpinitwds (f1, mpnw1)
call mpinitwds (f2, mpnw1)
call mpinitwds (f3, mpnw1)
call mpinitwds (f4, mpnw1)
call mpinitwds (f5, mpnw1)
call mpinitwds (sum1, mpnw1)
call mpinitwds (sum2, mpnw1)
call mpinitwds (sum3, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (rra, mpnw1)
call mpinitwds (td, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)
call mpabs (rr, rra, mpnw1)
call mpmdc (rra, d1, n1, 4)
d1 = d1 * 2.d0 ** n1

if (d1 < dfrac * mpnw1 * mpdpw) then
  call mpmul (rra, rra, t2, mpnw1)
  call mpmuld (t2, 0.25d0, t1, mpnw1)
  call mpdmc (1.d0, 0, f1, mpnw1)
  call mpdmc (1.d0, 0, f2, mpnw1)
  call mpdmc (1.d0, 0, f3, mpnw1)
  call mpdmc (0.d0, 0,  sum1, mpnw1)

  do k = 1, nua - 1
    call mpmuld (f1, dble (k), t2, mpnw1)
    call mpeq (t2, f1, mpnw1)
  enddo

  do k = 0, nua - 1
    if (k > 0) then
      call mpdivd (f1, dble (nua - k), t2, mpnw1)
      call mpeq (t2, f1, mpnw1)
      call mpmul (t1, f2, t2, mpnw1)
      call mpneg (t2, f2, mpnw1)
      call mpmuld (f3, dble (k), t2, mpnw1)
      call mpeq (t2, f3, mpnw1)
    endif
    call mpmul (f1, f2, t3, mpnw1)
    call mpdiv (t3, f3, t2, mpnw1)
    call mpadd (sum1, t2, t3, mpnw1)
    call mpeq (t3, sum1, mpnw1)
  enddo

  call mpmuld (sum1, 0.5d0, t2, mpnw1)
  call mpmuld (rra, 0.5d0, t3, mpnw1)
  call mpnpwr (t3, nua, t4, mpnw1)
  call mpdiv (t2, t4, sum1, mpnw1)

  call mpmuld (rra, 0.5d0, t2, mpnw1)
  call mplog (t2, t3, mpnw1)
  d1 = (-1.d0) ** (nua + 1)
  call mpmuld (t3, d1, t2, mpnw1)
  call mpbesselinr (nua, rra, t3, mpnw1)
  call mpmul (t2, t3, sum2, mpnw1)

  call mpneg (mpegammacon, f1, mpnw1)
  call mpeq (f1, f2, mpnw1)
  call mpdmc (1.d0, 0, f3, mpnw1)
  call mpdmc (1.d0, 0, f4, mpnw1)
  call mpdmc (1.d0, 0, f5, mpnw1)

  do k = 1, nua
    call mpdmc (1.d0, 0, t2, mpnw1)
    call mpdivd (t2, dble (k), t3, mpnw1)
    call mpadd (f2, t3, t4, mpnw1)
    call mpeq (t4, f2, mpnw1)
    call mpmuld (f5, dble (k), t2, mpnw1)
    call mpeq (t2, f5, mpnw1)
  enddo

  call mpadd (f1, f2, t2, mpnw1)
  call mpmul (t2, f3, t3, mpnw1)
  call mpmul (f4, f5, t4, mpnw1)
  call mpdiv (t3, t4, sum3, mpnw1)

  do k = 1, itrmax
    call mpdmc (1.d0, 0, t2, mpnw1)
    call mpdivd (t2, dble (k), t3, mpnw1)
    call mpadd (f1, t3, t4, mpnw1)
    call mpeq (t4, f1, mpnw1)
    call mpdivd (t2, dble (nua + k), t3, mpnw1)
    call mpadd (f2, t3, t4, mpnw1)
    call mpeq (t4, f2, mpnw1)
    call mpmul (t1, f3, t2, mpnw1)
    call mpeq (t2, f3, mpnw1)
    call mpmuld (f4, dble (k), t2, mpnw1)
    call mpeq (t2, f4, mpnw1)
    call mpmuld (f5, dble (nua + k), t2, mpnw1)
    call mpeq (t2, f5, mpnw1)
    call mpadd (f1, f2, t2, mpnw1)
    call mpmul (t2, f3, t3, mpnw1)
    call mpmul (f4, f5, t4, mpnw1)
    call mpdiv (t3, t4, t2, mpnw1)
    call mpadd (sum3, t2, t3, mpnw1)
    call mpeq (t3, sum3, mpnw1)

    call mpabs (t2, tc1, 4)
    call mpmul (eps, sum3, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 100
  enddo

  write (mpldb, 5)
5 format ('*** MPBESSELKNR: Loop end error 1')
  call mpabrt ( 528)

100 continue

  call mpmuld (rra, 0.5d0, t2, mpnw1)
  call mpnpwr (t2, nua, t3, mpnw1)
  d1 = (-1.d0)**nua * 0.5d0
  call mpmuld (t3, d1, t4, mpnw1)
  call mpmul (t4, sum3, t2, mpnw1)
  call mpeq (t2, sum3, mpnw1)
  call mpadd (sum1, sum2, t2, mpnw1)
  call mpadd (t2, sum3, t3, mpnw1)
else
!  sum1 = mpreal (1.d0, mpnw)
!  t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
!  t2 = mpreal (1.d0, mpnw)
!  t3 = mpreal (1.d0, mpnw)

  call mpdmc (1.d0, 0, sum1, mpnw1)
  d1 = 4.d0 * dble (nua)**2
  call mpdmc (d1, 0, t1, mpnw1)
  call mpdmc (1.d0, 0, tn, mpnw1)
  call mpdmc (1.d0, 0, td, mpnw1)

  do k = 1, itrmax
!  t2 = t2 * (t1 - (2.d0*k - 1.d0)**2)
!  t3 = t3 * dble (k) * 8.d0 * xa
!  t4 = t2 / t3
!  sum1 = sum1 + t4

    d1 = 2.d0 * k - 1.d0
    call mpdmc (d1, 0, t2, mpnw1)
    call mpmul (t2, t2, t3, mpnw1)
    call mpsub (t1, t3, t2, mpnw1)
    call mpmul (tn, t2, t3, mpnw1)
    call mpeq (t3, tn, mpnw1)
    call mpmuld (rra, 8.d0 * k, t2, mpnw1)
    call mpmul (td, t2, t3,  mpnw1)
    call mpeq (t3, td, mpnw1)
    call mpdiv (tn, td, t4, mpnw1)
    call mpadd (sum1, t4, t3, mpnw1)
    call mpeq (t3, sum1, mpnw1)

!   if (abs (t4) / abs (sum1) < eps) goto 110

    call mpabs (t4, tc1, 4)
    call mpmul (eps, sum1, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 110
  enddo

write (mpldb, 6)
6 format ('*** MPBESSELKNR: Loop end error 2')
call mpabrt ( 529)

110 continue

! t1 = sqrt (mppi (mpnw) / (2.d0 * xa)) / exp (xa)
! besseli = t1 * sum1

  call mpexp (rra, t1, mpnw1)
  call mpmuld (rra, 2.d0, t2, mpnw1)
  call mpdiv (mppicon, t2, t3, mpnw1)
  call mpsqrt (t3, t4, mpnw1)
  call mpdiv (t4, t1, t2, mpnw1)
  call mpmul (t2, sum1, t3, mpnw1)
endif

! if (x < 0.d0 .and. mod (nu, 2) /= 0) besselk = - besselk

if (mpsgn (rr) < 0 .and. mod (nu, 2) /= 0) then
  call mpneg (t3, t4, mpnw1)
  call mpeq (t4, t3, mpnw1)
endif
call mproun (t3, mpnw)
call mpeq (t3, ss, mpnw)
return
end subroutine mpbesselknr

subroutine mpbesselkr (qq, rr, ss, mpnw)

!   This evaluates the Bessel function BesselK (QQ,RR) for QQ and RR
!   both MPR. This uses DLMF formula 10.27.4.

implicit none
integer, intent(in):: mpnw
integer i1, i2, mpnw1, n1
real (mprknd) d1
integer (mpiknd), intent(in):: qq(0:), rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (qq) < mpnw + 4 .or. mpspacer (rr) < mpnw + 4 &
  .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELKR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 530)
endif

mpnw1 = mpnw + 1
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)

!   If QQ is integer, call mpbesselknr; if qq < 0 and rr <= 0, then error.

call mpinfr (qq, t1, t2, mpnw)
i1 = mpsgn (qq)
i2 = mpsgn (rr)
if (mpsgn (t2) == 0) then
  call mpmdc (qq, d1, n1, mpnw)
  d1 = d1 * 2.d0**n1
  n1 = nint (d1)
  call mpbesselknr (n1, rr, t1, mpnw)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (mpldb, 2)
2 format ('*** MPBESSELKR: First argument < 0 and second argument <= 0')
  call mpabrt ( 531)
endif

call mpneg (qq, t1, mpnw1)
call mpbesselir (t1, rr, t2, mpnw1)
call mpbesselir (qq, rr, t3, mpnw1)
call mpsub (t2, t3, t4, mpnw1)
call mpmul (qq, mppicon, t1, mpnw1)
call mpcssnr (t1, t2, t3, mpnw1)
call mpdiv (t4, t3, t2, mpnw1)
call mpmul (mppicon, t2, t3, mpnw1)
call mpmuld (t3, 0.5d0, t1, mpnw1)

120 continue

call mproun (t1, mpnw)
call mpeq (t1, ss, mpnw)
return
end subroutine mpbesselkr

subroutine mpbesselynr (nu, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselY (NU,RR).
!   NU is an integer. The algorithm is DLMF formula 10.8.1 for modest RR,
!   and DLMF 10.17.4 for large RR, relative to precision.

implicit none
integer, intent(in):: nu, mpnw
integer ic1, ic2, itrmax, k, mpnw1, nua, n1
real (mprknd) dfrac, d1, d2, egam, pi
parameter (itrmax = 1000000, dfrac = 1.5d0, egam = 0.5772156649015328606d0, &
  pi = 3.1415926535897932385d0)
integer (mpiknd), intent(in):: rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer (mpiknd) f1(0:2*mpnw+6), f2(0:2*mpnw+6), f3(0:2*mpnw+6), f4(0:2*mpnw+6), &
  f5(0:2*mpnw+6), rra(0:2*mpnw+6), rr2(0:2*mpnw+6), sum1(0:2*mpnw+6), &
  sum2(0:2*mpnw+6), sum3(0:2*mpnw+6), td1(0:2*mpnw+6), td2(0:2*mpnw+6), &
  tn1(0:2*mpnw+6), tn2(0:2*mpnw+6), t1(0:2*mpnw+6), t2(0:2*mpnw+6), &
  t3(0:2*mpnw+6), t4(0:2*mpnw+6), t41(0:2*mpnw+6), t42(0:2*mpnw+6), &
  t5(0:2*mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (rr) < mpnw + 4 .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELYNR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 532)
endif

!   Check for RR = 0.

if (mpsgn (rr) == 0) then
  write (mpldb, 2)
2 format ('*** MPBESSELYNR: argument is negative or too large')
  call mpabrt ( 533)
endif

!   Check if EGAMMA and PI have been precomputed.

mpnw1 = 2 * mpnw + 1
call mpmdc (mpegammacon, d1, n1, mpnw1)
if (n1 /= -1 .or. abs (d1 * 2.d0**n1 - egam) > mprdfz &
  .or. mpwprecr (mpegammacon) < mpnw1) then
  write (mpldb, 3) mpnw
3 format ('*** MPBESSELYNR: EGAMMA must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 534)
endif
call mpmdc (mppicon, d1, n1, mpnw1)
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw1) then
  write (mpldb, 4) mpnw1
4 format ('*** MPBESSELYNR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 535)
endif

nua = abs (nu)
call mpinitwds (f1, mpnw1)
call mpinitwds (f2, mpnw1)
call mpinitwds (f3, mpnw1)
call mpinitwds (f4, mpnw1)
call mpinitwds (f5, mpnw1)
call mpinitwds (rra, mpnw1)
call mpinitwds (rr2, mpnw1)
call mpinitwds (sum1, mpnw1)
call mpinitwds (sum2, mpnw1)
call mpinitwds (sum3, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t41, mpnw1)
call mpinitwds (t42, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (td1, mpnw1)
call mpinitwds (td2, mpnw1)
call mpinitwds (tn1, mpnw1)
call mpinitwds (tn2, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw*mpnbt, eps, 4)
call mpmdc (rr, d1, n1, 4)
d1 = abs (d1) * 2.d0 ** n1

if (d1 < dfrac * mpnw1 * mpdpw) then
  mpnw1 = min (mpnw + nint (d1 / (dfrac * mpdpw)), 2 * mpnw + 1)
  call mpabs (rr, rra, mpnw1)
  call mpmul (rra, rra, t2, mpnw1)
  call mpmuld (t2, 0.25d0, t1, mpnw1)
  call mpdmc (1.d0, 0, f1, mpnw1)
  call mpdmc (1.d0, 0, f2, mpnw1)
  call mpdmc (1.d0, 0, f3, mpnw1)
  call mpdmc (0.d0, 0,  sum1, mpnw1)

  do k = 1, nua - 1
    call mpmuld (f1, dble (k), t2, mpnw1)
    call mpeq (t2, f1, mpnw1)
  enddo

  do k = 0, nua - 1
    if (k > 0) then
      call mpdivd (f1, dble (nua - k), t2, mpnw1)
      call mpeq (t2, f1, mpnw1)
      call mpmul (t1, f2, t2, mpnw1)
      call mpeq (t2, f2, mpnw1)
      call mpmuld (f3, dble (k), t2, mpnw1)
      call mpeq (t2, f3, mpnw1)
    endif
    call mpmul (f1, f2, t3, mpnw1)
    call mpdiv (t3, f3, t2, mpnw1)
    call mpadd (sum1, t2, t3, mpnw1)
    call mpeq (t3, sum1, mpnw1)
  enddo

  call mpmuld (rra, 0.5d0, t3, mpnw1)
  call mpnpwr (t3, nua, t4, mpnw1)
  call mpdiv (sum1, t4, t3, mpnw1)
  call mpneg (t3, sum1, mpnw1)

  call mpmuld (rra, 0.5d0, t2, mpnw1)
  call mplog (t2, t3, mpnw1)
  call mpmuld (t3, 2.d0, t2, mpnw1)
  call mpbesseljnr (nua, rra, t3, mpnw1)
  call mpmul (t2, t3, sum2, mpnw1)

  call mpneg (mpegammacon, f1, mpnw1)
  call mpeq (f1, f2, mpnw1)
  call mpdmc (1.d0, 0, f3, mpnw1)
  call mpdmc (1.d0, 0, f4, mpnw1)
  call mpdmc (1.d0, 0, f5, mpnw1)

  do k = 1, nua
    call mpdmc (1.d0, 0, t2, mpnw1)
    call mpdivd (t2, dble (k), t3, mpnw1)
    call mpadd (f2, t3, t4, mpnw1)
    call mpeq (t4, f2, mpnw1)
    call mpmuld (f5, dble (k), t2, mpnw1)
    call mpeq (t2, f5, mpnw1)
  enddo

  call mpadd (f1, f2, t2, mpnw1)
  call mpmul (t2, f3, t3, mpnw1)
  call mpmul (f4, f5, t4, mpnw1)
  call mpdiv (t3, t4, sum3, mpnw1)

  do k = 1, itrmax
    call mpdmc (1.d0, 0, t2, mpnw1)
    call mpdivd (t2, dble (k), t3, mpnw1)
    call mpadd (f1, t3, t4, mpnw1)
    call mpeq (t4, f1, mpnw1)
    call mpdivd (t2, dble (nua + k), t3, mpnw1)
    call mpadd (f2, t3, t4, mpnw1)
    call mpeq (t4, f2, mpnw1)
    call mpmul (t1, f3, t2, mpnw1)
    call mpneg (t2, f3, mpnw1)
    call mpmuld (f4, dble (k), t2, mpnw1)
    call mpeq (t2, f4, mpnw1)
    call mpmuld (f5, dble (nua + k), t2, mpnw1)
    call mpeq (t2, f5, mpnw1)
    call mpadd (f1, f2, t2, mpnw1)
    call mpmul (t2, f3, t3, mpnw1)
    call mpmul (f4, f5, t4, mpnw1)
    call mpdiv (t3, t4, t2, mpnw1)
    call mpadd (sum3, t2, t3, mpnw1)
    call mpeq (t3, sum3, mpnw1)

    call mpabs (t2, tc1, 4)
    call mpmul (eps, sum3, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 100
  enddo

  write (mpldb, 6)
6 format ('*** MPBESSELYNR: Loop end error 1')
  call mpabrt ( 536)

100 continue

  call mpmuld (rra, 0.5d0, t2, mpnw1)
  call mpnpwr (t2, nua, t3, mpnw1)
  call mpmul (t3, sum3, t2, mpnw1)
  call mpneg (t2, sum3, mpnw1)

  call mpadd (sum1, sum2, t2, mpnw1)
  call mpadd (t2, sum3, t4, mpnw1)
  call mpeq (mppicon, t2, mpnw1)
  call mpdiv (t4, t2, t3, mpnw1)
else
! xa2 = xa**2
! t1 = mpreal (4.d0 * dble (nua)**2, mpnw)
! tn1 = mpreal (1.d0, mpnw)
! tn2 = (t1 - 1.d0) / 8.d0
! td1 = mpreal (1.d0, mpnw)
! td2 = xa
! sum1 = tn1 / td1
! sum2 = tn2 / t32

  mpnw1 = mpnw + 1
  call mpabs (rr, rra, mpnw1)
  call mpmul (rra, rra, rr2, mpnw1)
  d1 = 4.d0 * dble (nua)**2
  call mpdmc (d1, 0, t1, mpnw1)
  call mpdmc (1.d0, 0, tn1, mpnw1)
  call mpsub (t1, tn1, t2, mpnw1)
  call mpdivd (t2, 8.d0, tn2, mpnw1)
  call mpdmc (1.d0, 0, td1, mpnw1)
  call mpeq (rra, td2, mpnw1)
  call mpdiv (tn1, td1, sum1, mpnw1)
  call mpdiv (tn2, td2, sum2, mpnw1)

  do k = 1, itrmax
!   tn1 = -tn1 * (t1 - (2.d0*(2.d0*k-1.d0) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k-1.d0) - 1.d0)**2
    d2 = (2.d0*(2.d0*k) - 1.d0)**2
    call mpdmc (d1, 0, t2, mpnw1)
    call mpsub (t1, t2, t3, mpnw1)
    call mpdmc (d2, 0, t2, mpnw1)
    call mpsub (t1, t2, t5, mpnw1)
    call mpmul (t3, t5, t2, mpnw1)
    call mpmul (tn1, t2, t3, mpnw1)
    call mpneg (t3, tn1, mpnw1)

!   td1 = td1 * dble (2*k-1) * dble (2*k) * 64.d0 * xa2

    d1 = dble (2*k-1) * dble (2*k) * 64.d0
    call mpmuld (td1, d1, t2, mpnw1)
    call mpmul (t2, rr2, td1, mpnw1)

!   t41 = tn1 / td1
!   sum1 = sum1 + t41

    call mpdiv (tn1, td1, t41, mpnw1)
    call mpadd (sum1, t41, t2, mpnw1)
    call mpeq (t2, sum1, mpnw1)

!   tn2 = -tn2 * (t1 - (2.d0*(2.d0*k) - 1.d0)**2) * (t1 - (2.d0*(2.d0*k+1.d0) - 1.d0)**2)

    d1 = (2.d0*(2.d0*k) - 1.d0)**2
    d2 = (2.d0*(2.d0*k+1.d0) - 1.d0)**2
    call mpdmc (d1, 0, t2, mpnw1)
    call mpsub (t1, t2, t3, mpnw1)
    call mpdmc (d2, 0, t2, mpnw1)
    call mpsub (t1, t2, t5, mpnw1)
    call mpmul (t3, t5, t2, mpnw1)
    call mpmul (tn2, t2, t3, mpnw1)
    call mpneg (t3, tn2, mpnw1)

!   td2 = td2 * dble (2*k) * dble (2*k+1) * 64.d0 * xa2

    d1 = dble (2*k) * dble (2*k+1) * 64.d0
    call mpmuld (td2, d1, t2, mpnw1)
    call mpmul (t2, rr2, td2, mpnw1)

!  t42 = tn2 / td2
!  sum2 = sum2 + t42

    call mpdiv (tn2, td2, t42, mpnw1)
    call mpadd (sum2, t42, t2, mpnw1)
    call mpeq (t2, sum2, mpnw1)

!  if (abs (t41) / abs (sum1) < eps .and. abs (t42) / abs (sum2) < eps ) goto 110

    call mpabs (t41, tc1, 4)
    call mpmul (eps, sum1, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    call mpabs (t42, tc1, 4)
    call mpmul (eps, sum2, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic2, 4)
    if (ic1 <= 0 .and. ic2 <= 0) goto 110
  enddo

  write (mpldb, 5)
5 format ('*** MPBESSELYNR: Loop end error 2')
  call mpabrt ( 537)

110 continue

! t1 = xa - 0.5d0 * nua * pi - 0.25d0 * pi
! besselj = sqrt (2.d0 / (pi * xa)) * (cos (t1) * sum1 - sin (t1) * sum2)

  call mpmuld (mppicon, 0.5d0 * nua, t1, mpnw1)
  call mpsub (rra, t1, t2, mpnw1)
  call mpmuld (mppicon, 0.25d0, t1, mpnw1)
  call mpsub (t2, t1, t3, mpnw1)
  call mpcssnr (t3, t41, t42, mpnw1)
  call mpmul (t42, sum1, t1, mpnw1)
  call mpmul (t41, sum2, t2, mpnw1)
  call mpadd (t1, t2, t5, mpnw1)
  call mpmul (mppicon, rra, t1, mpnw1)
  call mpdmc (2.d0, 0, t2, mpnw1)
  call mpdiv (t2, t1, t3, mpnw1)
  call mpsqrt (t3, t1, mpnw1)
  call mpmul (t1, t5, t3, mpnw1)
endif

if (mod (nu, 2) /= 0) then
!   if (nu < 0 .and. x > 0.d0 .or. nu > 0 .and. x < 0.d0) bessely = - bessely

  ic1 = mpsgn (rr)
  if (nu < 0 .and. ic1 > 0 .or. nu > 0 .and. ic1 < 0) then
    call mpneg (t3, t4, mpnw1)
    call mpeq (t4, t3, mpnw1)
  endif
endif

call mproun (t3, mpnw)
call mpeq (t3, ss, mpnw)
return
end subroutine mpbesselynr

subroutine mpbesselyr (qq, rr, ss, mpnw)

!   This evaluates the modified Bessel function BesselY (QQ,RR).
!   NU is an integer. The algorithm is DLMF formula 10.2.2.

implicit none
integer, intent(in):: mpnw
integer i1, i2, mpnw1, n1
real (mprknd) d1
integer (mpiknd), intent(in):: qq(0:), rr(0:)
integer (mpiknd), intent(out):: ss(0:)
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (qq) < mpnw + 4 .or. mpspacer (rr) < mpnw + 4 &
  .or. mpspacer (ss) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPBESSELYR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 538)
endif

mpnw1 = mpnw + 1
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)

!   If QQ is integer, call mpbesselynr; if qq < 0 and rr <= 0, then error.

call mpinfr (qq, t1, t2, mpnw)
i1 = mpsgn (qq)
i2 = mpsgn (rr)
if (mpsgn (t2) == 0) then
  call mpmdc (qq, d1, n1, mpnw)
  d1 = d1 * 2.d0**n1
  n1 = nint (d1)
  call mpbesselynr (n1, rr, t1, mpnw)
  goto 120
elseif (i1 < 0 .and. i2 <= 0) then
  write (mpldb, 2)
2 format ('*** MPBESSELYR: First argument < 0 and second argument <= 0')
  call mpabrt ( 539)
endif

call mpmul (qq, mppicon, t1, mpnw1)
call mpcssnr (t1, t2, t3, mpnw1)
call mpbesseljr (qq, rr, t4, mpnw1)
call mpmul (t4, t2, t1, mpnw1)
call mpneg (qq, t2, mpnw1)
call mpbesseljr (t2, rr, t4, mpnw1)
call mpsub (t1, t4, t2, mpnw1)
call mpdiv (t2, t3, t1, mpnw1)

120 continue

call mproun (t1, mpnw)
call mpeq (t1, ss, mpnw)
return
end subroutine mpbesselyr

subroutine mpdigammabe (nb1, nb2, berne, x, y, mpnw)

!  This evaluates the digamma function, using asymptotic formula DLMF 5.11.2:
!  dig(x) ~ log(x) - 1/(2*x) - Sum_{k=1}^inf B[2k] / (2*k*x^(2*k)).
!  Before using this formula, the recursion dig(x+1) = dig(x) + 1/x is used
!  to shift the argument up by IQQ, where IQQ is set based on MPNW below.
!  The array berne contains precomputed even Bernoulli numbers (see MPBERNER
!  above). Its dimensions must be as shown below. NB2 must be greater than
!  1.4 x precision in decimal digits.

implicit none
integer, intent (in):: nb1, nb2, mpnw
integer (mpiknd), intent(in):: berne(0:nb1+5,nb2), x(0:)
integer (mpiknd), intent(out):: y(0:)
integer k, i1, i2, ic1, iqq, mpnw1, n1
real (mprknd) dber, dfrac, d1
parameter (dber = 1.4d0, dfrac = 0.4d0)
integer (mpiknd) f1(0:mpnw+6), sum1(0:mpnw+6), sum2(0:mpnw+6), &
  t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  t5(0:mpnw+6), xq(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (x) < mpnw + 4 .or. mpspacer (y) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPDIGAMMABE: Uninitialized or inadequately sized arrays')
  call mpabrt ( 540)
endif

mpnw1 = mpnw + 1
iqq = dfrac * mpnw1 * mpdpw
call mpinitwds (f1, mpnw1)
call mpinitwds (sum1, mpnw1)
call mpinitwds (sum2, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (xq, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw*mpnbt, eps, 4)
call mpdmc (1.d0, 0, f1, mpnw1)

!   Check if argument is less than or equal to 0 -- undefined.

if (mpsgn (x) <= 0) then
  write (mpldb, 2)
2 format ('*** MPDIGAMMABE: Argument is less than or equal to 0')
  call mpabrt ( 541)
endif

!   Check if berne array has been initialized.

i1 = 1000000000
i2 = 1000000000

do k = 1, nb2
  i1 = min (i1, mpspacer (berne(0:nb1+5,k)))
  i2 = min (i2, mpwprecr (berne(0:nb1+5,k)))
enddo

call mpmdc (berne(0:nb1+5,1), d1, n1, mpnw)
d1 = d1 * 2.d0 ** n1

if (i1 < mpnw + 6 .or. i2 < mpnw .or. abs (d1 - 1.d0 / 6.d0) > mprdfz .or. &
  nb2 < int (dber * mpdpw * (mpnw - 2))) then
  write (mpldb, 3) int (dber * mpdpw * (mpnw - 2))
3 format ('*** MPDIGAMMABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries using MPBERNE or MPBERNER.')
  call mpabrt ( 542)
endif

! sum1 = mpreal (0.d0, nwds)
! sum2 = mpreal (0.d0, nwds)
! xq = x + dble (iqq)

call mpdmc (0.d0, 0, sum1, mpnw1)
call mpdmc (0.d0, 0, sum2, mpnw1)
call mpdmc (dble (iqq), 0, t1, mpnw1)
call mpadd (x, t1, t2, mpnw1)
call mpeq (t2, xq, mpnw1)

do k = 0, iqq - 1
!   sum1 = sum1 + f1 / (x + dble (k))

  call mpdmc (dble (k), 0, t1, mpnw1)
  call mpadd (x, t1, t2, mpnw1)
  call mpdiv (f1, t2, t3, mpnw1)
  call mpadd (sum1, t3, t1, mpnw1)
  call mpeq (t1, sum1, mpnw1)
enddo

! t1 = mpreal (1.d0, nwds)
! t2 = xq ** 2

call mpdmc (1.d0, 0, t1, mpnw1)
call mpmul (xq, xq, t2, mpnw1)

do k = 1, nb2
!  t1 = t1 * t2
!  t3 = bb(k) / (2.d0 * dble (k) * t1)
!  sum2 = sum2 + t3

  call mpmul (t1, t2, t3, mpnw1)
  call mpeq (t3, t1, mpnw1)
  call mpmuld (t1, 2.d0 * dble (k), t4, mpnw1)
  call mpdiv (berne(0:nb1+5,k), t4, t3, mpnw1)
  call mpadd (sum2, t3, t4, mpnw1)
  call mpeq (t4, sum2, mpnw1)

!  if (abs (t3 / sum2) < eps) goto 100

  call mpabs (t3, tc1, 4)
  call mpmul (eps, sum2, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic1, 4)
  if (ic1 <= 0) goto 110
enddo

write (mpldb, 4)
4 format ('*** MPDIGAMMABE: Loop end error: Increase NB2')
call mpabrt ( 543)

110 continue

! digammax = -sum1 + log (xq) - 1.d0 / (2.d0 * xq) - sum2

call mpneg (sum1, t1, mpnw1)
call mplog (xq, t2, mpnw1)
call mpadd (t1, t2, t3, mpnw1)
call mpmuld (xq, 2.d0, t4, mpnw1)
call mpdiv (f1, t4, t5, mpnw1)
call mpsub (t3, t5, t2, mpnw1)
call mpsub (t2, sum2, t1, mpnw1)
call mproun (t1, mpnw)
call mpeq (t1, y, mpnw)
return
end subroutine mpdigammabe

subroutine mperfr (z, terf, mpnw)

!   This evaluates the erf function, using a combination of two series.
!   In particular, the algorithm is (where B = (mpnw + 1) * mpnbt, and
!   dcon is a constant defined below):

!   if (t == 0) then
!     erf = 0
!   elseif (z > sqrt(B*log(2))) then
!     erf = 1
!   elseif (z < -sqrt(B*log(2))) then
!     erf = -1
!   elseif (abs(z) < B/dcon + 8) then
!     erf = 2 / (sqrt(pi)*exp(z^2)) * Sum_{k>=0} 2^k * z^(2*k+1)
!             / (1.3....(2*k+1))
!   else
!     erf = 1 - 1 / (sqrt(pi)*exp(z^2))
!             * Sum_{k>=0} (-1)^k * (1.3...(2*k-1)) / (2^k * z^(2*k+1))
!   endif

implicit none
integer, intent(in):: mpnw
integer ic1, ic2, ic3, itrmx, k, mpnw1, nbt, n1
real (mprknd) dcon, d1, d2, pi
parameter (dcon = 100.d0, itrmx = 100000, pi = 3.1415926535897932385d0)
integer (mpiknd), intent(in):: z(0:)
integer (mpiknd), intent(out):: terf(0:)
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  t5(0:mpnw+6), t6(0:mpnw+6), t7(0:mpnw+6), z2(0:mpnw+6), tc1(0:9), &
  tc2(0:9), tc3(0:9), eps(0:9)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (z) < mpnw + 4 .or. mpspacer (terf) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPERFR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 544)
endif

mpnw1 = mpnw + 1

!   Check if Pi has been precomputed.

call mpmdc (mppicon, d1, n1, mpnw1)
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw1) then
  write (mpldb, 4) mpnw1
4 format ('*** MPERFR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 545)
endif

call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (t6, mpnw1)
call mpinitwds (t7, mpnw1)
call mpinitwds (z2, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

nbt = mpnw * mpnbt
d1 = aint (1.d0 + sqrt (nbt * log (2.d0)))
d2 = aint (nbt / dcon + 8.d0)
call mpdmc (d1, 0, t1, mpnw1)
call mpdmc (d2, 0, t2, mpnw1)
call mpcpr (z, t1, ic1, mpnw1)
! t1(2) = - t1(2)
call mpneg (t1, t3, mpnw1)
call mpeq (t3, t1, mpnw1)
call mpcpr (z, t1, ic2, mpnw1)
call mpcpr (z, t2, ic3, mpnw1)

if (mpsgn (z) == 0) then
  call mpdmc (0.d0, 0, terf, mpnw)
elseif (ic1 > 0) then
  call mpdmc (1.d0, 0, terf, mpnw)
elseif (ic2 < 0) then
  call mpdmc (-1.d0, 0, terf, mpnw)
elseif (ic3 < 0) then
  call mpmul (z, z, z2, mpnw1)
  call mpdmc (0.d0, 0, t1, mpnw1)
  call mpeq (z, t2, mpnw1)
  call mpdmc (1.d0, 0, t3, mpnw1)
  call mpdmc (1.d10, 0, t5, 4)

  do k = 0, itrmx
    if (k > 0) then
      call mpmuld (z2, 2.d0, t6, mpnw1)
      call mpmul (t6, t2, t7, mpnw1)
      call mpeq (t7, t2, mpnw1)
      d1 = 2.d0 * k + 1.d0
      call mpmuld (t3, d1, t6, mpnw1)
      call mpeq (t6, t3, mpnw1)
    endif

    call mpdiv (t2, t3, t4, mpnw1)
    call mpadd (t1, t4, t6, mpnw1)
    call mpeq (t6, t1, mpnw1)
    call mpdiv (t4, t1, t6, 4)
    call mpcpr (t6, eps, ic1, 4)
    call mpcpr (t6, t5, ic2, 4)
    if (ic1 <= 0 .or. ic2 >= 0) goto 120
    call mpeq (t6, t5, 4)
  enddo

write (mpldb, 3) 1, itrmx
3 format ('*** MPERFR: iteration limit exceeded',2i10)
call mpabrt ( 546)

120 continue

  call mpmuld (t1, 2.d0, t3, mpnw1)
  call mpsqrt (mppicon, t4, mpnw1)
  call mpexp (z2, t5, mpnw1)
  call mpmul (t4, t5, t6, mpnw1)
  call mpdiv (t3, t6, t7, mpnw1)
  call mproun (t7, mpnw)
  call mpeq (t7, terf, mpnw)
else
  call mpmul (z, z, z2, mpnw1)
  call mpdmc (0.d0, 0, t1, mpnw1)
  call mpdmc (1.d0, 0, t2, mpnw1)
  call mpabs (z, t3, mpnw1)
  call mpdmc (1.d10, 0, t5, 4)

  do k = 0, itrmx
    if (k > 0) then
      d1 = -(2.d0 * k - 1.d0)
      call mpmuld (t2, d1, t6, mpnw1)
      call mpeq (t6, t2, mpnw1)
      call mpmul (t2, t3, t6, mpnw1)
      call mpeq (t6, t3, mpnw1)
    endif

    call mpdiv (t2, t3, t4, mpnw1)
    call mpadd (t1, t4, t6, mpnw1)
    call mpeq (t6, t1, mpnw1)
    call mpdiv (t4, t1, t6, 4)
    call mpcpr (t6, eps, ic1, 4)
    call mpcpr (t6, t5, ic2, 4)
    if (ic1 <= 0 .or. ic2 >= 0) goto 130
    call mpeq (t6, t5, 4)
  enddo

write (mpldb, 3) 2, itrmx
call mpabrt ( 547)

130 continue

  call mpdmc (1.d0, 0, t2, mpnw1)
  call mpsqrt (mppicon, t3, mpnw1)
  call mpexp (z2, t4, mpnw1)
  call mpmul (t3, t4, t5, mpnw1)
  call mpdiv (t1, t5, t6, mpnw1)
  call mpsub (t2, t6, t7, mpnw1)
  call mproun (t7, mpnw)
  call mpeq (t7, terf, mpnw)
  if (mpsgn (z) < 0) then
    call mpneg (terf, t6, mpnw)
    call mpeq (t6, terf, mpnw)
  endif
endif

return
end subroutine mperfr

subroutine mperfcr (z, terfc, mpnw)

!   This evaluates the erfc function, using a combination of two series.
!   In particular, the algorithm is (where B = (mpnw + 1) * mpnbt, and
!   dcon is a constant defined below):

!   if (t == 0) then
!     erfc = 1
!   elseif (z > sqrt(B*log(2))) then
!     erfc = 0
!   elseif (z < -sqrt(B*log(2))) then
!     erfc = 2
!   elseif (abs(z) < B/dcon + 8) then
!     erfc = 1 - 2 / (sqrt(pi)*exp(z^2)) * Sum_{k>=0} 2^k * z^(2*k+1)
!               / (1.3....(2*k+1))
!   else
!     erfc = 1 / (sqrt(pi)*exp(z^2))
!             * Sum_{k>=0} (-1)^k * (1.3...(2*k-1)) / (2^k * z^(2*k+1))
!   endif

implicit none
integer, intent(in):: mpnw
integer ic1, ic2, ic3, itrmx, k, mpnw1, nbt, n1
real (mprknd) dcon, d1, d2, pi
parameter (dcon = 100.d0, itrmx = 100000, pi = 3.1415926535897932385d0)
integer (mpiknd), intent(in):: z(0:)
integer (mpiknd), intent(out):: terfc(0:)
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  t5(0:mpnw+6), t6(0:mpnw+6), t7(0:mpnw+6), z2(0:mpnw+6), tc1(0:9), &
  tc2(0:9), tc3(0:9), eps(0:9)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (z) < mpnw + 4 .or. mpspacer (terfc) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPERFCR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 548)
endif

mpnw1 = mpnw + 1

!   Check if Pi has been precomputed.

call mpmdc (mppicon, d1, n1, mpnw1)
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw1) then
  write (mpldb, 4) mpnw1
4 format ('*** MPERFCR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 549)
endif

call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (t6, mpnw1)
call mpinitwds (t7, mpnw1)
call mpinitwds (z2, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

nbt = mpnw * mpnbt
d1 = aint (1.d0 + sqrt (nbt * log (2.d0)))
d2 = aint (nbt / dcon + 8.d0)
call mpdmc (d1, 0, t1, mpnw1)
call mpdmc (d2, 0, t2, mpnw1)
call mpcpr (z, t1, ic1, mpnw1)
call mpneg (t1, t3, mpnw1)
call mpeq (t3, t1,  mpnw1)
call mpcpr (z, t1, ic2, mpnw1)
call mpcpr (z, t2, ic3, mpnw1)

if (mpsgn (z) == 0) then
  call mpdmc (1.d0, 0, terfc, mpnw)
elseif (ic1 > 0) then
  call mpdmc (0.d0, 0, terfc, mpnw)
elseif (ic2 < 0) then
  call mpdmc (2.d0, 0, terfc, mpnw)
elseif (ic3 < 0) then
  call mpmul (z, z, z2, mpnw1)
  call mpdmc (0.d0, 0, t1, mpnw1)
  call mpeq (z, t2, mpnw1)
  call mpdmc (1.d0, 0, t3, mpnw1)
  call mpdmc (1.d10, 0, t5, 4)

  do k = 0, itrmx
    if (k > 0) then
      call mpmuld (z2, 2.d0, t6, mpnw1)
      call mpmul (t6, t2, t7, mpnw1)
      call mpeq (t7, t2, mpnw1)
      d1 = 2.d0 * k + 1.d0
      call mpmuld (t3, d1, t6, mpnw1)
      call mpeq (t6, t3, mpnw1)
    endif

    call mpdiv (t2, t3, t4, mpnw1)
    call mpadd (t1, t4, t6, mpnw1)
    call mpeq (t6, t1, mpnw1)
    call mpdiv (t4, t1, t6, 4)
    call mpcpr (t6, eps, ic1, 4)
    call mpcpr (t6, t5, ic2, 4)
    if (ic1 <= 0 .or. ic2 >= 0) goto 120
    call mpeq (t6, t5, 4)
  enddo

write (mpldb, 3) 1, itrmx
3 format ('*** MPERFCR: iteration limit exceeded',2i10)
call mpabrt ( 550)

120 continue

  call mpdmc (1.d0, 0, t2, mpnw1)
  call mpmuld (t1, 2.d0, t3, mpnw1)
  call mpsqrt (mppicon, t4, mpnw1)
  call mpexp (z2, t5, mpnw1)
  call mpmul (t4, t5, t6, mpnw1)
  call mpdiv (t3, t6, t7, mpnw1)
  call mpsub (t2, t7, t6, mpnw1)
  call mproun (t6, mpnw)
  call mpeq (t6, terfc, mpnw)
else
  call mpmul (z, z, z2, mpnw1)
  call mpdmc (0.d0, 0, t1, mpnw1)
  call mpdmc (1.d0, 0, t2, mpnw1)
  call mpabs (z, t3, mpnw1)
  call mpdmc (1.d10, 0, t5, 4)

  do k = 0, itrmx
    if (k > 0) then
      d1 = -(2.d0 * k - 1.d0)
      call mpmuld (t2, d1, t6, mpnw1)
      call mpeq (t6, t2, mpnw1)
      call mpmul (t2, t3, t6, mpnw1)
      call mpeq (t6, t3, mpnw1)
    endif

    call mpdiv (t2, t3, t4, mpnw1)
    call mpadd (t1, t4, t6, mpnw1)
    call mpeq (t6, t1, mpnw1)
    call mpdiv (t4, t1, t6, 4)
    call mpcpr (t6, eps, ic1, 4)
    call mpcpr (t6, t5, ic2, 4)
    if (ic1 <= 0 .or. ic2 >= 0) goto 130
    call mpeq (t6, t5, 4)
  enddo

write (mpldb, 3) 2, itrmx
call mpabrt ( 551)

130 continue

  call mpsqrt (mppicon, t3, mpnw1)
  call mpexp (z2, t4, mpnw1)
  call mpmul (t3, t4, t5, mpnw1)
  call mpdiv (t1, t5, t6, mpnw1)
  if (mpsgn (z) < 0) then
    call mpdmc (2.d0, 0, t2, mpnw1)
    call mpsub (t2, t6, t7, mpnw1)
    call mpeq (t7, t6, mpnw1)
  endif

  call mproun (t6, mpnw)
  call mpeq (t6, terfc, mpnw)
endif

return
end subroutine mperfcr

subroutine mpexpint (x, y, mpnw)

!   This evaluates the exponential integral function Ei(x):
!   Ei(x) = - incgamma (0, -x)

implicit none
integer, intent(in):: mpnw
integer (mpiknd), intent(in):: x(0:)
integer (mpiknd), intent(out):: y(0:)
integer (mpiknd) t1(0:mpnw+5), t2(0:mpnw+5), t3(0:mpnw+5)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (x) < mpnw + 4 .or. mpspacer (y) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPEXPINT: Uninitialized or inadequately sized arrays')
  call mpabrt ( 552)
endif

if (mpsgn (x) == 0) then
  write (mpldb, 2)
2 format ('*** MPEXPINT: argument is zero')
  call mpabrt ( 553)
endif

call mpinitwds (t1, mpnw)
call mpinitwds (t2, mpnw)
call mpinitwds (t3, mpnw)
call mpdmc (0.d0, 0, t1, mpnw)
call mpneg (x, t2, mpnw)
call mpincgammar (t1, t2, t3, mpnw)
call mpneg (t3, y, mpnw)
return
end

subroutine mpgammar (t, z, mpnw)

!   This evaluates the gamma function, using an algorithm of R. W. Potter.
!   The argument t must not exceed 10^8 in size (this limit is set below),
!   must not be zero, and if negative must not be integer.

!   In the parameter statement below:
!     itrmx = limit of number of iterations in series; default = 100000.
!     con1 = 1/2 * log (10) to DP accuracy.
!     dmax = maximum size of input argument.

implicit none
integer, intent(in):: mpnw
integer i, i1, ic1, itrmx, j, mpnw1, nt, n1, n2, n3
real (mprknd) alpha, al2, dmax, d1, d2, d3, pi
parameter (al2 = 0.69314718055994530942d0, dmax = 1d8, itrmx = 100000, &
  pi = 3.1415926535897932385d0)
integer (mpiknd), intent(in):: t(0:)
integer (mpiknd), intent(out):: z(0:)
integer (mpiknd) f1(0:mpnw+6), sum1(0:mpnw+6), sum2(0:mpnw+6), tn(0:mpnw+6), &
  t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6), &
  t6(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (t) < mpnw + 4 .or. mpspacer (z) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPGAMMAR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 554)
endif

mpnw1 = mpnw + 1

!   Check if Pi has been precomputed.

call mpmdc (mppicon, d1, n1, mpnw1)
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw1) then
  write (mpldb, 4) mpnw1
4 format ('*** MPGAMMAR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 555)
endif

call mpinitwds (f1, mpnw1)
call mpinitwds (sum1, mpnw1)
call mpinitwds (sum2, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (t6, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

call mpmdc (t, d1, n1, mpnw)
d1 = d1 * 2.d0**n1
call mpnint (t, t1, mpnw)
call mpcpr (t, t1, ic1, mpnw)
i1 = mpsgn (t)
if (i1 == 0 .or. d1 > dmax .or. (i1 < 0 .and. ic1 == 0)) then
  write (mpldb, 2) dmax
2 format ('*** MPGAMMAR: input argument must have absolute value <=',f10.0,','/ &
  'must not be zero, and if negative must not be an integer.')
  call mpabrt ( 556)
endif

call mpdmc (1.d0, 0, f1, mpnw1)

!   Find the integer and fractional parts of t.

call mpinfr (t, t2, t3, mpnw1)

if (mpsgn (t3) == 0) then

!   If t is a positive integer, then apply the usual factorial recursion.

  call mpmdc (t2, d2, n2, mpnw1)
  nt = d2 * 2.d0 ** n2
  call mpeq (f1, t1, mpnw1)

  do i = 2, nt - 1
    call mpmuld (t1, dble (i), t2, mpnw1)
    call mpeq (t2, t1, mpnw1)
  enddo

  call mproun (t1, mpnw)
  call mpeq (t1, z, mpnw)
  goto 120
elseif (mpsgn (t) > 0) then

!   Apply the identity Gamma[t+1] = t * Gamma[t] to reduce the input argument
!   to the unit interval.

  call mpmdc (t2, d2, n2, mpnw1)
  nt = d2 * 2.d0 ** n2
  call mpeq (f1, t1, mpnw1)
  call mpeq (t3, tn, mpnw1)

  do i = 1, nt
    call mpdmc (dble (i), 0, t4, mpnw1)
    call mpsub (t, t4, t5, mpnw1)
    call mpmul (t1, t5, t6, mpnw1)
    call mpeq (t6, t1, mpnw1)
  enddo
else

!   Apply the gamma identity to reduce a negative argument to the unit interval.

  call mpsub (f1, t, t4, mpnw1)
  call mpinfr (t4, t3, t5, mpnw1)
  call mpmdc (t3, d3, n3, mpnw1)
  nt = d3 * 2.d0 ** n3

  call mpeq (f1, t1, mpnw1)
  call mpsub (f1, t5, t2, mpnw1)
  call mpeq (t2, tn, mpnw1)

  do i = 0, nt - 1
    call mpdmc (dble (i), 0, t4, mpnw1)
    call mpadd (t, t4, t5, mpnw1)
    call mpdiv (t1, t5, t6, mpnw1)
    call mpeq (t6, t1, mpnw1)
  enddo
endif

!   Calculate alpha = bits of precision * log(2) / 2, then take the next integer
!   value mod 4, so that d2 = 0.25 * alpha^2 can be calculated exactly in DP.

alpha = 4.d0 * aint ((0.5d0 * mpnbt * al2 * (mpnw1 + 1)) / 4.d0 + 1.d0)
d2 = 0.25d0 * alpha**2

call mpeq (tn, t2, mpnw1)
call mpdiv (f1, t2, t3, mpnw1)
call mpeq (t3, sum1, mpnw1)

!   Evaluate the series with t.

do j = 1, itrmx
  call mpdmc (dble (j), 0, t6, mpnw1)
  call mpadd (t2, t6, t4, mpnw1)
  call mpmuld (t4, dble (j), t5, mpnw1)
  call mpdiv (t3, t5, t6, mpnw1)
  call mpmuld (t6, d2, t3, mpnw1)
  call mpadd (sum1, t3, t4, mpnw1)
  call mpeq (t4, sum1, mpnw1)

  call mpabs (t3, tc1, 4)
  call mpmul (eps, sum1, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic1, 4)
  if (ic1 <= 0) goto 100
enddo

write (mpldb, 3) 1, itrmx
3 format ('*** MPGAMMAR: iteration limit exceeded',2i10)
call mpabrt ( 557)

100 continue

call mpneg (tn, t2, mpnw1)
call mpdiv (f1, t2, t3, mpnw1)
call mpeq (t3, sum2, mpnw1)

!   Evaluate the same series with -t.

do j = 1, itrmx
  call mpdmc (dble (j), 0, t6, mpnw1)
  call mpadd (t2, t6, t4, mpnw1)
  call mpmuld (t4, dble (j), t5, mpnw1)
  call mpdiv (t3, t5, t6, mpnw1)
  call mpmuld (t6, d2, t3, mpnw1)
  call mpadd (sum2, t3, t4, mpnw1)
  call mpeq (t4, sum2, mpnw1)

  call mpabs (t3, tc1, 4)
  call mpmul (eps, sum2, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic1, 4)
  if (ic1 <= 0) goto 110
enddo

write (mpldb, 3) 2, itrmx
call mpabrt ( 558)

110 continue

!   Compute sqrt (mppic * sum1 / (tn * sin (mppic * tn) * sum2))
!   and (alpha/2)^tn terms.

call mpeq (mppicon, t2, mpnw1)
call mpmul (t2, tn, t3, mpnw1)
call mpcssnr (t3, t4, t5, mpnw1)

call mpmul (t5, sum2, t6, mpnw1)
call mpmul (tn, t6, t5, mpnw1)
call mpmul (t2, sum1, t3, mpnw1)
call mpdiv (t3, t5, t6, mpnw1)
call mpneg (t6, t4, mpnw1)
call mpeq (t4, t6, mpnw1)
call mpsqrt (t6, t2, mpnw1)
call mpdmc (0.5d0 * alpha, 0, t3, mpnw1)
call mplog (t3, t4, mpnw1)
call mpmul (tn, t4, t5, mpnw1)
call mpexp (t5, t6, mpnw1)
call mpmul (t2, t6, t3, mpnw1)
call mpmul (t1, t3, t4, mpnw1)

!   Round to mpnw words precision.

call mproun (t4, mpnw)
call mpeq (t4, z, mpnw)

120 continue

return
end subroutine mpgammar

subroutine mphurwitzzetan (is, aa, zz, mpnw)

!   This returns the Hurwitz zeta function of IS and AA, using an algorithm from:
!   David H. Bailey and Jonathan M. Borwein, "Crandall's computation of the
!   incomplete gamma function and the Hurwitz zeta function with applications to
!   Dirichlet L-series," Applied Mathematics and Computation, vol. 268C (Oct 2015),
!   pg. 462-477, preprint at:
!   https://www.davidhbailey.com/dhbpapers/lerch.pdf
!   This is limited to IS >= 2 and 0 < AA < 1.

implicit none
integer i1, ic1, ic2, ic3, itrmax, k, n1, mpnw, mpnw1
real (mprknd) d1, dk, pi
parameter (itrmax = 1000000, pi = 3.1415926535897932385d0)
integer, intent(in):: is
integer (mpiknd), intent(in):: aa(0:)
integer (mpiknd), intent(out):: zz(0:)
integer (mpiknd) gs1(0:mpnw+6), gs2(0:mpnw+6), ss(0:mpnw+6), &
  sum1(0:mpnw+6), sum2(0:mpnw+6), sum3(0:mpnw+6), ss1(0:mpnw+6), ss2(0:mpnw+6), &
  ss3(0:mpnw+6), ss4(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), s3(0:mpnw+6), &
  t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6), &
  t6(0:mpnw+6), t7(0:mpnw+6), t8(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (aa) < mpnw + 4 .or. mpspacer (zz) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPHURWITZZETAN: Uninitialized or inadequately sized arrays')
  call mpabrt ( 559)
endif

mpnw1 = mpnw + 1

!   Check if Pi has been precomputed.

call mpmdc (mppicon, d1, n1, mpnw1)
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw1) then
  write (mpldb, 2) mpnw1
2 format ('*** MPHURWITZZETAN: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 560)
endif

call mpinitwds (gs1, mpnw1)
call mpinitwds (gs2, mpnw1)
call mpinitwds (ss, mpnw1)
call mpinitwds (sum1, mpnw1)
call mpinitwds (sum2, mpnw1)
call mpinitwds (sum3, mpnw1)
call mpinitwds (ss1, mpnw1)
call mpinitwds (ss2, mpnw1)
call mpinitwds (ss3, mpnw1)
call mpinitwds (ss4, mpnw1)
call mpinitwds (s1, mpnw1)
call mpinitwds (s2, mpnw1)
call mpinitwds (s3, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (t6, mpnw1)
call mpinitwds (t7, mpnw1)
call mpinitwds (t8, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

if (is <= 0) then
  write (mpldb, 3)
3 format ('*** MPHURWITZZETAN: IS less than or equal to 0:')
  call mpabrt ( 561)
endif

call mpdmc (0.d0, 0, t1, mpnw1)
call mpdmc (1.d0, 0, t2, mpnw1)
call mpcpr (aa, t1, ic1, mpnw1)
call mpcpr (aa, t2, ic2, mpnw1)
if (ic1 <= 0 .or. ic2 >= 0) then
  write (mpldb, 4)
4 format ('*** MPHURWITZZETAN: AA must be in the range (0, 1)')
  call mpabrt ( 562)
endif

! ss = mpreal (dble (is), mpnw)
! ss1 = 0.5d0 * ss
! ss2 = 0.5d0 * (ss + 1.d0)
! ss3 = 0.5d0 * (1.d0 - ss)
! ss4 = 1.d0 - 0.5d0 * ss

call mpdmc (dble(is), 0, ss, mpnw1)
call mpmuld (ss, 0.5d0, ss1, mpnw1)
call mpdmc (1.d0, 0, t1, mpnw1)
call mpadd (t1, ss, t2, mpnw1)
call mpmuld (t2, 0.5d0, ss2, mpnw1)
call mpsub (t1, ss, t2, mpnw1)
call mpmuld (t2, 0.5d0, ss3, mpnw1)
call mpsub (t1, ss1, ss4, mpnw1)

! gs1 = gamma (ss1)
! gs2 = gamma (ss2)
! t1 = pi * aa ** 2

call mpgammar (ss1, gs1, mpnw1)
call mpgammar (ss2, gs2, mpnw1)
call mpmul (aa, aa, t2, mpnw1)
call mpmul (mppicon, t2, t1, mpnw1)

! sum1 = (incgamma (ss1, t1) / gs1 + incgamma (ss2, t1) / gs2) / abs (aa)**is
! sum2 = mpreal (0.d0, mpnw)
! sum3 = mpreal (0.d0, mpnw)

call mpincgammar (ss1, t1, t2, mpnw1)
call mpdiv (t2, gs1, t3, mpnw1)
call mpincgammar (ss2, t1, t2, mpnw1)
call mpdiv (t2, gs2, t4, mpnw1)
call mpadd (t3, t4, t2, mpnw1)
call mpabs (aa, t3, mpnw1)
call mpnpwr (t3, is, t4, mpnw1)
call mpdiv (t2, t4, sum1, mpnw1)
call mpdmc (0.d0, 0, sum2, mpnw1)
call mpdmc (0.d0, 0, sum3, mpnw1)

do k = 1, itrmax
  dk = dble (k)

!  t1 = pi * (dk + aa)**2
!  t2 = pi * (-dk + aa)**2
!  t3 = dk**2 * pi
!  t4 = 2.d0 * pi * dk * aa

  call mpdmc (dk, 0, t5, mpnw1)
  call mpadd (t5, aa, t6, mpnw1)
  call mpmul (t6, t6, t5, mpnw1)
  call mpmul (mppicon, t5, t1, mpnw1)
  call mpdmc (-dk, 0, t5, mpnw1)
  call mpadd (t5, aa, t6, mpnw1)
  call mpmul (t6, t6, t7, mpnw1)
  call mpmul (mppicon, t7, t2, mpnw1)
  call mpmul (t5, t5, t6, mpnw1)
  call mpmul (mppicon, t6, t3, mpnw1)
  call mpmuld (mppicon, 2.d0 * dk, t5, mpnw1)
  call mpmul (t5, aa, t4, mpnw1)

!  s1 = (incgamma (ss1, t1) / gs1 + incgamma (ss2, t1) / gs2) / abs (dk + aa)**is

  call mpincgammar (ss1, t1, t5, mpnw1)
  call mpdiv (t5, gs1, t6, mpnw1)
  call mpincgammar (ss2, t1, t5, mpnw1)
  call mpdiv (t5, gs2, t7, mpnw1)
  call mpadd (t6, t7, t5, mpnw1)
  call mpdmc (dk, 0, t6, mpnw1)
  call mpadd (t6, aa, t7, mpnw1)
  call mpabs (t7, t6, mpnw1)
  call mpnpwr (t6, is, t7, mpnw1)
  call mpdiv (t5, t7, s1, mpnw1)

!  s2 = (incgamma (ss1, t2) / gs1 - incgamma (ss2, t2) / gs2) / abs (-dk + aa)**is

  call mpincgammar (ss1, t2, t5, mpnw1)
  call mpdiv (t5, gs1, t6, mpnw1)
  call mpincgammar (ss2, t2, t5, mpnw1)
  call mpdiv (t5, gs2, t7, mpnw1)
  call mpsub (t6, t7, t5, mpnw1)
  call mpdmc (-dk, 0, t6, mpnw1)
  call mpadd (t6, aa, t7, mpnw1)
  call mpabs (t7, t6, mpnw1)
  call mpnpwr (t6, is, t7, mpnw1)
  call mpdiv (t5, t7, s2, mpnw1)

!  sum1 = sum1 + s1
!  sum2 = sum2 + s2

  call mpadd (sum1, s1, t5, mpnw1)
  call mpeq (t5, sum1, mpnw1)
  call mpadd (sum2, s2, t5, mpnw1)
  call mpeq (t5, sum2, mpnw1)

!  s3 = (incgamma (ss3, t3) * cos (t4)/ gs1 + incgamma (ss4, t3) * sin (t4) / gs2) &
!    / mpreal (dk, mpnw)**(1-is)

  call mpincgammar (ss3, t3, t5, mpnw1)
  call mpcssnr (t4, t6, t7, mpnw1)
  call mpmul (t5, t6, t8, mpnw1)
  call mpdiv (t8, gs1, t5, mpnw1)
  call mpincgammar (ss4, t3, t6, mpnw1)
  call mpmul (t6, t7, t8, mpnw1)
  call mpdiv (t8, gs2, t6, mpnw1)
  call mpadd (t5, t6, t7, mpnw1)
  call mpdmc (dk, 0, t5, mpnw1)
  i1 = 1 - is
  call mpnpwr (t5, i1, t6, mpnw1)
  call mpdiv (t7, t6, s3, mpnw1)

!  sum3 = sum3 + s3

  call mpadd (sum3, s3, t5, mpnw1)
  call mpeq (t5, sum3, mpnw1)

!  if (abs (s1) < eps * abs (sum1) .and. abs (s2) < eps * abs (sum2) .and. &
!    abs (s3) < eps * abs (sum3)) goto 100

  call mpabs (s1, tc1, 4)
  call mpmul (eps, sum1, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic1, 4)
  call mpabs (s2, tc1, 4)
  call mpmul (eps, sum2, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic2, 4)
  call mpabs (s3, tc1, 4)
  call mpmul (eps, sum3, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic3, 4)
  if (ic1 <= 0 .and. ic2 <= 0 .and. ic3 <= 0) goto 100
enddo

write (mpldb, 5)
5 format ('*** MPHURWITZZETAN: Loop end error')
call mpabrt ( 563)

100 continue

if (mod (is, 2) == 0) then
!  t1 = pi ** (is / 2) / ((ss - 1.d0) * gamma (ss1))

  i1 = is / 2
  call mpnpwr (mppicon, i1, t2, mpnw1)
  call mpdmc (1.d0, 0, t3, mpnw1)
  call mpsub (ss, t3, t4, mpnw1)
  call mpgammar (ss1, t5, mpnw1)
  call mpmul (t4, t5, t6, mpnw1)
  call mpdiv (t2, t6, t1, mpnw1)
else
!  t1 = sqrt (pi) * pi ** ((is - 1) / 2) / ((ss - 1.d0) * gamma (ss1))

  i1 = (is - 1) / 2
  call mpnpwr (mppicon, i1, t2, mpnw1)
  call mpsqrt (mppicon, t3, mpnw1)
  call mpmul (t2, t3, t4, mpnw1)
  call mpdmc (1.d0, 0, t2, mpnw1)
  call mpsub (ss, t2, t3, mpnw1)
  call mpgammar (ss1, t5, mpnw1)
  call mpmul (t3, t5, t6, mpnw1)
  call mpdiv (t4, t6, t1, mpnw1)
endif

!t2 = pi ** is / sqrt (pi)

call mpnpwr (mppicon, is, t3, mpnw1)
call mpsqrt (mppicon, t4, mpnw1)
call mpdiv (t3, t4, t2, mpnw1)

! hurwitzzetan = t1 + 0.5d0 * sum1 + 0.5d0 * sum2 + t2 * sum3

call mpmuld (sum1, 0.5d0, t3, mpnw1)
call mpmuld (sum2, 0.5d0, t4, mpnw1)
call mpmul (sum3, t2, t5, mpnw1)
call mpadd (t1, t3, t6, mpnw1)
call mpadd (t6, t4, t7, mpnw1)
call mpadd (t7, t5, t1, mpnw1)
call mproun (t1, mpnw)
call mpeq (t1, zz, mpnw)

return
end subroutine mphurwitzzetan

subroutine mphurwitzzetanbe (nb1, nb2, berne, iss, aa, zz, mpnw)

!  This evaluates the Hurwitz zeta function, using the combination of
!  the definition formula (for large iss), and an Euler-Maclaurin scheme
!  (see formula 25.2.9 of the DLMF). The array berne contains precomputed
!  even Bernoulli numbers (see MPBERNER above). Its dimensions must be as
!  shown below. NB2 must be greater than 1.4 x precision in decimal digits.

implicit none
integer i1, i2, ic1, iqq, itrmax, k, n1, mpnw, mpnw1
real (mprknd) d1, dber, dcon, dp
parameter (itrmax = 1000000, dber = 1.4d0, dcon = 0.6d0)
integer, intent(in):: nb1, nb2, iss
integer (mpiknd), intent(in):: berne(0:nb1+5,nb2), aa(0:)
integer (mpiknd), intent(out):: zz(0:)
integer (mpiknd) aq(0:mpnw+6), aq2(0:mpnw+6), s1(0:mpnw+6), s2(0:mpnw+6), &
  s3(0:mpnw+6), s4(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), &
  t4(0:mpnw+6), t5(0:mpnw+6), t6(0:mpnw+6), eps(0:9), f1(0:9), tc1(0:9), &
  tc2(0:9), tc3(0:9)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (aa) < mpnw + 4 .or. mpspacer (zz) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPHURWITZZETANBE: Uninitialized or inadequately sized arrays')
  call mpabrt ( 564)
endif

mpnw1 = mpnw + 1

!   Check if berne array has been initialized.

i1 = 1000000000
i2 = 1000000000

do k = 1, nb2
  i1 = min (i1, mpspacer (berne(0:nb1+5,k)))
  i2 = min (i2, mpwprecr (berne(0:nb1+5,k)))
enddo

call mpmdc (berne(0:nb1+5,1), d1, n1, mpnw)
d1 = d1 * 2.d0 ** n1

if (i1 < mpnw + 6 .or. i2 < mpnw .or. abs (d1 - 1.d0 / 6.d0) > mprdfz .or. &
  nb2 < int (dber * mpdpw * (mpnw - 2))) then
  write (mpldb, 3) int (dber * mpdpw * (mpnw - 2))
3 format ('*** MPHURWITZZETANBE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries by calling MPBERNE or MPBERNER.')
  call mpabrt ( 565)
endif

call mpinitwds (aq, mpnw1)
call mpinitwds (aq2, mpnw1)
call mpinitwds (s1, mpnw1)
call mpinitwds (s2, mpnw1)
call mpinitwds (s3, mpnw1)
call mpinitwds (s4, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (t6, mpnw1)
call mpinitwds (eps, 4)
call mpinitwds (f1, 4)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)
call mpdmc (1.d0, 0, f1, 4)
call mpdmc (0.d0, 0, s1, mpnw1)
call mpdmc (0.d0, 0, s2, mpnw1)
call mpdmc (0.d0, 0, s3, mpnw1)
call mpdmc (0.d0, 0, s4, mpnw1)

if (iss <= 0) then
  write (mpldb, 4)
4 format ('*** MPHURWITZZETANBE: ISS <= 0')
  call mpabrt ( 566)
endif

if (mpsgn (aa) < 0) then
  write (mpldb, 5)
5 format ('*** MPHURWITZZETANBE: AA < 0')
  call mpabrt ( 567)
endif

dp = anint (mpnw1 * mpdpw)

!   If iss > a certain value, then use definition formula.

if (iss > 2.303d0 * dp / log (2.515d0 * dp)) then
  do k = 0, itrmax
!    t1 = 1.d0 / (aa + dble (k))**iss
!    s1 = s1 + t1

    call mpdmc (dble (k), 0, t1, mpnw1)
    call mpadd (aa, t1, t2, mpnw1)
    call mpnpwr (t2, iss, t3, mpnw1)
    call mpdiv (f1, t3, t1, mpnw1)
    call mpadd (s1, t1, t2,  mpnw1)
    call mpeq (t2, s1, mpnw1)

!    if (abs (t1 / s1) < eps) goto 110

    call mpabs (t1, tc1, 4)
    call mpmul (eps, s1, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 110
  enddo

  write (6, 6)
6 format ('*** MPHURWITZZETANBE: Loop end error 1')
  call mpabrt ( 568)
endif

call mpmdc (aa, d1, n1, mpnw1)
d1 = d1 * 2.d0**n1
iqq = max (dcon * mpnw1 * mpdpw - d1, 0.d0)

do k = 0, iqq - 1
!  s1 = s1 + 1.d0 / (aa + dble (k))**iss

  call mpdmc (dble (k), 0, t1, mpnw1)
  call mpadd (aa, t1, t2, mpnw1)
  call mpnpwr (t2, iss, t3, mpnw1)
  call mpdiv (f1, t3, t1, mpnw1)
  call mpadd (s1, t1, t2, mpnw1)
  call mpeq (t2, s1, mpnw1)
enddo

! aq = aa + dble (iqq)

call mpdmc (dble (iqq), 0, t1, mpnw1)
call mpadd (aa, t1, aq, mpnw1)

! s2 = 1.d0 / (dble (iss - 1) * aq**(iss -  1))

call mpdmc (dble (iss - 1), 0, t1, mpnw1)
call mpnpwr (aq, iss - 1, t2, mpnw1)
call mpmul (t1, t2, t3,  mpnw1)
call mpdiv (f1, t3, s2, mpnw1)

! s3 = 1.d0 / (2.d0 * aq**iss)

call mpnpwr (aq, iss, t1, mpnw1)
call mpmuld (t1, 2.d0, t2, mpnw1)
call mpdiv (f1, t2, s3, mpnw1)

! t1 = mpreal (dble (iss), nwds)
! t2 = mpreal (1.d0, nwds)
! t3 = aq**(iss - 1)
! aq2 = aq**2

call mpdmc (dble (iss), 0, t1, mpnw1)
call mpdmc (1.d0, 0, t2, mpnw1)
call mpnpwr (aq, iss - 1, t3, mpnw1)
call mpmul (aq, aq, aq2, mpnw1)

do k = 1, nb2
!  if (k > 1) t1 = t1 * dble (iss + 2*k - 3) * dble (iss + 2*k - 2)

  if (k > 1) then
    call mpmuld (t1, dble (iss + 2*k - 3), t5, mpnw1)
    call mpmuld (t5, dble (iss + 2*k - 2), t1, mpnw1)
  endif

!  t2 = t2 * dble (2 * k - 1) * dble (2 * k)

  call mpmuld (t2, dble (2 * k - 1), t5, mpnw1)
  call mpmuld (t5, dble (2 * k), t2, mpnw1)

!  t3 = t3 * aq2
!  t4 = rb(k) * t1 / (t2 * t3)
!  s4 = s4 + t4

  call mpmul (t3, aq2, t5, mpnw1)
  call mpeq (t5, t3, mpnw1)
  call mpmul (berne(0:nb1+5,k), t1, t5, mpnw1)
  call mpmul (t2, t3, t6, mpnw1)
  call mpdiv (t5, t6, t4, mpnw1)
  call mpadd (s4, t4, t5, mpnw1)
  call mpeq (t5, s4, mpnw1)

!  if (abs (t4) < eps) goto 110

  call mpabs (t4, tc1, 4)
  call mpmul (eps, s4, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic1, 4)
  if (ic1 <= 0) goto 110
enddo

write (6, 7)
7 format ('*** MPHURWITZZETANBE: End loop error 2; call MPBERNE with larger NB.')
call mpabrt ( 569)

110 continue

! hurwitz_be = s1 + s2 + s3 + s4

call mpadd (s1, s2, t5, mpnw1)
call mpadd (t5, s3, t6, mpnw1)
call mpadd (t6, s4, s1, mpnw1)
call mproun (s1, mpnw)
call mpeq (s1, zz, mpnw)

return
end subroutine mphurwitzzetanbe

subroutine mphypergeompfq (np, nq, nw, aa, bb, xx, yy, mpnw)

!  This returns the HypergeometricPFQ function, namely the sum of the infinite series

!  Sum_0^infinity poch(aa(1),n)*poch(aa(2),n)*...*poch(aa(np),n) /
!      poch(bb(1),n)*poch(bb(2),n)*...*poch(bb(nq),n) * xx^n / n!

!  This subroutine evaluates the HypergeometricPFQ function directly according to
!  this definitional formula. The arrays aa and bb must be dimensioned as shown below.
!  NP and NQ are limited to [1,10].

implicit none
integer, intent(in):: np, nq, nw, mpnw
integer i1, i2, ic1, itrmax, j, k, mpnw1, npq
parameter (itrmax = 1000000, npq = 10)
integer (mpiknd), intent(in):: aa(0:nw+5,np), bb(0:nw+5,nq), xx(0:)
integer (mpiknd), intent(out):: yy(0:)
integer (mpiknd) sum(0:mpnw+6), td(0:mpnw+6), tn(0:mpnw+6), t1(0:mpnw+6), &
  t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
integer (mpiknd) ix8

!  End of declaration.

i1 = 1000000
i2 = 1000000

do k = 1, np
  i1 = min (i1, mpspacer(aa(0:nw+5,k)))
enddo

do k = 1, nq
  i2 = min (i2, mpspacer(bb(0:nw+5,k)))
enddo

if (mpnw < 4 .or. min (i1, i2) < mpnw + 4 .or. mpspacer (xx) < mpnw + 4 .or. &
  mpspacer (yy) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPHYPERGEOMPFQ: Uninitialized or inadequately sized arrays')
  call mpabrt ( 570)
endif

mpnw1 = mpnw + 1
call mpinitwds (sum, mpnw1)
call mpinitwds (td, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

if (np < 1 .or. np > npq .or. nq < 1 .or. nq > npq) then
  write (mpldb, 2) npq
2 format ('*** MPHYPERGEOMPFQ: NP and NQ must be between 1 and',i4)
  call mpabrt ( 571)
endif

call mpdmc (1.d0, 0, sum, mpnw1)
call mpdmc (1.d0, 0, td, mpnw1)
call mpdmc (1.d0, 0, tn, mpnw1)

do k = 1, itrmax
  call mpdmc (dble (k - 1), 0, t1, mpnw1)

  do j = 1, np
    call mpadd (t1, aa(0:nw+5,j), t2, mpnw1)
    call mpmul (tn, t2, t3, mpnw1)
    call mpeq (t3, tn, mpnw1)
  enddo

  do j = 1, nq
    call mpadd (t1, bb(0:nw+5,j), t2, mpnw1)
    call mpmul (td, t2, t3, mpnw1)
    call mpeq (t3, td, mpnw1)
  enddo

  call mpmul (tn, xx, t2, mpnw1)
  call mpeq (t2, tn, mpnw1)
  call mpmuld (td, dble (k), t3,  mpnw1)
  call mpeq (t3, td, mpnw1)
  call mpdiv (tn, td, t1, mpnw1)
  call mpadd (sum, t1, t2, mpnw1)
  call mpeq (t2, sum, mpnw1)

  call mpabs (t1, tc1, 4)
  call mpmul (eps, sum, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic1, 4)
  if (ic1 <= 0) goto 100
enddo

    write (mpldb, 3) itrmax
3   format ('*** MPHYPERGEOMPFQ: Loop end error',i10)
    call mpabrt ( 572)

100  continue

call mproun (sum, mpnw)
call mpeq (sum, yy, mpnw)
return
end subroutine mphypergeompfq

subroutine mpincgammar (s, z, g, mpnw)

!  This returns the incomplete gamma function, using a combination of formula
!  8.7.3 of the DLMF (for modest-sized z), formula 8.11.2 (for large z),
!  a formula from the Wikipedia page for the case S = 0, and another formula
!  from the Wikipedia page for the case S = negative integer. The formula
!  for the case S = 0 requires increased working precision, up to 2.5X normal,
!  depending on the size of Z.

implicit none
integer, intent(in):: mpnw
integer ic1, itrmax, k, mpnw1, mpnw2, nn, n1, n2
real (mprknd) d1, d2, bits, dmax, egam
parameter (dmax = 0.833d0, itrmax = 1000000, egam = 0.5772156649015328606d0)
integer (mpiknd), intent(in):: s(0:), z(0:)
integer (mpiknd), intent(out):: g(0:)
integer (mpiknd) t0(0:5*mpnw/2+6), t1(0:5*mpnw/2+6), t2(0:5*mpnw/2+6), &
  t3(0:5*mpnw/2+6), t4(0:5*mpnw/2+6), t5(0:5*mpnw/2+6), f1(0:5*mpnw/2+6), &
  tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (s) < mpnw + 4 .or. mpspacer (z) < mpnw + 4 &
  .or. mpspacer (g) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPINCGAMMAR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 573)
endif

n1 = mpsgn (s)
n2 = mpsgn (z)
if (n2 == 0 .or. n1 /= 0 .and. n2 < 0) then
  write (mpldb, 2)
2 format ('*** MPINCGAMMAR: The second argument must not be zero,'/ &
    'and must not be negative unless the first is zero.')
  call mpabrt ( 574)
endif

!   Check if EGAMMA has been precomputed.

mpnw1 = mpnw + 1
call mpmdc (mpegammacon, d1, n1, mpnw1)
if (n1 /= -1 .or. abs (d1 * 2.d0**n1 - egam) > mprdfz &
  .or. mpwprecr (mpegammacon) < mpnw1) then
  write (mpldb, 3) mpnw
3 format ('*** MPINCGAMMAR: EGAMMA must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 575)
endif

call mpinitwds (t0, 5*mpnw/2+1)
call mpinitwds (t1, 5*mpnw/2+1)
call mpinitwds (t2, 5*mpnw/2+1)
call mpinitwds (t3, 5*mpnw/2+1)
call mpinitwds (t4, 5*mpnw/2+1)
call mpinitwds (t5, 5*mpnw/2+1)
call mpinitwds (f1, 5*mpnw/2+1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

call mpdmc (1.d0, 0, f1, mpnw1)
call mpmdc (z, d1, n1, mpnw1)
d1 = d1 * 2.d0 ** n1
bits = mpnw1 * mpnbt

if (abs (d1) < dmax * bits) then

!   This is for modest-sized z.

  call mpinfr (s, t1, t2, mpnw1)
  call mpcpr (s, t1, ic1, mpnw1)
  call mpmdc (s, d2, n2, mpnw1)
  nn = d2 * 2.d0**n2

  if (ic1 == 0 .and. nn == 1) then

!   S = 1; result is exp (-z).

    call mpneg (z, t0,  mpnw1)
    call mpexp (t0, t1, mpnw1)
    goto 200
  elseif (ic1 == 0 .and. nn <= 0) then

!    S is zero or a negative integer -- use a different algorithm. In
!    either event, first compute incgamma for S = 0. For large Z, the
!    working precision must be increased, up to 2.5X times normal.

!    mpnw2 = min (mpnw1 + 1.5d0 * d1 / (dmax * bits) * mpnw, 5*mpnw/2+1.d0)
    mpnw2 = mpnw1
    call mpeq (z, t0, mpnw2)
    call mpeq (z, t1, mpnw2)
    call mpdmc (1.d0, 0, t2, mpnw2)

    do k = 2, itrmax
      if (mod (k, 2) == 1) then
        d1 = dble (k)
        call mpdivd (f1, d1, t3, mpnw2)
        call mpadd (t2, t3, t4, mpnw2)
        call mpeq (t4, t2, mpnw2)
      endif
      call mpmul (z, t1, t3, mpnw2)
      d1 = 2.d0 * dble (k)
      call mpdivd (t3, d1, t1, mpnw2)
      call mpmul (t1, t2, t3, mpnw2)
      call mpadd (t0, t3, t4, mpnw2)
      call mpeq (t4, t0, mpnw2)

      call mpabs (t3, tc1, 4)
      call mpmul (eps, t0, tc3, 4)
      call mpabs (tc3, tc2, 4)
      call mpcpr (tc1, tc2, ic1, 4)
      if (ic1 <= 0) goto 100
    enddo

    write (mpldb, 4)
4   format ('*** MPINCGAMMAR: Loop end error 1')
    call mpabrt ( 576)

100  continue

    call mpneg (mpegammacon, t1, mpnw1)
    call mpabs (z, t3, mpnw1)
    call mplog (t3, t2, mpnw1)
    call mpsub (t1, t2, t3,  mpnw1)
    call mpmuld (z, -0.5d0, t4, mpnw1)
    call mpexp (t4, t5, mpnw1)
    call mpmul (t5, t0, t4, mpnw1)
    call mpadd (t3, t4, t1, mpnw1)
    if (nn == 0) goto 200

!   S is negative integer (not zero).

    nn = abs (nn)
    call mpdmc (1.d0, 0, t0, mpnw1)
    call mpeq (t0, t2, mpnw1)

    do k = 1, nn - 1
      call mpmuld (t0, dble (k), t2, mpnw1)
      call mpeq (t2, t0, mpnw1)
    enddo

    call mpmuld (t0, dble (nn), t5, mpnw1)

    do k = 1, nn - 1
      call mpmul (t2, z, t3, mpnw1)
      call mpdivd (t3, dble (nn - k), t4, mpnw1)
      call mpneg (t4, t2, mpnw1)
      call mpadd (t0, t2, t3, mpnw1)
      call mpeq (t3, t0, mpnw1)
    enddo

    call mpexp (z, t2, mpnw1)
    call mpdiv (t0, t2, t3, mpnw1)
    call mpnpwr (z, nn, t4, mpnw1)
    call mpdiv (t3, t4, t2, mpnw1)

    if (mod (nn, 2) == 0) then
      call mpadd (t2, t1, t3, mpnw1)
    else
      call mpsub (t2, t1, t3, mpnw1)
    endif
    call mpdiv (t3, t5, t1, mpnw1)
    goto 200
  endif

  call mpgammar (s, t1, mpnw1)
  call mpmul (s, t1, t3, mpnw1)
  call mpdiv (f1, t3, t2, mpnw1)
  call mpeq (t2, t0, mpnw1)

  do k = 1, itrmax
    call mpmul (t2, z, t5, mpnw1)
    call mpdmc (dble (k), 0, t3, mpnw1)
    call mpadd (s, t3, t4, mpnw1)
    call mpdiv (t5, t4, t2, mpnw1)
    call mpadd (t0, t2, t3, mpnw1)
    call mpeq (t3, t0, mpnw1)

    call mpabs (t2, tc1, 4)
    call mpmul (eps, t0, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 110
  enddo

  write (mpldb, 3) itrmax
  call mpabrt ( 577)

110 continue

  call mppower (z, s, t2, mpnw1)
  call mpexp (z, t3, mpnw1)
  call mpdiv (t2, t3, t4, mpnw1)
  call mpmul (t4, t0, t5, mpnw1)
  call mpsub (f1, t5, t2, mpnw1)
  call mpmul (t1, t2, t3, mpnw1)
  call mpeq (t3, t1, mpnw1)
  goto 200
else

!   This is for large z. Note that if S is a positive integer, this loop
!   is finite.

  call mpdmc (1.d0, 0, t0, mpnw1)
  call mpdmc (1.d0, 0, t1, mpnw1)

  do k = 1, itrmax
    call mpdmc (dble (k), 0, t2, mpnw1)
    call mpsub (s, t2, t3, mpnw1)
    call mpmul (t1, t3, t4, mpnw1)
    call mpdiv (t4, z, t1, mpnw1)
    call mpadd (t0, t1, t2, mpnw1)
    call mpeq (t2, t0, mpnw1)

    call mpabs (t1, tc1, 4)
    call mpmul (eps, t0, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 120
  enddo

  write (mpldb, 5)
5 format ('*** MPINCGAMMAR: Loop end error 2')
  call mpabrt ( 578)

120 continue

   call mpsub (s, f1, t2, mpnw1)
   call mppower (z, t2, t3, mpnw1)
   call mpexp (z, t4, mpnw1)
   call mpdiv (t3, t4, t2, mpnw1)
   call mpmul (t2, t0, t1, mpnw1)
   goto 200
endif

200 continue

call mproun (t1, mpnw)
call mpeq (t1, g, mpnw)

return
end subroutine mpincgammar

subroutine mppolygamma (nn, x, y, mpnw)

!   This returns polygamma (nn, x) for nn >= 0 and 0 < x < 1, by calling
!   mphurwitzzetan.

implicit none
integer, intent(in):: nn, mpnw
integer ic1, ic2, k, mpnw1
integer (mpiknd), intent(in):: x(0:)
integer (mpiknd), intent(out):: y(0:)
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (x) < mpnw + 4 .or. mpspacer (y) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPPOLYGAMMA: Uninitialized or inadequately sized arrays')
  call mpabrt ( 579)
endif

mpnw1 = mpnw
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)

if (nn <= 0) then
  write (mpldb, 2)
2 format ('*** MPPOLYGAMMA: NN <= 0')
  call mpabrt ( 580)
endif

call mpdmc (0.d0, 0, t1, mpnw1)
call mpdmc (1.d0, 0, t2, mpnw1)
call mpcpr (x, t1, ic1, mpnw1)
call mpcpr (x, t2, ic2, mpnw1)
if (ic1 <= 0 .or. ic2 >= 0) then
  write (mpldb, 3)
3 format ('*** MPPOLYGAMMA: X must be in the range (0, 1)')
  call mpabrt ( 581)
endif

call mpdmc (1.d0, 0, t1, mpnw1)

do k = 1, nn
  call mpmuld (t1, dble(k), t2, mpnw1)
  call mpeq (t2, t1, mpnw1)
enddo

if (mod (nn + 1, 2) == 1) then
  call mpneg (t1, t2, mpnw1)
  call mpeq (t2, t1, mpnw1)
endif
call mphurwitzzetan (nn + 1, x, t2, mpnw1)
call mpmul (t1, t2, t3, mpnw1)
call mproun (t3, mpnw)
call mpeq (t3, y, mpnw)

return
end subroutine mppolygamma

subroutine mppolygammabe (nb1, nb2, berne, nn, x, y, mpnw)

!  This returns polygamma (nn, x) for nn >= 0, by calling mphurwitzzetanbe.
!  The array berne contains precomputed even Bernoulli numbers (see MPBERNER
!  above). Its dimensions must be as shown below. NB2 must be greater than
!  1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb1, nb2, nn, mpnw
integer i1, i2, k, mpnw1, n1
real (mprknd) dber, d1
parameter (dber = 1.4d0)
integer (mpiknd), intent(in):: berne(0:nb1+5,nb2), x(0:)
integer (mpiknd), intent(out):: y(0:)
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (x) < mpnw + 4 .or. mpspacer (y) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPPOLYGAMMABE: Uninitialized or inadequately sized arrays')
  call mpabrt ( 582)
endif

mpnw1 = mpnw

!   Check if berne array has been initialized.

i1 = 1000000000
i2 = 1000000000

do k = 1, nb2
  i1 = min (i1, mpspacer (berne(0:nb1+5,k)))
  i2 = min (i2, mpwprecr (berne(0:nb1+5,k)))
enddo

call mpmdc (berne(0:nb1+5,1), d1, n1, mpnw)
d1 = d1 * 2.d0 ** n1

if (i1 < mpnw + 6 .or. i2 < mpnw .or. abs (d1 - 1.d0 / 6.d0) > mprdfz .or. &
  nb2 < int (dber * mpdpw * (mpnw - 2))) then
  write (mpldb, 2) int (dber * mpdpw * (mpnw - 2))
2 format ('*** MPPOLYGAMMABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries by calling MPBERNE or MPBERER.')
  call mpabrt ( 583)
endif

call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)

if (nn <= 0) then
  write (mpldb, 3)
3 format ('*** MPPOLYGAMMABE: NN <= 0')
  call mpabrt ( 584)
endif

if (mpsgn (x) < 0) then
  write (mpldb, 4)
4 format ('*** MPPOLYGAMMABE: X < 0')
  call mpabrt ( 585)
endif

call mpdmc (1.d0, 0, t1, mpnw1)

do k = 1, nn
  call mpmuld (t1, dble(k), t2, mpnw1)
  call mpeq (t2, t1, mpnw1)
enddo

if (mod (nn + 1, 2) == 1) then
  call mpneg (t1, t2, mpnw1)
  call mpeq (t2, t1, mpnw1)
endif
call mphurwitzzetanbe (nb1, nb2, berne, nn + 1, x, t2, mpnw1)
call mpmul (t1, t2, t3, mpnw1)
call mproun (t3, mpnw)
call mpeq (t3, y, mpnw)

return
end subroutine mppolygammabe

subroutine mppolylogini (na, nn, arr, mpnw)

!   Initializes the MP array arr with data for mppolylogneg.
!   NN must be in the range (-nmax, -1).

implicit none
integer, intent(in):: na, nn, mpnw
integer i1, i2, k, n, nmax, nna, mpnw1
parameter (nmax = 1000)
integer (mpiknd), intent(out):: arr(0:na+5,1:abs(nn))
integer (mpiknd) aa(0:mpnw+6,2,abs(nn)), t1(0:mpnw+6), t2(0:mpnw+6)
integer (mpiknd) ix8

!  End of declaration.

if (nn == 0 .or. abs (nn) > nmax) then
  write (mpldb, 1) -nmax
1 format ('*** MPPOLYLOGINI: N = 0 or N > ',i6)
  call mpabrt ( 586)
endif

nna = abs (nn)
i1 = 1000000000

do k = 1, nna
  i1 = min (i1, mpspacer (arr(0:na+5,k)))
enddo

if (mpnw < 4 .or. i1 < mpnw + 6) then
  write (mpldb, 2)
2 format ('*** MPPOLYLOGINI: Uninitialized or inadequately sized arrays')
  call mpabrt ( 587)
endif

mpnw1 = mpnw + 1
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
i1 = 2
i2 = 1
call mpinitwds (aa(0:mpnw+6,1,1), mpnw1)
call mpinitwds (aa(0:mpnw+6,2,1), mpnw1)
call mpdmc (1.d0, 0, aa(0:mpnw+6,1,1), mpnw1)
call mpdmc (1.d0, 0, aa(0:mpnw+6,2,1), mpnw1)

do k = 2, nna
  call mpinitwds (aa(0:mpnw+6,1,k), mpnw1)
  call mpinitwds (aa(0:mpnw+6,2,k), mpnw1)
  call mpdmc (0.d0, 0, aa(0:mpnw+6,1,k), mpnw1)
  call mpdmc (0.d0, 0, aa(0:mpnw+6,2,k), mpnw1)
enddo

do n = 2, nna
  i1 = 3 - i1
  i2 = 3 - i1

  do k = 2, n
    call mpmuld (aa(0:mpnw+6,i1,k-1), dble (n + 1 - k), t1, mpnw1)
    call mpmuld (aa(0:mpnw+6,i1,k), dble (k), t2, mpnw1)
    call mpadd (t1, t2, aa(0:mpnw+6,i2,k), mpnw1)
  enddo
enddo

do k = 1, nna
  call mpeq (aa(0:mpnw+6,i2,k), arr(0:na+5,k), mpnw)
enddo

return
end subroutine mppolylogini

subroutine mppolylogneg (na, nn, arr, x, y, mpnw)

!   This returns polylog (nn, x) for the case nn < 0. Before calling this,
!   one must call mppolylognini to initialize the array arr for this NN.
!   The dimensions of arr must be as shown below.
!   NN must be in the range (-nmax, -1).
!   The parameter nmxa is the maximum number of additional words of
!   precision needed to overcome cancelation errors when x is negative,
!   for nmax = 1000.

implicit none
integer, intent(in):: na, nn, mpnw
integer i1, i2, k, mpnw1, n1, n2, nna, nmax, nmxa
parameter (nmax = 1000, nmxa = 8525 / mpnbt + 1)
real (mprknd) d1, d2
integer (mpiknd), intent(in):: arr(0:na+5,1:abs(nn)), x(0:)
integer (mpiknd), intent(out):: y(0:)
integer (mpiknd) t1(0:mpnw+6+nmxa), t2(0:mpnw+6+nmxa), t3(0:mpnw+6+nmxa), &
  t4(0:mpnw+6+nmxa)
integer (mpiknd) ix8

!  End of declaration.

if (nn < -nmax .or. nn >= 0) then
  write (mpldb, 1) -nmax
1 format ('*** MPPOLYLOGNEG: N is <',i6,' or n >= 0.'/ &
  'For n >= 0, call mppolylogpos or polylog_pos.')
  call mpabrt ( 588)
endif

nna = abs (nn)
i1 = 1000000000
i2 = 1000000000

do k = 1, nna
  i1 = min (i1, mpspacer (arr(0:na+5,k)))
  i2 = min (i2, mpwprecr (arr(0:na+5,k)))
enddo

call mpmdc (arr(0:na+5,1), d1, n1, mpnw)
d1 = d1 * 2.d0 ** n1
call mpmdc (arr(0:na+5,nna), d2, n2, mpnw)
d2 = d2 * 2.d0 ** n2

if (mpnw < 4 .or. i1 < mpnw + 6 .or. i2 < mpnw .or. d1 /= 1.d0 .or. d2 /= 1.d0 &
  .or. mpspacer (x) < mpnw + 4 .or. mpspacer (y) < mpnw + 6) then
  write (mpldb, 2)
2 format ('*** MPPOLYLOGNEG: Uninitialized or inadequately sized arrays'/ &
  'Call mppolylogini or polylog_ini to initialize array. See documentation.')
  call mpabrt ( 589)
endif

mpnw1 = mpnw + 1
call mpinitwds (t1, mpnw + nmxa + 1)
call mpinitwds (t2, mpnw + nmxa + 1)
call mpinitwds (t3, mpnw + nmxa + 1)
call mpinitwds (t4, mpnw + nmxa + 1)

if (mpsgn (x) < 0) then
  i1 = (nna + 1) / 2
  call mpmdc (arr(0:mpnw+5,i1), d1, n1, mpnw1)
  mpnw1 = mpnw1 + (n1 + 1) / mpnbt + 1
endif

call mpeq (x, t1, mpnw1)
call mpeq (t1, t2, mpnw1)

do k = 2, nna
  call mpmul (x, t1, t3, mpnw1)
  call mpeq (t3, t1, mpnw1)
  call mpmul (arr(0:mpnw+5,k), t1, t3, mpnw1)
  call mpadd (t2, t3, t4, mpnw1)
  call mpeq (t4, t2, mpnw1)
enddo

call mpdmc (1.d0, 0, t3, mpnw1)
call mpsub (t3, x, t4, mpnw1)
call mpnpwr (t4, nna + 1, t3, mpnw1)
call mpdiv (t2, t3, t4, mpnw1)
call mproun (t4, mpnw)
call mpeq (t4, y, mpnw)

return
end subroutine mppolylogneg

subroutine mppolylogpos (nn, x, y, mpnw)

!   This returns polylog (nn, x) for the case nn >= 0.

implicit none
integer, intent(in):: nn, mpnw
integer ic1, itrmax, k, mpnw1
parameter (itrmax = 1000000)
integer (mpiknd), intent(in):: x(0:)
integer (mpiknd), intent(out):: y(0:)
integer (mpiknd) t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), t4(0:mpnw+6), &
  t5(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps (0:9)
integer (mpiknd) ix8

!  End of declaration.

if (nn < 0) then
  write (mpldb, 1)
1 format ('*** MPPOLYLOGPOS: N is less than zero.'/ &
  'For negative n, call mppolylogneg or polylog_neg. See documentation.')
  call mpabrt ( 590)
endif

if (mpnw < 4 .or. mpspacer (x) < mpnw + 4 .or. mpspacer (y) < mpnw + 6) then
  write (mpldb, 2)
2 format ('*** MPPOLYLOGPOS: Uninitialized or inadequately sized arrays')
  call mpabrt ( 591)
endif

mpnw1 = mpnw + 1
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

call mpabs (x, t1, mpnw)
call mpdmc (1.d0, 0, t2, mpnw)
call mpcpr (t1, t2, ic1, mpnw)
if (ic1 >= 0) then
  write (mpldb, 3)
3 format ('*** MPPOLYLOGPOS: |X| must be less than one.')
  call mpabrt ( 592)
endif

if (nn == 0) then
  call mpdmc (1.d0, 0, t1, mpnw1)
  call mpsub (t1, x, t2, mpnw1)
  call mpdiv (x, t2, t3, mpnw1)
  call mproun (t3, mpnw)
  call mpeq (t3, y, mpnw)
else
  call mpeq (x, t1, mpnw1)
  call mpeq (x, t2, mpnw1)

  do k = 2, itrmax
    call mpmul (x, t2, t3, mpnw1)
    call mpeq (t3, t2, mpnw1)
    call mpdmc (dble (k), 0, t3, mpnw1)
    call mpnpwr (t3, nn, t4, mpnw1)
    call mpdiv (t2, t4, t3, mpnw1)
    call mpadd (t1, t3, t4, mpnw1)
    call mpeq (t4, t1, mpnw1)

    call mpabs (t3, tc1, 4)
    call mpmul (eps, t1, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 100
  enddo

  write (mpldb, 4)
4 format ('*** MPPOLYLOGPOS: Loop end error')
  call mpabrt ( 593)

100 continue

  call mproun (t1, mpnw)
  call mpeq (t1, y, mpnw)
endif

return
end subroutine mppolylogpos

subroutine mpstruvehn (nu, ss, zz, mpnw)

!   This returns the StruveH function with integer arg NU and MPFR argument SS.

implicit none
integer, intent(in):: mpnw
integer ic1, itrmax, mpnw1, k, n1
real (mprknd) dmax, d1, pi
parameter (itrmax = 1000000, dmax = 1000.d0, pi = 3.1415926535897932385d0)
integer, intent(in):: nu
integer (mpiknd), intent(in):: ss(0:)
integer (mpiknd), intent(out):: zz(0:)
integer (mpiknd) sum(0:2*mpnw+6), td1(0:2*mpnw+6), td2(0:2*mpnw+6), tn1(0:2*mpnw+6), &
  tnm1(0:2*mpnw+6), t1(0:2*mpnw+6), t2(0:2*mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (ss) < mpnw + 4 .or. mpspacer (zz) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPSTRUVEHN: Uninitialized or inadequately sized arrays')
  call mpabrt ( 594)
endif

if (nu < 0) then
  write (mpldb, 2)
2 format ('*** MPSTRUVEHN: NU < 0')
  call mpabrt ( 595)
endif

call mpmdc (ss, d1, n1, mpnw)
d1 = abs (d1) * 2.d0**n1
if (d1 > dmax) then
  write (mpldb, 3)
3 format ('*** MPSTRUVEHN: ABS(SS) >',f8.2)
  call mpabrt ( 596)
endif

mpnw1 = mpnw * (1.d0 + d1 / dmax)

!   Check if Pi has been precomputed.

call mpmdc (mppicon, d1, n1, mpnw1)
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw1) then
  write (mpldb, 4) mpnw1
4 format ('*** MPSTRUVEHN: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 597)
endif

call mpinitwds (sum, 2*mpnw+1)
call mpinitwds (td1, 2*mpnw+1)
call mpinitwds (td2, 2*mpnw+1)
call mpinitwds (tn1, 2*mpnw+1)
call mpinitwds (tnm1, 2*mpnw+1)
call mpinitwds (t1, 2*mpnw+1)
call mpinitwds (t2, 2*mpnw+1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw*mpnbt, eps, 4)

! tn1 = mpreal (1.d0, nwds1)
! tnm1 = -0.25d0 * mpreal (ss, nwds1)**2
! td1 = 0.5d0 * sqrt (mppi (nwds1))
! td2 = td1

call mpdmc (1.d0, 0, tn1, mpnw1)
call mpmul (ss, ss, t1, mpnw1)
call mpmuld (t1, -0.25d0, tnm1, mpnw1)
call mpsqrt (mppicon, t1, mpnw1)
call mpmuld (t1, 0.5d0, td1, mpnw1)
call mpeq (td1, td2, mpnw1)

! do k = 1, nu
!  td2 = (k + 0.5d0) * td2
! enddo

do k = 1, nu
  d1 = k + 0.5d0
  call mpmuld (td2, d1, t1, mpnw1)
  call mpeq (t1, td2, mpnw1)
enddo

! sum = tn1 / (td1 * td2)

call mpmul (td1, td2, t2, mpnw1)
call mpdiv (tn1, t2, sum, mpnw1)

do k = 1, itrmax

!  tn1 = tnm1 * tn1
!  td1 = (k + 0.5d0) * td1
!  td2 = (nu + k + 0.5d0) * td2
!  t1 = tn1 / (td1 * td2)
!  sum = sum + t1

  call mpmul (tnm1, tn1, t1, mpnw1)
  call mpeq (t1, tn1, mpnw1)
  d1 = k + 0.5d0
  call mpmuld (td1, d1, t1, mpnw1)
  call mpeq (t1, td1, mpnw1)
  d1 = nu + k + 0.5d0
  call mpmuld (td2, d1, t1, mpnw1)
  call mpeq (t1, td2, mpnw1)
  call mpmul (td1, td2, t2, mpnw1)
  call mpdiv (tn1, t2, t1, mpnw1)
  call mpadd (sum, t1, t2, mpnw1)
  call mpeq (t2, sum, mpnw1)

!  if (abs (t1) < eps) goto 100

  call mpabs (t1, tc1, 4)
  call mpmul (eps, sum, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic1, 4)
  if (ic1 <= 0) goto 100
enddo

write (mpldb, 5)
5 format ('*** MPSTRUVEHN: Loop end error')
call mpabrt ( 598)

100 continue

! struvehn = (0.5d0 * ss)**(nu + 1) * sum

call mpmuld (ss, 0.5d0, t1, mpnw1)
n1 = nu + 1
call mpnpwr (t1, n1, t2, mpnw1)
call mpmul (t2, sum, t1, mpnw1)
call mproun (t1, mpnw)
call mpeq (t1, zz, mpnw)
return
end

subroutine mpzetar (ss, zz, mpnw)

!   This returns the zeta function of an MPR argument SS using an algorithm
!   due to Peter Borwein.

implicit none
integer, intent(in):: mpnw
integer i, ic1, iss, itrmax, j, mpnw1, n, n1, n2
real (mprknd) dfrac, d1, d2, pi
parameter (itrmax = 1000000, dfrac = 1.d0+ceiling(mpdpw), pi = 3.1415926535897932385d0)
integer (mpiknd), intent(in):: ss(0:)
integer (mpiknd), intent(out):: zz(0:)
integer (mpiknd) f1(0:mpnw+6), s(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), &
  t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6), tn(0:mpnw+6), tt(0:mpnw+6), &
  tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
real (mprknd) sgn
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (ss) < mpnw + 4 .or. mpspacer (zz) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPZETAR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 599)
endif

mpnw1 = mpnw + 1

!   Check if Pi has been precomputed.

call mpmdc (mppicon, d1, n1, mpnw1)
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw1) then
  write (mpldb, 4) mpnw1
4 format ('*** MPZETAR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 600)
endif

call mpinitwds (f1, mpnw1)
call mpinitwds (s, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (tt, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

call mpdmc (1.d0, 0, f1, mpnw1)
call mpcpr (ss, f1, ic1, mpnw1)
call mpinfr (ss, t1, t2, mpnw1)

if (ic1 == 0) then
  write (mpldb, 2)
2 format ('*** MPZETAR: argument is 1')
  call mpabrt ( 601)
elseif (mpsgn (t2) == 0) then

!   The argument is an integer value. Call mpzetaintr instead.

  call mpmdc (ss, d1, n1, mpnw)
  iss = d1 * 2.d0**n1
  call mpzetaintr (iss, t1, mpnw)
  goto 200
elseif (mpsgn (ss) < 0) then

!   If arg < 0, compute zeta(1-ss), and later apply Riemann's formula.

  call mpsub (f1, ss, tt, mpnw1)
else
  call mpeq (ss, tt, mpnw1)
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = mpnbt * mpnw * log (2.d0) / log (2.d0 * mpnbt * mpnw / 3.d0)
call mpmdc (tt, d2, n2, mpnw1)
d2 = d2 * 2.d0 ** n2

if (d2 > d1) then

!   Evaluate the infinite series.

  call mpdmc (1.d0, 0, t1, mpnw1)

  do i = 2, itrmax
    call mpdmc (dble (i), 0, t4, mpnw1)
    call mppower (t4, tt, t2, mpnw1)
    call mpdiv (f1, t2, t3, mpnw1)
    call mpadd (t1, t3, t2, mpnw1)
    call mpeq (t2, t1, mpnw1)

    call mpabs (t3, tc1, 4)
    call mpmul (eps, t1, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 200
  enddo

  write (mpldb, 3) 1, itrmax
3 format ('*** MPZETAR: iteration limit exceeded',2i10)
  call mpabrt ( 602)
endif

n = dfrac * mpnw1
call mpdmc (2.d0, 0, t1, mpnw1)
call mpnpwr (t1, n, tn, mpnw1)
call mpneg (tn, t1, mpnw1)
call mpdmc (0.d0, 0, t2, mpnw1)
call mpdmc (0.d0, 0, s, mpnw1)

sgn = 1.d0

do j = 0, 2 * n - 1
  call mpdmc (dble (j + 1), 0, t4, mpnw1)
  call mppower (t4, tt, t3, mpnw1)
  call mpdiv (t1, t3, t4, mpnw1)
  call mpmuld (t4, sgn, t5, mpnw1)
  call mpadd (s, t5, t4, mpnw1)
  call mpeq (t4, s, mpnw1)
  sgn = - sgn

  if (j .lt. n - 1) then
    call mpdmc (0.d0, 0, t2, mpnw1)
  elseif (j .eq. n - 1) then
    call mpdmc (1.d0, 0, t2, mpnw1)
  else
    call mpmuld (t2, dble (2 * n - j), t3, mpnw1)
    call mpdivd (t3, dble (j + 1 - n), t2, mpnw1)
  endif
  call mpadd (t1, t2, t3, mpnw1)
  call mpeq (t3, t1, mpnw1)
enddo

call mpsub (f1, tt, t3, mpnw1)
call mpdmc (2.d0, 0, t2, mpnw1)
call mppower (t2, t3, t4, mpnw1)
call mpsub (f1, t4, t2, mpnw1)
call mpmul (tn, t2, t3, mpnw1)
call mpdiv (s, t3, t1, mpnw1)
call mpneg (t1, t2, mpnw1)
call mpeq (t2, t1, mpnw1)

!   If original argument was negative, apply Riemann's formula.

if (mpsgn (ss) < 0) then
  call mpgammar (tt, t3, mpnw1)
  call mpmul (t1, t3, t2, mpnw1)
  call mpmul (mppicon, tt, t1, mpnw1)
  call mpmuld (t1, 0.5d0, t3, mpnw1)
  call mpcssnr (t3, t4, t5, mpnw1)
  call mpmul (t2, t4, t1, mpnw1)
  call mpmuld (mppicon, 2.d0, t2, mpnw1)
  call mppower (t2, tt, t3, mpnw1)
  call mpdiv (t1, t3, t2, mpnw1)
  call mpmuld (t2, 2.d0, t1, mpnw1)
endif

200 continue

call mproun (t1, mpnw)
call mpeq (t1, zz, mpnw)
return
end subroutine mpzetar

subroutine mpzetaintr (iss, zz, mpnw)

!   This returns the zeta function of the integer argument ISS using an algorithm
!   due to Peter Borwein.

implicit none
integer, intent(in):: mpnw
integer i, ic1, itrmax, j, mpnw1, n, n1, itt
real (mprknd) dfrac, d1, pi
parameter (itrmax = 1000000, dfrac = 1.d0+ceiling(mpdpw), pi = 3.1415926535897932385d0)
integer, intent(in):: iss
integer (mpiknd), intent(out):: zz(0:)
integer (mpiknd) f1(0:mpnw+6), s(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), &
  t3(0:mpnw+6), t4(0:mpnw+6), t5(0:mpnw+6), tn(0:mpnw+6), &
  tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
real (mprknd) sgn
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (zz) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPZETAINTR: Uninitialized or inadequately sized arrays')
  call mpabrt ( 603)
endif

mpnw1 = mpnw + 1

!   Check if Pi has been precomputed.

call mpmdc (mppicon, d1, n1, mpnw1)
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw1) then
  write (mpldb, 4) mpnw1
4 format ('*** MPZETAINTR: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 604)
endif

call mpinitwds (f1, mpnw1)
call mpinitwds (s, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (tn, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

call mpdmc (1.d0, 0, f1, mpnw1)

if (iss == 1) then
  write (mpldb, 2)
2 format ('*** MPZETAINTR: argument is 1')
  call mpabrt ( 605)
elseif (iss == 0) then

!   Argument is zero -- result is -1/2.

  call mpdmc (-0.5d0, 0, t1, mpnw1)
  goto 200
elseif (iss < 0) then

!   If argument is a negative even integer, the result is zero.

  if (mod (iss, 2) == 0) then
    call mpdmc (0.d0, 0, t1, mpnw1)
    goto 200
  endif

!   Otherwise if arg < 0, compute zeta(1-is), and later apply Riemann's formula.

  itt = 1 - iss
else
  itt = iss
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = mpnbt * mpnw * log (2.d0) / log (2.d0 * mpnbt * mpnw / 3.d0)

if (itt > d1) then

!   Evaluate the infinite series.

  call mpdmc (1.d0, 0, t1, mpnw1)

  do i = 2, itrmax
    call mpdmc (dble (i), 0, t4, mpnw1)
    call mpnpwr (t4, itt, t2, mpnw1)
    call mpdiv (f1, t2, t3, mpnw1)
    call mpadd (t1, t3, t2, mpnw1)
    call mpeq (t2, t1, mpnw1)

    call mpabs (t3, tc1, 4)
    call mpmul (eps, t1, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 200
  enddo

  write (mpldb, 3) 1, itrmax
3 format ('*** MPZETAINTR: iteration limit exceeded',2i10)
  call mpabrt ( 606)
endif

n = dfrac * mpnw1
call mpdmc (2.d0, 0, t1, mpnw1)
call mpnpwr (t1, n, tn, mpnw1)
call mpneg (tn, t1, mpnw1)
call mpdmc (0.d0, 0, t2, mpnw1)
call mpdmc (0.d0, 0, s, mpnw1)

sgn = 1.d0

do j = 0, 2 * n - 1
  call mpdmc (dble (j + 1), 0, t4, mpnw1)
  call mpnpwr (t4, itt, t3, mpnw1)
  call mpdiv (t1, t3, t4, mpnw1)
  call mpmuld (t4, sgn, t5, mpnw1)
  call mpadd (s, t5, t4, mpnw1)
  call mpeq (t4, s, mpnw1)
  sgn = - sgn

  if (j .lt. n - 1) then
    call mpdmc (0.d0, 0, t2, mpnw1)
  elseif (j .eq. n - 1) then
    call mpdmc (1.d0, 0, t2, mpnw1)
  else
    call mpmuld (t2, dble (2 * n - j), t3, mpnw1)
    call mpdivd (t3, dble (j + 1 - n), t2, mpnw1)
  endif

  call mpadd (t1, t2, t3, mpnw1)
  call mpeq (t3, t1, mpnw1)
enddo

call mpdmc (2.d0, 0, t2, mpnw1)
call mpnpwr (t2, 1 - itt, t4, mpnw1)
call mpsub (f1, t4, t2, mpnw1)
call mpmul (tn, t2, t3, mpnw1)
call mpdiv (s, t3, t1, mpnw1)
call mpneg (t1, t2, mpnw1)
call mpeq (t2, t1, mpnw1)

!   If original argument was negative, apply Riemann's formula.

if (iss < 0) then
  call mpdmc (1.d0, 0, t3, mpnw1)
  do i = 1, itt - 1
    call mpmuld (t3, dble (i), t4, mpnw1)
    call mpeq (t4, t3, mpnw1)
  enddo

  call mpmul (t1, t3, t2, mpnw1)
  call mpmuld (mppicon, dble (itt), t1, mpnw1)
  call mpmuld (t1, 0.5d0, t3, mpnw1)
  call mpcssnr (t3, t4, t5, mpnw1)
  call mpmul (t2, t4, t1, mpnw1)
  call mpmuld (mppicon, 2.d0, t2, mpnw1)
  call mpnpwr (t2, itt, t3, mpnw1)
  call mpdiv (t1, t3, t2, mpnw1)
  call mpmuld (t2, 2.d0, t1, mpnw1)
endif

200 continue

call mproun (t1, mpnw)
call mpeq (t1, zz, mpnw)
return
end subroutine mpzetaintr

subroutine mpzetabe (nb1, nb2, berne, s, z, mpnw)

!  This evaluates the Riemann zeta function, using the combination of
!  the definition formula (for large s), and an Euler-Maclaurin scheme
!  (see formula 25.2.9 of the DLMF). The array berne contains precomputed
!  even Bernoulli numbers (see MPBERNER above). Its dimensions must be as
!  shown below. NB2 must be greater than 1.4 x precision in decimal digits.

implicit none
integer, intent(in):: nb1, nb2, mpnw
integer i, i1, i2, ic1, itrmax, k, mpnw1, n1, n2, nn
real (mprknd) dber, dfrac, dlogb, d1, d2, pi
parameter (itrmax = 1000000, dber = 1.4d0, dfrac = 8.5d0, &
  dlogb = 33.27106466687737d0, pi = 3.1415926535897932385d0)
integer (mpiknd), intent(in):: berne(0:nb1+5,nb2), s(0:)
integer (mpiknd), intent(out):: z(0:)
integer (mpiknd) t0(0:mpnw+6), t1(0:mpnw+6), t2(0:mpnw+6), t3(0:mpnw+6), &
  t4(0:mpnw+6), t5(0:mpnw+6), t6(0:mpnw+6), t7(0:mpnw+6), t8(0:mpnw+6), &
  t9(0:mpnw+6), tt(0:mpnw+6), f1(0:mpnw+6), tc1(0:9), tc2(0:9), tc3(0:9), eps(0:9)
integer (mpiknd) ix8

!  End of declaration.

if (mpnw < 4 .or. mpspacer (s) < mpnw + 4 .or. mpspacer (z) < mpnw + 6) then
  write (mpldb, 1)
1 format ('*** MPZETABE: Uninitialized or inadequately sized arrays')
  call mpabrt ( 607)
endif

i = 0
k = 0
mpnw1 = mpnw + 1

!   Check if Pi has been precomputed.

call mpmdc (mppicon, d1, n1, mpnw1)
if (n1 /= 1 .or. abs (d1 * 2.d0**n1 - pi) > mprdfz &
  .or. mpwprecr (mppicon) < mpnw1) then
  write (mpldb, 5) mpnw1
5 format ('*** MPZETABE: Pi must be precomputed to precision',i9,' words.'/ &
  'See documentation for details.')
  call mpabrt ( 608)
endif

call mpinitwds (t0, mpnw1)
call mpinitwds (t1, mpnw1)
call mpinitwds (t2, mpnw1)
call mpinitwds (t3, mpnw1)
call mpinitwds (t4, mpnw1)
call mpinitwds (t5, mpnw1)
call mpinitwds (t6, mpnw1)
call mpinitwds (t7, mpnw1)
call mpinitwds (t8, mpnw1)
call mpinitwds (t9, mpnw1)
call mpinitwds (tt, mpnw1)
call mpinitwds (f1, mpnw1)

call mpinitwds (tc1, 4)
call mpinitwds (tc2, 4)
call mpinitwds (tc3, 4)
call mpinitwds (eps, 4)
call mpdmc (2.d0, 0, tc1, 4)
call mpnpwr (tc1, -mpnw1*mpnbt, eps, 4)

!   Check if argument is 1 -- undefined.

call mpdmc (1.d0, 0, t0, mpnw1)
call mpcpr (s, t0, ic1, mpnw1)
if (ic1 == 0) then
  write (mpldb, 2)
2 format ('*** MPZETABE: argument is 1')
  call mpabrt ( 609)
endif

!   Check if berne array has been initialized.

i1 = 1000000000
i2 = 1000000000

do k = 1, nb2
  i1 = min (i1, mpspacer (berne(0:nb1+5,k)))
  i2 = min (i2, mpwprecr (berne(0:nb1+5,k)))
enddo

call mpmdc (berne(0:nb1+5,1), d1, n1, mpnw)
d1 = d1 * 2.d0 ** n1

if (i1 < mpnw + 6 .or. i2 < mpnw .or. abs (d1 - 1.d0 / 6.d0) > mprdfz .or. &
  nb2 < int (dber * mpdpw * (mpnw - 2))) then
  write (mpldb, 3) int (dber * mpdpw * (mpnw - 2))
3 format ('*** MPZETABE: Array of even Bernoulli coefficients must be initialized'/ &
   'with at least',i8,' entries.')
  call mpabrt ( 610)
endif

call mpdmc (1.d0, 0, f1, mpnw1)

!   Check if argument is zero.  If so, result is - 1/2.

if (mpsgn (s) == 0) then
  call mpdmc (-0.5d0, 0, t1, mpnw)
  goto 200
endif

!   Check if argument is negative.

if (mpsgn (s) < 0) then

!   Check if argument is a negative even integer.  If so, the result is zero.

  call mpmuld (s, 0.5d0, t1, mpnw1)
  call mpinfr (t1, t2, t3, mpnw1)
  if (mpsgn (t3) == 0) then
    call mpdmc (0.d0, 0, t1, mpnw1)
    goto 200
  endif

!   Otherwise compute zeta(1-s), and later apply the reflection formula.

  call mpsub (f1, s, tt, mpnw1)
else
  call mpeq (s, tt, mpnw1)
endif

!  Check if argument is large enough that computing with definition is faster.

d1 = dlogb * mpnw1 / log (32.d0 * mpnw1)
call mpmdc (tt, d2, n2, mpnw1)
d2 = d2 * 2.d0**n2
if (d2 > d1) then
  call mpdmc (1.d0, 0, t1, mpnw1)

  do i = 2, itrmax
    call mpdmc (dble (i), 0, t4, mpnw1)
    call mppower (t4, tt, t2, mpnw1)
    call mpdiv (f1, t2, t3, mpnw1)
    call mpadd (t1, t3, t2, mpnw1)
    call mpeq (t2, t1, mpnw1)

    call mpabs (t3, tc1, 4)
    call mpmul (eps, t1, tc3, 4)
    call mpabs (tc3, tc2, 4)
    call mpcpr (tc1, tc2, ic1, 4)
    if (ic1 <= 0) goto 200
  enddo

  write (mpldb, 4) 1, itrmax
4 format ('*** MPZETABE: iteration limit exceeded',2i10)
  call mpabrt ( 611)
endif

call mpdmc (1.d0, 0, t0, mpnw1)
nn = dfrac * mpnw1

do k = 2, nn
  call mpdmc (dble (k), 0, t2, mpnw1)
  call mppower (t2, tt, t1, mpnw1)
  call mpdiv (f1, t1, t2, mpnw1)
  call mpadd (t0, t2, t3, mpnw1)
  call mpeq (t3, t0, mpnw1)
enddo

call mpdmc (dble (nn), 0, t2, mpnw1)
call mpsub (tt, f1, t3, mpnw1)
call mpmul (t1, t3, t4, mpnw1)
call mpdiv (t2, t4, t3, mpnw1)
call mpadd (t0, t3, t2, mpnw1)
call mpdmc (0.5d0, 0, t3, mpnw1)
call mpdiv (t3, t1, t4, mpnw1)
call mpsub (t2, t4, t0, mpnw1)

call mpeq (tt, t3, mpnw1)
d1 = 12.d0 * dble (nn)
call mpmuld (t1, d1, t4, mpnw1)
call mpdiv (t3, t4, t2, mpnw1)
call mpmuld (t1, dble (nn), t5, mpnw1)
call mpdmc (dble (nn), 0, t6, mpnw1)
call mpmul (t6, t6, t9, mpnw1)

do k = 2, min (nb2, itrmax)
  call mpdmc (dble (2 * k - 2), 0, t4, mpnw1)
  call mpadd (tt, t4, t6, mpnw1)
  call mpdmc (dble (2 * k - 3), 0, t7, mpnw1)
  call mpadd (tt, t7, t8, mpnw1)
  call mpmul (t6, t8, t7, mpnw1)
  call mpmul (t3, t7, t4, mpnw1)
  call mpdmc (dble (2 * k - 1), 0, t6, mpnw1)
  call mpdmc (dble (2 * k - 2), 0, t7, mpnw1)
  call mpmul (t6, t7, t8, mpnw1)
  call mpdiv (t4, t8, t3, mpnw1)
  call mpmul (t5, t9, t6, mpnw1)
  call mpeq (t6, t5, mpnw1)
  call mpmul (t3, berne(0:nb1+5,k), t4, mpnw1)
  call mpmuld (t5, dble (2 * k), t6, mpnw1)
  call mpdiv (t4, t6, t7, mpnw1)
  call mpadd (t2, t7, t4, mpnw1)
  call mpeq (t4, t2, mpnw1)

  call mpabs (t7, tc1, 4)
  call mpmul (eps, t2, tc3, 4)
  call mpabs (tc3, tc2, 4)
  call mpcpr (tc1, tc2, ic1, 4)
  if (ic1 <= 0) goto 110
enddo

write (mpldb, 3) 2, min (nb2, itrmax)
call mpabrt ( 612)

110 continue

call mpadd (t0, t2, t1, mpnw1)

!   If original argument was negative, apply the reflection formula.

if (mpsgn (s) < 0) then
  call mpgammar (tt, t3, mpnw1)
  call mpmul (t1, t3, t2, mpnw1)
  call mpmul (mppicon, tt, t1, mpnw1)
  call mpmuld (t1, 0.5d0, t3, mpnw1)
  call mpcssnr (t3, t4, t5, mpnw1)
  call mpmul (t2, t4, t1, mpnw1)
  call mpmuld (mppicon, 2.d0, t2, mpnw1)
  call mppower (t2, tt, t3, mpnw1)
  call mpdiv (t1, t3, t2, mpnw1)
  call mpmuld (t2, 2.d0, t1, mpnw1)
endif

200 continue

call mproun (t1, mpnw)
call mpeq (t1, z, mpnw)

return
end subroutine mpzetabe

end module mpfune
