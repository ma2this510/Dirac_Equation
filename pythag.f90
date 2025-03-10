      function pythag(a,b)
USE mpmodule
      implicit none
      type (mp_real)    pythag,a,b
!
!     finds sqrt(a**2+b**2) without overflow or destructive underflow
!
      type (mp_real)    p,r,s,t,u,eps
      type (mp_real)    zero, one, epsilonn
!

      zero = '0.d0'
      one = '1.d0'
      eps = epsilonn(one)


      p = max(abs(a),abs(b))
      if (p == zero) then
         pythag = zero
      else
         r = min(abs(a),abs(b))/p
         if (r < eps) then
            pythag = p
         else
            pythag = p*sqrt(one+r**2)
         end if
      end if

      return
      end
 
