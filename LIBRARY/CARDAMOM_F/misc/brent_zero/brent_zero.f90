module brent_zero
  contains
function zbrent ( called_from, f, a, b,  t_2, ftol )

!*****************************************************************************80
!
!! ZERO seeks the root of a function F(X) in an interval [A, B].
!
!  Discussion:
!
!    The interval [A, B] must be a change of sign interval for F.
!    That is, F(A) and F(B) must be of opposite signs.  Then
!    assuming that F is continuous implies the existence of at least
!    one value C between A and B for which F(C) = 0.
!
!    The location of the zero is determined to within an accuracy
!    of 6*MACHEPS*abs ( C ) + 2*T.
!
!    Thanks to Thomas Secretin for pointing out a transcription error in the
!    setting of the value of P, 11 February 2013.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    11 February 2013
!
!  Author:
!
!    Original FORTRAN77 version by Richard Brent.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Richard Brent, 
!    Algorithms for Minimization Without Derivatives, 
!    Dover, 2002, 
!    ISBN: 0-486-41998-3, 
!    LC: QA402.5.B74.
!
!  Parameters:
!
!    Input, real ( kind = dp ) A, B, the endpoints of the change of 
!    sign interval.
!
!
!    Input, real ( kind = dp ) T_2, a positive error tolerance, 2*T.

!    Input, real ( kind = dp ) ftol, tolerance on function magnitude.
!
!    Input, external real ( kind = dp ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose zero is being sought.
!
!    Output, real ( kind = dp ) zbrent, the estimated value of a zero of
!    the function F.
!
!    2024 jklebes modified for CARDAMOM: 
!       - pulled this function from brent.f90 and wrapped it in a module
!       - renamed from ZERO to zbrent
!       - kind 8 -> dp
!       - wrap function argument in interface.  Not pure.
!       - arguments number and order to match calls in CARDAMOM
!       - added loop counter and limit ITMAX = 8 after cardamom
!       - adjusted tolerances setting to match CARDAMOM
!         - added ftol
!         - MACHEP-brent had input argument, cardamom had hard-coded 6d-8, 
!                    I use intrinsic
!       - argument t halved on entry to match cardamom zbrent
  implicit none
  integer, parameter:: dp = kind(1.d0)
  real ( kind = dp )  :: zbrent
  character(len=*), intent(in):: called_from  ! name of procedure calling (used to pass through for errors)
  real ( kind = dp ), intent(in):: a, b
  real ( kind = dp ) c
  real ( kind = dp ) d
  real ( kind = dp ) e
  real ( kind = dp ) fa
  real ( kind = dp ) fb
  real ( kind = dp ) fc
  real ( kind = dp ) m
  real ( kind = dp ) machep
  real ( kind = dp ) p
  real ( kind = dp ) q
  real ( kind = dp ) r
  real ( kind = dp ) s
  real ( kind = dp ) sa
  real ( kind = dp ) sb
  real ( kind = dp ), intent(in) ::  ftol   ! tolerance on magnitude of f
  real ( kind = dp ), intent(in) ::  t_2  ! input, = 2*t 
  real ( kind = dp )  ::  t  
  real ( kind = dp ) tol      ! for iteratively updated tolerance
  integer            :: iter
  integer, parameter:: ITMAX = 8

  
  interface
     function f( val )
      integer, parameter:: dp = selected_real_kind(15, 9)
      real ( kind = dp ), intent(in):: val
      real ( kind = dp )            :: f
    end function f
  end interface

  machep = epsilon(0d0)
  t = 0.5_dp*t_2

!
!  Make local copies of A and B.
!
  sa = a
  sb = b
  fa = f ( sa )
  fb = f ( sb )

  c = sa
  fc = fa
  e = sb-sa
  d = e

  !t = 0.5d0*tol

  do iter = 1, ITMAX

    if ( abs ( fc ) < abs ( fb ) ) then

      sa = sb
      sb = c
      c = sa
      fa = fb
      fb = fc
      fc = fa

    end if

    tol = 2.0D+00*machep*abs ( sb ) + t
    m = 0.5D+00 * ( c-sb )

    if ( abs ( m ) <= tol .or. abs(fb) < ftol ) then
      zbrent = sb
      return
    end if

    if ( abs ( e ) < tol .or. abs ( fa ) <= abs ( fb ) ) then

      e = m
      d = e

    else

      s = fb/fa

      if ( sa == c ) then

        p = 2.0D+00*m * s
        q = 1.0D+00-s

      else

        q = fa/fc
        r = fb/fc
        p = s * ( 2.0D+00*m * q * ( q-r ) - ( sb-sa ) * ( r-1.0D+00 ) )
        q = ( q-1.0D+00 ) * ( r-1.0D+00 ) * ( s-1.0D+00 )

      end if

      if ( 0.0D+00 < p ) then
        q = - q
      else
        p = - p
      end if

      s = e
      e = d

      if ( 2.0D+00*p < 3.0D+00*m * q-abs ( tol*q ) .and. &
        p < abs ( 0.5D+00*s * q ) ) then
        d = p/q
      else
        e = m
        d = e
      end if

    end if

    sa = sb
    fa = fb

    if ( tol < abs ( d ) ) then
      sb = sb+d
    else if ( 0.0D+00 < m ) then
      sb = sb+tol
    else
      sb = sb-tol
    end if

    fb = f ( sb )

    if ( ( 0.0D+00 < fb .and. 0.0D+00 < fc ) .or. &
         ( fb <= 0.0D+00 .and. fc <= 0.0D+00 ) ) then
      c = sa
      fc = fa
      e = sb-sa
      d = e
    end if

  end do

  zbrent = sb

  return
end function
end module
