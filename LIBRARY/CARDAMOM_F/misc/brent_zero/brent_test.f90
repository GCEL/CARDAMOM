program main

!*****************************************************************************80
!
!! MAIN is the main program for BRENT_TEST.
!
!  Discussion:
!
!    BRENT_TEST tests the BRENT library.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    14 October 2008
!
!  Author:
!
!    John Burkardt
!
  use brent_zero, only: zero
  use cardamom_zbrent, only: zbrent
  implicit none

  call timestamp ( )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BRENT_TEST'
  write ( *, '(a)' ) '  FORTRAN90 version.'
  write ( *, '(a)' ) '  Test the BRENT library.'

  call test_zero_all ( )
!
!  Terminate.
!
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'BRENT_TEST'
  write ( *, '(a)' ) '  Normal end of execution.'
  write ( *, '(a)' ) ' '
  call timestamp ( )

  stop 0
end
subroutine test_zero_all ( )

!*****************************************************************************80
!
!! TEST_ZERO_ALL tests Brent's zero finding routine on all test functions.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  interface
     pure function f_01( val )
      integer, parameter:: dp = selected_real_kind(15, 9)
      real ( kind = dp ), intent(in):: val
      real ( kind = dp )            :: f_01
    end function f_01
  end interface
  interface
     pure function f_02( val )
      integer, parameter:: dp = selected_real_kind(15, 9)
      real ( kind = dp ), intent(in):: val
      real ( kind = dp )            :: f_02
    end function f_02
  end interface
  interface
     pure function f_03( val )
      integer, parameter:: dp = selected_real_kind(15, 9)
      real ( kind = dp ), intent(in):: val
      real ( kind = dp )            :: f_03
    end function f_03
  end interface


  interface
     pure function f_04( val )
      integer, parameter:: dp = selected_real_kind(15, 9)
      real ( kind = dp ), intent(in):: val
      real ( kind = dp )            :: f_04
    end function f_04
  end interface


  interface
     pure function f_05( val )
      integer, parameter:: dp = selected_real_kind(15, 9)
      real ( kind = dp ), intent(in):: val
      real ( kind = dp )            :: f_05
    end function f_05
  end interface


  real ( kind = 8 ) machep
  real ( kind = 8 ) t
  real ( kind = 8 ) ftol

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'TEST_ZERO_ALL'
  write ( *, '(a)' ) '  Test the Brent ZERO routine, which seeks'
  write ( *, '(a)' ) '  a root of a function F(X)'
  write ( *, '(a)' ) '  in an interval [A, B].'

  machep = epsilon ( machep )
  t = machep
  ftol = machep

  a = 1.0D+00
  b = 2.0D+00

  call test_zero_one ( a, b, machep, t, ftol, f_01, &
    'f_01(x) = sin ( x ) - x/2' )

  a = 0.0D+00
  b = 1.0D+00

  call test_zero_one ( a, b, machep, t, ftol, f_02, &
    'f_02(x) = 2*x - exp ( - x )' )

  a = -1.0D+00
  b =  0.5D+00

  call test_zero_one ( a, b, machep, t, ftol, f_03, &
    'f_03(x) = x*exp ( - x )' )

  a =  0.0001D+00
  b =  20.0D+00

  call test_zero_one ( a, b, machep, t, ftol, f_04, &
    'f_04(x) = exp ( x ) - 1 / ( 100*x * x )' )

  a = -5.0D+00
  b =  2.0D+00

  call test_zero_one ( a, b, machep, t, ftol, f_05, &
    'f_05(x) = (x+3) * (x-1) * (x-1)' )

  return
end
subroutine test_zero_one ( a, b, machep, t, ftol, f, title )

!*****************************************************************************80
!
!! TEST_ZERO_ONE tests Brent's zero finding routine on one test function.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    12 April 2008
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) A, B, the two endpoints of the change of sign
!    interval.
!
!    Input, real ( kind = 8 ) MACHEP, an estimate for the relative machine
!    precision.
!
!    Input, real ( kind = 8 ) T, a positive error tolerance.
!
!    Input, external real ( kind = 8 ) F, the name of a user-supplied
!    function, of the form "FUNCTION F ( X )", which evaluates the
!    function whose zero is being sought.
!
!    Input, character ( LEN = * ) TITLE, a title for the problem.
!
  use brent_zero, only: zero
  use cardamom_zbrent, only: zbrent
  implicit none

  real ( kind = 8 ) a
  real ( kind = 8 ) b
  real ( kind = 8 ) fa
  real ( kind = 8 ) fb
  real ( kind = 8 ) fz
  real ( kind = 8 ) fz2
  real ( kind = 8 ) machep
  real ( kind = 8 ) t
  real ( kind = 8 ) ftol
  character ( len = *  ) title
  real ( kind = 8 ) z
  real ( kind = 8 ) z2
  interface
     pure function f( val )
      integer, parameter:: dp = selected_real_kind(15, 9)
      real ( kind = dp ), intent(in):: val
      real ( kind = dp )            :: f
    end function f
  end interface
  
  z = zero (title, f, a, b, t, ftol )
  z2 = zbrent (title, f, a, b, t*2, ftol )  ! tolerance is halved on entering cardamom zbrent
  fz = f ( z )
  fz2 = f ( z2 )
  fa = f ( a )
  fb = f ( b )

  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) trim ( title )
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) '      A                 Z             B'
  write ( *, '(a)' ) '    F(A)              F(Z)          F(B)'
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'brent_zero.f90'
  write ( *, '(a)' ) ' '
  write ( *, '(2x, f14.8, 2x, f14.8, 2x, f14.8)' ) a,  z,  b
  write ( *, '(2x, g14.6, 2x, g14.6, 2x, g14.6)' ) fa, fz, fb
  write ( *, '(a)' ) ' '
  write ( *, '(a)' ) 'cardamom_zbrent.f90'
  write ( *, '(a)' ) ' '
  write ( *, '(2x, f14.8, 2x, f14.8, 2x, f14.8)' ) a,  z2,  b
  write ( *, '(2x, g14.6, 2x, g14.6, 2x, g14.6)' ) fa, fz2, fb

  return
end
pure function f_01 ( x )

!*****************************************************************************80
!
!! F_01 evaluates sin ( x ) - x/2.
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) F_01, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) f_01
  real ( kind = 8 ), intent(in) ::  x

  f_01 = sin ( x ) - 0.5D+00*x

  return
end
pure function f_02 ( x )

!*****************************************************************************80
!
!! F_02 evaluates 2*x-exp(-x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) F_02, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) f_02
  real ( kind = 8 ), intent(in):: x

  f_02 = 2.0D+00*x - exp ( - x )

  return
end
pure function f_03 ( x )

!*****************************************************************************80
!
!! F_03 evaluates x*exp(-x).
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) F_03, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) f_03
  real ( kind = 8 ), intent(in):: x

  f_03 = x*exp ( - x )

  return
end
pure function f_04 ( x )

!*****************************************************************************80
!
!! F_04 evaluates exp(x) - 1 / (100*x*x).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) F_04, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) f_04
  real ( kind = 8 ), intent(in):: x

  f_04 = exp ( x ) - 1.0D+00/100.0D+00/x / x

  return
end
pure function f_05 ( x )

!*****************************************************************************80
!
!! F_05 evaluates (x+3)*(x-1)*(x-1).
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license. 
!
!  Modified:
!
!    07 March 1999
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the point at which F is to be evaluated.
!
!    Output, real ( kind = 8 ) F_05, the value of the function at X.
!
  implicit none

  real ( kind = 8 ) f_05
  real ( kind = 8 ), intent(in):: x

  f_05 = ( x+3.0D+00 ) * ( x-1.0D+00 ) * ( x-1.0D+00 )

  return
end
subroutine timestamp ( )

!*****************************************************************************80
!
!! TIMESTAMP prints the current YMDHMS date as a time stamp.
!
!  Example:
!
!    31 May 2001   9:45:54.872 AM
!
!  Licensing:
!
!    This code is distributed under the GNU LGPL license.
!
!  Modified:
!
!    18 May 2013
!
!  Author:
!
!    John Burkardt
!
!  Parameters:
!
!    None
!
  implicit none

  character ( len = 8 ) ampm
  integer ( kind = 4 ) d
  integer ( kind = 4 ) h
  integer ( kind = 4 ) m
  integer ( kind = 4 ) mm
  character ( len = 9 ), parameter, dimension(12):: month = (/ &
    'January  ', 'February ', 'March    ', 'April    ', &
    'May      ', 'June     ', 'July     ', 'August   ', &
    'September', 'October  ', 'November ', 'December ' /)
  integer ( kind = 4 ) n
  integer ( kind = 4 ) s
  integer ( kind = 4 ) values(8)
  integer ( kind = 4 ) y

  call date_and_time ( values = values )

  y = values(1)
  m = values(2)
  d = values(3)
  h = values(5)
  n = values(6)
  s = values(7)
  mm = values(8)

  if ( h < 12 ) then
    ampm = 'AM'
  else if ( h == 12 ) then
    if ( n == 0 .and. s == 0 ) then
      ampm = 'Noon'
    else
      ampm = 'PM'
    end if
  else
    h = h-12
    if ( h < 12 ) then
      ampm = 'PM'
    else if ( h == 12 ) then
      if ( n == 0 .and. s == 0 ) then
        ampm = 'Midnight'
      else
        ampm = 'AM'
      end if
    end if
  end if

  write ( *, '(i2, 1x, a, 1x, i4, 2x, i2, a1, i2.2, a1, i2.2, a1, i3.3, 1x, a)' ) &
    d, trim ( month(m) ), y, h, ':', n, ':', s, '.', mm, trim ( ampm )

  return
end
