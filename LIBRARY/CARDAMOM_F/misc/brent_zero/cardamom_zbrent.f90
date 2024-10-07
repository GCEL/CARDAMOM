module cardamom_zbrent
  contains
  
  ! The Cardamom-native zbrent function, formerly included in every model file. 
  ! Also ultimately descended form Brent "Algorithms for minimization without derivatives" (1973)
  ! A copy is kept here for reference and testing
  double precision function zbrent( called_from, func, x1, x2, tol, toltol)

    ! This is a bisection routine. When ZBRENT is called, we provide a    !
    ! reference to a particular function and also two values which bound  !
    ! the arguments for the function of interest. ZBRENT finds a root of  !
    ! the function (i.e. the point where the function equals zero), that  !
    ! lies between the two bounds.                                        !
    ! There are five exit conditions:                                     !
    ! 1) The first proposal for the root of the function equals zero      !
    ! 2) The proposal range has been reduced to less then tol             !
    ! 3) The magnitude of the function is less than toltol                !
    ! 4) Maximum number of iterations has been reached                    !
    ! 5) The root of the function does now lie between supplied bounds    !
    ! For a full description see Press et al. (1986).                     !

    implicit none

    ! arguments..
    character(len=*), intent(in):: called_from    ! name of procedure calling (used to pass through for errors)
    double precision, intent(in):: tol, toltol, x1, x2

    ! Interfaces are the correct way to pass procedures as arguments.
    interface
      double precision pure function func( xval )
        double precision, intent(in):: xval
      end function func
    end interface

    ! local variables..
    integer            :: iter
    integer, parameter:: ITMAX = 8
    double precision   :: a, b, c, d, e, fa, fb, fc, p, q, r, s, tol1, tol0, xm
    double precision, parameter:: EPS = 6d-8

    ! calculations...
    a  = x1
    b  = x2
    fa = func( a )
    fb = func( b )
    tol0 = tol*0.5d0

    ! Check that we haven't (by fluke) already started with the root...
    if ( abs(fa) < toltol ) then
        zbrent = a
        return
    elseif ( abs(fb) < toltol ) then
        zbrent = b
        return
    end if
    ! Ensure the supplied x-values give y-values that lie either
    ! side of the root and if not flag an error message...
!    if (fa*fb > 0d0) then
        ! tell me otherwise what is going on
!!       print*,"Supplied values must bracket the root of the function.",new_line('x'),  &
!!         "     ","You supplied x1:",x1, new_line('x'),                     &
!!         "     "," and x2:",x2, new_line('x'),                             &
!!         "     "," which give function values of fa :",fa, new_line('x'),  &
!!         "     "," and fb:",fb, " .",new_line('x'),                        &
!!         " zbrent was called by: ",trim(called_from)
!    end if
    c = b
    fc = fb

    do iter = 1, ITMAX

      ! If the new value (f(c)) doesn't bracket
      ! the root with f(b) then adjust it..
!      if ( sign(1d0, fb) .eq. sign(1d0, fc) ) then
      if (fb*fc > 0d0) then
        c  = a
        fc = fa
        d  = b-a
        e  = d
      end if
      if ( abs(fc) .lt. abs(fb) ) then
        a  = b
        b  = c
        c  = a
        fa = fb
        fb = fc
        fc = fa
      end if
      tol1 = EPS*abs(b) + tol0
      xm   = 0.5d0 * ( c-b )
!      if ( ( abs(xm) .le. tol1 ) .or. ( fb .eq. 0d0 ) ) then
      if ( ( abs(xm) .le. tol1 ) .or. ( abs(fb) < toltol ) ) then
        zbrent = b
        return
      end if
      if ( ( abs(e) .ge. tol1 ) .and. ( abs(fa) .gt. abs(fb) ) ) then
        s = fb/fa
        if ( a .eq. c ) then
          p = 2d0*xm*s
          q = 1d0-s
        else
          q = fa/fc
          r = fb/fc
          p = s * ( 2d0*xm*q * ( q-r ) - ( b-a ) * ( r-1d0 ) )
          q = ( q-1d0 ) * ( r-1d0 ) * ( s-1d0 )
        end if
        if ( p .gt. 0d0 ) q = -q
        p = abs( p )
        if ( (2d0*p) .lt. min( 3d0*xm*q-abs(tol1*q), abs(e*q) ) ) then
          e = d
          d = p/q
        else
          d = xm
          e = d
        end if
      else
        d = xm
        e = d
      end if
      a  = b
      fa = fb
      if ( abs(d) .gt. tol1 ) then
        b = b+d
      else
        b = b+sign( tol1, xm )
      end if
      fb = func(b)
    enddo

    zbrent = b

  end function zbrent
end module
