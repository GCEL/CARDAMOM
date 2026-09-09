module test_functions
  use samplers_shared, only: PARINFO
  ! A struct with parameter limits
  type(PARINFO):: PI_xy
  double precision:: x_ideal = 5.1
  double precision:: y_ideal = 5.0
  double precision:: x_lower = 0d0
  double precision:: x_upper = 6d0
  double precision:: y_lower = 1d0
  double precision:: y_upper = 9d0

contains

subroutine ll_normal(pars, npars, res, id) bind(C)
!! A test function E = A*(x-x_0)^2+B*(y-y_0)^2
  integer, intent(in):: npars
double precision, dimension(npars), intent(inout):: pars  ! has to be inout because of C compatibilty
double precision, intent(out):: res
integer, intent(in):: id
double precision:: x, y  ! the pars to fit
double precision:: x_0, y_0  ! The correct, energy/loglikelihood-minimizing answer will be x = x0, y = y0
double precision:: A, B
x = pars(1)
y = pars(2)
x_0 = x_ideal
y_0 = y_ideal
A = 1.0
B = 1.6  ! covariance matrix expected to have inversely proportional entries on diagonal
! and zeros on off-diagonal
res = -(A*(x-x_0)**2+B*(y-y_0)**2)
end subroutine ll_normal

subroutine ll_step(pars, npars, res, id) bind(C)
  !! A test function that is a stepped rectangular well
  !! It's "correct" with value ll = 0 for x=(0, 5) and y=(1, 9)
  !! and has penalty of-5 in steps for values outside of the target domain
  integer, intent(in):: npars
double precision, dimension(npars), intent(inout):: pars
double precision, intent(out):: res
  !! result : loglikelihood penalty
integer, intent(in) :: id
double precision:: x, y  ! x, y values recieved as pars vector
double precision:: x_1, y_1, x_2, y_2  ! the bounds of the target domain
x_1 = 0d0
x_2 = 6d0
y_1 = 1d0
y_2 = 9d0
x = pars(1)
y = pars(2)
res = 0d0
if (x < x_1) then
  res = res-5 * ceiling((x_1-x)/5)
else if (x > x_2) then
  res = res-5 * ceiling((x-x_2)/5)
endif
if (y < y_1) then
  res = res-5 * ceiling((y_1-y)/5)
else if (y > y_2) then
  res = res-5 * ceiling((y-y_2)/5)
endif
end subroutine ll_step

subroutine ll_bounded(pars, npars, res, id) bind(C)
  !! A test function that is quadratic potential,
  !! plus hard boundaries resticting it to the domain
  !! that can be found by ll_step
  integer, intent(in):: npars
double precision, dimension(npars), intent(inout):: pars
double precision, intent(out):: res
  !! result : loglikelihood penalty
integer, intent(in) :: id
double precision:: x, y  ! x, y values recieved as pars vector
double precision:: x_1, y_1, x_2, y_2  ! the bounds of the target domain
double precision:: P = 0d0
x_1 = 0d0
x_2 = 6d0
y_1 = 1d0
y_2 = 9d0
x = pars(1)
y = pars(2)
res = 0d0
A = 1.0
B = 1.6
if (x >= x_1 .and. x <= x_2 .and. y >= y_1 .and. y <= y_2 ) then
  res = -(A*(x-x_ideal)**2+B*(y-y_ideal)**2)
else
  res = log(P)  ! loglikelihood = -Infinity as signal to always reject
endif
end subroutine ll_bounded

subroutine init_PI()
  PI_xy%npars = 2
  if (.not. allocated(PI_xy%parmin)) allocate(PI_xy%parmin(PI_xy%npars))
  if (.not. allocated(PI_xy%parmax)) allocate(PI_xy%parmax(PI_xy%npars))
  PI_xy%parmin(1) = -3.0
  PI_xy%parmax(1) = 110
  PI_xy%parmin(2) = 2.5
  PI_xy%parmax(2) = 7.5

  if (.not. allocated(PI_xy%paradj)) allocate(PI_xy%paradj(PI_xy%npars))
  PI_xy%paradj(1) = 4.0
  PI_xy%paradj(2) = 0.0

  if (.not. allocated(PI_xy%fix_pars)) allocate(PI_xy%fix_pars(PI_xy%npars))
  PI_xy%fix_pars = .false.
end subroutine init_PI

end module test_functions
