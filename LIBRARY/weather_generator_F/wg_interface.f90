
  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! Created: O. Browne (University of Edinburgh)
  ! Subsequent contributions by:
  ! T. L. Smallman (University of Edinburgh, t.l.smallman@ed.ac.uk)
  ! For exceptions and further references see code within

  !! This provides an interface b/wn the wrapper and whatever !!
  !! weather generator we wish to call.                       !!

  !****************************************

  subroutine weathergeneratorinterface( latitude   , nos_days   , & ! in
                                        days_in_yr , time       , & ! in
                                        sat_avg_in , sat_max_in , & ! in
                                        sat_min_in , ppt_in     , & ! in
                                        swrad_in   , coa_in     , & ! in
                                        rh_avg_in  , rh_max_in  , & ! in
                                        rh_min_in  , wind_in    , & ! in
                                        sat_out    , ppt_out    , & ! out
                                        swrad_out  , coa_out    , & ! out
                                        rh_out     , wind_out    )  ! out

    use weather_generator

    implicit none

    ! arguments..
    ! location and timing information
    integer,intent(in)              :: nos_days       !
    real,dimension(nos_days),intent(in) :: latitude & ! degrees
                                          ,time       ! day of time series
    real,intent(in)                 :: days_in_yr

    ! declare daily input variables
    real,dimension(nos_days),intent(in) :: sat_avg_in & ! daily avg surface air temperature (oC)
                                             ,sat_max_in & ! daily max surface air temperature (oC)
                                             ,sat_min_in & ! daily min surface air temperature (oC)
                                             ,ppt_in     & ! mean (precip) (kg.m-2.s-1)
                                             ,swrad_in   & ! mean short wave in (W.m-2)
                                             ,coa_in     & ! mean CO2 (ppm)
                                             ,rh_avg_in  & ! daily avg rel humidity (frac)
                                             ,rh_max_in  & ! daily max rel humidity (frac)
                                             ,rh_min_in  & ! daily min rel humidity (frac)
                                             ,wind_in      ! mean wind (ms.-1)

    real,dimension(nos_days*nint(hours_per_day)),intent(inout) :: sat_out   & ! hourly surface air temperature (oC)
                                                                 ,ppt_out   & ! hourly precip (kg.m-2.s-1)
                                                                 ,swrad_out & ! hourly sw radiation (W.m-2)
                                                                 ,coa_out   & ! hourly CO2 (ppm)
                                                                 ,rh_out    & ! hourly rel humidity (frac)
                                                                 ,wind_out    !

    ! local variables..
    integer :: i, a
    real :: day_number
    real, dimension(nint(hours_per_day)) :: coa_hrly, ppt_hrly, rh_hrly, sat_hrly, sw_hrly, wind_hrly

    ! calculate needed timing information
    a = nint(hours_per_day)

    do i = 1 , nos_days

       ! calculate day-of-year..
       if (time(i) == days_in_yr) then
           day_number = days_in_yr
       else
           day_number = real(ceiling(amod(time(i), days_in_yr)))
       endif

       ! Call the weather-generator with the daily-means..
       call daily_to_hourly( latitude(i),   day_number,    days_in_yr,  &  ! in
                             coa_in(i),     sat_avg_in(i),       &         ! in
                             sat_max_in(i), sat_min_in(i),       &         ! in
                             ppt_in(i),     swrad_in(i),         &         ! in
                             rh_avg_in(i),  rh_max_in(i),        &         ! in
                             rh_min_in(i),  wind_in(i),          &         ! in
                             coa_hrly,      sat_hrly,            &         ! out
                             ppt_hrly,      sw_hrly,             &         ! out
                             rh_hrly,       wind_hrly            )         ! out

       ! Update the output variables..
         coa_out( (((i-1)*a)+1) : i*a ) =  coa_hrly
         ppt_out( (((i-1)*a)+1) : i*a ) =  ppt_hrly
          rh_out( (((i-1)*a)+1) : i*a ) =   rh_hrly
         sat_out( (((i-1)*a)+1) : i*a ) =  sat_hrly
       swrad_out( (((i-1)*a)+1) : i*a ) =   sw_hrly
        wind_out( (((i-1)*a)+1) : i*a ) = wind_hrly

    enddo

  end subroutine weathergeneratorinterface
