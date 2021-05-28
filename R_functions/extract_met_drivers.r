
###
## Function to extract meteorology data from global gridded data
###

# These functions are by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
# Exceptions are given within specific functions.

  est_potential_radiation <- function(doy,latitude,diurnal_range,mean_diuranal_range,surf_pressure_Pa,water_vapour_Pa,hardcode_Tt) {

    # Function to calculate an estimate of potential short wave radation (MJ.m-2.day-1).
    # Day of year and latitude (degrees) are used to first estimate potential clear skys radiation.
    # Cloudiness fraction is then estimated from dirnal range and used to scale down the estimate.
    # Ref: Thornton & Running (1999) Agri Forest Met 93, 211-228

    # The approach used here assumes vectoriation to generate 24 hour estimate

    # parameters
    So = 1360 ; deg_to_rad = pi/180.0 ; hours_of_day = seq(1,24,1) ; hours_in_day = 24 ; seconds_per_hour = 60*60
    sea_surface_Pa = 101325
    alpha = -6.1e-5 # Pa-1 impact on transmittance per Pa of water vapour
    Tnadir = 0.87 # maximum clear sky transmittance on a dry day

    # convert latitude in radians
    latitude_radians = latitude * deg_to_rad

    # calculate declination
    declination = - asin ( sin ( 23.45 * deg_to_rad ) * cos ( 2.0 * pi * ( doy + 10.0 ) / 365.0 ) )
    #declination = max(0.0,declination)
    # calculate angle of sun each hour
    hourangle = deg_to_rad * 15. * ( hours_of_day - 0.5 * hours_in_day ) * 24. / hours_in_day
    # estimate cosine of solar zenith angle
    cos_solar_zenith_angle = sin(latitude_radians) * sin(declination) + cos(latitude_radians) * cos(declination) * cos(hourangle)
    # calculate potential radiation (MJ.m-2.s-1)
    est_potential_radiation = So * cos_solar_zenith_angle
    est_potential_radiation = pmax(0.0,est_potential_radiation)
    if (hardcode_Tt == 0) {
        b0 = 0.031 ; b1 = 0.201 ; b2 = 0.185 ; B = b0 + b1 * exp(-b2 * mean_diuranal_range) ; C = 2.4
        # estimate optical air mass for given solar angle
        m0 = cos_solar_zenith_angle**-1
        # estimate Maximum total transmittance fraction
        Ttmax = Tnadir**m0
        if (surf_pressure_Pa != 0) { Ttmax = Ttmax ** (surf_pressure_Pa/sea_surface_Pa) }
        if (water_vapour_Pa != 0) { Ttmax = Ttmax + alpha*water_vapour_Pa }
        # estimate Total transmittance fraction
        Tt = Ttmax*(1-exp(-B * diurnal_range**C))
    } else {
        if (hardcode_Tt > Tnadir | hardcode_Tt < 0.1) {return(print("Inputted Transmittance fraction is outside of allowable bounds"))}
        Tt = hardcode_Tt
    }
    # upscale to MJ.m-2.day-1
    est_potential_radiation = sum(est_potential_radiation*Tt*seconds_per_hour)*1e-6

    # explicit return
    return(est_potential_radiation)

  } # function est_potential_radiation

daily_mean <-function(var, interval, missing_allowed) {

   # work out how many intervals fit
   # i.e. number of days possible
   nos_days=ceiling(length(var)/interval)
   output=array(NaN, dim=c(nos_days))
   b=0
   for (i in seq(1, nos_days)) {
        if (length(which(is.na(var[b:(b+interval)]))) < missing_allowed) {
          output[i] = mean(var[b:(b+interval)], na.rm=T)
        } else {
          output[i] = NaN
        }
        b = b + interval
   }
   # clean up
   rm(nos_days,b,i) ; gc()
   return(output)

} # end function daily_mean

daily_sum <-function(var, interval, missing_allowed) {

   # work out how many intervals fit
   # i.e. number of days possible
   nos_days=ceiling(length(var)/interval)
   output=array(NaN, dim=c(nos_days))
   b=0
   for (i in seq(1, nos_days)) {
       if (length(which(is.na(var[b:(b+interval)]))) < missing_allowed) {
           output[i]=sum(var[b:(b+interval)], na.rm=T)
       } else {
         output[i]=NaN
       }
       b=b+interval
   }
   # clean up
   rm(nos_days,b,i) ; gc()
   return(output)

} # end function daily sum

calc_photoperiod_sec<-function(lat,days){

   # Function calculates the day length in hours based on day of year and latitude (degrees).
   # the output is daylength converted to seconds
   # Ref: NEEDS TO BE ADDED

   declin    = - asin ( sin ( 23.45 * ( pi / 180 ) ) * cos ( 2. * pi * ( days + 10. ) / 365. ) )
   sinld     = sin ( lat*(pi/180.) ) * sin ( declin )
   cosld     = cos ( lat*(pi/180.) ) * cos ( declin )
   aob       = sinld / cosld
   aob       = pmax(-1.0,pmin(1.0,sinld / cosld))
   daylength = 12.0 * ( 1. + 2. * asin ( aob ) / pi )
   # convert hours to seconds
   daylength=daylength*3600
   # clean up
   rm(declin,sinld,cosld,aob) ; gc()
   # now return
   return(daylength)

} # end function calc_photoperiod_sec

sp_humidity_to_vpd<-function(sp_moist,atmos_press,air_temperature) {

   # Converts specific humidity (kg/kg) into VPD (Pa)
   # Determine vapour pressure (Pa); based on specific humidity, air temperature (oC),
   # air pressure (Pa->(hPa) and ratio of gas constants (p95 McIlveen 1986, Basic Meteorology -
   # a physical outline)

   # calculate vapour pressure of the air
   vpair=((sp_moist*(atmos_press*1.0e-2))/0.62197)*1.0e2
   # Saturation vapour pressure (Pa) calculation from Jones p110; uses
   # absolute air temperature (oC)
   vpsat = (0.061375*exp((17.502*air_temperature)/(240.97+air_temperature)))*1.0e4
   # Difference between the vapour pressure of saturation and air, i.e. the
   # VPD (Pa)
   vpd_pa = vpsat - vpair
   # clean up
   rm(vpair,vpsat) ; gc()
   # return to user
   return(vpd_pa)

} # end function sp_humidity_to_vpd

vpd_to_rh<-function(vpd_in,air_temperature) {

   # Converts VPD (Pa) to rel humidity (frac)
   # Determine vapour pressure (Pa)"," based on specific humidity, air
   # pressure (Pa input)
   # (hPa) and ratio of gas constants (p95 McIlveen 1986, Basic Meteorology -
   # a physical outline)

   # Saturation vapour pressure (Pa) calculation from Jones p110"," uses
   # absolute air temperature (oC)
   vpsat = (0.061375*exp((17.502*air_temperature)/(240.97+air_temperature)))*1.0e4
   # Difference between the vapour pressure of saturation and air, i.e. the
   # VPD (Pa)
   rh = vpd_in/vpsat
   rh[rh > 1] <- 1
   rh[rh < 0] <- 0
   # clean up
   rm(vpsat) ; gc()
   # return to user
   return(rh)

} # end function vpd_to_rh

dew_temp_to_sp_humidity<-function(dew_airt,airt,pressure) {

  # Ref: p95 McIlveen 1986, Basic Meteorology -
  # a physical outline

  # dew_airt (oC)
  # airt (oC)
  # pressure (Pa->hPa)
  # vapour pressues (hPa->Pa)

  ## Specific humidity (kg/kg)
  vapour_pressure = 6.11*10**((7.5*dew_airt)/(237.3+dew_airt))
  dew_temp_to_sp_humidity = 0.622*vapour_pressure/(pressure*1e-2)

} # end function dew_temp_to_sp_humidity

rh_to_vpd<-function(rh_in,air_temperature) {

   # Converts rel humidity (frac) into VPD (Pa)
   # Determine vapour pressure (Pa)"," based on specific humidity, air
   # pressure (Pa input)
   # (hPa) and ratio of gas constants (p95 McIlveen 1986, Basic Meteorology -
   # a physical outline)

   # Saturation vapour pressure (Pa) calculation from Jones p110"," uses
   # absolute air temperature (oC)
   vpsat = (0.061375*exp((17.502*air_temperature)/(240.97+air_temperature)))*1.0e4
   # Difference between the vapour pressure of saturation and air, i.e. the
   # VPD (Pa)
   vpair = vpsat*rh_in
   vpd_pa = vpsat-vpair
   # clean up
   rm(vpsat) ; gc()
   # return to user
   return(vpd_pa)

} # end function rh_to_vpd

sp_humidity_to_rh <- function(qair, temp, press = 1013.25){

   # Ref: p95 McIlveen 1986, Basic Meteorology -
   # a physical outline

   # temperature in oC, specific humidity pressure in (hPa)
   es <-  6.112 * exp((17.67 * temp)/(temp + 243.5))
   e <- qair * press / (0.378 * qair + 0.622)
   rh <- e / es
   rh[rh > 1] <- 1
   rh[rh < 0] <- 0
   # clean up
   rm(es,e) ; gc()
   return(rh)
}

nos_days_in_year<-function(year) {

    # is current year a leap or not
    nos_days = 365
    mod=as.numeric(year)-round((as.numeric(year)/4))*4
    if (mod == 0) {
        nos_days = 366
        mod=as.numeric(year)-round((as.numeric(year)/100))*100
        if (mod == 0) {
            nos_days  = 365
            mod=as.numeric(year)-round((as.numeric(year)/400))*400
            if (mod == 0) {
                nos_days  = 366
            }
        }
    }

    # clean up
    rm(mod) ; gc()

    # return to user
    return(nos_days)

} # function to determine the number of days in year

extract_met_drivers<-function(n,timestep_days,start_year,end_year,latlon_wanted,met_in,met_source,site_name) {

  if (met_source == "site_specific") {

      # load met driver file for site specific information
      infile=paste(path_to_site_obs,site_name,"_timeseries_met.csv",sep="")
      # time information
      doy = read_site_specific_obs("doy",infile) # decimal day

      if (doy[1] != -9999) {
          # how many time steps before we change days?
          i = floor(doy[1]) ; j = 1
          while (i == floor(doy[j])) {
             j = j + 1
          }
          steps_in_day = j-1
          input_step_size = 24 / steps_in_day # hours
          # check whether this is covering just a single day or multiple...
          deltaDay = doy[(steps_in_day+1)] - doy[1]
          if (deltaDay > 1) {
              # each time step it actually greater than a day,
              # we must now check that this aligns with our expectations
              # and make suitable correction
              if ((deltaDay - timestep_days[1]) > 1) {
                  # mis-match between expected time stpe and actual...
                  stop("mis-match in meteorology drivers between expected time step and actual...")
              }
              # if we have made it to this point then we must be ok,
              # now apply correction to determine the correct timing variables
              steps_in_day = 1/mean(timestep_days)
              input_step_size = 24 / steps_in_day # hours
          } # deltaDay > 1
      } else {
          # currently assumed defaults
          steps_in_day = 24   # steps per day
          input_step_size = 1 # hours
      } # doy[1] != -9999

      # Max, min and average timestep air temperatures (oC)
      maxt = read_site_specific_obs("maxt_C",infile)
      if (maxt[1] == -9999) {stop("No max temperature provided: maxt in C must be provided")}
      mint = read_site_specific_obs("mint_C",infile)
      if (mint[1] == -9999) {stop("No min temperature provided: mint in C must be provided")}
      airt = read_site_specific_obs("airt_C",infile)       # if no mean air temperature available assume mean of max / min
      if (airt[1] == -9999) {airt = (maxt + mint) * 0.5}

      swrad = read_site_specific_obs("swrad_Wm2",infile) # W.m-2
      if (swrad[1] == -9999) {
          # try and look for shortwave in MJ/m2/day
          swrad=read_site_specific_obs("swrad_MJm2day",infile) # MJ/m2/day
          if (swrad[1] == -9999) {stop("No short wave radiation found: either swrad_Wm2 (in W/m2) or swrad_MJm2day must be provided")}
      } else {
          # assume we have got the W/m2, which we need to conver to MJ/m2/day
          swrad = swrad * input_step_size * 3600 * 1e-6
      }

      # atmospheric CO2 concentration (ppm)
      co2 = read_site_specific_obs("co2_ppm",infile)
      # if no co2 provided assume global mean value
      if (co2[1] == -9999) {co2=rep(400,length.out=length(doy))}

      # liquid + ice precipitation (kgH2O.m-2.s-1)
      precip = read_site_specific_obs("precip_kgm2s",infile)

      # vapour pressure deficit kPa, units converted below
      vpd = read_site_specific_obs("vpd_kPa",infile)
      if (vpd[1] == -9999) {
          # if no VPD search for relative humidity
          vpd = read_site_specific_obs("rh_fraction",infile)
          if (vpd[1] == -9999) {stop('No vapour pressure deficit (vpd; kPa) vpd or relative humidity (rh; 0-1) information has been provided')}
          if (max(vpd) >= 1 | min(vpd) <= 0 ) {
              stop('Relative humidity is out of range (0-1)')
          } else {
              # assume humidity data is good and conver to VPD
              vpd = rh_to_vpd(vpd,airt)
          }
      } else {
          # assume VPD has been found successfully
          vpd = vpd * 1000 # kPa->Pa
      } # vpd or not?

      # Wind speed (m/s)
      wind_spd=read_site_specific_obs("wind_spd_ms",infile)
      # if now wind speed data use global mean
      if (wind_spd[1] == -9999) {wind_spd=rep(3.23,length.out=length(doy))} # CRU global mean wind speed (m.s-1)

      # Construct day of year time series
      years_to_load=as.numeric(start_year):as.numeric(end_year)
      for (yr in seq(1,length(years_to_load))) {
           if (yr  == 1) {
               doy = 1:nos_days_in_year(years_to_load[yr])
           } else {
               doy = append(doy,1:nos_days_in_year(years_to_load[yr]))
           }
      }
      # timestep_days is of length == 1 when using a daily time step, for many of the calculations below
      # required timestep_days to be a vector of length equal to the analysis...make it so...
      if (length(timestep_days) == 1) {timestep_days = rep(timestep_days, length.out = length(doy))}

      # declare output variables
      maxt_out = 0 ; mint_out = 0 ; swrad_out = 0 ; co2_out = 0 ; precip_out = 0 ; vpd_out = 0 ; avgTemp_out = 0 ; wind_spd_out = 0
      vpd_lagged_out = 0 ; photoperiod_out = 0 ; avgTmax_out = 0

      if (steps_in_day > 1) {
          # loop through days to generate daily mean values first
          # lagged variables for GSI calculated afterwards
          for (daily in seq(1,length(swrad),steps_in_day)) {
               if (maxt[1] != -9999 & mint[1] != -9999) {
                   maxt_out=append(maxt_out,max(maxt[daily:(daily+steps_in_day-1)]))
                   mint_out=append(mint_out,min(mint[daily:(daily+steps_in_day-1)]))
                   avgTemp_out=append(avgTemp_out,(mint[daily:(daily+steps_in_day-1)]+maxt[daily:(daily+steps_in_day-1)])*0.5)
                   avgTmax_out=append(avgTmax_out,max(maxt[daily:(daily+steps_in_day-1)]))
               } else {
                   maxt_out=append(maxt_out,max(airt[daily:(daily+steps_in_day-1)]))
                   mint_out=append(mint_out,min(airt[daily:(daily+steps_in_day-1)]))
                   avgTemp_out=append(avgTemp_out,mean(airt[daily:(daily+steps_in_day-1)]))
                   avgTmax_out=append(avgTmax_out,max(airt[daily:(daily+steps_in_day-1)]))
               }
               # Short wave radiation (W.m-2)
               swrad_out=append(swrad_out,sum(swrad[daily:(daily+steps_in_day-1)]))
               # precipitation mean over time period (kgH2O.m-2.s-1)
               precip_out=append(precip_out,mean(precip[daily:(daily+steps_in_day-1)]))
               # wind speed mean over time period (m/s)
               wind_spd_out=append(wind_spd_out,mean(wind_spd[daily:(daily+steps_in_day-1)]))
               # cumulative precip lagged over a given number of days, in this case 42
               co2_out=append(co2_out,mean(co2[daily:(daily+steps_in_day-1)]))
               vpd_out=append(vpd_out,mean(vpd[daily:(daily+steps_in_day-1)]))
         } # looping within days

         # remove initial values from datasets
         swrad_out=swrad_out[-1] ; maxt_out=maxt_out[-1]
         mint_out=mint_out[-1]   ; co2_out=co2_out[-1]
         precip_out=precip_out[-1]
         avgTemp_out=avgTemp_out[-1]
         avgTmax_out=avgTmax_out[-1]
         vpd_out=vpd_out[-1] ; wind_spd_out=wind_spd_out[-1]

      } else {

         # currently provided drivers cover greater than a day, so just pass drivers directly
         maxt_out = maxt ; mint_out = mint ; avgTemp_out = airt ; avgTmax_out = maxt
         swrad_out = swrad ; precip_out = precip ; wind_spd_out = wind_spd
         co2_out = co2 ; vpd_out = vpd

      } # if there are more than 1 time step per day...

      # determine the actual daily positions
      run_day_selector = cumsum(timestep_days)

      # to allow for consistency of the rolling mean calculations we need to expand each value
      # by the number of days each represents...
      avgTmax_out = rep(avgTmax_out, times = timestep_days)
      vpd_lagged_out = rep(vpd_out, times = timestep_days)

      # rolling averaged for GSI
      avg_days = 21 # assume that the first 21 days are just the actual values, We expect this should result in a small error only
      # create photoperiod information; add 21 days to the output
      photoperiod_out = calc_photoperiod_sec(latlon_wanted[1],c(seq(365,(365-(avg_days-2)),-1),doy))
      # now take the daily values and turn them into rolling 21 day averages
      photoperiod_out = rollapply(photoperiod_out,avg_days,mean,na.rm=FALSE)
      avgTmax_out = rollapply(avgTmax_out,avg_days,mean,na.rm=FALSE)
      vpd_lagged_out = rollapply(vpd_lagged_out,avg_days,mean,na.rm=FALSE)
      # GSI adjustment
      avgTmax_out = append(avgTmax_out[1:(avg_days-1)],avgTmax_out)
      vpd_lagged_out = append(vpd_lagged_out[1:(avg_days-1)],vpd_lagged_out)
      # construct output
      met = list(run_day=run_day_selector,mint=mint_out,maxt=maxt_out,swrad=swrad_out,co2=co2_out
                ,doy=doy[run_day_selector],precip=precip_out,avgTmax=avgTmax_out[run_day_selector]
                ,photoperiod=photoperiod_out[run_day_selector],vpd_lagged=vpd_lagged_out[run_day_selector]
                ,avgTemp=avgTemp_out,vpd=vpd_out,wind_spd=wind_spd_out)

  } else { # site specific or from global databases

      #
      # Extraction from global databases
      #

      # calculate approximate offset for time zone
      #offset = round(latlon_wanted[2] * 24 / 360, digits=0)
      ## should I be applying the offset here?

      # sub-select for sites
      swrad_out = met_in$swrad[n,] ; maxt_out = met_in$maxt[n,] ; precip_out = met_in$precip[n,]
      vpd_out = met_in$vpd[n,] ; pressure_out = met_in$pressure[n,]
      mint_out = met_in$mint[n,] ; wind_spd_out = met_in$wind_spd[n,]
      co2_out = met_in$co2 ; avgTmax_out = maxt_out
      avgTemp_out = (mint_out + maxt_out) * 0.5

      # user update
      print(paste("Met data extracted for current location ",Sys.time(),sep=""))

      avg_days = 21 # assume that the first 21 days are just the actual values
      # create photoperiod information; add 21 days to the output
      photoperiod_out=calc_photoperiod_sec(latlon_wanted[1],c(seq(365,(365-(avg_days-2)),-1),met_in$doy))

      # now take the daily values and turn them into rolling 21 day averages
      photoperiod_out=rollapply(photoperiod_out,avg_days,mean,na.rm=FALSE)
      avgTmax_out=rollapply(avgTmax_out,avg_days,mean,na.rm=FALSE)
      vpd_lagged_out=rollapply(vpd_out,avg_days,mean,na.rm=FALSE)

      # make adjustments to the time series should we have read in an extra year for GSI calculation
      if (met_in$extra_year) {

          # we now need to remove the additional portion of the datasets from the front of them
          adjustment_end=length(swrad_out) ; adjustment_begin=length(swrad_out)-(length(met_in$doy)-1)
          swrad_out=swrad_out[adjustment_begin:adjustment_end]
          maxt_out=maxt_out[adjustment_begin:adjustment_end]
          mint_out=mint_out[adjustment_begin:adjustment_end]
          avgTemp_out=avgTemp_out[adjustment_begin:adjustment_end]
          wind_spd_out=wind_spd_out[adjustment_begin:adjustment_end]
          precip_out=precip_out[adjustment_begin:adjustment_end]
          co2_out=co2_out[adjustment_begin:adjustment_end]
          adjustment_end=length(avgTmax_out) ; adjustment_begin=length(avgTmax_out)-(length(met_in$doy)-1)
          avgTmax_out=avgTmax_out[adjustment_begin:adjustment_end]
          vpd_lagged_out=vpd_lagged_out[adjustment_begin:adjustment_end]
          vpd_out=vpd_out[adjustment_begin:adjustment_end]

      } else {

          # if we do not have all the needed for simplisity we will assume that the first 21 days are repeated
          avgTmax_out=append(avgTmax_out[1:(avg_days-1)],avgTmax_out)
          vpd_lagged_out=append(vpd_lagged_out[1:(avg_days-1)],vpd_lagged_out)

      } # extra year?

      if (length(timestep_days) == 1 & timestep_days[1] == 1) {

          # well actually we do nothing
          run_day_selector=seq(1,length(met_in$run_day),timestep_days)

      } else {

          # generally this now deals with time steps which are not daily.
          # However if not monthly special case
          if (length(timestep_days) == 1) {
              run_day_selector = seq(1,length(met_in$run_day),timestep_days)
              timestep_days = rep(timestep_days, length.out=length(met_in$run_day))
          }

          print("...calculating monthly or weekly averages for met drivers")
          # determine the actual daily positions
          run_day_selector=cumsum(timestep_days)
          # create needed variables
          swrad_agg=array(NA,dim=length(run_day_selector)) ; maxt_agg=array(NA,dim=length(run_day_selector))
          mint_agg=array(NA,dim=length(run_day_selector))
          precip_agg=array(NA,dim=length(run_day_selector)) ; wind_spd_agg=array(NA,dim=length(run_day_selector))
          co2_agg=array(NA,dim=length(run_day_selector)) ; avgTmax_agg=array(NA,dim=length(run_day_selector))
          vpd_agg=array(NA,dim=length(run_day_selector)) ; photoperiod_agg=array(NA,dim=length(run_day_selector))
          vpd_lagged_agg=array(NA,dim=length(run_day_selector))
          avgTemp_agg=array(NA,dim=length(run_day_selector))
          for (y in seq(1,length(run_day_selector))) {
               swrad_agg[y]=mean(swrad_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]])
               maxt_agg[y]=mean(maxt_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]])
               mint_agg[y]=mean(mint_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]])
               avgTemp_agg[y]=mean(avgTemp_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]])
               wind_spd_agg[y]=mean(wind_spd_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]])
               precip_agg[y]=mean(precip_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]])
               co2_agg[y]=mean(co2_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]])
               avgTmax_agg[y]=mean(avgTmax_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]])
               vpd_agg[y]=mean(vpd_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]])
               vpd_lagged_agg[y]=mean(vpd_lagged_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]])
               photoperiod_agg[y]=mean(photoperiod_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]])
          }
          # update with new output information
          swrad_out=swrad_agg ; maxt_out=maxt_agg ; mint_out=mint_agg
          precip_out=precip_agg ; vpd_lagged_out=vpd_lagged_agg
          co2_out=co2_agg ; avgTmax_out=avgTmax_agg ; vpd_out=vpd_agg ; photoperiod_out=photoperiod_agg
          avgTemp_out=avgTemp_agg ; wind_spd_out = wind_spd_agg

          # clean up
          rm(y,swrad_agg,maxt_agg,mint_agg,precip_agg,co2_agg,avgTmax_agg,
             vpd_agg,photoperiod_agg,avgTemp_agg,wind_spd_agg)
          gc(reset=TRUE,verbose=FALSE)

      } # monthly aggregation etc

      # Temperature information for all analyses provide temperature in Kelvin.
      # At this final point convert all (except rolling mean) to oC for CARDAMOM
      mint_out = mint_out - 273.15 ; maxt_out = maxt_out - 273.15
      avgTemp_out = avgTemp_out - 273.15 ; avgTmax_out = avgTmax_out - 273.15

      # output variables
      met=list(run_day=met_in$run_day[run_day_selector],mint=mint_out,maxt=maxt_out,
               swrad=swrad_out,co2=co2_out,doy=met_in$doy[run_day_selector],precip=precip_out,
               avgTmax=avgTmax_out,photoperiod=photoperiod_out,vpd_lagged=vpd_lagged_out,
               avgTemp=avgTemp_out,vpd=vpd_out,wind_spd=wind_spd_out)

      # clean up
      rm(swrad_out,mint_out,maxt_out,precip_out,co2_out,vpd_lagged_out,
         photoperiod_out,avgTmax_out,avgTemp_out,vpd_out)

  } # conditional for site_specific or global datasets

  # return met drivers now combined for a single site
  return(met)

} # function end
