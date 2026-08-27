#########################################################################################
# CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
# CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to 
# assimilate observations and ecological theory to retrieve parameters for the 
# DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
# used as a fully integrated component of CARDAMOM or independently. 
# Copyright (C) 2024  University of Edinburgh,
#                     Mathew Williams (mat.williams@ed.ac.uk), 
#                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
# UoE = University of Edinburgh

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ########## File specific description ##########
# This function is to create binary input files for CARDMOM. 
# The original Matlab function was development by A. A. Bloom 
# (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman 
# (t.l.smallman@ed.ac.uk, UoE) & J. F. Exbrayat (formally UoE).
#
#########################################################################################

# TEMPLATE FOR ALL DALEC MCMC DATA files
# Static Elements: 1-50
# Parameter Priors: 51-150
# Parameter prior uncertainty: 151-250
# Parameter prior weightings: 251-350
# Other priors: 351-400
# Other prior uncertainties: 401-450
# Other prior weights: 451-500
# TEMPORAL DRIVERS & DATA: 501-end

# Define global variables to contain the variable names for the 'met' and 'obs' arrays
met_array_names <<- c("Total number of days ran (starts = 31 for January if monthly)",
                      "Daily minimum temperatures (C)",
                      "Daily maximum temperatures (C)",
                      "Incoming daily shortwave radiation (MJ/m2/d)",
                      "Atmospheric CO2 concentration (ppm)",
                      "Julian day of year for the middle of the timestep",
                      "Precipitation (kgH2O/m2/s)",
                      "Biomass removal (fraction)",
                      "Burned area (fraction)",
                      "21 day rolling mean of daily max temperature (C)",
                      "21 day rolling mean of daily photoperiod (s)",
                      "21 day rolling mean of daily VPD (Pa)",
                      "Biomass removal - management type",
                      "Daily mean temperature (C)",
                      "Wind speed (m/s)",
                      "Daily vapour pressure deficit (Pa)")

obs_array_names <<- c("GPP (gC/m2/day)",
                      "GPP variance (gC/m2/day)",
                      "Lag period over which to average (steps)",
                      "Leaf area index (m2/m2)",
                      "Leaf area index variance (m2/m2)",
                      "Lag period over which to average (steps)",
                      "Net Ecosystem Exchange of CO2 (gC/m2/day)",
                      "Net Ecosystem Exchange of CO2 variance (gC/m2/day)",
                      "Lag period over which to average (steps)",
                      "Fire C emission (gC/m2/day)",
                      "Fire C emission variance (gC/m2/day)",
                      "Lag period over which to average (steps)",
                      "Ecosystem respiration (Ra + Rh gC/m2/day)",
                      "Ecosystem respiration (Ra + Rh) variance (gC/m2/day)",
                      "Lag period over which to average (steps)",
                      "Foliar stock (gC/m2)",
                      "Foliar stock variance (gC/m2)",
                      "Lag period over which to average (steps)",
                      "Wood stock (above + below gC/m2)",
                      "Wood stock (above + below) variance (gC/m2)",
                      "Lag period over which to average (steps)",
                      "Fine root stock (gC/m2)",
                      "Fine root stock variance (gC/m2)",
                      "Lag period over which to average (steps)",
                      "Foliar + fine root litter stock (gC/m2)",
                      "Foliar + fine root litter stock variance (gC/m2)" ,
                      "Lag period over which to average (steps)",
                      "Soil organic matter stock (gC/m2)",
                      "Soil organic matter stock variance (gC/m2)",
                      "Lag period over which to average (steps)",
                      "Above ground biomass stock (gC/m2)",
                      "Above ground biomass stock variance (gC/m2)",
                      "Lag period over which to average (steps)",
                      "Fraction absorbed PAR (0-1)",
                      "Fraction absorbed PAR variance (0-1)",
                      "Lag period over which to average (steps)",
                      "Coarse root stock (gC/m2)",
                      "Coarse root stock variance (gC/m2)",
                      "Lag period over which to average (steps)",
                      "Evapotranspiration (kgH2O/m2/day)",
                      "Evapotranspiration variance (kgH2O/m2/day)",
                      "Lag period over which to average (steps)",
                      "Snow water equivalent (kgH2O/m2)",
                      "Snow water equivalent variance (kgH2O/m2)",
                      "Lag period over which to average (steps)",
                      "Net Biome Exchange (Reco + Fire - GPP) of CO2 (gC/m2/day)",
                      "Net Biome Exchange variance (gC/m2/day)",
                      "Lag period over which to average (steps)",
                      "Mean woody gross productivity (gC/m2/day)",
                      "Mean woody gross productivity variance (gC/m2/day)",
                      "Lag period over which to average (steps)",
                      "Mean woody loss (gC/m2/day)",
                      "Mean woody loss variance (gC/m2/day)",
                      "Lag period over which to average (steps)",
                      "Mean leaf litter flux (gC/m2/day)",
                      "Mean leaf litter flux variance (gC/m2/day)",
                      "Lag period over which to average (steps)",
                      "Mean woody net change (gC/m2/day)",
                      "Mean woody net change variance (gC/m2/day)",
                      "Lag period over which to average (steps)",
                      "Mean AGB gross productivity (gC/m2/day)",
                      "Mean AGB gross productivity variance (gC/m2/day)",
                      "Lag period over which to average (steps)",
                      "Mean AGB loss (gC/m2/day)",
                      "Mean AGB loss variance (gC/m2/day)",
                      "Lag period over which to average (steps)",
                      "Mean leaf litter flux (gC/m2/day)",
                      "Mean leaf litter flux variance (gC/m2/day)",
                      "Lag period over which to average (steps)",
                      "Mean AGB net change (gC/m2/day)",
                      "Mean AGB net change variance (gC/m2/day)",
                      "Lag period over which to average (steps)",
                      "Extracted C due to harvest (gC/m2/day)",
                      "Extracted C due to harvest variance (gC/m2/day)",
                      "Lag period over which to average (steps)",
                      "Surface soil moisture (0-30cm m3/m3)",
                      "Surface soil moisture variance (m3/m3)",
                      "Lag period over which to average (steps)")

# Logistic model
logistic_model<-function(conc, gradient, mid_point, ref_val) {

  # Gradient = maximum gradient at midpoint
  # ref_val = reference value, 

  # Increasing logistic function
  logistic_model = 1 - (1+exp(gradient*(conc-mid_point)))**(-1)
  logistic_model = ref_val * logistic_model
  return(logistic_model)
  
}           

binary_data<-function(met,OBS,file,EDC,lat_degrees,ctessel_pft,modelname,parameter_type,nopars,noyears) {

  # Inform the user
  if (use_parallel == FALSE) {print(paste("writing out binary...",Sys.time(),sep=""))}

  # set model ID
  if (modelname == "ACM") {
      modelid = 0
  } else if (modelname == "DALEC.D1.F2.001") {
      modelid = 1
  } else if (modelname == "DALEC.C1.D1.F2.P1.002"){
      modelid = 2
  } else if (modelname == "DALEC.A1.C1.D2.F2.H1.P1.003"){
      modelid = 3
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P1.004"){
      modelid = 4
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P1.R1.005") {
      modelid = 5
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P1.R1.006") {
      modelid = 6
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P2.R1.007") {
      modelid = 7
  } else if (modelname == "DALEC.A1.C2.D2.F2.H1.P3.R1.008") {
      modelid = 8
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P3.R1.009") {
      modelid = 9
  } else if (modelname == "DALEC.A1.C2.D2.F2.H1.P4.R2.010") {
      modelid = 10
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P4.R2.011") {
      modelid = 11
  } else if (modelname == "DALEC.C4.D1.F2.012") {
      modelid = 12
  } else if (modelname == "DALEC.C5.D1.F2.P1.013") {
      modelid = 13
  } else if (modelname == "DALEC.C3.M1.014") {
      modelid = 14
  } else if (modelname == "DALEC.A3.C3.H2.M1.015") {
      modelid = 15
  } else if (modelname == "DALEC.A3.H1.M2.016") {
      modelid = 16
  } else if (modelname == "DALEC.A3.H2.M2.017") {
      modelid = 17
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P2.018") {
      modelid = 18
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P2.R3.019"){
      modelid = 19
  } else if (modelname == "DALEC.A2.C1.D2.F2.H2.P1.020"){
      modelid = 20
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P5.021") {
      modelid = 21
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P6.022") {
      modelid = 22
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P7.R2.023") {
      modelid = 23
  } else if (modelname == "DALEC.A1.C1.D2.F2.H4.P1.024") {
      modelid = 24
  } else if (modelname == "DALEC...025") {
      modelid = 25
  } else if (modelname == "DALEC.A4.C6.D2.F2.H3.P10.026") {
      modelid = 26
  } else if (modelname == "DALEC_1005") {
      modelid = 27
  } else if (modelname == "DALEC_1005a") {
      modelid = 28
  } else if (modelname == "DALEC.A1.C1.D2.F2.H3.P1.029") {
      modelid = 29
  } else if (modelname == "DALEC.A3.C1.D2.F2.H2.P1.030") {
      modelid = 30
  } else if (modelname == "DALEC.A4.C6.D2.F2.H2.P11.031") {
      modelid = 31
  } else if (modelname == "DALEC.A1.C1.D2.F2.H6.P1.R5.032") {
      modelid = 32
  } else if (modelname == "DALEC.A4.C6.D2.F2.H3.P12.033") {
      modelid = 33
  } else if (modelname == "") {
      modelid = 34
  } else if (modelname == "") {
      modelid = 35
  } else if (modelname == "DALEC.A1.C8.D2.F2.H2.P1.R4.036") {
      modelid = 36
  } else if (modelname == "DALEC.A1.C1.D2.F2.H5.P1.037"){
      modelid = 37
  } else if (modelname == "") {
      modelid = 38
  } else if (modelname == "") {
      modelid = 39
  } # model ID assignment

  pass = TRUE
  if (min(met$mint) < -200) {pass = FALSE ; print(summary(met$mint)) ; print('mint error in binary_data')} # Celcius
  if (min(met$maxt) < -200) {pass = FALSE ; print(summary(met$maxt)) ; print('maxt error in binary_data')} # Celcius
  if (min(met$vpd) < 0) {pass = FALSE ; print(summary(met$vpd)) ; print('VPD error in binary_data')} # Pa
  if (min(met$swrad) < 0 | max(met$swrad) > 36) { pass = FALSE ; print(summary(met$swrad)) ; print('RAD error in binary_data')} # MJ/m2/day
  if (lat_degrees < -90 | lat_degrees > 90) { pass = FALSE ; print(lat_degrees) ; stop('Latitude passed to binary_data is not -90/90 degrees')} # degrees only

  # Assuming forcings and observations pass criterior we will generate the files
  if (pass) {

      # some drivers may be passed as single values assuming this will apply across the whole time series
      # if this is so we need to now adjust that
      if (length(OBS$forest_management) != length(met$run_day)) {
          # we will assume that this is a constant value and will now repeast it
          OBS$forest_management = array(OBS$forest_management, dim=c(length(met$run_day)))
      }

      steps_per_year = floor(length(met$run_day)[1] / noyears)

      # extract information from list to array
      if (modelname == "ACM") {
          MET = array(NA,dim=c(length(met$run_day),(length(met)+2)))
      } else {
          MET = array(NA,dim=c(length(met$run_day),(length(met)+3)))
      }

      # Assign forcings into the array for saving
      MET[,1]  = met$run_day
      MET[,2]  = met$mint           # (Average) daily minimum temperature (oC) 
      MET[,3]  = met$maxt           # (Average) daily maximum temperature (oC)
      MET[,4]  = met$swrad          # Average daily sum shortwave radiation (MJ/m2/d)
      MET[,5]  = met$co2#+200       # CO2 ppm
      MET[,6]  = met$doy
      MET[,7]  = pmax(0,met$precip) # kgH2O/m2/s
      MET[,8]  = OBS$deforestation  # fraction
      MET[,9]  = OBS$burnt_area     # fraction
      MET[,10] = met$avgTmax       # C
      MET[,11] = met$photoperiod   # Seconds
      MET[,12] = met$vpd_lagged    # Pa
      MET[,13] = OBS$forest_management # type
      if (modelname == "ACM") {
          MET[,14] = met$avgN        # gN/m2
          MET[,15] = met$lai         # m2/m2
          MET[,16] = met$lat         # degree (-90/90)
          MET[,17] = met$wind_spd    # m/s
          MET[,18] = met$vpd         # Pa
          MET[,19] = met$Rtot        # MPa
          MET[,20] = met$top_sand    # %
          MET[,21] = met$bot_sand    # %
          MET[,22] = met$top_clay    # %
          MET[,23] = met$bot_clay    # %
      } else {
          MET[,14] = met$avgTemp          # Celcius
          MET[,15] = pmax(0,met$wind_spd) # m/s
          MET[,16] = pmax(0,met$vpd)      # Pa
      } # acm or not

      # Create time series observation matrix, i.e. things with uncertainty associated.
      OBSMAT = array(-9999.0,dim=c(length(met$run_day),66))
      # Line makes the correct array size but with -9999 in place of all
      OBSMAT[,1]  = OBS$GPP                   # GPP (gC/m2/day)
      OBSMAT[,2]  = OBS$GPP_unc               # GPP variance (gC/m2/day)
      OBSMAT[,3]  = OBS$GPP_lag               # GPP lag (steps)
      OBSMAT[,4]  = OBS$LAI                   # Leaf area index (m2/m2)
      OBSMAT[,5]  = OBS$LAI_unc               # Leaf area index variance
      OBSMAT[,6]  = OBS$LAI_lag               # Leaf area index lag (steps)
      OBSMAT[,7]  = OBS$NEE                   # Net Ecosystem Exchange of CO2 (gC/m2/day)
      OBSMAT[,8]  = OBS$NEE_unc               # Net Ecosystem Exchange of CO2 variance
      OBSMAT[,9]  = OBS$NEE_lag               # Net Ecosystem Exchange of CO2 lag (steps)
      OBSMAT[,10] = OBS$Fire                  # Fire C emission (gC/m2day)
      OBSMAT[,11] = OBS$Fire_unc              # Fire C emission variance
      OBSMAT[,12] = OBS$Fire_lag              # Fire C emission lag (steps)
      OBSMAT[,13] = OBS$Reco                  # Ecosystem respiration (Ra + Rh; gC/m2/day)
      OBSMAT[,14] = OBS$Reco_unc              # Ecosystem respiration (Ra + Rh) variance
      OBSMAT[,15] = OBS$Reco_lag              # Ecosystem respiration (Ra + Rh) lag
      OBSMAT[,16] = OBS$Cfol_stock            # Foliar stock (gC/m2)
      OBSMAT[,17] = OBS$Cfol_stock_unc        # Foliar stock variance
      OBSMAT[,18] = OBS$Cfol_stock_lag        # Foliar stock lag (steps)
      OBSMAT[,19] = OBS$Cwood_stock           # Wood stock (above + below; gC/m2)
      OBSMAT[,20] = OBS$Cwood_stock_unc       # Wood stock (above + below) variance
      OBSMAT[,21] = OBS$Cwood_stock_lag       # Wood stock (above + below) lag (steps)
      OBSMAT[,22] = OBS$Croots_stock          # Fine root stock (gC/m2)
      OBSMAT[,23] = OBS$Croots_stock_unc      # Fine root stock variance
      OBSMAT[,24] = OBS$Croots_stock_lag      # Fine root stock lag (steps)
      OBSMAT[,25] = OBS$Clit_stock            # Foliar + fine root litter stock (gC/m2)
      OBSMAT[,26] = OBS$Clit_stock_unc        # Foliar + fine root litter stock variance
      OBSMAT[,27] = OBS$Clit_stock_lag        # Foliar + fine root litter stock lag (steps)
      OBSMAT[,28] = OBS$Csom_stock            # Soil organic matter stock (gC/m2)
      OBSMAT[,29] = OBS$Csom_stock_unc        # Soil organic matter stock variance
      OBSMAT[,30] = OBS$Csom_stock_lag        # Soil organic matter stock lag
      OBSMAT[,31] = OBS$Cagb_stock            # Above ground biomass stock (gC/m2)
      OBSMAT[,32] = OBS$Cagb_stock_unc        # Above ground biomass stock variance
      OBSMAT[,33] = OBS$Cagb_stock_lag        # Above ground biomass stock lag (steps)
      OBSMAT[,34] = OBS$fAPAR                 # Fraction absorbed PAR
      OBSMAT[,35] = OBS$fAPAR_unc             # Fraction absorbed PAR variance
      OBSMAT[,36] = OBS$fAPAR_lag             # Fraction absorbed PAR lag (steps)
      OBSMAT[,37] = OBS$Ccoarseroot_stock     # Coarse root stock (gC/m2)
      OBSMAT[,38] = OBS$Ccoarseroot_stock_unc # Coarse root stock variance
      OBSMAT[,39] = OBS$Ccoarseroot_stock_unc # Coarse root stock lag (steps)
      OBSMAT[,40] = OBS$ET                    # Evapotranspiration (kgH2O/m2/day)
      OBSMAT[,41] = OBS$ET_unc                # Evapotranspiration variance
      OBSMAT[,42] = OBS$ET_lag                # Evapotranspiration lag (steps)
      OBSMAT[,43] = OBS$SWE                   # Snow water equivalent (kgH2O/m2)
      OBSMAT[,44] = OBS$SWE_unc               # Snow water equivalent variance
      OBSMAT[,45] = OBS$SWE_lag               # Snow water equivalent lag (steps)
      OBSMAT[,46] = OBS$nbe                   # Net Biome Exchange (Reco + Fire - GPP) of CO2 (gC/m2/day)
      OBSMAT[,47] = OBS$nbe_unc               # Net Biome Exchange variance
      OBSMAT[,48] = OBS$nbe_lag               # Net Biome Exchange lag (steps)
      OBSMAT[,49] = OBS$Cwood_growth          # Mean woody productivity over lag period (gC/m2/day)
      OBSMAT[,50] = OBS$Cwood_growth_unc      # Mean woody productivity varince
      OBSMAT[,51] = OBS$Cwood_growth_lag      # Lag period over which to average  (steps)
      OBSMAT[,52] = OBS$Cwood_loss            # Mean woody loss over lag period (gC/m2/day)
      OBSMAT[,53] = OBS$Cwood_loss_unc        # Mean woody loss varince
      OBSMAT[,54] = OBS$Cwood_loss_lag        # Lag period over which to average  (steps)
      OBSMAT[,55] = OBS$foliage_to_litter     # Mean litter flux over lag period (gC/m2/day)
      OBSMAT[,56] = OBS$foliage_to_litter_unc # Mean litter flux varince
      OBSMAT[,57] = OBS$foliage_to_litter_lag # Lag period over which to average (steps)
      OBSMAT[,58] = OBS$Cwood_change          # Mean woody net change over lag period (gC/m2/day)
      OBSMAT[,59] = OBS$Cwood_change_unc      # Mean woody net change varince
      OBSMAT[,60] = OBS$Cwood_change_lag      # Lag period over which to average  (steps)
      OBSMAT[,61] = OBS$harvest               # Extracted C due to harvest over lag period (gC/m2/day)
      OBSMAT[,62] = OBS$harvest_unc           # Extracted C due to harvest varince
      OBSMAT[,63] = OBS$harvest_lag           # Lag period over which to average (steps)
      OBSMAT[,64] = OBS$soilwater             # Surface (0-30cm) soil water content (m3/m3)
      OBSMAT[,65] = OBS$soilwater_unc         # Surface (0-30cm) soil water content variance (m3/m3)
      OBSMAT[,66] = OBS$soilwater_lag         # Surface (0-30cm) soil water content lag (step)
      OBSMAT[,67] = OBS$Cagb_change          # Mean AGB net change over lag period (gC/m2/day)
      OBSMAT[,68] = OBS$Cagb_change_unc      # Mean AGB net change varince
      OBSMAT[,69] = OBS$Cagb_change_lag      # Lag period over which to average  (steps)
      OBSMAT[,70] = OBS$Cagb_growth          # Mean AGB productivity over lag period (gC/m2/day)
      OBSMAT[,71] = OBS$Cagb_growth_unc      # Mean AGB productivity varince
      OBSMAT[,72] = OBS$Cagb_growth_lag      # Lag period over which to average  (steps)
      OBSMAT[,73] = OBS$Cagb_loss            # Mean AGB loss over lag period (gC/m2/day)
      OBSMAT[,74] = OBS$Cagb_loss_unc        # Mean AGB loss varince
      OBSMAT[,75] = OBS$Cagb_loss_lag        # Lag period over which to average  (steps)

      # STATIC DATA (1-50)
      # Model ID      = static_data[1]; DALEC_CDEA, DALEC.A1.C2.D2.F2.H2.P4.R2. etc
      # LAT           = static_data[2]; Latitude of site(Degrees)
      # nodays        = static_data[3]; Number of days (or time steps) in simulation
      # nomet         = static_data[4]; Number of met variables
      # noobs         = static_data[5]; Number of observation streams
      # EDC           = static_data[6]; EDCs on (1) or off (0)
      # pft           = static_data[7]; CTESSEL plant functional type, only used by ACM_TESSEL
      # yield class   = static_data[8]; NOT IN USE
      # age           = static_data[9]; NOT IN USE
      # nos. pars     = static_data[10]; number of parameters to be optimised for model
      # random search = static_data[11]; force random starting points for all parameters
      # top_sand      = static_data[12]; top soil (0-30cm) sand % content
      # bot_sand      = static_data[13]; bottom soil (31cm-maxdepth) sand % content
      # top_clay      = static_data[14]; top soil (0-30cm) clay % content
      # bot_clay      = static_data[15]; bottom soil (31cm-maxdepth) clay % content

      # if force_random_search == 1 then CARDAMOM ignores parameter priors even if present in the file during the EDC initialisation
      force_random_search = -9999 #; OBS$age = -9999
      # pass static information
      static_data = rep(-9999.0,length.out=50)
      tmp = c(modelid,lat_degrees,dim(MET)[1],dim(MET)[2],dim(OBSMAT)[2],
              EDC,ctessel_pft,OBS$yield_class,OBS$age,nopars,force_random_search,
              OBS$top_sand[1],OBS$bot_sand[1],OBS$top_clay[1],OBS$bot_clay[1])
      static_data[1:length(tmp)] = tmp

      # Define model parameter prior information
      PARPRIORS = rep(-9999.0,length.out=100)   # Prior estimate
      PARPRIORUNC = rep(-9999.0,length.out=100) # Gaussian uncertainty estimate
      PARPRIORWEIGHT = rep(1,length.out=100) # Weighting factor, e.g. RaGPP is assumed to be applicable on annual basis, so if 12 years = 12
      #For all other multiparameter user-defined priors
      OTHERPRIORS = rep(-9999.0,length.out=50)
      OTHERPRIORUNC = rep(-9999.0,length.out=50)
      OTHERPRIORWEIGHT = rep(1,length.out=50) # Weighting factor, e.g. RaGPP is assumed to be applicable on annual basis, so if 12 years = 12

      ###
      # Potentially useful information from now defunct and removed model structures
      ###
      # PARPRIORS[2]   = 51.70631  ; PARPRIORUNC[2]   = 124.6938  # C:N root (gC/gN) Kattge et al., (2011)
      # OTHERPRIORS[3] = 27.295    ; OTHERPRIORUNC[3] = 11.03755  # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
      # PARPRIORS[11]  = 0.2764618 ; PARPRIORUNC[11]  = 0.2014871 # log10 avg foliar N (gN.m-2)
      # PARPRIORS[15]  = 416.6667  ; PARPRIORUNC[15]  = 326.3762  # C:N wood (gC/gN) Kattge et al., (2011)
      # PARPRIORS[26]  = 11.197440 ; PARPRIORUNC[26]  = 1.32313*3 # NUE prior derived from Kattge et al., (2011)
      # PARPRIORS[43] = 275.1452   ; PARPRIORUNC[43] = 296.2767   # Leaf lifespan prior form Kattge et al., 2011, based on log10 gauusian distribution
      # OTHERPRIORS[6] = 0.5       ; OTHERPRIORUNC[6] = 0.05      # Rleaf:Rauto ratio See Atkin et al., various for approximate ratio
      # PARPRIORS[36]  = 1.639     ; PARPRIORUNC[36]  = 0.125     # Rm_leaf N**exponent (gC/gN) Reich et al., (2008)
      # PARPRIORS[37]  = 0.778     ; PARPRIORUNC[37]  = 0.133     # Rm_leaf intercept (gC/gN) Reich et al., (2008)
      # PARPRIORS[38]  = 1.352     ; PARPRIORUNC[38]  = 0.150     # Rm_root N**exponent (gC/gN) Reich et al., (2008)
      # PARPRIORS[39]  = 0.997     ; PARPRIORUNC[39]  = 0.082     # Rm_root intercept (gC/gN) Reich et al., (2008)
      # PARPRIORS[40]  = 1.344     ; PARPRIORUNC[40]  = 0.150     # Rm_wood N**exponent (gC/gN) Reich et al., (2008)
      # PARPRIORS[41]  = 0.946     ; PARPRIORUNC[41]  = 0.107     # Rm_wood intercept (gC/gN) Reich et al., (2008)
      # OTHERPRIORS[1] = 0.56      ; OTHERPRIORUNC[1] = 0.12      # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      ###

      # Assign model specific parameter priors
     if (modelname == "DALEC.D1.F2.001") {
          PARPRIORS[2]  = 0.54                 ; PARPRIORUNC[2] = 0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[10] = OBS$RhetQ10          ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)          
          PARPRIORS[11] = 16.9                 ; PARPRIORUNC[11] = 7.502147 # Ceff: derived from multiple trait values from Kattge et al., (2011)
          PARPRIORS[12] = OBS$lca              ; PARPRIORUNC[12] = OBS$lca_unc #; PARPRIORWEIGHT[12] = noyears
          PARPRIORS[13] = OBS$Cfol_initial     ; PARPRIORUNC[13] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[14] = OBS$Croots_initial   ; PARPRIORUNC[14] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[15] = OBS$Cwood_initial    ; PARPRIORUNC[15] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[16] = OBS$Clit_initial     ; PARPRIORUNC[16] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[17] = OBS$Csom_initial     ; PARPRIORUNC[17] = OBS$Csom_initial_unc # Csom prior
          # Other priors
#          OTHERPRIORS[5] =                    ; OTHERPRIORUNC[5] = # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[12] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4] = -9999 ; OBSMAT[filter,5] = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }             
      } else if (modelname == "DALEC.C1.D1.F2.P1.002") {

          PARPRIORS[2]  = 0.54                 ; PARPRIORUNC[2] = 0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[10] = OBS$RhetQ10          ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[11] = 16.9                 ; PARPRIORUNC[11] = 7.502147 # Ceff: derived from multiple trait values from Kattge et al., (2011)
          PARPRIORS[17] = OBS$lca              ; PARPRIORUNC[17] = OBS$lca_unc
          PARPRIORS[19] = OBS$Cfol_initial     ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial   ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial    ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial     ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial     ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          #PARPRIORS[24] = 0.5                  ; PARPRIORUNC[24] = 0.25 # Resilience factor
          #PARPRIORS[25] = 0.5                  ; PARPRIORUNC[25] = 0.25 # Foliar combustion completeness
          #PARPRIORS[26] = 0.1                  ; PARPRIORUNC[26] = 0.25 # Root / wood combustion completeness
          PARPRIORS[27] = 0.01                 ; PARPRIORUNC[27] = 0.05 # Soil combustion completeness
          # Other priors
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }               
      } else if (modelname == "DALEC.A1.C1.D2.F2.H1.P1.003") {
          PARPRIORS[2] = 0.54                              ; PARPRIORUNC[2]  = 0.12 #; PARPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[11] = 1.89*14.77735                    ; PARPRIORUNC[11] = 1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                               # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                               # These observational constraints are not the same and would lead to
                                                                               # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
#          PARPRIORS[11] = 65.0                             ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
          PARPRIORS[11] = 21.1491                           ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values 
                                                                                        # from Kattge et al., (2011)
                                                                                        # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                                        # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[24] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[24] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[26] = OBS$MaxRootDepth                 ; PARPRIORUNC[26] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[26] = 1.0                              ; PARPRIORUNC[26] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          #PARPRIORS[27] = 0.87                             ; PARPRIORUNC[27] = 0.41 # Resilience factor
          #PARPRIORS[28] = 0.5                              ; PARPRIORUNC[28] = 0.25 # Foliar combustion completeness
          #PARPRIORS[29] = 0.1                              ; PARPRIORUNC[29] = 0.25 # Root / wood combustion completeness
          PARPRIORS[30] = 0.01                             ; PARPRIORUNC[30] = 0.05 # Soil combustion completeness
          #PARPRIORS[31] = 0.25                             ; PARPRIORUNC[31] = 0.25 # Foliage + root litter combustion completeness
          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }                    
      } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P1.004") {
          PARPRIORS[2] = 0.54                              ; PARPRIORUNC[2]  = 0.12 #; PARPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Global Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[10] = 0.0334798                        ; PARPRIORUNC[10] = 0.015 #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[11] = 1.89*14.77735                    ; PARPRIORUNC[11] = 1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                               # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                               # These observational constraints are not the same and would lead to
                                                                               # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
#          PARPRIORS[11] = 65.0                             ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
          PARPRIORS[11] = 21.1491                            ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values 
                                                                                        # from Kattge et al., (2011)
                                                                                        # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                                        # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = min(19.9,max(0.36,OBS$MaxRootDepth))  ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
          PARPRIORS[12] = OBS$labile_release_timing         ; PARPRIORUNC[12] = OBS$labile_release_timing_unc # timing of max bud burst
          PARPRIORS[14] = OBS$labile_release_period         ; PARPRIORUNC[14] = OBS$labile_release_period_unc # timing of bud burst period
#          PARPRIORS[27] = 1.0                              ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          #PARPRIORS[28] = 0.87                             ; PARPRIORUNC[28] = 0.41 # Resilience factor
          #PARPRIORS[29] = 0.5                              ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
          #PARPRIORS[30] = 0.1                              ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
          PARPRIORS[31] = 0.01                             ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          #PARPRIORS[32] = 0.25                             ; PARPRIORUNC[32] = 0.25 # Foliage + root litter combustion completeness
          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          OTHERPRIORS[6] = OBS$MTTsom          ; OTHERPRIORUNC[6] = OBS$MTTsom_unc # Prior on the MTT of soil
          OTHERPRIORS[9] = OBS$leaflifespan    ; OTHERPRIORUNC[9] = OBS$leaflifespan_unc # Prior on the leaf lifespan
          # Hack to remove loss terms in pixels with disturbance
          if (max(MET[,8]) > 0 | max(MET[,9]) > 0) {
              OBSMAT[,52] = -9999            # Mean woody loss over lag period (gC/m2/day)
              OBSMAT[,53] = -9999            # Mean woody loss varince
              OBSMAT[,54] = -9999            # Lag period over which to average  (steps)          
          }
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }                    
      } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P1.R1.005") {
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[11] = 21.1491                          ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values 
                                                                                        # from Kattge et al., (2011)
                                                                                        # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                                        # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = OBS$MaxRootDepth                 ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                              ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          #PARPRIORS[28] = 0.87                             ; PARPRIORUNC[28] = 0.41 # Resilience factor
          #PARPRIORS[29] = 0.5                              ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
          #PARPRIORS[30] = 0.1                              ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
          PARPRIORS[31] = 0.01                             ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          #PARPRIORS[32] = 0.25                             ; PARPRIORUNC[32] = 0.25 # Foliage + root litter combustion completeness
          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[2] = 0.54                ; OTHERPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }                    
      } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P1.R1.006") {
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[11] = 21.1491                          ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values 
                                                                                        # from Kattge et al., (2011)
                                                                                        # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                                        # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = OBS$MaxRootDepth                 ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                              ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          PARPRIORS[31] = 0.01                             ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[2] = 0.54                ; OTHERPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#NOTUSE          OTHERPRIORS[3] = 27.295        ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }                    
      } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P2.R1.007") {
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[11] = 21.1491                          ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values 
                                                                                        # from Kattge et al., (2011)
                                                                                        # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                                        # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = OBS$MaxRootDepth                 ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                              ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          PARPRIORS[31] = 0.01                             ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[2] = 0.54                ; OTHERPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#NOTUSE          OTHERPRIORS[3] = 27.295        ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }                    
      } else if (modelname == "DALEC.A1.C2.D2.F2.H1.P3.R1.008") {
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[11] = 21.1491                          ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values 
                                                                                        # from Kattge et al., (2011)
                                                                                        # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                                        # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = OBS$MaxRootDepth                 ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                              ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          #PARPRIORS[30] = 0.87                             ; PARPRIORUNC[30] = 0.41 # Resilience factor
          #PARPRIORS[31] = 0.5                              ; PARPRIORUNC[31] = 0.25 # Foliar combustion completeness
          #PARPRIORS[32] = 0.1                              ; PARPRIORUNC[32] = 0.25 # Root / wood combustion completeness
          PARPRIORS[33] = 0.01                             ; PARPRIORUNC[33] = 0.05 # Soil combustion completeness
          #PARPRIORS[34] = 0.25                             ; PARPRIORUNC[34]] = 0.25 # Foliage + root litter combustion completeness
          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[2] = 0.54                ; OTHERPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#NOTUSE          OTHERPRIORS[3] = 27.295        ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }                    
      } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P3.R1.009") {
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[11] = 21.1491                          ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values 
                                                                                        # from Kattge et al., (2011)
                                                                                        # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                                        # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = OBS$MaxRootDepth                 ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                              ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          PARPRIORS[31] = 0.01                             ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[2] = 0.54                ; OTHERPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#NOTUSE          OTHERPRIORS[3] = 27.295        ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }                    
      } else if (modelname == "DALEC.A1.C2.D2.F2.H1.P4.R2.010") {
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[11] = 0.2764618                        ; PARPRIORUNC[11] = 0.2014871 # log10 avg foliar N (gN.m-2)
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = OBS$MaxRootDepth                 ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                              ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          PARPRIORS[31] = 0.01                             ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          PARPRIORS[40] = 11.197440                        ; PARPRIORUNC[40] = 9.3  # NUE prior derived from Kattge et al., (2011), based on log10 gaussian distribution
          PARPRIORS[43] = 275.1452                         ; PARPRIORUNC[43] = 296.2767 # Leaf lifespan prior form Kattge et al., 2011, based on log10 gauusian distribution
          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[2] = 0.54                ; OTHERPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#NOTUSE          OTHERPRIORS[3] = 27.295        ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }                    
      } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P4.R2.011") {
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[11] = 0.2764618                        ; PARPRIORUNC[11] = 0.2014871 # log10 avg foliar N (gN.m-2)
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = OBS$MaxRootDepth                 ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                              ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          PARPRIORS[31] = 0.01                             ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          PARPRIORS[41] = 11.197440                        ; PARPRIORUNC[41] = 9.3  # NUE prior derived from Kattge et al., (2011), based on log10 gaussian distribution
          PARPRIORS[44] = 275.1452                         ; PARPRIORUNC[44] = 296.2767 # Leaf lifespan prior form Kattge et al., 2011, based on log10 gauusian distribution
          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[2] = 0.54                            ; OTHERPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#NOTUSE          OTHERPRIORS[3] = 27.295        ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
          OTHERPRIORS[4] = 0.66                            ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential             ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }                    
      } else if (modelname == "DALEC.C4.D1.F2.012") {
          PARPRIORS[1]  = 0.54                 ; PARPRIORUNC[1]  = 0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[6]  = OBS$RhetQ10          ; PARPRIORUNC[6]  = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[7]  = 16.9                 ; PARPRIORUNC[7]  = 7.502147 # Ceff: derived from multiple trait values from Kattge et al., (2011)
          PARPRIORS[8]  = OBS$lca              ; PARPRIORUNC[8]  = OBS$lca_unc #; PARPRIORWEIGHT[8] = noyears
          PARPRIORS[9]  = OBS$Cfol_initial     ; PARPRIORUNC[9]  = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[10] = OBS$Cwood_initial    ; PARPRIORUNC[10] = OBS$Cwood_initial_unc # Croot + Cwood prior
          PARPRIORS[11] = OBS$Csom_initial     ; PARPRIORUNC[11] = OBS$Csom_initial_unc # Csom + Clitter prior
      } else if (modelname == "DALEC.C5.D1.F2.P1.013") {
          PARPRIORS[1]  = 0.54                 ; PARPRIORUNC[1]  = 0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[6]  = OBS$RhetQ10          ; PARPRIORUNC[6]  = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[7]  = 16.9                 ; PARPRIORUNC[7]  = 7.502147 # Ceff: derived from multiple trait values from Kattge et al., (2011)
          PARPRIORS[15] = OBS$Cfol_initial     ; PARPRIORUNC[15] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[16] = OBS$Cwood_initial    ; PARPRIORUNC[16] = OBS$Cwood_initial_unc # Croot + Cwood prior
          PARPRIORS[17] = OBS$Csom_initial     ; PARPRIORUNC[17] = OBS$Csom_initial_unc # Csom + Clitter prior
      } else if (modelname == "DALEC.C3.M1.014") {
       
          # Crop model specific adjustment - disallow observational constraints during the first 12 months
          OBSMAT[1:steps_per_year,1:dim(OBSMAT)[2]] = -9999
          # Remove early and late growing season LAI which are often pretty poor
          filter = which(MET[,6] < 95 | MET[,6] > 250)
          OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
          OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates

          PARPRIORS[11] = 11.197440              ; PARPRIORUNC[11] = 9.3 # NUE prior derived from Kattge et al., (2011), based on log10 gaussian distribution      
          PARPRIORS[13] = 0.21875                ; PARPRIORUNC[13] = 0.01 # Respiratory costs of labile transfer
          PARPRIORS[12] = OBS$planting_doy       ; PARPRIORUNC[12] = OBS$planting_doy_unc
          PARPRIORS[15] = OBS$growing_season_doy ; PARPRIORUNC[15] = OBS$growing_season_doy_unc
          PARPRIORS[17] = OBS$lca                ; PARPRIORUNC[17] = OBS$lca_unc
          PARPRIORS[19] = OBS$Cfol_initial       ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial     ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial      ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial       ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial       ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[35] = 0.99                   ; PARPRIORUNC[35] = 0.1 # labile turnover rate (/day)
          PARPRIORS[36] = 4.088                  ; PARPRIORUNC[36] = 0.6052851 # Constant for canopy N dilution model (gN/m2leaf)
          PARPRIORS[37] = -0.0252                ; PARPRIORUNC[37] = 0.00092   # Coefficient relating foliar C to N dilution
          # Other priors
          OTHERPRIORS[1] = 0.54                  ; OTHERPRIORUNC[1] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          # Yield:GPP Winter Wheat ATEC experiment plus He et al., (2018), doi: 10.3390/rs10030372
          # Values from He et al., Spring Wheat 0.24, Barley 0.42, Duram Wheat 0.22, Alfalfa 0.55, Pea 0.28, Maize 0.44
          OTHERPRIORS[8] = 0.38                ; OTHERPRIORUNC[8] = 0.087 #; OTHERPRIORWEIGHT[8] = noyears 
      } else if (modelname == "DALEC.A3.C3.H2.M1.015") {
    
          # Crop model specific adjustment - disallow observational constraints during the first 12 months
          OBSMAT[1:steps_per_year,1:dim(OBSMAT)[2]] = -9999
          # Remove early and late growing season LAI which are often pretty poor
          filter = which(MET[,6] < 95 | MET[,6] > 250)
          OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
          OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates

#          # Parameter priors for Winter Wheat (yes something better needs to be done for the storing of these)
#          # derived from the ATEC experiment field (but not assimilated) or Sus et al., (2010)
#          PARPRIORS[2] = 0.44           ; PARPRIORUNC[2]  = 0.08         # Fraction of GPP allocated to autotrophic pool
#          PARPRIORS[3] = 0.04           ; PARPRIORUNC[3]  = 0.02         # Development rate coefficient DS 0-1
#          PARPRIORS[4] = 0.035          ; PARPRIORUNC[4]  = 0.02         # Development rate coefficient DS 1-2
#          PARPRIORS[7] = 0.03           ; PARPRIORUNC[7]  = 0.03         # Potential turnover rate of foliage due to self-shading (fraction/day)
#          PARPRIORS[8] = 22.5           ; PARPRIORUNC[8]  = 5.0          # No. of vernalisation days for plants to be 50 % vernalised
#          PARPRIORS[11] = 21.1491       ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = noyears # Ceff: derived from multiple trait values from Kattge et al., (2011)
#          PARPRIORS[13] = 125.0         ; PARPRIORUNC[13] = 10.0         # Phenological heat units for seed emergence (aka growing degree days)
#          #PARPRIORS[12]=OBS$planting_doy       ; PARPRIORUNC[12]=OBS$planting_doy_unc # Sow day of year, applied as p12%%365.25
#          #PARPRIORS[14]=OBS$growing_season_doy ; PARPRIORUNC[14]=OBS$growing_season_doy_unc  # Growing season length sowing->harvest days
#          if (max(OBS$LAI) > 0) {
#              # Prior on canopy N derived from ATEC experiment assuming max LAI is related to canopy N
#              PARPRIORS[15] = min(7.0,max(OBS$LAI) * 0.488 + 2.5131)
#              PARPRIORUNC[15] = 1.2 # mean confidence interval of linear regression for LAI ranges 1-6
#              #PARPRIORWEIGHT[15] = noyears
#          } else {
#              PARPRIORS[15] = 4.088        ; PARPRIORUNC[15] = 0.6052851    # Constant for canopy N dilution model (gN/m2leaf)
#          }
#          PARPRIORS[16] = -0.0252      ; PARPRIORUNC[16] = 0.00092      # Coefficient relating foliar C to N dilution
#          PARPRIORS[17] = 20.40        ; PARPRIORUNC[17] = 5            # Leaf carbon per unit area (gC/m2)
#                                                                        # Prior values from Penning de Vries (1982), Simulation of Ecophysiological Processes of Growth in Several Annual Crops
#                                                                        # Barley = 15.60, Maize = 21.60, Potato = 14.40, Rice = 21.12, Sorghum = 19.20, Soybean = 19.20, SugarBeet = 24.00, 
#                                                                        # Sugarcane = 33.60, Sunflower = 25.92, WinterWheat = 20.40 , SprintWheat = 24.00
#          #PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
#          PARPRIORS[19] = OBS$Cfol_initial       ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
#          PARPRIORS[20] = OBS$Croots_initial     ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
#          PARPRIORS[21] = OBS$Cwood_initial      ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
#          PARPRIORS[22] = OBS$Clit_initial       ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
#          PARPRIORS[23] = OBS$Csom_initial       ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
#          PARPRIORS[26] = 273.15       ; PARPRIORUNC[26] = 5.0          # min temperature for development (K)
#          PARPRIORS[27] = 308.15       ; PARPRIORUNC[27] = 5.0          # max temperature for development (K)
#          PARPRIORS[28] = 298.15       ; PARPRIORUNC[28] = 5.0          # optimum temperature for development (K)
#          PARPRIORS[29] = 271.85       ; PARPRIORUNC[29] = 5.0          # min temperature for vernalisation (K)
#          PARPRIORS[30] = 288.85       ; PARPRIORUNC[30] = 5.0          # max temperature for vernalisation (K)
#          PARPRIORS[31] = 278.05       ; PARPRIORUNC[31] = 5.0          # optimum temperature for vernalisation (K)
#          PARPRIORS[32] = 8.25         ; PARPRIORUNC[32] = 2.0          # critical photoperiod for development to begin (days)
#          PARPRIORS[33] = 0.25         ; PARPRIORUNC[33] = 0.2          # sensitivity to photoperiod above critical threshold
#          PARPRIORS[35] = 0.99         ; PARPRIORUNC[35] = 0.1          # autotrophic turnover rate (/day)
#          PARPRIORS[36] = 50.0         ; PARPRIORUNC[36] = 25.0         # fine root biomass (g/m2) needed to reach 50 % of maximum depth
#          PARPRIORS[37] = 1.5          ; PARPRIORUNC[37] = 0.25         # maximum rooting depth (m)

          # ATEC LAI and N fertiliser addition is a highly uncertain saturating function.
          # i.e. maxLAI = N addition (kgN/ha) * 0.011840 + 2.216000

          # Parameter priors for Winter Wheat (yes something better needs to be done for the storing of these)
          # derived from the ATEC experiment field, Sus et al., (2010), or updated based on daily CARDAMOM-DALEC.15 analysis
          if (max(OBS$LAI) > 0) {
              # Linear fit between allocation to Ra and LAI
              # R2 = 0.94      Estimate Std. Error t value Pr(>|t|)    
              #(Intercept)     0.50433    0.01257   40.13 1.55e-09 ***
              #max_lai_yield  -0.03998    0.00364  -10.98 1.15e-05 ***
              # Fraction of GPP allocated to autotrophic pool
              PARPRIORS[2] = max(OBS$LAI) * -0.03998 + 0.50433
              PARPRIORUNC[2] = 0.044 # mean confidence interval of linear regression for LAI ranges 1-8              
          } else {
              PARPRIORS[2] = 0.44           ; PARPRIORUNC[2]  = 0.08         # Fraction of GPP allocated to autotrophic pool
          }

          PARPRIORS[3] = 0.04           ; PARPRIORUNC[3]  = 0.02         # Development rate coefficient DS 0-1
          PARPRIORS[4] = 0.023          ; PARPRIORUNC[4]  = 0.02         # Development rate coefficient DS 1-2
          PARPRIORS[5] = 0.008          ; PARPRIORUNC[5]  = 0.03         # turnover rate foliage (frac/day)
          PARPRIORS[6] = 0.010          ; PARPRIORUNC[6]  = 0.03         # TOR stem* - 1% loss per year value (day-1)
          PARPRIORS[7] = 0.03           ; PARPRIORUNC[7]  = 0.03         # Potential turnover rate of foliage due to self-shading (fraction/day)
          PARPRIORS[8] = 22.5           ; PARPRIORUNC[8]  = 5.0          # No. of vernalisation days for plants to be 50 % vernalised
          if (max(OBS$LAI) > 0) {
              # Linear fit between NUE ~ max LAI
              # R2 = 0.8689   Estimate Std. Error t value Pr(>|t|)    
              #(Intercept)    9.5774     2.2492   4.258 0.003755 ** 
              #max_lai_yield  4.7810     0.6514   7.340 0.000157 ***
              PARPRIORS[11] = max(OBS$LAI) * 4.7810 + 9.5774
              PARPRIORUNC[11] = 4.0 # mean confidence interval of linear regression for LAI ranges 1-8
          } else {
              PARPRIORS[11] = 21.1491       ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = noyears # NUE: derived from multiple trait values from Kattge et al., (2011)
          }
          PARPRIORS[13] = 125.0         ; PARPRIORUNC[13] = 10.0         # Phenological heat units for seed emergence (aka growing degree days)
          #PARPRIORS[12]=OBS$planting_doy       ; PARPRIORUNC[12]=OBS$planting_doy_unc # Sow day of year, applied as p12%%365.25
          #PARPRIORS[14]=OBS$growing_season_doy ; PARPRIORUNC[14]=OBS$growing_season_doy_unc  # Growing season length sowing->harvest days
          if (max(OBS$LAI) > 0) {
              # Linear fit between initial fol N ~ max LAI
              #R = 0.88       Estimate Std. Error t value Pr(>|t|)    
              #(Intercept)    1.9607     0.4154   4.720  0.002157 ***
              #max_lai_yield  0.9473     0.1203   7.875  0.000101 ***              
              #PARPRIORS[15] = min(8.0,max(OBS$LAI) * 0.9473 + 1.9607)
              #PARPRIORUNC[15] = 1.5  # mean confidence interval of linear regression for LAI ranges 1-8
              # Logistic fit between initial fol N ~ max LAI
              # R2  = 0.94  
              # gradient mid_point   ref_val 
              #   0.7865    1.9120    7.0713 
              # residual sum-of-squares: 0.5801
              PARPRIORS[15] = logistic_model(max(OBS$LAI), gradient = 0.7865, mid_point = 1.912, ref_val = 7.0713)
              PARPRIORUNC[15] = 0.58  # mean confidence interval of linear regression for LAI ranges 1-8
              #PARPRIORWEIGHT[15] = noyears
          } else {
              # Prior on canopy N derived from ATEC experiment assuming max LAI is related to canopy N          
              PARPRIORS[15] = 4.088        ; PARPRIORUNC[15] = 0.6052851    # Constant for canopy N dilution model (gN/m2leaf)
          }
          PARPRIORS[16] = -0.0252      ; PARPRIORUNC[16] = 0.00092      # Coefficient relating foliar C to N dilution
          PARPRIORS[17] = 20.40        ; PARPRIORUNC[17] = 5            # Leaf carbon per unit area (gC/m2)
                                                                        # Prior values from Penning de Vries (1982), Simulation of Ecophysiological Processes of Growth in Several Annual Crops
                                                                        # Barley = 15.60, Maize = 21.60, Potato = 14.40, Rice = 21.12, Sorghum = 19.20, Soybean = 19.20, SugarBeet = 24.00, 
                                                                       # Sugarcane = 33.60, Sunflower = 25.92, WinterWheat = 20.40 , SprintWheat = 24.00
          #PARPRIORS[17] = OBS$lca                ; PARPRIORUNC[17] = OBS$lca_unc
          #PARPRIORS[19] = OBS$Cfol_initial       ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          #PARPRIORS[20] = OBS$Croots_initial     ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          #PARPRIORS[21] = OBS$Cwood_initial      ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          #PARPRIORS[22] = OBS$Clit_initial       ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial       ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[26] = 273.15       ; PARPRIORUNC[26] = 5.0          # min temperature for development (K)
          PARPRIORS[27] = 308.15       ; PARPRIORUNC[27] = 5.0          # max temperature for development (K)
          PARPRIORS[28] = 291.76       ; PARPRIORUNC[28] = 2.3          # optimum temperature for development (K)
          PARPRIORS[29] = 270.47       ; PARPRIORUNC[29] = 1.5          # min temperature for vernalisation (K)
          PARPRIORS[30] = 288.85       ; PARPRIORUNC[30] = 5.0          # max temperature for vernalisation (K)
          PARPRIORS[31] = 278.05       ; PARPRIORUNC[31] = 5.0          # optimum temperature for vernalisation (K)
          PARPRIORS[32] = 8.25         ; PARPRIORUNC[32] = 2.0          # critical photoperiod for development to begin (days)
          PARPRIORS[33] = 0.25         ; PARPRIORUNC[33] = 0.2          # sensitivity to photoperiod above critical threshold          
          PARPRIORS[34] = 0.03         ; PARPRIORUNC[34] = 0.03         # rate of labile turnover (/day)
          PARPRIORS[35] = 0.99         ; PARPRIORUNC[35] = 0.1          # autotrophic turnover rate (/day)
          PARPRIORS[36] = 50.0         ; PARPRIORUNC[36] = 25.0         # fine root biomass (g/m2) needed to reach 50 % of maximum depth
          PARPRIORS[37] = 1.5          ; PARPRIORUNC[37] = 0.25         # maximum rooting depth (m)
          # Other priors
          OTHERPRIORS[1] = 0.54        ; OTHERPRIORUNC[1] = 0.12 #; OTHERPRIORWEIGHT[1] = noyears  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          # Yield:GPP Winter Wheat ATEC experiment plus He et al., (2018), doi: 10.3390/rs10030372
          # Values from He et al., Winter Wheat 0.38, Spring Wheat 0.24, Barley 0.42, Duram Wheat 0.22, Alfalfa 0.55, Pea 0.28, Maize 0.44
          OTHERPRIORS[8] = 0.38        ; OTHERPRIORUNC[8] = 0.087*0.5 #; OTHERPRIORWEIGHT[8] = noyears
      } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P2.018") {
          PARPRIORS[2]  = 0.54                 ; PARPRIORUNC[2] = 0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[10] = OBS$RhetQ10          ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[11] = 21.1491              ; PARPRIORUNC[11] = 8.534234 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                            #       Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                            # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca              ; PARPRIORUNC[17] = OBS$lca_unc
          PARPRIORS[19] = OBS$Cfol_initial     ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial   ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial    ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial     ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial     ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[27] = OBS$MaxRootDepth     ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
          #PARPRIORS[28] = 0.5                  ; PARPRIORUNC[28] = 0.25 # Resilience factor
          #PARPRIORS[29] = 0.5                  ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
          #PARPRIORS[30] = 0.1                  ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
          PARPRIORS[31] = 0.01                 ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          # Other priors
          #OTHERPRIORS[1] =       ; OTHERPRIORUNC[1] =  # Initial soil water fraction 
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }         
      } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P2.R3.019") {
          PARPRIORS[1]  = 0.5                  ; PARPRIORUNC[1] = 0.125 # fraction of litter decomposition to Csom
          PARPRIORS[10] = OBS$RhetQ10          ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[11] = 21.1491              ; PARPRIORUNC[11] = 8.534234 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                          # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                          # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca              ; PARPRIORUNC[17] = OBS$lca_unc
          PARPRIORS[19] = OBS$Cfol_initial     ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial   ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial    ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial     ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial     ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[27] = OBS$MaxRootDepth     ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
          PARPRIORS[31] = 0.01                 ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          # Other priors
          #OTHERPRIORS[1] =       ; OTHERPRIORUNC[1] =  # Initial soil water fraction 
          OTHERPRIORS[2] = 0.54                ; OTHERPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#NOTUSE          OTHERPRIORS[3] = 27.295        ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }         
  } else if (modelname == "DALEC.A2.C1.D2.F2.H2.P1.020") {
          PARPRIORS[2] = 0.54                              ; PARPRIORUNC[2]  = 0.12 #; PARPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[11] = 65.0                             ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
          PARPRIORS[11] = 24.43                            ; PARPRIORUNC[11]= 10.4 #; PARPRIORWEIGHT[11] = 1 # Vcmax: Median of CARDAMOM analysis assimilating 4 GPP products
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = OBS$MaxRootDepth                 ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                              ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          PARPRIORS[31] = 0.01                             ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }                    
      } else if (modelname == "DALEC_1005") {
          PARPRIORS[2]  = 0.54                   ; PARPRIORUNC[2]  = 0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[10] = OBS$RhetQ10            ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[11] = 16.9                   ; PARPRIORUNC[11] = 7.502147 # Ceff: derived from multiple trait values from Kattge et al., (2011)
          PARPRIORS[17] = OBS$lca                ; PARPRIORUNC[17] = OBS$lca_unc
          PARPRIORS[19] = OBS$Cfol_initial       ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial     ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial      ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial       ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial       ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          # Other priors
          OTHERPRIORS[5] = OBS$Cwood_potential   ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
      } else if (modelname == "DALEC_1005a") {
          PARPRIORS[2]  = 0.54                   ; PARPRIORUNC[2]  = 0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[10] = OBS$RhetQ10            ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[11] = 16.9                   ; PARPRIORUNC[11] = 7.502147 # Ceff: derived from multiple trait values from Kattge et al., (2011)
          PARPRIORS[17] = OBS$lca                ; PARPRIORUNC[17] = OBS$lca_unc
          PARPRIORS[19] = OBS$Cfol_initial       ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial     ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial      ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial       ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial       ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          # Other priors
          OTHERPRIORS[5] = OBS$Cwood_potential   ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A1.C1.D2.F2.H3.P1.029") {
          PARPRIORS[2] = 0.54                              ; PARPRIORUNC[2] = 0.12 #; PARPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#          PARPRIORS[11] = 1.89*14.77735                    ; PARPRIORUNC[11] = 1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                                # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                                # These observational constraints are not the same and would lead to
                                                                                # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
#          PARPRIORS[11] = 65.0                             ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
          PARPRIORS[11] = 21.1491                          ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                              # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                              # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17]=OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = OBS$MaxRootDepth                 ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
          PARPRIORS[33] = OBS$minLWP                       ; PARPRIORUNC[33] = OBS$minLWP_unc # minLWP prior
          PARPRIORS[31] = 0.01                             ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          # Other priors
          #OTHERPRIORS[1] =       ; OTHERPRIORUNC[1] =  # Initial soil water fraction 
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
      } else if (modelname == "DALEC.A4.C6.D2.F2.H2.P11.031") {
          PARPRIORS[10] = OBS$RhetQ10         ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[11] = 65.0               ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
#          PARPRIORS[11] = 60.0               ; PARPRIORUNC[11]= 20.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax: derived from multiple trait values from Kattge et al., (2011)
                                                                          # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                          # due to the different temperature response functions used in ACM2 vs ACM 1
#          PARPRIORS[11] = 24.43                  ; PARPRIORUNC[11]= 10.4 #; PARPRIORWEIGHT[11] = 1 # Vcmax: Median of CARDAMOM analysis assimilating 4 GPP products
          PARPRIORS[11] = 54.165                  ; PARPRIORUNC[11]= 20.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax: Mean of reported PFT values from Oliver et al., (2022)
                                                                          # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                          # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca                ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial       ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial     ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial      ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial       ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial       ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[27] = OBS$MaxRootDepth       ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                    ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, based on median from Fan et al., (2017) https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          PARPRIORS[31] = 0.01                 ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
#          PARPRIORS[33] = 0.05                 ; PARPRIORUNC[33] = 0.05 # labile:biomass at which growth limited by 50 %
          PARPRIORS[36] = 1.0                 ; PARPRIORUNC[36] = 5.0 # temperature at which foliage and root growth totally suppressed (oC)
          PARPRIORS[37] = 5.0                 ; PARPRIORUNC[37] = 1.0 # temperature at which wood growth totally suppressed (oC)
          PARPRIORS[43] = -2.0                ; PARPRIORUNC[43] = 0.5 # minimum leaf water potential (MPa)
          # Other priors
          #OTHERPRIORS[1] =       ; OTHERPRIORUNC[1] =  # Initial soil water fraction 
          OTHERPRIORS[2] = 0.54                ; OTHERPRIORUNC[2] = 0.12 #; OTHERPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4] = -9999 ; OBSMAT[filter,5] = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }           
      } else if (modelname == "DALEC.A1.C1.D2.F2.H6.P1.R5.032") {
          PARPRIORS[2] = 0.54                              ; PARPRIORUNC[2]  = 0.12 #; PARPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Global Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[10] = 0.0334798                        ; PARPRIORUNC[10] = 0.015 #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[11] = 1.89*14.77735                    ; PARPRIORUNC[11] = 1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                               # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                               # These observational constraints are not the same and would lead to
                                                                               # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
#          PARPRIORS[11] = 65.0                             ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
          PARPRIORS[11] = 21.1491                            ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values 
                                                                                        # from Kattge et al., (2011)
                                                                                        # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                                        # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = OBS$MaxRootDepth                 ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                              ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          #PARPRIORS[28] = 0.87                             ; PARPRIORUNC[28] = 0.41 # Resilience factor
          #PARPRIORS[29] = 0.5                              ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
          #PARPRIORS[30] = 0.1                              ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
          PARPRIORS[31] = 0.01                             ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          #PARPRIORS[32] = 0.25                             ; PARPRIORUNC[32] = 0.25 # Foliage + root litter combustion completeness



          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          OTHERPRIORS[6] = OBS$MTTsom          ; OTHERPRIORUNC[6] = OBS$MTTsom_unc # Prior on the MTT of soil
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }                                       
      } else if (modelname == "DALEC.A4.C6.D2.F2.H3.P12.033") {
          PARPRIORS[10] = OBS$RhetQ10            ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[11] = 65.0               ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
#          PARPRIORS[11] = 60.0               ; PARPRIORUNC[11]= 20.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax: derived from multiple trait values from Kattge et al., (2011)
                                                                          # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                          # due to the different temperature response functions used in ACM2 vs ACM 1
#          PARPRIORS[11] = 24.43                  ; PARPRIORUNC[11]= 10.4 #; PARPRIORWEIGHT[11] = 1 # Vcmax: Median of CARDAMOM analysis assimilating 4 GPP products
          PARPRIORS[11] = 54.165                  ; PARPRIORUNC[11]= 20.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax: Mean of reported PFT values from Oliver et al., (2022)
                                                                          # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                          # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[12] = -2.0                   ; PARPRIORUNC[12] = 0.5 # minimum leaf water potential (MPa)
          PARPRIORS[17] = OBS$lca                ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial       ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial     ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial      ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial       ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial       ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[27] = OBS$MaxRootDepth       ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                    ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, based on median from Fan et al., (2017) https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          PARPRIORS[31] = 0.01                   ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
#          PARPRIORS[33] = 0.01                 ; PARPRIORUNC[33] = 0.05 # labile:biomass at which growth limited by 50 %
          PARPRIORS[36] = 0.0                 ; PARPRIORUNC[36] = 5.0 # temperature at which foliage and root growth totally suppressed (oC)
          PARPRIORS[37] = 5.0                 ; PARPRIORUNC[37] = 1.0 # temperature at which wood growth totally suppressed (oC)
          # Other priors
          #OTHERPRIORS[1] =       ; OTHERPRIORUNC[1] =  # Initial soil water fraction 
          OTHERPRIORS[2] = 0.54                ; OTHERPRIORUNC[2] = 0.12 #; OTHERPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
              #if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4] = -9999 ; OBSMAT[filter,5] = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
              #}
          }                    
      } else if (modelname == "DALEC.A3.C1.D2.F2.H2.P1.030") {
          PARPRIORS[2] = 0.54                              ; PARPRIORUNC[2]  = 0.12 #; PARPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Global Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[10] = 0.0334798                        ; PARPRIORUNC[10] = 0.015 #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[11] = 1.89*14.77735                    ; PARPRIORUNC[11] = 1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                               # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                               # These observational constraints are not the same and would lead to
                                                                               # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
#          PARPRIORS[11] = 65.0                             ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
          PARPRIORS[11] = 21.1491                            ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values 
                                                                                        # from Kattge et al., (2011)
                                                                                        # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                                        # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = OBS$MaxRootDepth                 ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                              ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          #PARPRIORS[28] = 0.87                             ; PARPRIORUNC[28] = 0.41 # Resilience factor
          #PARPRIORS[29] = 0.5                              ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
          #PARPRIORS[30] = 0.1                              ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
          PARPRIORS[31] = 0.01                             ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          #PARPRIORS[32] = 0.25                             ; PARPRIORUNC[32] = 0.25 # Foliage + root litter combustion completeness
          PARPRIORS[33] = 0.43                 ; PARPRIORUNC[33] = 0.05 # Canopy NIR reflectance
          PARPRIORS[34] = 0.16                 ; PARPRIORUNC[34] = 0.05 # Canopy PAR reflectance
          PARPRIORS[35] = 0.26                 ; PARPRIORUNC[35] = 0.05 # Canopy NIR transmittance
          PARPRIORS[36] = 0.16                 ; PARPRIORUNC[36] = 0.05 # Canopy PAR transmittance
          PARPRIORS[37] = 0.023                ; PARPRIORUNC[37] = 0.05 # Soil reflectance to near infrared radiation
          PARPRIORS[38] = 0.033                ; PARPRIORUNC[38] = 0.05 # Foliage + root litter combustion completeness
          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          OTHERPRIORS[6] = OBS$MTTsom          ; OTHERPRIORUNC[6] = OBS$MTTsom_unc # Prior on the MTT of soil
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }                    
      } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P5.021") {
          PARPRIORS[2] = 0.54                              ; PARPRIORUNC[2]  = 0.12 #; PARPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Global Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[10] = 0.0334798                        ; PARPRIORUNC[10] = 0.015 #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[11] = 1.89*14.77735                    ; PARPRIORUNC[11] = 1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                               # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                               # These observational constraints are not the same and would lead to
                                                                               # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
#          PARPRIORS[11] = 65.0                             ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
          PARPRIORS[11] = 21.1491                            ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values 
                                                                                        # from Kattge et al., (2011)
                                                                                        # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                                        # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = OBS$MaxRootDepth                 ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                              ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          #PARPRIORS[28] = 0.87                             ; PARPRIORUNC[28] = 0.41 # Resilience factor
          #PARPRIORS[29] = 0.5                              ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
          #PARPRIORS[30] = 0.1                              ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
          PARPRIORS[31] = 0.01                             ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          #PARPRIORS[32] = 0.25                             ; PARPRIORUNC[32] = 0.25 # Foliage + root litter combustion completeness
          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          OTHERPRIORS[6] = OBS$MTTsom          ; OTHERPRIORUNC[6] = OBS$MTTsom_unc # Prior on the MTT of soil
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }                    
      } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P6.022") {
          PARPRIORS[2]  = 0.54                 ; PARPRIORUNC[2] = 0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[10] = OBS$RhetQ10          ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[11] = 21.1491              ; PARPRIORUNC[11] = 8.534234 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                            #       Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                            # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca              ; PARPRIORUNC[17] = OBS$lca_unc
          PARPRIORS[19] = OBS$Cfol_initial     ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial   ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial    ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial     ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial     ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[27] = OBS$MaxRootDepth     ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
          #PARPRIORS[28] = 0.5                  ; PARPRIORUNC[28] = 0.25 # Resilience factor
          #PARPRIORS[29] = 0.5                  ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
          #PARPRIORS[30] = 0.1                  ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
          PARPRIORS[31] = 0.01                 ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          # Other priors
          #OTHERPRIORS[1] =       ; OTHERPRIORUNC[1] =  # Initial soil water fraction 
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }         
      } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P7.R2.023") {
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[11] = 0.2764618                        ; PARPRIORUNC[11] = 0.2014871 # log10 avg foliar N (gN.m-2)
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = OBS$MaxRootDepth                 ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                              ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          #PARPRIORS[28] = 0.87                             ; PARPRIORUNC[28] = 0.41 # Resilience factor
          #PARPRIORS[29] = 0.5                              ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
          #PARPRIORS[30] = 0.1                              ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
          PARPRIORS[31] = 0.01                             ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          PARPRIORS[32] = 0.25                             ; PARPRIORUNC[32] = 0.25 # Foliage + root litter combustion completeness
          PARPRIORS[35] = 30                               ; PARPRIORUNC[35] = 2.5   # CGI optimum temperature curve
          PARPRIORS[46] = 0                                ; PARPRIORUNC[46] = 1.0   # Absolute min temperature for leaf growth
          PARPRIORS[41] = 11.197440                        ; PARPRIORUNC[41] = 9.3   # NUE prior derived from Kattge et al., (2011), based on log10 gaussian distribution
          PARPRIORS[42] = 1.639                            ; PARPRIORUNC[42] = 0.125 # Rm_leaf N**exponent (gC/gN) Reich et al., (2008)
          PARPRIORS[43] = 0.778                            ; PARPRIORUNC[43] = 0.133 # Rm_leaf intercept (gC/gN) Reich et al., (2008)          
          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[2] = 0.54                            ; OTHERPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#NOTUSE          OTHERPRIORS[3] = 27.295        ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
          OTHERPRIORS[4] = 0.66                            ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential             ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }                    
      } else if (modelname == "DALEC.A1.C1.D2.F2.H4.P1.024") {
          PARPRIORS[2] = 0.54                              ; PARPRIORUNC[2]  = 0.12 #; PARPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Global Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[10] = 0.0334798                        ; PARPRIORUNC[10] = 0.015 #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[11] = 1.89*14.77735                    ; PARPRIORUNC[11] = 1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                               # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                               # These observational constraints are not the same and would lead to
                                                                               # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
#          PARPRIORS[11] = 65.0                             ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
          PARPRIORS[11] = 21.1491                            ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values 
                                                                                        # from Kattge et al., (2011)
                                                                                        # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                                        # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = OBS$MaxRootDepth                 ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                              ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          #PARPRIORS[28] = 0.87                             ; PARPRIORUNC[28] = 0.41 # Resilience factor
          #PARPRIORS[29] = 0.5                              ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
          #PARPRIORS[30] = 0.1                              ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
          PARPRIORS[31] = 0.01                             ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          #PARPRIORS[32] = 0.25                             ; PARPRIORUNC[32] = 0.25 # Foliage + root litter combustion completeness



          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          OTHERPRIORS[6] = OBS$MTTsom          ; OTHERPRIORUNC[6] = OBS$MTTsom_unc # Prior on the MTT of soil
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }                                       
      } else if (modelname == "DALEC...025") {

#NOT IN USE

      } else if (modelname == "DALEC.A4.C6.D2.F2.H3.P10.026") {
          PARPRIORS[10] = OBS$RhetQ10            ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[11] = 65.0               ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
#          PARPRIORS[11] = 60.0               ; PARPRIORUNC[11]= 20.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax: derived from multiple trait values from Kattge et al., (2011)
                                                                          # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                          # due to the different temperature response functions used in ACM2 vs ACM 1
#          PARPRIORS[11] = 24.43                  ; PARPRIORUNC[11]= 10.4 #; PARPRIORWEIGHT[11] = 1 # Vcmax: Median of CARDAMOM analysis assimilating 4 GPP products
          PARPRIORS[11] = 54.165                  ; PARPRIORUNC[11]= 20.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax: Mean of reported PFT values from Oliver et al., (2022)
                                                                          # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                          # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[12] = -2.0                   ; PARPRIORUNC[12] = 0.5 # minimum leaf water potential (MPa)
          PARPRIORS[17] = OBS$lca                ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial       ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial     ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial      ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial       ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial       ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[27] = OBS$MaxRootDepth       ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                    ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, based on median from Fan et al., (2017) https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          PARPRIORS[31] = 0.01                   ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
#          PARPRIORS[33] = 0.01                 ; PARPRIORUNC[33] = 0.05 # labile:biomass at which growth limited by 50 %
          PARPRIORS[36] = 0.0                 ; PARPRIORUNC[36] = 5.0 # temperature at which foliage and root growth totally suppressed (oC)
          PARPRIORS[37] = 5.0                 ; PARPRIORUNC[37] = 1.0 # temperature at which wood growth totally suppressed (oC)
          # Other priors
          #OTHERPRIORS[1] =       ; OTHERPRIORUNC[1] =  # Initial soil water fraction 
          OTHERPRIORS[2] = 0.54                ; OTHERPRIORUNC[2] = 0.12 #; OTHERPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
              #if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4] = -9999 ; OBSMAT[filter,5] = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
              #}
          }                    
#          PARPRIORS[11] = 0.2764618            ; PARPRIORUNC[11] = 0.2014871 # log10 avg foliar N (gN.m-2)
#          PARPRIORS[42] = 11.197440            ; PARPRIORUNC[42] = 9.3  # NUE prior derived from Kattge et al., (2011), based on log10 gaussian distribution
#          PARPRIORS[43] = 275.1452             ; PARPRIORUNC[43] = 296.2767 # Leaf lifespan prior form Kattge et al., 2011, based on log10 gauusian distribution
#          OTHERPRIORS[3] = 27.295              ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
      } else if (modelname == "DALEC.A3.H1.M2.016") {
          # Override the default deforestation fraction forcing with the grassland LAI change variable.
          # Note the positive values, i.e. implied growth, will be ignored by the model
          MET[,8]  = OBS$lai_change  # m2/m2

          PARPRIORS[2]  = 0.54                 ; PARPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[11] = 21.1491              ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                            # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                            # due to the different temperature response functions used in ACM2 vs ACM 1
          #PARPRIORS[15] = OBS$lca              ; PARPRIORUNC[15] = OBS$lca_unc
          PARPRIORS[15] = 32                   ; PARPRIORUNC[15]=13 # Grass prior for UK purpose
          PARPRIORS[17] = OBS$Cfol_initial     ; PARPRIORUNC[17] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[18] = OBS$Croots_initial   ; PARPRIORUNC[18] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[19] = OBS$Clit_initial     ; PARPRIORUNC[19] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial     ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$MaxRootDepth     ; PARPRIORUNC[25] = OBS$MaxRootDepth_unc # Maximum rooting depth prior,           
          # other priors
#          OTHERPRIORS[2] =        ; OTHERPRIORUNC[2] =  # Initial soil water fraction 
#          OTHERPRIORS[3] = 27.295              ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
           OTHERPRIORS[4] = 0.66               ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      } else if (modelname == "DALEC.A3.H2.M2.017") {
          # Override the default deforestation fraction forcing with the grassland LAI change variable.
          # Note the positive values, i.e. implied growth, will be ignored by the model
          MET[,8]  = OBS$lai_change  # m2/m2

          PARPRIORS[2]  = 0.54                 ; PARPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[11] = 21.1491              ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                            # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                            # due to the different temperature response functions used in ACM2 vs ACM 1
          #PARPRIORS[15] = OBS$lca              ; PARPRIORUNC[15] = OBS$lca_unc
          PARPRIORS[15] = 32                   ; PARPRIORUNC[15]=13 # Grass prior for UK purpose
          PARPRIORS[17] = OBS$Cfol_initial     ; PARPRIORUNC[17] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[18] = OBS$Croots_initial   ; PARPRIORUNC[18] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[19] = OBS$Clit_initial     ; PARPRIORUNC[19] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial     ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$MaxRootDepth     ; PARPRIORUNC[25] = OBS$MaxRootDepth_unc # Maximum rooting depth prior,           
          # other priors
#          OTHERPRIORS[2] =        ; OTHERPRIORUNC[2] =  # Initial soil water fraction 
#          OTHERPRIORS[3] = 27.295              ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
           OTHERPRIORS[4] = 0.66               ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      } else if (modelname == "DALEC.A1.C8.D2.F2.H2.P1.R4.036") {
          PARPRIORS[1] = 0.5                               ; PARPRIORUNC[1]  = 0.1 #; PARPRIORWEIGHT[10] = 1 # Efficiency of Litter decomposition to soil (0-1)
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
          PARPRIORS[11] = 21.1491                          ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values 
                                                                                        # from Kattge et al., (2011)
                                                                                        # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                                        # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = OBS$MaxRootDepth                 ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                              ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          #PARPRIORS[28] = 0.87                             ; PARPRIORUNC[28] = 0.41 # Resilience factor
          #PARPRIORS[29] = 0.5                              ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
          #PARPRIORS[30] = 0.1                              ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
          PARPRIORS[31] = 0.01                             ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          #PARPRIORS[32] = 0.25                             ; PARPRIORUNC[32] = 0.25 # Foliage + root litter combustion completeness
          #PARPRIORS[33] = 0.25                             ; PARPRIORUNC[33] = 0.25 # wood litter combustion completeness
          PARPRIORS[41] = 0.13                             ; PARPRIORUNC[41] = 0.03 # Foliar lignin fraction (Mediavilla et al., (2008), Functional Ecology, 22: 787-793. https://doi.org/10.1111/j.1365-2435.2008.01453.x)
          PARPRIORS[42] = 0.20                             ; PARPRIORUNC[42] = 0.08 # Fine root lignin fraction (Rough number from grasslands)
          PARPRIORS[43] = 0.28                             ; PARPRIORUNC[43] = 0.15 # Wood lignin fraction (Cornwell et al., 2009, doi: https://doi.org/10.1111/j.1365-2486.2009.01916.x)
          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[2] = 0.54                ; OTHERPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#NOTUSE          OTHERPRIORS[3] = 27.295        ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }           
      } else if (modelname == "DALEC.A1.C1.D2.F2.H5.P1.037") {
          PARPRIORS[2] = 0.54                              ; PARPRIORUNC[2]  = 0.12 #; PARPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
          PARPRIORS[10] = OBS$RhetQ10                      ; PARPRIORUNC[10] = OBS$RhetQ10_unc #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[10] = 0.0334798                        ; PARPRIORUNC[10] = 0.015 #; PARPRIORWEIGHT[10] = 1 # Heterotrophic exponential temperature response (Q10 = 1.4, Hashimoto et al., 2015; doi:10.5194/bg-12-4121-2015)
#          PARPRIORS[11] = 1.89*14.77735                    ; PARPRIORUNC[11] = 1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                               # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                               # These observational constraints are not the same and would lead to
                                                                               # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
#          PARPRIORS[11] = 65.0                             ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
          PARPRIORS[11] = 21.1491                            ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values 
                                                                                        # from Kattge et al., (2011)
                                                                                        # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                                        # due to the different temperature response functions used in ACM2 vs ACM 1
          PARPRIORS[17] = OBS$lca                          ; PARPRIORUNC[17] = OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
          PARPRIORS[19] = OBS$Cfol_initial                 ; PARPRIORUNC[19] = OBS$Cfol_initial_unc # Cfoliar prior
          PARPRIORS[20] = OBS$Croots_initial               ; PARPRIORUNC[20] = OBS$Croots_initial_unc # Croots prior
          PARPRIORS[21] = OBS$Cwood_initial                ; PARPRIORUNC[21] = OBS$Cwood_initial_unc # Cwood prior
          PARPRIORS[22] = OBS$Clit_initial                 ; PARPRIORUNC[22] = OBS$Clit_initial_unc # Clitter prior
          PARPRIORS[23] = OBS$Csom_initial                 ; PARPRIORUNC[23] = OBS$Csom_initial_unc # Csom prior
          PARPRIORS[25] = OBS$frac_Cwood_coarse_root_prior ; PARPRIORUNC[25] = OBS$frac_Cwood_coarse_root_prior_unc # Fraction Cwood coarse root prior
          PARPRIORS[27] = OBS$MaxRootDepth                 ; PARPRIORUNC[27] = OBS$MaxRootDepth_unc # Maximum rooting depth prior, 
#          PARPRIORS[27] = 1.0                              ; PARPRIORUNC[27] = 0.5 # Maximum rooting depth prior, 
#                                                                                   # based on median from Fan et al., (2017) 
#                                                                                   # https://www.pnas.org/doi/epdf/10.1073/pnas.1712381114
          #PARPRIORS[28] = 0.87                             ; PARPRIORUNC[28] = 0.41 # Resilience factor
          #PARPRIORS[29] = 0.5                              ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
          #PARPRIORS[30] = 0.1                              ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
          PARPRIORS[31] = 0.01                             ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
          #PARPRIORS[32] = 0.25                             ; PARPRIORUNC[32] = 0.25 # Foliage + root litter combustion completeness
          PARPRIORS[33] = -2.0                             ; PARPRIORUNC[33] = 0.5 # minimum leaf water potential (MPa)          
          # Other priors
          #OTHERPRIORS[1] =        ; OTHERPRIORUNC[1] =  # Initial soil water fraction (GLEAM v3.1a)
          OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
          OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
          OTHERPRIORS[6] = OBS$MTTsom          ; OTHERPRIORUNC[6] = OBS$MTTsom_unc # Prior on the MTT of soil
          OTHERPRIORS[7] = 0.7                 ; OTHERPRIORUNC[7] = 0.1 ; OTHERPRIORWEIGHT[7] = noyears # Prior on the unstressed CiCa ratio
          # Hack to remove LAI observations out of growing season for high LCA areas
          if (PARPRIORS[17] > 100) {
#              if (lat_degrees > 50) {
                  filter = which(MET[,6] < 175 | MET[,6] > 250)
                  OBSMAT[filter,4]  = -9999 ; OBSMAT[filter,5]  = -9999 # Filter LAI estimates
                  OBSMAT[filter,34] = -9999 ; OBSMAT[filter,35] = -9999 # Filter fAPAR estimates
#              }
          }                                       
      } else if (modelname == "ACM") {

          # For ACM_GPP_ET
          # p(1) = nitrogen use efficiency at optimum temperature (oC)
          #,unlimited by CO2, light and photoperiod (34gC/gN)
#          PARPRIORS[1] = 35.65 ; PARPRIORUNC[1] = 35.65 * 0.1 #26.19
          # p(2) = maximum temperature at which photosynthesis occurs (oC) (57.05oC)
#          PARPRIORS[2] = 57.05 ; PARPRIORUNC[2] = 57.05 * 0.1
          # p(3) = optimum temperature for photosynthesis (oC) (30oC)
#          PARPRIORS[3] = 30.0 ; PARPRIORUNC[3] = 30.0 * 0.1
          # p(4) = kurtosis for temperature response of photosynthesis (0.185912)
          PARPRIORS[4] = 0.185912 ; PARPRIORUNC[4] = 0.185912 * 0.01
          # p(5) = SPA apparent quantium yield gC/MJ PAR = ~3.2
#          PARPRIORS[5] = 3.2 ; PARPRIORUNC[5] = 0.1
          # p(6) = min leaf water potential (MPa)
#          PARPRIORS[6] = -2.0 ; PARPRIORUNC[6] = abs(-2.0*0.1)
          # p(7) = Coefficient linking soil isothermal->net adjustment and LAI
          #        SPA based prior is -2.7108547 W/m2 (SE +/- 0.0222038)
                   PARPRIORS[7] = -2.7108547 ; PARPRIORUNC[7] = 0.0222038*1.98
          # p(8) = iWUE (gC/m2leaf/dayl/mmolH2Ogs/s)
          # Actual value used in SPA is 8.4e-8 (Williams et al., 1996)
          # Other reported values are 9.0e-8 -> 1.8e-7 (Bonan et al., 2014)
          # NOTE: As this is applied at daily time step and the
          #       hourly time step activities in SPA operate across
          #       a non-linear diurnal cycle the true equivalent value is effectively unknown.
          #PARPRIORS[8] = 8.4e-8 ; PARPRIORUNC[8] = 0.1
          # p(9) = Soil SW absorption (fraction)
          PARPRIORS[9] = 0.972 ; PARPRIORUNC[9] = 0.972 * 0.1
          # p(10) = Max fractional reduction of longwave release from canopy
          #         Prior from offline SPA calibration = 0.9517081 +/- 0.0001011 SE
          PARPRIORS[10] = 0.9517081 ; PARPRIORUNC[10] = 0.0001011*1.98
          # p(11) = LAI adjustment for long wave release from canopy
          #         Prior from offline SPA calibration = 4.6917871 +/- 0.0013296 SE
          PARPRIORS[11] = 4.6917871 ; PARPRIORUNC[11] = 0.0013296*1.98
          # p(12) = Coefficient linking soil isothermal->net adjustment and absorbed SW
          #         SPA based prior is -0.0357603 W/m2 (SE +/- 0.0004877)
          PARPRIORS[12] = -0.0357603 ; PARPRIORUNC[12] = 0.0004877*1.98
          # p(13) = Constant relating soil isothermal->net adjustment soil radiation (W/m2)
          #         SPA based prior is 3.4806352 W/m2 (SE +/- 0.0636849)
          PARPRIORS[13] = 3.4806352 ; PARPRIORUNC[13] = 0.0636849*1.98
          # p(14) = Canopy PAR transmittance (fraction)
          PARPRIORS[14] = 0.16 ; PARPRIORUNC[14] = 0.16 * 0.01
          # p(15) = Canopy NIR transmittance (fraction)
          PARPRIORS[15] = 0.26 ; PARPRIORUNC[15] = 0.26 * 0.01
          # p(16) = Canopy PAR reflectance (fraction)
          PARPRIORS[16] = 0.16 ; PARPRIORUNC[16] = 0.16 * 0.01
          # p(17) = Canopy NIR reflectance (fraction)
          PARPRIORS[17] = 0.43 ; PARPRIORUNC[17] = 0.43 * 0.01
          # p(18) = Coefficient linking canopy isothermal->net adjustment and absorbed SW
          #         SPA based prior is 0.0154225 (SE +/- 0.0005798)
          PARPRIORS[18] = 0.0154225 ; PARPRIORUNC[18] = 0.0005798*1.98
          # p(19) = Constant relating canopy isothermal->net adjustment canopy radiation (W/m2)
          #         SPA based prior is 0.0577857 W/m2 (SE +/- 0.0217731)
          PARPRIORS[19] = 0.0577857 ; PARPRIORUNC[19] = 0.0217731*1.98
          # p(20) = Coefficient linking canopy isothermal->net adjustment and LAI
          #         SPA based prior is 2.4526437 (SE +/- 0.0229691)
          PARPRIORS[20] = 2.4526437 ; PARPRIORUNC[20] = 0.0229691*1.98
      } # model specific priora

      # Combine the forcings and observations matrices together
      DATA_TEMP = t(cbind(MET,OBSMAT))
      # combine the static data
      DATA_STAT = c(PARPRIORS,PARPRIORUNC,PARPRIORWEIGHT,OTHERPRIORS,OTHERPRIORUNC,OTHERPRIORWEIGHT)

      # open the binary file
      zz <- file(file, "wb")
      # write with 8 bit precision (i.e. double)
      writeBin(as.double(static_data), zz)
      writeBin(as.double(DATA_STAT), zz)
      writeBin(as.double(as.vector(DATA_TEMP)), zz)
      close(zz)

  } # forcings pass condition

} # end function binary_daya

## Use byte compile
binary_data<-cmpfun(binary_data)
