
###
## Function to create binary input files
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE) & J. F. Exbrayat (UoE).

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
                      "Leaf area index (m2/m2)",
                      "Leaf area index variance",
                      "Net Ecosystem Exchange of CO2 (gC/m2/day)",
                      "Net Ecosystem Exchange of CO2 variance",
                      "Fire C emission (gC/m2/day)",
                      "Fire C emission variance",
                      "Ecosystem respiration (Ra + Rh gC/m2/day)",
                      "Ecosystem respiration (Ra + Rh) variance",
                      "Foliar stock (gC/m2)",
                      "Foliar stock variance",
                      "Wood stock (above + below gC/m2)",
                      "Wood stock (above + below) variance",
                      "Fine root stock (gC/m2)",
                      "Fine root stock variance",
                      "Foliar + fine root litter stock (gC/m2)",
                      "Foliar + fine root litter stock variance",
                      "Soil organic matter stock (gC/m2)",
                      "Soil organic matter stock variance",
                      "Above ground biomass stock (gC/m2)",
                      "Above ground biomass stock variance",
                      "Fraction absorbed PAR",
                      "Fraction absorbed PAR variance",
                      "Empty",
                      "Empty",
                      "Coarse root stock (gC/m2)",
                      "Coarse root stock variance",
                      "Annual foliar maximum (gC/m2)",
                      "Annual foliar maximum variance",
                      "Evapotranspiration (kgH2O/m2/day)",
                      "Evapotranspiration variance",
                      "Snow water equivalent (kgH2O/m2)",
                      "Snow water equivalent variance",
                      "Net Biome Exchange (Reco + Fire - GPP) of CO2 (gC/m2/day)",
                      "Net Biome Exchange variance",
                      "Mean woody productivity over lag period (gC/m2/day)",
                      "Mean woody productivity varince",
                      "Lag period over which to average (steps)",
                      "Mean woody natural mortality over lag period (gC/m2/day)",
                      "Mean woody natural mortality varince",
                      "Lag period over which to average (steps)",
                      "Mean litter flux over lag period (gC/m2/day)",
                      "Mean litter flux varince",
                      "Lag period over which to average (steps)",
                      "Mean woody net increment over lag period (gC/m2/day)",
                      "Mean woody net increment varince",
                      "Lag period over which to average (steps)",
                      "Extracted C due to harvest over lag period (gC/m2/day)",
                      "Extracted C due to harvest varince",
                      "Lag period over which to average (steps)")

binary_data<-function(met,OBS,file,EDC,latlon_in,ctessel_pft,modelname,parameter_type,nopars,noyears) {

  # Inform the user
  if (use_parallel == FALSE) {print(paste("writing out binary...",Sys.time(),sep=""))}

  # set model ID
  if (modelname == "ACM") {
    modelid = 0
  } else if (modelname == "DALEC.D1.F2.#") {
    modelid = 1
  } else if (modelname == "DALEC.C1.D1.F2.P1.#"){
    modelid = 2
  } else if (modelname == "DALEC.A1.C1.D2.F2.H1.P1.#"){
    modelid = 3
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P1.#"){
    modelid = 4
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P1.R1.#") {
    modelid = 5
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P1.R1.#") {
    modelid = 6
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P2.R1.#") {
    modelid = 7
  } else if (modelname == "DALEC.A1.C2.D2.F2.H1.P3.R1.#") {
    modelid = 8
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P3.R1.#") {
    modelid = 9
  } else if (modelname == "DALEC.A1.C2.D2.F2.H1.P4.R2.#") {
    modelid = 10
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P4.R2.#") {
    modelid = 11
  } else if (modelname == "DALEC.C4.D1.F2.#") {
    modelid = 12
  } else if (modelname == "DALEC.C5.D1.F2.P1.#") {
    modelid = 13
  } else if (modelname == "DALEC.C3.M1.#") {
    modelid = 14
  } else if (modelname == "DALEC.A3.C3.H2.M1.#") {
    modelid = 15
  } else if (modelname == "DALEC.M2.#") {
    modelid = 16
  } else if (modelname == "DALEC.A3.H2.M2.#"){
    modelid = 17
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P2.#"){
    modelid = 18
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P2.R3.#"){
    modelid = 19
  } else if (modelname == "DALEC.A2.C1.D2.F2.H2.P1.#"){
    modelid = 20
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P5.#") {
    modelid = 21
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P6.#") {
    modelid = 22
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P7.R2.#") {
    modelid = 23
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P8.R2.#") {
    modelid = 24
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P9.R2.#") {
    modelid = 25
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P10.R2.#") {
    modelid = 26
  } else if (modelname == "DALEC_1005") {
    modelid = 27
  } else if (modelname == "DALEC_1005a") {
    modelid = 28
  } else if (modelname == "DALEC.A1.C1.D2.F2.H3.P1.#") {
    modelid = 29
  } else if (modelname == "DALEC.A3.C1.D2.F2.H2.P1.#") {
    modelid = 30
  } else if (modelname == "DALEC.A4.C6.D2.F2.H2.P11.#"){
    modelid = 31
  } else if (modelname == ""){
    modelid = 32
  } else if (modelname == ""){
    modelid = 33
  } else if (modelname == ""){
    modelid = 34
  } else if (modelname == ""){
    modelid = 35
  } else if (modelname == "") {
    modelid = 36
  } else if (modelname == "") {
    modelid = 37
  } else if (modelname == "") {
    modelid = 38
  } else if (modelname == "") {
    modelid = 39
  } # model ID assignment

  # some drivers may be passed as single values assuming this will apply across the whole time series
  # if this is so we need to now adjust that
  if (length(OBS$forest_management) != length(met$run_day)) {
      # we will assume that this is a constant value and will now repeast it
      OBS$forest_management=array(OBS$forest_management, dim=c(length(met$run_day)))
  }

  steps_per_year = floor(length(met$run_day)[1] / noyears)

  # extract information from list to array
  if (modelname == "ACM") {
      MET = array(NA,dim=c(length(met$run_day),(length(met)+2)))
  } else {
      MET = array(NA,dim=c(length(met$run_day),(length(met)+3)))
  }

  MET[,1]  = met$run_day
  MET[,2]  = met$mint  ; if (min(met$mint) < -200) {print(summary(met$mint)) ; stop('mint error in binary_data')} # Celcius
  MET[,3]  = met$maxt  ; if (min(met$maxt) < -200) {print(summary(met$maxt)) ; stop('maxt error in binary_data')} # Celcius
  MET[,4]  = met$swrad ; if (min(met$swrad) < 0 | max(met$swrad) > 36) {print(summary(met$swrad)) ; stop('RAD error in binary_data')} # MJ/m2/day
  MET[,5]  = met$co2#+200 # ppm
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
  }

  # Create time series observation matrix, i.e. things with uncertainty associated.
  # Currently space for 18 time series of observation and its uncertainty.
  # Uncertainty is assumed to be the Gaussian variance in same units as the observation itself.
  # NOTE: that not all models are currently coded to be compatible with all observation streams.
  OBSMAT = array(-9999.0,dim=c(length(met$run_day),51))
  # Line makes the correct array size but with -9999 in place of all
  OBSMAT[,1] = OBS$GPP                    # GPP (gC/m2/day)
  OBSMAT[,2] = OBS$GPP_unc                # GPP variance (gC/m2/day)
  OBSMAT[,3] = OBS$LAI                    # Leaf area index (m2/m2)
  OBSMAT[,4] = OBS$LAI_unc                # Leaf area index variance
  OBSMAT[,5] = OBS$NEE                    # Net Ecosystem Exchange of CO2 (gC/m2/day)
  OBSMAT[,6] = OBS$NEE_unc                # Net Ecosystem Exchange of CO2 variance
  OBSMAT[,7] = OBS$Fire                   # Fire C emission (gC/m2day)
  OBSMAT[,8] = OBS$Fire_unc               # Fire C emission variance
  OBSMAT[,9] = OBS$Reco                   # Ecosystem respiration (Ra + Rh; gC/m2/day)
  OBSMAT[,10] = OBS$Reco_unc              # Ecosystem respiration (Ra + Rh) variance
  OBSMAT[,11] = OBS$Cfol_stock            # Foliar stock (gC/m2)
  OBSMAT[,12] = OBS$Cfol_stock_unc        # Foliar stock variance
  OBSMAT[,13] = OBS$Cwood_stock           # Wood stock (above + below; gC/m2)
  #OBSMAT[which(OBSMAT[,13] > 0),14] = 250
  OBSMAT[,14] = OBS$Cwood_stock_unc       # Wood stock (above + below) variance
  OBSMAT[,15] = OBS$Croots_stock          # Fine root stock (gC/m2)
  OBSMAT[,16] = OBS$Croots_stock_unc      # Fine root stock variance
  OBSMAT[,17] = OBS$Clit_stock            # Foliar + fine root litter stock (gC/m2)
  OBSMAT[,18] = OBS$Clit_stock_unc        # Foliar + fine root litter stock variance
  OBSMAT[,19] = OBS$Csom_stock            # Soil organic matter stock (gC/m2)
  OBSMAT[,20] = OBS$Csom_stock_unc        # Soil organic matter stock variance
  OBSMAT[,21] = OBS$Cagb_stock            # Above ground biomass stock (gC/m2)
  OBSMAT[,22] = OBS$Cagb_stock_unc        # Above ground biomass stock variance
  OBSMAT[,23] = OBS$fAPAR                 # Fraction absorbed PAR
  OBSMAT[,24] = OBS$fAPAR_unc             # Fraction absorbed PAR varianceq
  OBSMAT[,25] = -9999                     # Empty
  OBSMAT[,26] = -9999                     # Empty
  OBSMAT[,27] = OBS$Ccoarseroot_stock     # Coarse root stock (gC/m2)
  OBSMAT[,28] = OBS$Ccoarseroot_stock_unc # Coarse root stock variance
  OBSMAT[,29] = OBS$Cfolmax_stock         # Annual foliar maximum (gC/m2)
  OBSMAT[,30] = OBS$Cfolmax_stock_unc     # Annual foliar maximum variance
  OBSMAT[,31] = OBS$Evap                  # Evapotranspiration (kgH2O/m2/day)
  OBSMAT[,32] = OBS$Evap_unc              # Evapotranspiration variance
  OBSMAT[,33] = OBS$SWE                   # Snow water equivalent (kgH2O/m2)
  OBSMAT[,34] = OBS$SWE_unc               # Snow water equivalent variance
  OBSMAT[,35] = OBS$nbe                   # Net Biome Exchange (Reco + Fire - GPP) of CO2 (gC/m2/day)
  OBSMAT[,36] = OBS$nbe_unc               # Net Biome Exchange variance
  OBSMAT[,37] = OBS$Cwood_inc             # Mean woody productivity over lag period (gC/m2/day)
  OBSMAT[,38] = OBS$Cwood_inc_unc         # Mean woody productivity varince
  OBSMAT[,39] = OBS$Cwood_inc_lag         # Lag period over which to average  (steps)
  OBSMAT[,40] = OBS$Cwood_mortality       # Mean woody natural mortality over lag period (gC/m2/day)
  OBSMAT[,41] = OBS$Cwood_mortality_unc   # Mean woody natural mortality varince
  OBSMAT[,42] = OBS$Cwood_mortality_lag   # Lag period over which to average  (steps)
  OBSMAT[,43] = OBS$foliage_to_litter     # Mean litter flux over lag period (gC/m2/day)
  OBSMAT[,44] = OBS$foliage_to_litter_unc # Mean litter flux varince
  OBSMAT[,45] = OBS$foliage_to_litter_lag # Lag period over which to average (steps)
  #OBSMAT[,46] = OBS$Cwood_net_inc         # Mean woody net increment over lag period (gC/m2/day)
  #OBSMAT[,47] = OBS$Cwood_net_inc_unc     # Mean woody net increment varince
  #OBSMAT[,48] = OBS$Cwood_net_inc_lag     # Lag period over which to average  (steps)
  OBSMAT[,49] = OBS$harvest               # Extracted C due to harvest over lag period (gC/m2/day)
  OBSMAT[,50] = OBS$harvest_unc           # Extracted C due to harvest varince
  OBSMAT[,51] = OBS$harvest_lag           # Lag period over which to average  (steps)
  
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
  tmp = c(modelid,latlon_in[1],dim(MET)[1],dim(MET)[2],dim(OBSMAT)[2],
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
  #PARPRIORS[11]=0.2764618	; PARPRIORUNC[11]=0.2014871 # log10 avg foliar N (gN.m-2)
  # PARPRIORS[2]=51.70631   ; PARPRIORUNC[2]=124.6938 # C:N root (gC/gN) Kattge et al., (2011)
  # PARPRIORS[15]=416.6667  ; PARPRIORUNC[15]=326.3762 # C:N wood (gC/gN) Kattge et al., (2011)
  # PARPRIORS[26]=11.197440 ; PARPRIORUNC[26]=1.32313*3 # NUE prior derived from Kattge et al., (2011)
  # PARPRIORS[36]=1.639     ; PARPRIORUNC[36]=0.125   # Rm_leaf N**exponent (gC/gN) Reich et al., (2008)
  # PARPRIORS[37]=0.778     ; PARPRIORUNC[37]=0.133   # Rm_leaf intercept (gC/gN) Reich et al., (2008)
  # PARPRIORS[38]=1.352     ; PARPRIORUNC[38]=0.150   # Rm_root N**exponent (gC/gN) Reich et al., (2008)
  # PARPRIORS[39]=0.997     ; PARPRIORUNC[39]=0.082   # Rm_root intercept (gC/gN) Reich et al., (2008)
  # PARPRIORS[40]=1.344     ; PARPRIORUNC[40]=0.150   # Rm_wood N**exponent (gC/gN) Reich et al., (2008)
  # PARPRIORS[41]=0.946     ; PARPRIORUNC[41]=0.107   # Rm_wood intercept (gC/gN) Reich et al., (2008)
  # OTHERPRIORS[1] = 0.46   ; OTHERPRIORUNC[1] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
  ###

  # Assign model specific parameter priors
  if (modelname == "DALEC.C1.D1.F2.P1.#") {
      PARPRIORS[2] =0.54                   ; PARPRIORUNC[2]=0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      PARPRIORS[11]=16.9                   ; PARPRIORUNC[11]=7.502147 # Ceff: derived from multiple trait values from Kattge et al., (2011)
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_un
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      #PARPRIORS[24] = 0.5                  ; PARPRIORUNC[24] = 0.25 # Resilience factor
      #PARPRIORS[25] = 0.5                  ; PARPRIORUNC[25] = 0.25 # Foliar combustion completeness
      #PARPRIORS[26] = 0.1                  ; PARPRIORUNC[26] = 0.25 # Root / wood combustion completeness
      PARPRIORS[27] = 0.01                 ; PARPRIORUNC[27] = 0.05 # Soil combustion completeness
      # Other priors
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.C3.M1.#") {
   
      # Crop model specific adjustment - disallow observational constraints during the first 12 months
      OBSMAT[1:steps_per_year,1:dim(OBSMAT)[2]] = -9999

      PARPRIORS[11]=11.197440              ; PARPRIORUNC[11]=9.3 # NUE prior derived from Kattge et al., (2011), based on log10 gaussian distribution      
      PARPRIORS[13]=0.21875                ; PARPRIORUNC[13]=0.01 # Respiratory costs of labile transfer
      PARPRIORS[12]=OBS$planting_doy       ; PARPRIORUNC[12]=OBS$planting_doy_unc
      PARPRIORS[15]=OBS$harvest_doy        ; PARPRIORUNC[15]=OBS$harvest_doy_unc
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      PARPRIORS[35]=0.99                   ; PARPRIORUNC[35]=0.1 # labile turnover rate (/day)
      PARPRIORS[36]=4.088                  ; PARPRIORUNC[36]=0.6052851 # Constant for canopy N dilution model (gN/m2leaf)
      PARPRIORS[37]=-0.0252                ; PARPRIORUNC[37]=0.00092   # Coefficient relating foliar C to N dilution
      # Other priors
      OTHERPRIORS[1] = 0.54                ; OTHERPRIORUNC[1]=0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      # Yield:GPP Winter Wheat ATEC experiment plus He et al., (2018), doi: 10.3390/rs10030372
      # Values from He et al., Spring Wheat 0.24, Barley 0.42, Duram Wheat 0.22, Alfalfa 0.55, Pea 0.28, Maize 0.44
      OTHERPRIORS[8] = 0.38                ; OTHERPRIORUNC[8]=0.087 ; OTHERPRIORWEIGHT[8] = noyears 
  } else if (modelname == "DALEC.A3.C3.H2.M1.#") {

      # Crop model specific adjustment - disallow observational constraints during the first 12 months
      OBSMAT[1:steps_per_year,1:dim(OBSMAT)[2]] = -9999

 #     # Parameter priors for Winter Wheat (yes something better needs to be done for the storing of these)
 #     # derived from the ATEC experiment field (but not assimilated) or Sus et al., (2010)
 #     PARPRIORS[2] = 0.44           ; PARPRIORUNC[2]  = 0.08         # Fraction of GPP allocated to autotrophic pool
 #     PARPRIORS[3] = 0.04           ; PARPRIORUNC[3]  = 0.02         # Development rate coefficient DS 0-1
 #     PARPRIORS[4] = 0.035          ; PARPRIORUNC[4]  = 0.02         # Development rate coefficient DS 1-2
 #     PARPRIORS[7] = 0.03           ; PARPRIORUNC[7]  = 0.03         # Potential turnover rate of foliage due to self-shading (fraction/day)
 #     PARPRIORS[8] = 22.5           ; PARPRIORUNC[8]  = 5.0          # No. of vernalisation days for plants to be 50 % vernalised
 #     PARPRIORS[11] = 21.1491       ; PARPRIORUNC[11] = 8.534234 #; PARPRIORWEIGHT[11] = noyears # Ceff: derived from multiple trait values from Kattge et al., (2011)
 #     PARPRIORS[13] = 125.0         ; PARPRIORUNC[13] = 10.0         # Phenological heat units for seed emergence (aka growing degree days)
 #     #PARPRIORS[12]=OBS$planting_doy       ; PARPRIORUNC[12]=OBS$planting_doy_unc # Sow day of year, applied as p12%%365.25
 #     #PARPRIORS[14]=OBS$harvest_doy        ; PARPRIORUNC[14]=OBS$harvest_doy_unc  # Growing season length sowing->harvest days
 #     if (max(OBS$LAI) > 0) {
 #         # Prior on canopy N derived from ATEC experiment assuming max LAI is related to canopy N
 #         PARPRIORS[15] = min(7.0,max(OBS$LAI) * 0.488 + 2.5131)
 #         PARPRIORUNC[15] = 1.2 # mean confidence interval of linear regression for LAI ranges 1-6
 #         #PARPRIORWEIGHT[15] = noyears
 #     } else {
 #         PARPRIORS[15] = 4.088        ; PARPRIORUNC[15] = 0.6052851    # Constant for canopy N dilution model (gN/m2leaf)
 #     }
 #     PARPRIORS[16] = -0.0252      ; PARPRIORUNC[16] = 0.00092      # Coefficient relating foliar C to N dilution
 #     PARPRIORS[17] = 20.40        ; PARPRIORUNC[17] = 5            # Leaf carbon per unit area (gC/m2)
 #                                                                   # Prior values from Penning de Vries (1982), Simulation of Ecophysiological Processes of Growth in Several Annual Crops
 #                                                                   # Barley = 15.60, Maize = 21.60, Potato = 14.40, Rice = 21.12, Sorghum = 19.20, Soybean = 19.20, SugarBeet = 24.00, 
 #                                                                   # Sugarcane = 33.60, Sunflower = 25.92, WinterWheat = 20.40 , SprintWheat = 24.00
 #     #PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
 #     #PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
 #     #PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
 #     #PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
 #     #PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
 #     PARPRIORS[23] = OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
 #     PARPRIORS[26] = 273.15       ; PARPRIORUNC[26] = 5.0          # min temperature for development (K)
 #     PARPRIORS[27] = 308.15       ; PARPRIORUNC[27] = 5.0          # max temperature for development (K)
 #     PARPRIORS[28] = 298.15       ; PARPRIORUNC[28] = 5.0          # optimum temperature for development (K)
#      PARPRIORS[29] = 271.85       ; PARPRIORUNC[29] = 5.0          # min temperature for vernalisation (K)
#      PARPRIORS[30] = 288.85       ; PARPRIORUNC[30] = 5.0          # max temperature for vernalisation (K)
#      PARPRIORS[31] = 278.05       ; PARPRIORUNC[31] = 5.0          # optimum temperature for vernalisation (K)
#      PARPRIORS[32] = 8.25         ; PARPRIORUNC[32] = 2.0          # critical photoperiod for development to begin (days)
#      PARPRIORS[33] = 0.25         ; PARPRIORUNC[33] = 0.2          # sensitivity to photoperiod above critical threshold
#      PARPRIORS[35] = 0.99         ; PARPRIORUNC[35] = 0.1          # autotrophic turnover rate (/day)
#      PARPRIORS[36] = 50.0         ; PARPRIORUNC[36] = 25.0         # fine root biomass (g/m2) needed to reach 50 % of maximum depth
#      PARPRIORS[37] = 1.5          ; PARPRIORUNC[37] = 0.25         # maximum rooting depth (m)

# ATEC LAI and N fertiliser addition is a highly uncertain saturating function.
# i.e. maxLAI = N addition (kgN/ha) * 0.011840 + 2.216000

      # Parameter priors for Winter Wheat (yes something better needs to be done for the storing of these)
      # derived from the ATEC experiment field (but not assimilated) or Sus et al., (2010)
      if (max(OBS$LAI) > 0) {
          # Prior on allocation of GPP to the autotrophic respiration pool
          # derived from ATEC experiment assuming max LAI is related to canopy N
          PARPRIORS[2] = min(4,max(OBS$LAI)) * -0.020653 + 0.478990       
          PARPRIORUNC[2] = 0.003040314*2 # mean confidence interval of linear regression for LAI ranges 1-6
      } else {
          PARPRIORS[2] = 0.405598979   ; PARPRIORUNC[2]  = 0.06497        # Fraction of GPP allocated to autotrophic pool
      }
      PARPRIORS[3] = 0.04           ; PARPRIORUNC[3]  = 0.02         # Development rate coefficient DS 0-1
      PARPRIORS[4] = 0.02           ; PARPRIORUNC[4]  = 0.02         # Development rate coefficient DS 1-
      if (max(OBS$LAI) > 0) {
          # Prior on canopy turnover due to DS derived from ATEC experiment assimilated into DALEC.15 assuming max LAI is related
          PARPRIORS[5] = min(4,max(OBS$LAI)) * -0.0019073 + 0.0155355
          PARPRIORUNC[5] = 0.001491601*2 # mean confidence interval of linear regression for LAI ranges 1-6
      } else {
          PARPRIORS[5] = 0.0086         ; PARPRIORUNC[5]  = 0.05     # Potential turnover rate of foliage due to DS (fraction/day)
      }
      PARPRIORS[7] = 0.038          ; PARPRIORUNC[7]  = 0.03         # Potential turnover rate of foliage due to self-shading (fraction/day)
      PARPRIORS[8] = 21.4108        ; PARPRIORUNC[8]  = 5.0          # No. of vernalisation days for plants to be 50 % vernalised
      PARPRIORS[11] = 19.410719948  ; PARPRIORUNC[11] = 7.049778  #; PARPRIORWEIGHT[11] = noyears # Ceff: derived from multiple trait values from Kattge et al., (2011)
      if (max(OBS$LAI) > 0) {
          # Prior on canopynitrogen use efficiency derived from ATEC experiment assimilated into DALEC.15 assuming max LAI is related to NUE
          PARPRIORS[11] = max(10,-9.9513 + 15.3898*min(4,max(OBS$LAI)) + -1.8213*min(4,max(OBS$LAI))**2)         
          PARPRIORUNC[11] = 1.008786*2 # mean confidence interval of linear regression for LAI ranges 1-6
      } else {
          PARPRIORS[11] = 19.410719948       ; PARPRIORUNC[11] = 7.049778  #; PARPRIORWEIGHT[11] = noyears # Ceff: derived from multiple trait values from Kattge et al., (2011)
      }
      PARPRIORS[13] = 125.0         ; PARPRIORUNC[13] = 10.0         # Phenological heat units for seed emergence (aka growing degree days)
      #PARPRIORS[12]=OBS$planting_doy       ; PARPRIORUNC[12]=OBS$planting_doy_unc # Sow day of year, applied as p12%%365.25
      #PARPRIORS[14]=OBS$harvest_doy        ; PARPRIORUNC[14]=OBS$harvest_doy_unc  # Growing season length sowing->harvest days
      PARPRIORS[14] = 225.0         ; PARPRIORUNC[14] = 32.0         # Growing season length sowing->harvest days
      if (max(OBS$LAI) > 0) {
          # Prior on canopy N derived from ATEC experiment assuming max LAI is related to canopy N
          #PARPRIORS[15] = min(7.0,max(OBS$LAI) * 0.488 + 2.5131)
          #PARPRIORUNC[15] = 0.6052851 # mean confidence interval of linear regression for LAI ranges 1-6
          # Prior on canopy N derived from ATEC experiment assimilated into DALEC.15 assuming max LAI is related to canopy N
          PARPRIORS[15] = min(7.0,max(OBS$LAI) * 0.61421 + 2.23757)          
          PARPRIORUNC[15] = 0.7182558*0.5 # mean confidence interval of linear regression for LAI ranges 1-6#
          #PARPRIORWEIGHT[15] = noyears
      } else {
          PARPRIORS[15] = 4.33        ; PARPRIORUNC[15] = 0.6052851    # Constant for canopy N dilution model (gN/m2leaf)
      }
      PARPRIORS[16] = -0.0252      ; PARPRIORUNC[16] = 0.00092      # Coefficient relating foliar C to N dilution
      PARPRIORS[17] = 20.40        ; PARPRIORUNC[17] = 3.8          # Leaf carbon per unit area (gC/m2)
                                                                    # Prior values from Penning de Vries (1982), Simulation of Ecophysiological Processes of Growth in Several Annual Crops
                                                                    # Barley = 15.60, Maize = 21.60, Potato = 14.40, Rice = 21.12, Sorghum = 19.20, Soybean = 19.20, SugarBeet = 24.00, 
                                                                    # Sugarcane = 33.60, Sunflower = 25.92, WinterWheat = 20.40 , SprintWheat = 24.00
      #PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      #PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      #PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      #PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      #PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23] = OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      PARPRIORS[26] = 274.07       ; PARPRIORUNC[26] = 1.4          # min temperature for development (K)
      PARPRIORS[27] = 306.35       ; PARPRIORUNC[27] = 2.3          # max temperature for development (K)
      PARPRIORS[28] = 291.55       ; PARPRIORUNC[28] = 2.5          # optimum temperature for development (K)
      PARPRIORS[29] = 270.46       ; PARPRIORUNC[29] = 1.6          # min temperature for vernalisation (K)
      PARPRIORS[30] = 288.73       ; PARPRIORUNC[30] = 2.2          # max temperature for vernalisation (K)
      PARPRIORS[31] = 277.66       ; PARPRIORUNC[31] = 1.6          # optimum temperature for vernalisation (K)
      PARPRIORS[32] = 7.85         ; PARPRIORUNC[32] = 1.28         # critical photoperiod for development to begin (days)
      PARPRIORS[33] = 0.22         ; PARPRIORUNC[33] = 0.2          # sensitivity to photoperiod above critical threshold
      PARPRIORS[34] = 0.005386     ; PARPRIORUNC[34] = 0.01         # labile turnover rate (fraction/day)
      PARPRIORS[35] = 0.99         ; PARPRIORUNC[35] = 0.1          # autotrophic turnover rate (/day)
      PARPRIORS[36] = 50.0         ; PARPRIORUNC[36] = 25.0         # fine root biomass (g/m2) needed to reach 50 % of maximum depth
      PARPRIORS[37] = 1.5          ; PARPRIORUNC[37] = 0.25         # maximum rooting depth (m)
      # Other priors
      OTHERPRIORS[1] = 0.54        ; OTHERPRIORUNC[1] = 0.12 #; OTHERPRIORWEIGHT[1] = noyears  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      # Yield:GPP Winter Wheat ATEC experiment plus He et al., (2018), doi: 10.3390/rs10030372
      # Values from He et al., Spring Wheat 0.24, Barley 0.42, Duram Wheat 0.22, Alfalfa 0.55, Pea 0.28, Maize 0.44
      OTHERPRIORS[8] = 0.38        ; OTHERPRIORUNC[8] = 0.087 #; OTHERPRIORWEIGHT[8] = noyears
  } else if (modelname == "DALEC_1005") {
      PARPRIORS[2] =0.54                   ; PARPRIORUNC[2]=0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      PARPRIORS[11]=16.9                   ; PARPRIORUNC[11]=7.502147 # Ceff: derived from multiple trait values from Kattge et al., (2011)
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      # Other priors
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC_1005a") {
      PARPRIORS[2] =0.54                   ; PARPRIORUNC[2]=0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      PARPRIORS[11]=16.9                   ; PARPRIORUNC[11]=7.502147 # Ceff: derived from multiple trait values from Kattge et al., (2011)
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      # Other priors
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A1.C1.D2.F2.H1.P1.#") {
      PARPRIORS[2] =0.54                   ; PARPRIORUNC[2]=0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#      PARPRIORS[11]=1.89*14.77735          ; PARPRIORUNC[11]=1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                             # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                             # These observational constraints are not the same and would lead to
                                                                             # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
      PARPRIORS[11]=21.1491                ; PARPRIORUNC[11]=8.534234 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                      # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                      # due to the different temperature response functions used in ACM2 vs ACM 1
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      #PARPRIORS[24] = 0.5                  ; PARPRIORUNC[24] = 0.25 # Resilience factor
      #PARPRIORS[25] = 0.5                  ; PARPRIORUNC[25] = 0.25 # Foliar combustion completeness
      #PARPRIORS[26] = 0.1                  ; PARPRIORUNC[26] = 0.25 # Root / wood combustion completeness
      PARPRIORS[27] = 0.01                 ; PARPRIORUNC[27] = 0.05 # Soil combustion completeness
      # Other priors
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P1.#") {
      PARPRIORS[2] = 0.54                  ; PARPRIORUNC[2] = 0.12 #; PARPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#      PARPRIORS[11]=1.89*14.77735          ; PARPRIORUNC[11]=1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                           # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                           # These observational constraints are not the same and would lead to
                                                                           # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
#      PARPRIORS[11] = 65.0               ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
      PARPRIORS[11]=21.1491                ; PARPRIORUNC[11]=8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                      # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                      # due to the different temperature response functions used in ACM2 vs ACM 1
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      PARPRIORS[25]=OBS$frac_Cwood_coarse_root_prior       ; if (OBS$frac_Cwood_coarse_root_prior != -9999) {PARPRIORUNC[25]=OBS$frac_Cwood_coarse_root_prior_unc} # Fraction Cwood coarse root prior
      #PARPRIORS[28] = 0.87                 ; PARPRIORUNC[28] = 0.41 # Resilience factor
      #PARPRIORS[29] = 0.5                  ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
      #PARPRIORS[30] = 0.1                  ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
      PARPRIORS[31] = 0.01                 ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
      #PARPRIORS[32] = 0.25                 ; PARPRIORUNC[32] = 0.25 # Foliage + root litter combustion completeness
      # Other priors
      OTHERPRIORS[1] = OBS$soilwater       ; OTHERPRIORUNC[1] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 ; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A1.C1.D2.F2.H3.P1.#") {
      PARPRIORS[2] = 0.54                  ; PARPRIORUNC[2] = 0.12 #; PARPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#      PARPRIORS[11]=1.89*14.77735          ; PARPRIORUNC[11]=1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                           # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                           # These observational constraints are not the same and would lead to
                                                                           # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
#      PARPRIORS[11] = 65.0               ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
      PARPRIORS[11]=21.1491                ; PARPRIORUNC[11]=8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                      # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                      # due to the different temperature response functions used in ACM2 vs ACM 1
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      PARPRIORS[25]=OBS$frac_Cwood_coarse_root_prior       ; if (OBS$frac_Cwood_coarse_root_prior != -9999) {PARPRIORUNC[25]=OBS$frac_Cwood_coarse_root_prior_unc} # Fraction Cwood coarse root prior
      PARPRIORS[33]=OBS$minLWP             ; ; if (OBS$minLWP != -9999) {PARPRIORUNC[33]=OBS$minLWP_unc} # minLWP prior
      #PARPRIORS[28] = 0.87                ; PARPRIORUNC[28] = 0.41 # Resilience factor
      #PARPRIORS[29] = 0.5                  ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
      #PARPRIORS[30] = 0.1                  ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
      PARPRIORS[31] = 0.01                 ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
      #PARPRIORS[32] = 0.25                 ; PARPRIORUNC[32] = 0.25 # Foliage + root litter combustion completeness
      #PARPRIORS[33] = -1.8                 ; PARPRIORUNC[33] = 0.5 # minimum leaf water potential
      # Other priors
      OTHERPRIORS[1] = OBS$soilwater       ; OTHERPRIORUNC[1] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 ; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A2.C1.D2.F2.H2.P1.#") {
      PARPRIORS[2] = 0.54                  ; PARPRIORUNC[2] = 0.12 #; PARPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#      PARPRIORS[11]=1.89*14.77735          ; PARPRIORUNC[11]=1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                           # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                           # These observational constraints are not the same and would lead to
                                                                           # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
#      PARPRIORS[11] = 65.0               ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
#      PARPRIORS[11] = 60.0                  ; PARPRIORUNC[11]= 20.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax: derived from multiple trait values from Kattge et al., (2011)
                                                                      # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                      # due to the different temperature response functions used in ACM2 vs ACM 1
      PARPRIORS[11] = 24.43                  ; PARPRIORUNC[11]= 10.4 #; PARPRIORWEIGHT[11] = 1 # Vcmax: Median of CARDAMOM analysis assimilating 4 GPP products
#      PARPRIORS[11] = 54.165                  ; PARPRIORUNC[11]= 20.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax: Mean of reported PFT values from Oliver et al., (2022)
#                                                                      # Note that this prior is difference from DALEC.C1.D1.F2.P1.
#                                                                      # due to the different temperature response functions used in ACM2 vs ACM 1
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      #PARPRIORS[28] = 0.87                ; PARPRIORUNC[28] = 0.41 # Resilience factor
      #PARPRIORS[29] = 0.5                  ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
      #PARPRIORS[30] = 0.1                  ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
      PARPRIORS[31] = 0.01                 ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
      #PARPRIORS[32] = 0.25                 ; PARPRIORUNC[32] = 0.25 # Foliage + root litter combustion completeness
      # Other priors
      OTHERPRIORS[1] = OBS$soilwater       ; OTHERPRIORUNC[1] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A4.C6.D2.F2.H2.P11.#") {
      PARPRIORS[2] = 0.54                  ; PARPRIORUNC[2] = 0.12 #; PARPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#      PARPRIORS[11] = 65.0               ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
#      PARPRIORS[11] = 60.0                  ; PARPRIORUNC[11]= 20.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax: derived from multiple trait values from Kattge et al., (2011)
                                                                      # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                      # due to the different temperature response functions used in ACM2 vs ACM 1
#      PARPRIORS[11] = 24.43                  ; PARPRIORUNC[11]= 10.4 #; PARPRIORWEIGHT[11] = 1 # Vcmax: Median of CARDAMOM analysis assimilating 4 GPP products
      PARPRIORS[11] = 54.165                  ; PARPRIORUNC[11]= 20.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax: Mean of reported PFT values from Oliver et al., (2022)
                                                                      # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                      # due to the different temperature response functions used in ACM2 vs ACM 1
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      #PARPRIORS[28] = 0.87                ; PARPRIORUNC[28] = 0.41 # Resilience factor
      #PARPRIORS[29] = 0.5                  ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
      #PARPRIORS[30] = 0.1                  ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
      PARPRIORS[31] = 0.01                 ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
      #PARPRIORS[32] = 0.25                 ; PARPRIORUNC[32] = 0.25 # Foliage + root litter combustion completeness
      PARPRIORS[33] = 0.05                 ; PARPRIORUNC[33] = 0.05 # labile:biomass at which growth limited by 50 %
      PARPRIORS[36] = 1.0                 ; PARPRIORUNC[36] = 5.0 # temperature at which root growth totally suppressed (oC)
      PARPRIORS[37] = 5.0                 ; PARPRIORUNC[37] = 1.0 # temperature at which wood growth totally suppressed (oC)
      # Other priors
      OTHERPRIORS[1] = OBS$soilwater       ; OTHERPRIORUNC[1] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 #; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A3.C1.D2.F2.H2.P1.#") {
      PARPRIORS[2] = 0.54                  ; PARPRIORUNC[2] = 0.12 #; PARPRIORWEIGHT[2] = noyears # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#      PARPRIORS[11]=1.89*14.77735          ; PARPRIORUNC[11]=1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                           # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                           # These observational constraints are not the same and would lead to
                                                                           # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
#      PARPRIORS[11] = 65.0               ; PARPRIORUNC[11] = 30.0 #; PARPRIORWEIGHT[11] = 1 # Vcmax (gC/m2/day): Wullscheller (1993)
      PARPRIORS[11]=21.1491                ; PARPRIORUNC[11]=8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                      # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                      # due to the different temperature response functions used in ACM2 vs ACM 1
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc #; PARPRIORWEIGHT[17] = noyears
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      #PARPRIORS[28] = 0.87                 ; PARPRIORUNC[28] = 0.41 # Resilience factor
      #PARPRIORS[29] = 0.5                  ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
      #PARPRIORS[30] = 0.1                  ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
      PARPRIORS[31] = 0.01                 ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
      #PARPRIORS[32] = 0.25                 ; PARPRIORUNC[32] = 0.25 # Foliage + root litter combustion completeness
      PARPRIORS[33] = 0.43                 ; PARPRIORUNC[33] = 0.05 # Canopy NIR reflectance
      PARPRIORS[34] = 0.16                 ; PARPRIORUNC[34] = 0.05 # Canopy PAR reflectance
      PARPRIORS[35] = 0.26                 ; PARPRIORUNC[35] = 0.05 # Canopy NIR transmittance
      PARPRIORS[36] = 0.16                 ; PARPRIORUNC[36] = 0.05 # Canopy PAR transmittance
      PARPRIORS[37] = 0.023                ; PARPRIORUNC[37] = 0.05 # Soil reflectance to near infrared radiation
      PARPRIORS[38] = 0.033                ; PARPRIORUNC[38] = 0.05 # Foliage + root litter combustion completeness
      # Other priors
      OTHERPRIORS[1] = OBS$soilwater       ; OTHERPRIORUNC[1] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 ; OTHERPRIORWEIGHT[4] = noyears # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P5.#") {
      PARPRIORS[2] =0.54                   ; PARPRIORUNC[2]=0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#      PARPRIORS[11]=1.89*14.77735          ; PARPRIORUNC[11]=1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                    # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                    # These observational constraints are not the same and would lead to
                                                                    # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
      PARPRIORS[11]=21.1491                ; PARPRIORUNC[11]=8.534234 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                      #       Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                      # due to the different temperature response functions used in ACM2 vs ACM 1
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      #PARPRIORS[28] = 0.5                  ; PARPRIORUNC[28] = 0.25 # Resilience factor
      #PARPRIORS[29] = 0.5                  ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
      #PARPRIORS[30] = 0.1                  ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
      PARPRIORS[31] = 0.01                 ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
      # Other priors
      OTHERPRIORS[1] = OBS$soilwater       ; OTHERPRIORUNC[1] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P6.#") {
      PARPRIORS[2] =0.54                   ; PARPRIORUNC[2]=0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#      PARPRIORS[11]=1.89*14.77735         ; PARPRIORUNC[11]=1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                    # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                    # These observational constraints are not the same and would lead to
                                                                    # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
      PARPRIORS[11]=21.1491                ; PARPRIORUNC[11]=8.534234 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                      #       8.534234 is equal to the limb for 95%CI, therefore half that makes an approximate SD
                                                                      #       Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                      # due to the different temperature response functions used in ACM2 vs ACM 1
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      #PARPRIORS[28] = 0.5                  ; PARPRIORUNC[28] = 0.25 # Resilience factor
      #PARPRIORS[29] = 0.5                  ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
      #PARPRIORS[30] = 0.1                  ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
      PARPRIORS[31] = 0.01                 ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
      # Other priors
      OTHERPRIORS[1] = OBS$soilwater       ; OTHERPRIORUNC[1] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P2.#") {
      PARPRIORS[2] =0.54                   ; PARPRIORUNC[2]=0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#      PARPRIORS[11]=1.89*14.77735          ; PARPRIORUNC[11]=1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                    # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                    # These observational constraints are not the same and would lead to
                                                                    # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
      PARPRIORS[11]=21.1491                ; PARPRIORUNC[11]=8.534234 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                      #       Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                      # due to the different temperature response functions used in ACM2 vs ACM 1
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      #PARPRIORS[28] = 0.5                  ; PARPRIORUNC[28] = 0.25 # Resilience factor
      #PARPRIORS[29] = 0.5                  ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
      #PARPRIORS[30] = 0.1                  ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
      PARPRIORS[31] = 0.01                 ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
      # Other priors
      OTHERPRIORS[1] = OBS$soilwater       ; OTHERPRIORUNC[1] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P1.R1.#") {
      PARPRIORS[1]=0.5                     ; PARPRIORUNC[1]=0.125 # fraction of litter decomposition to Csom
#      PARPRIORS[11]=1.89*14.77735          ; PARPRIORUNC[11]=1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                    # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                    # These observational constraints are not the same and would lead to
                                                                    # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
      PARPRIORS[11]=21.1491                ; PARPRIORUNC[11]=8.534234 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                      #       Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                      # due to the different temperature response functions used in ACM2 vs ACM 1
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      #PARPRIORS[28] = 0.5                  ; PARPRIORUNC[28] = 0.25 # Resilience factor
      #PARPRIORS[29] = 0.5                  ; PARPRIORUNC[29] = 0.25 # Foliar combustion completeness
      #PARPRIORS[30] = 0.1                  ; PARPRIORUNC[30] = 0.25 # Root / wood combustion completeness
      PARPRIORS[31] = 0.01                 ; PARPRIORUNC[31] = 0.05 # Soil combustion completeness
      # Other priors
      OTHERPRIORS[1] = OBS$soilwater       ; OTHERPRIORUNC[1] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
      OTHERPRIORS[2] = 0.54                ; OTHERPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#NOTUSE      OTHERPRIORS[3] = 27.295        ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P1.R1.#") {
      PARPRIORS[1]=0.5                     ; PARPRIORUNC[1]=0.125 # fraction of litter decomposition to Csom
#      PARPRIORS[11]=1.89*14.77735          ; PARPRIORUNC[11]=1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                    # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                    # These observational constraints are not the same and would lead to
                                                                    # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
      PARPRIORS[11]=21.1491                ; PARPRIORUNC[11]=8.534234 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                      # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                      # due to the different temperature response functions used in ACM2 vs ACM 1
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      #PARPRIORS[30] = 0.5                  ; PARPRIORUNC[30] = 0.25 # Resilience factor
      #PARPRIORS[31] = 0.5                  ; PARPRIORUNC[31] = 0.25 # Foliar combustion completeness
      #PARPRIORS[32] = 0.1                  ; PARPRIORUNC[32] = 0.25 # Root / wood combustion completeness
      PARPRIORS[33] = 0.01                 ; PARPRIORUNC[33] = 0.05 # Soil combustion completeness
      # Other priors
      OTHERPRIORS[1] = OBS$soilwater       ; OTHERPRIORUNC[1] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
      OTHERPRIORS[2] = 0.54                ; OTHERPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#NOTUSE      OTHERPRIORS[3] = 27.295        ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P2.R1.#") {
      PARPRIORS[1]=0.5                     ; PARPRIORUNC[1]=0.125 # fraction of litter decomposition to Csom
#      PARPRIORS[11]=1.89*14.77735          ; PARPRIORUNC[11]=1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                    # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                    # These observational constraints are not the same and would lead to
                                                                    # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
      PARPRIORS[11]=21.1491                ; PARPRIORUNC[11]=8.534234 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                      # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                      # due to the different temperature response functions used in ACM2 vs ACM 1
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      #PARPRIORS[31] = 0.5                  ; PARPRIORUNC[31] = 0.25 # Resilience factor
      #PARPRIORS[32] = 0.5                  ; PARPRIORUNC[32] = 0.25 # Foliar combustion completeness
      #PARPRIORS[33] = 0.1                  ; PARPRIORUNC[33] = 0.25 # Root / wood combustion completeness
      PARPRIORS[34] = 0.01                 ; PARPRIORUNC[34] = 0.05 # Soil combustion completeness
      # Other priors
      OTHERPRIORS[1] = OBS$soilwater       ; OTHERPRIORUNC[1] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
      OTHERPRIORS[2] = 0.54                ; OTHERPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#NOTUSE      OTHERPRIORS[3] = 27.295        ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P2.R3.#") {
      PARPRIORS[1]=0.5                     ; PARPRIORUNC[1]=0.125 # fraction of litter decomposition to Csom
#      PARPRIORS[11]=1.89*14.77735          ; PARPRIORUNC[11]=1.89*0.4696238 # Derived from ACM2 recalibration.
                                                                    # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                    # These observational constraints are not the same and would lead to
                                                                    # overestimation of GPP (SPA = 34, ACM2 = 15), but here multiple by avN (1.89) to get Ceff
#      PARPRIORS[11]=21.1491                ; PARPRIORUNC[11]=8.534234 # Ceff: derived from multiple trait values from Kattge et al., (2011)
#                                                                      # Note that this prior is difference from DALEC.C1.D1.F2.P1.
#                                                                      # due to the different temperature response functions used in ACM2 vs ACM 1
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      #PARPRIORS[31] = 0.5                  ; PARPRIORUNC[31] = 0.25 # Resilience factor
      #PARPRIORS[32] = 0.5                  ; PARPRIORUNC[32] = 0.25 # Foliar combustion completeness
      #PARPRIORS[33] = 0.1                  ; PARPRIORUNC[33] = 0.25 # Root / wood combustion completeness
      PARPRIORS[34] = 0.01                 ; PARPRIORUNC[34] = 0.05 # Soil combustion completeness
      # Other priors
      OTHERPRIORS[1] = OBS$soilwater       ; OTHERPRIORUNC[1] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
      OTHERPRIORS[2] = 0.54                ; OTHERPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#NOTUSE      OTHERPRIORS[3] = 27.295        ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
 } else if (modelname == "DALEC.D1.F2.#") {
      PARPRIORS[2] = 0.54                  ; PARPRIORUNC[2]=0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      PARPRIORS[11] = 16.9                 ; PARPRIORUNC[11]=7.502147 # Ceff: derived from multiple trait values from Kattge et al., (2011)
      PARPRIORS[13] = OBS$Cfol_initial     ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[13] = OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[14] = OBS$Croots_initial   ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[14] = OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[15] = OBS$Cwood_initial    ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[15] = OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[16] = OBS$Clit_initial     ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[16] = OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[17] = OBS$Csom_initial     ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[17] = OBS$Csom_initial_unc} # Csom prior
      # Other priors
#      OTHERPRIORS[5] =                    ; OTHERPRIORUNC[5] = # Steady state attractor for wood
  } else if (modelname == "DALEC.C5.D1.F2.P1.#") {
      PARPRIORS[1] = 0.54                  ; PARPRIORUNC[1]=0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      PARPRIORS[7] = 16.9                  ; PARPRIORUNC[7]=7.502147 # Ceff: derived from multiple trait values from Kattge et al., (2011)
      PARPRIORS[15] = OBS$Cfol_initial     ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[15] = OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[16] = OBS$Cwood_initial    ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[16] = OBS$Cwood_initial_unc} # Croot + Cwood prior
      PARPRIORS[17] = OBS$Csom_initial     ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[17] = OBS$Csom_initial_unc} # Csom + Clitter prior
  } else if (modelname == "DALEC.C4.D1.F2.#") {
      PARPRIORS[1] = 0.54                  ; PARPRIORUNC[1]=0.12  # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      PARPRIORS[7] = 16.9                  ; PARPRIORUNC[7]=7.502147 # Ceff: derived from multiple trait values from Kattge et al., (2011)
      PARPRIORS[9] = OBS$Cfol_initial      ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[9] = OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[10] = OBS$Cwood_initial    ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[10] = OBS$Cwood_initial_unc} # Croot + Cwood prior
      PARPRIORS[11] = OBS$Csom_initial     ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[11] = OBS$Csom_initial_unc} # Csom + Clitter prior
  } else if (modelname == "DALEC.A1.C2.D2.F2.H1.P3.R1.#") {
      PARPRIORS[11]=0.2764618              ; PARPRIORUNC[11]=0.2014871 # log10 avg foliar N (gN.m-2)
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      PARPRIORS[1]=0.5              ; PARPRIORUNC[1]=0.125 # fraction of litter decomposition to Csom
#      PARPRIORS[36]=14.77735        ; PARPRIORUNC[36]=0.4696238*2 # Derived from ACM2 recalibration.
                                                                   # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                   # These observational constraints are not the same and would lead to
                                                                   # overestimation of GPP (SPA = 34, ACM2 = 15)
      PARPRIORS[36]=11.197440          ; PARPRIORUNC[36]=9.3         # NUE prior derived from Kattge et al., (2011), based on log10 gaussian distribution
      # other priors
      OTHERPRIORS[1] = 0.54   ; OTHERPRIORUNC[1] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      # POSITION 2 used for water which does not apply here
#      OTHERPRIORS[3] = 27.295 ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A1.C2.D2.F2.H1.P4.R2.#") {
      PARPRIORS[1]=0.5                     ; PARPRIORUNC[1]=0.125 # fraction of litter decomposition to Csom
      PARPRIORS[11]=0.2764618              ; PARPRIORUNC[11]=0.2014871 # log10 avg foliar N (gN.m-2)
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      PARPRIORS[42]=11.197440              ; PARPRIORUNC[42]=9.3  # NUE prior derived from Kattge et al., (2011), based on log10 gaussian distribution
#      PARPRIORS[42]=14.77735               ; PARPRIORUNC[42]=0.4696238*2 # Derived from ACM2 recalibration.
                                                                # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                # These observational constraints are not the same and would lead to
                                                                # overestimation of GPP (SPA = 34, ACM2 = 15)
      PARPRIORS[43]=275.1452               ; PARPRIORUNC[43]=296.2767 # Leaf lifespan prior form Kattge et al., 2011, based on log10 gauusian distribution
      # other priors
      OTHERPRIORS[1] = 0.54                ; OTHERPRIORUNC[1] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
#NOT IN USE      OTHERPRIORS[2] = OBS$soilwater ; OTHERPRIORUNC[2] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
#      OTHERPRIORS[3] = 27.295              ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P4.R2.#") {
      PARPRIORS[1]=0.5                     ; PARPRIORUNC[1]=0.125 # fraction of litter decomposition to Csom
      PARPRIORS[11]=0.2764618              ; PARPRIORUNC[11]=0.2014871 # log10 avg foliar N (gN.m-2)
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      PARPRIORS[42]=11.197440              ; PARPRIORUNC[42]=9.3  # NUE prior derived from Kattge et al., (2011), based on log10 gaussian distribution
#      PARPRIORS[42]=14.77735               ; PARPRIORUNC[42]=0.4696238*2 # Derived from ACM2 recalibration.
                                                                # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                # These observational constraints are not the same and would lead to
                                                                # overestimation of GPP (SPA = 34, ACM2 = 15)
#      PARPRIORS[43]=275.1452              ; PARPRIORUNC[43]=296.2767 # Leaf lifespan prior form Kattge et al., 2011, based on log10 gauusian distribution
      # other priors
      OTHERPRIORS[1] = 0.54                ; OTHERPRIORUNC[1] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      OTHERPRIORS[2] = OBS$soilwater       ; OTHERPRIORUNC[2] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
#      OTHERPRIORS[3] = 27.295              ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood

  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P7.R2.#") {
      PARPRIORS[1]=0.5                     ; PARPRIORUNC[1]=0.125 # fraction of litter decomposition to Csom
      PARPRIORS[11]=0.2764618              ; PARPRIORUNC[11]=0.2014871 # log10 avg foliar N (gN.m-2)
      PARPRIORS[14]=0                      ; PARPRIORUNC[14]=5  # minimum CGI temperature
      PARPRIORS[15]=30                     ; PARPRIORUNC[15]=5 # optimum CMI and CGI temperatures
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      PARPRIORS[32]=-1.8                   ; PARPRIORUNC[32]=1 # minLWP (MPa)
      PARPRIORS[42]=11.197440              ; PARPRIORUNC[42]=9.3  # NUE prior derived from Kattge et al., (2011), based on log10 gaussian distribution
#      PARPRIORS[42]=14.77735               ; PARPRIORUNC[42]=0.4696238 # Derived from ACM2 recalibration.
                                                                # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                # These observational constraints are not the same and would lead to
                                                                # overestimation of GPP (SPA = 34, ACM2 = 15)
#      PARPRIORS[43]=275.1452               ; PARPRIORUNC[43]=296.2767 # Leaf lifespan prior form Kattge et al., 2011, based on log10 gauusian distribution
      # other priors
      OTHERPRIORS[1] = 0.54                ; OTHERPRIORUNC[1] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      OTHERPRIORS[2] = OBS$soilwater       ; OTHERPRIORUNC[2] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
#      OTHERPRIORS[3] = 27.295              ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
      OTHERPRIORS[6] = 0.5                 ; OTHERPRIORUNC[6] = 0.05 # Rleaf:Rauto ratio See Atkin et al., various for approximate ratio
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P8.R2.#") {
      PARPRIORS[1]=0.5                     ; PARPRIORUNC[1]=0.125 # fraction of litter decomposition to Csom
      PARPRIORS[11]=0.2764618              ; PARPRIORUNC[11]=0.2014871 # log10 avg foliar N (gN.m-2)
      PARPRIORS[14]=0                      ; PARPRIORUNC[14]=5  # minimum CGI temperature
      PARPRIORS[15]=30                     ; PARPRIORUNC[15]=5 # optimum CMI and CGI temperatures
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      PARPRIORS[32]=-1.8                   ; PARPRIORUNC[32]=1 # minLWP (MPa)
      PARPRIORS[42]=11.197440              ; PARPRIORUNC[42]=9.3  # NUE prior derived from Kattge et al., (2011), based on log10 gaussian distribution
#      PARPRIORS[42]=14.77735               ; PARPRIORUNC[42]=0.4696238 # Derived from ACM2 recalibration.
                                                                # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                # These observational constraints are not the same and would lead to
                                                                # overestimation of GPP (SPA = 34, ACM2 = 15)
#      PARPRIORS[43]=275.1452               ; PARPRIORUNC[43]=296.2767 # Leaf lifespan prior form Kattge et al., 2011, based on log10 gauusian distribution
      # other priors
      OTHERPRIORS[1] = 0.54                ; OTHERPRIORUNC[1] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      OTHERPRIORS[2] = OBS$soilwater       ; OTHERPRIORUNC[2] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
#      OTHERPRIORS[3] = 27.295              ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
      OTHERPRIORS[6] = 0.5                 ; OTHERPRIORUNC[6] = 0.05 # Rleaf:Rauto ratio See Atkin et al., various for approximate ratio
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P9.R2.#") {
      PARPRIORS[1]=0.5                     ; PARPRIORUNC[1]=0.125 # fraction of litter decomposition to Csom
      PARPRIORS[11]=0.2764618              ; PARPRIORUNC[11]=0.2014871 # log10 avg foliar N (gN.m-2)
      PARPRIORS[14]=0                      ; PARPRIORUNC[14]=5  # minimum CGI temperature
      PARPRIORS[15]=30                     ; PARPRIORUNC[15]=5 # optimum CMI and CGI temperatures
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      PARPRIORS[32]=-1.8                   ; PARPRIORUNC[32]=1 # minLWP (MPa)
      PARPRIORS[42]=11.197440              ; PARPRIORUNC[42]=9.3  # NUE prior derived from Kattge et al., (2011), based on log10 gaussian distribution
#      PARPRIORS[42]=14.77735               ; PARPRIORUNC[42]=0.4696238 # Derived from ACM2 recalibration.
                                                                # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                # These observational constraints are not the same and would lead to
                                                                # overestimation of GPP (SPA = 34, ACM2 = 15)
#      PARPRIORS[43]=275.1452               ; PARPRIORUNC[43]=296.2767 # Leaf lifespan prior form Kattge et al., 2011, based on log10 gauusian distribution
      # other priors
      OTHERPRIORS[1] = 0.54                ; OTHERPRIORUNC[1] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      OTHERPRIORS[2] = OBS$soilwater       ; OTHERPRIORUNC[2] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
#      OTHERPRIORS[3] = 27.295              ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
      OTHERPRIORS[6] = 0.5                 ; OTHERPRIORUNC[6] = 0.05 # Rleaf:Rauto ratio See Atkin et al., various for approximate ratio
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P10.R2.#") {
      PARPRIORS[1]=0.5                     ; PARPRIORUNC[1]=0.125 # fraction of litter decomposition to Csom
      PARPRIORS[11]=0.2764618              ; PARPRIORUNC[11]=0.2014871 # log10 avg foliar N (gN.m-2)
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      PARPRIORS[42]=11.197440              ; PARPRIORUNC[42]=9.3  # NUE prior derived from Kattge et al., (2011), based on log10 gaussian distribution
#      PARPRIORS[42]=14.77735               ; PARPRIORUNC[42]=0.4696238 # Derived from ACM2 recalibration.
                                                                # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                # These observational constraints are not the same and would lead to
                                                                # overestimation of GPP (SPA = 34, ACM2 = 15)
#      PARPRIORS[43]=275.1452               ; PARPRIORUNC[43]=296.2767 # Leaf lifespan prior form Kattge et al., 2011, based on log10 gauusian distribution
      # other priors
      OTHERPRIORS[1] = 0.54                ; OTHERPRIORUNC[1] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      OTHERPRIORS[2] = OBS$soilwater       ; OTHERPRIORUNC[2] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
#      OTHERPRIORS[3] = 27.295              ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
      OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P3.R1.#") {
      PARPRIORS[11]=0.2764618	             ; PARPRIORUNC[11]=0.2014871 # log10 avg foliar N (gN.m-2)
      PARPRIORS[17]=OBS$lca                ; PARPRIORUNC[17]=OBS$lca_unc
#      PARPRIORS[17]=35.5                   ; PARPRIORUNC[17]=35.5*0.23 # Kiuic LCA prior
      PARPRIORS[19]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[19]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[20]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[20]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[21]=OBS$Cwood_initial      ; if (OBS$Cwood_initial != -9999) {PARPRIORUNC[21]=OBS$Cwood_initial_unc} # Cwood prior
      PARPRIORS[22]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[22]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      PARPRIORS[1]=0.5            ; PARPRIORUNC[1]=0.125 # fraction of litter decomposition to Csom
      PARPRIORS[34]=-0.0005602    ; PARPRIORUNC[34]= 6.9e-05 # GSI sensitivity for leaf scenescence
      PARPRIORS[36]=11.197440     ; PARPRIORUNC[36]=9.3 # NUE prior derived from Kattge et al., (2011), based on log10 gaussian distribution
#      PARPRIORS[36]=14.77735       ; PARPRIORUNC[36]=0.4696238*2 # Derived from ACM2 recalibration.
                                                                # Note despite having the same name as ecosystem property of Amax per gN or SPA's kappaC
                                                                # These observational constraints are not the same and would lead to
                                                                # overestimation of GPP (SPA = 34, ACM2 = 15)
      PARPRIORS[43] = 0.07        ; PARPRIORUNC[43]=0.036 # Combustion completeness for fine root and wood
      PARPRIORS[44] = 0.016       ; PARPRIORUNC[44]=0.005 # Combustion completeness for soil
      PARPRIORS[45] = 0.11        ; PARPRIORUNC[45]=0.025 # Combustion completeness for foliage and fine root litter
      PARPRIORS[46] = 0.11        ; PARPRIORUNC[46]=0.028 # Combustion completeness for wood litter
      # other priors
      OTHERPRIORS[1] = 0.54        ; OTHERPRIORUNC[1] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      OTHERPRIORS[2] = OBS$soilwater ; OTHERPRIORUNC[2] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
#      OTHERPRIORS[3] = 27.295      ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
      OTHERPRIORS[4] = 0.66        ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
      OTHERPRIORS[5] = OBS$Cwood_potential ; OTHERPRIORUNC[5] = OBS$Cwood_potential_unc # Steady state attractor for wood
  } else if (modelname == "DALEC.M2.#") {
      PARPRIORS[2] = 0.54                  ; PARPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      PARPRIORS[10]=16.9                   ; PARPRIORUNC[10]=7.502147 # Ceff: derived from multiple trait values from Kattge et al., (2011)
      #PARPRIORS[15]=OBS$lca                ; PARPRIORUNC[15]=OBS$lca_unc
      PARPRIORS[15]=32                     ; PARPRIORUNC[15]=13
      PARPRIORS[17]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[17]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[18]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[18]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[19]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[19]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      # other priors
#      OTHERPRIORS[2] = OBS$soilwater       ; OTHERPRIORUNC[2] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
#      OTHERPRIORS[3] = 27.295              ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
  } else if (modelname == "DALEC.A3.H2.M2.#") {
      PARPRIORS[2] = 0.54                  ; PARPRIORUNC[2] = 0.12 # Ra:GPP Collalti & Prentice (2019), Tree Physiology, 10.1093/treephys/tpz034
      PARPRIORS[11]=21.1491                ; PARPRIORUNC[11]=8.534234 #; PARPRIORWEIGHT[11] = 1 # Ceff: derived from multiple trait values from Kattge et al., (2011)
                                                                      # Note that this prior is difference from DALEC.C1.D1.F2.P1.
                                                                      # due to the different temperature response functions used in ACM2 vs ACM 1
      #PARPRIORS[15]=OBS$lca                ; PARPRIORUNC[15]=OBS$lca_unc
      PARPRIORS[15]=32                     ; PARPRIORUNC[15]=13
      PARPRIORS[17]=OBS$Cfol_initial       ; if (OBS$Cfol_initial != -9999) {PARPRIORUNC[17]=OBS$Cfol_initial_unc} # Cfoliar prior
      PARPRIORS[18]=OBS$Croots_initial     ; if (OBS$Croots_initial != -9999) {PARPRIORUNC[18]=OBS$Croots_initial_unc} # Croots prior
      PARPRIORS[19]=OBS$Clit_initial       ; if (OBS$Clit_initial != -9999) {PARPRIORUNC[19]=OBS$Clit_initial_unc} # Clitter prior
      PARPRIORS[23]=OBS$Csom_initial       ; if (OBS$Csom_initial != -9999) {PARPRIORUNC[23]=OBS$Csom_initial_unc} # Csom prior
      # other priors
#      OTHERPRIORS[2] = OBS$soilwater       ; OTHERPRIORUNC[2] = OBS$soilwater_unc # Initial soil water fraction (GLEAM v3.1a)
#      OTHERPRIORS[3] = 27.295              ; OTHERPRIORUNC[3] = 11.03755 # Foliar C:N (gC/gN) prior derived from Kattge et al., (2011)
       OTHERPRIORS[4] = 0.66                ; OTHERPRIORUNC[4] = 0.12 # Prior on mean annual ET/P See Zhang et al., (2018) doi:10.5194/hess-22-241-2018
  } else if (modelname == "ACM") {

      # For ACM_GPP_ET
      # p(1) = nitrogen use efficiency at optimum temperature (oC)
      #,unlimited by CO2, light and photoperiod (34gC/gN)
#      PARPRIORS[1] = 35.65 ; PARPRIORUNC[1] = 35.65 * 0.1 #26.19
      # p(2) = maximum temperature at which photosynthesis occurs (oC) (57.05oC)
#      PARPRIORS[2] = 57.05 ; PARPRIORUNC[2] = 57.05 * 0.1
      # p(3) = optimum temperature for photosynthesis (oC) (30oC)
#      PARPRIORS[3] = 30.0 ; PARPRIORUNC[3] = 30.0 * 0.1
      # p(4) = kurtosis for temperature response of photosynthesis (0.185912)
      PARPRIORS[4] = 0.185912 ; PARPRIORUNC[4] = 0.185912 * 0.01
      # p(5) = SPA apparent quantium yield gC/MJ PAR = ~3.2
#    PARPRIORS[5] = 3.2 ; PARPRIORUNC[5] = 0.1
      # p(6) = min leaf water potential (MPa)
#      PARPRIORS[6] = -2.0 ; PARPRIORUNC[6] = abs(-2.0*0.1)
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
  }

  # HACK: If file name contains specific COMPLEX experiment code
  # indicating the file should not have any observational constraint reset all
#  if (grepl("exp1f",file) | grepl("exp2f",file) | grepl("exp3f",file) | grepl("exp4f",file)) {
#      #ONLY USED FOR LOG NORMALLY PRESCRIBED PARAMETER PRIORS
#      PARPRIORS = rep(-9999.0,length.out=100)
#      PARPRIORUNC = rep(-9999.0,length.out=100)
#      #For all other multiparameter user-defined priors
#      OTHERPRIORS = rep(-9999.0,length.out=100)
#      OTHERPRIORUNC = rep(-9999.0,length.out=100)
#  }

  # Combine the forcings and observations matrices together
  DATA_TEMP = t(cbind(MET,OBSMAT))
  # combine the static data
  DATA_STAT = c(PARPRIORS,PARPRIORUNC,PARPRIORWEIGHT,OTHERPRIORS,OTHERPRIORUNC,OTHERPRIORWEIGHT)

  # Possible location to set all NA values to -9999?
  # This would allow bit by bit replacement of NA to -9999 to be removed
#  static_data[is.na(static_data)] = -9999
#  DATA_STAT[is.na(DATA_STAT)] -9999
#  DATA_TEMP[is.na(DATA_TEMP)] -9999

  # open the binary file
  zz <- file(file, "wb")
  # write with 8 bit precision (i.e. double)
  writeBin(as.double(static_data), zz)
  writeBin(as.double(DATA_STAT), zz)
  writeBin(as.double(as.vector(DATA_TEMP)), zz)
  close(zz)

}
## Use byte compile
binary_data<-cmpfun(binary_data)
