   
###
## Generic script to generate summary information on the simulated C budget
## and compare CARDAMOM estimates with its assimilated and independent observations
## NOTE: we assume this is a gridded run
###

## Sections of code you may want to change are flagged
## "PointsOfChange"

#### TO DO
# Add code to read in MODIS GPP products
# Add code to automatically aggregate each observation dataset to the resolution of the analysis where possible?

###
## Load needed libraries and framework functions
###

# Read library
library(fields)
library(compiler)
library(RColorBrewer)
library(plotrix)
library(zoo)
library(raster)
library(ncdf4)
library(abind)

# Load any CARDAMOM functions which might be useful
source("~/WORK/GREENHOUSE/models/CARDAMOM/R_functions/generate_wgs_grid.r")
source("~/WORK/GREENHOUSE/models/CARDAMOM/R_functions/calc_pixel_area.r")
source("~/WORK/GREENHOUSE/models/CARDAMOM/R_functions/read_binary_file_format.r")
source("~/WORK/GREENHOUSE/models/CARDAMOM/R_functions/function_closest2d.r")
source("~/WORK/GREENHOUSE/models/CARDAMOM/R_functions/plotconfidence.r")
source("~/WORK/GREENHOUSE/models/CARDAMOM/R_functions/read_src_model_priors.r")

# Function to determine the number of days in any given year
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

ensemble_within_range<-function(target,proposal) {

   # Determine what proportion of a proposed PDF is within a target range
   # Returned value 0-1

   t_range = range(target, na.rm=TRUE)
   in_range = length(which(proposal >= t_range[1] & proposal <= t_range[2]))
   return(in_range / length(proposal))

} # ensemble_within_range

fudgeit <- function(){
  # fudgeit.leg.lab, label for the colour scale must be added as a global variable
  # function to plot a legend to the smoothScatter plot
  xm <- get('xm', envir = parent.frame(1))
  ym <- get('ym', envir = parent.frame(1))
  z  <- get('dens', envir = parent.frame(1))
  colramp <- get('colramp', parent.frame(1))
  fields::image.plot(xm,ym,z, col = colramp(255), legend.only = T, legend.line = 2,
                     axis.args = list(hadj = 0.4), horizontal = FALSE,
                     legend.cex = 0.9, legend.lab=fudgeit.leg.lab, add = F,
#                     smallplot = c(.78,.81,0.28,0.85))
                     smallplot = c(0.97-0.12,1.0-0.12,0.28,0.85))
} # end function fudgeit

###
## Analysis specific information and generic creation
###

###
## Load analysis

# PointsOfChange

#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_RmHeskel_Rg_CWD_wMRT_MHMCMC/global_4deg_C15/infofile.RData")
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT_MHMCMC/Miombo_0.25deg_allWood/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT_MHMCMC/Trendyv9_historical/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT_MHMCMC/ODA_extension_Africa/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT_MHMCMC/ODA_extension_Africa_lca/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT_MHMCMC/ODA_extension_Africa_gpp/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT_MHMCMC/ODA_extension_Africa_lca_gpp/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT_MHMCMC/ODA_extension_Africa_5%CI_agb/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT_MHMCMC/ODA_extension_Africa_10%CI_agb/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT_MHMCMC/ODA_extension_Africa_20%CI_agb/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT_MHMCMC/ODA_extension_Africa_actualCI_agb/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT_MHMCMC/ODA_extension_Africa_one_AGB/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_RmRg_CWD_wMRT_MHMCMC/ODA_extension_Africa_250gCI_agb/infofile.RData")
load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))
load(paste(PROJECT$results_processedpath,PROJECT$name,"_parameter_maps.RData",sep=""))

# Set output path for figures and tables
#out_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/ESSD_update/figures_4deg_C15/"
#out_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/LTSS_CARBON_INTEGRATION/figures_africa/"
out_dir = "~/WORK/GREENHOUSE/models/CARDAMOM/SECO/figures/"

# Specify the position within the stored ensemble for the median estimate and the desired uncertainty bands
mid_quant = 4 ; low_quant = 2 ; high_quant = 6
wanted_quant = c(low_quant,3,mid_quant,5,high_quant)

# Extract timing information
run_years = as.numeric(PROJECT$start_year) : as.numeric(PROJECT$end_year)
nos_years = length(as.numeric(PROJECT$start_year) : as.numeric(PROJECT$end_year))
steps_per_year = length(PROJECT$model$timestep_days) / nos_years

###
## Determine needed spatial information

# generate the lat / long grid again
output = generate_wgs84_grid(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution)
grid_lat = array(output$lat, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
grid_long = array(output$long,dim=c(PROJECT$long_dim,PROJECT$lat_dim))
cardamom_ext = output$cardamom_ext
# then generate the area estimates for each pixel (m)
area = calc_pixel_area(output$lat,output$long,PROJECT$resolution)
# this output is in vector form and we need matching array shapes so...
area = array(area, dim=c(PROJECT$long_dim,PROJECT$lat_dim))

### 
## Create land mask / boundary overlays needed

# PointsOfChange

# This will be used to filter the analysis to include specific locations only
use_filter = TRUE
if (use_filter) {
    #  Design a user created / loaded filter 
    landfilter = array(NA,dim=dim(grid_parameters$obs_wood_gCm2))
    landfilter[which(grid_parameters$obs_wood_gCm2 > 0)] = 1
} else { 
    # Use this option if you don't want to filter
    landfilter = array(1,dim=dim(grid_parameters$obs_wood_gCm2)) 
}

# Load any map information you might want to overlay, such as biome or national boundaries
# if you don't want to have any set as false

# PointsOfChange

add_biomes = FALSE
if (add_biomes) {
    # Read in shape file for boundaries
    landmask = shapefile("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/SECO/analysis/ssa_wwf_dissolve/ssa_wwf_dissolve.shp")
    # just to be sure enforce the projection to WGS-84
    landmask = spTransform(landmask, crs(cardamom_ext))
    # Clip to the extent of the CARDAMOM analysis
    landmask = crop(landmask, cardamom_ext)
} else {
    # load global shape file for land sea mask
    landmask = shapefile("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/R_functions/global_map/national_boundaries/ne_10m_admin_0_countries.shx")
    # just to be sure enforce the projection to WGS-84
    landmask = spTransform(landmask,crs(cardamom_ext))
    # subset by continent 
#    landmask = subset(landmask, CONTINENT == "South America") # Change continent to target area or comment out if spanning zones
    # Clip to the extent of the CARDAMOM analysis
    landmask = crop(landmask, cardamom_ext)
} 

###
## Create plotting colour schemes needed

# set up colour scheme
#smoothScatter_colours=colorRampPalette(c("white",brewer.pal(9,"YlOrRd")))
#smoothScatter_colours=colorRampPalette(c("white",brewer.pal(9,"Blues")))
#smoothScatter_colours=colorRampPalette(c("white","grey92",brewer.pal(9,"Blues")[2:9]))
smoothScatter_colours=colorRampPalette(c("white",rep(rev(brewer.pal(11,"Spectral")),each=3)))
colour_choices_default = colorRampPalette(brewer.pal(11,"Spectral")) 
colour_choices_sign = colorRampPalette(brewer.pal(11,"PRGn"))
#colour_choices_gain = colorRampPalette(brewer.pal(9,"Greens"))
colour_choices_gain = colorRampPalette(brewer.pal(9,"YlGnBu"))
#colour_choices_loss = colorRampPalette(brewer.pal(9,"Reds"))
colour_choices_loss = colorRampPalette(brewer.pal(9,"YlOrRd"))
colour_choices_CI = colorRampPalette(brewer.pal(9,"Purples"))
# Model specific colour choices
scenario_colours = colorRampPalette(brewer.pal(8,"Dark2"))
scenario_colours = scenario_colours(4)
model_colours = colorRampPalette(brewer.pal(12,"Paired"))
model_colours = model_colours(5)
obs_colours = colorRampPalette(brewer.pal(8,"Dark2"))
obs_colours = obs_colours(4)

# array sizes are always the same so
colour_choices_default = colour_choices_default(100)
colour_choices_sign = colour_choices_sign(100)
colour_choices_gain = colour_choices_gain(100)
colour_choices_loss = colour_choices_loss(100)
colour_choices_CI = colour_choices_CI(100)

###
## Aggregate information needed for calibration data comparison

# Extract gridded information on the observations
dims = dim(grid_output$mean_lai_m2m2)
# Soil prior
SoilCPrior = array(NA, dim=c(dims[1], dims[2]))
# Mean annual LAI obs
LAIobs = array(NA, dim=c(dims[1],dims[2],nos_years))
# Disturbance
HarvestFraction = array(NA, dim=c(dims[1], dims[2]))
BurnedFraction = array(NA, dim=c(dims[1], dims[2]))
FireFreq = array(NA, dim=c(dims[1],dims[2]))
# Observed wood trends information
WoodCobs = array(NA, dim=c(dims[1], dims[2],length(PROJECT$model$timestep_days)))
WoodCobs_CI = array(NA, dim=c(dims[1], dims[2],length(PROJECT$model$timestep_days)))
WoodCobs_trend_map = array(NA, dim=c(dims[1], dims[2]))
WoodCobs_trend = rep(NA, PROJECT$nosites)
mean_obs_wood = rep(NA, PROJECT$nosites)
WoodCobs_mean_CI = rep(0, PROJECT$nosites)
# Modelled wood trends information
WoodC = array(NA, dim=c(dims[1], dims[2],length(PROJECT$model$timestep_days)))
wood_trend_map = array(NA, dim=c(dims[1], dims[2]))
wood_trend = rep(NA, PROJECT$nosites)
mean_wood = rep(NA, PROJECT$nosites)
# Loop through every site
for (n in seq(1, PROJECT$nosites)) {
     # Check that location has run
     if (is.na(grid_output$i_location[n]) == FALSE & is.na(grid_output$j_location[n]) == FALSE & is.na(landfilter[grid_output$i_location[n],grid_output$j_location[n]]) == FALSE) {
         # Read in pixel driving data
         drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep=""))
         # Determine forest harvest intensity
         HarvestFraction[grid_output$i_location[n],grid_output$j_location[n]] = sum(drivers$met[,8]) / nos_years
         # Determine mean annual fire intensity and frequency
         BurnedFraction[grid_output$i_location[n],grid_output$j_location[n]] = sum(drivers$met[,9]) / nos_years
         FireFreq[grid_output$i_location[n],grid_output$j_location[n]] = length(which(drivers$met[,9] > 0)) / nos_years
         # Load any priors
         SoilCPrior[grid_output$i_location[n],grid_output$j_location[n]] = drivers$parpriors[23]
         # Clear missing data from and extract observed LAI
         drivers$obs[which(drivers$obs[,3] == -9999),3] = NA
         LAIobs[grid_output$i_location[n],grid_output$j_location[n],] = rollapply(drivers$obs[,3], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)
         # If wood stock estimate available get that too
         tmp = which(drivers$obs[,13] > 0)
         if (length(tmp) > 0) {
             for (t in seq(1, length(tmp))) {
                  # Observational constraint
                  WoodCobs[grid_output$i_location[n],grid_output$j_location[n],tmp[t]] = drivers$obs[tmp[t],13]
                  WoodCobs_CI[grid_output$i_location[n],grid_output$j_location[n],tmp[t]] = drivers$obs[tmp[t],14]
                  # Corresponding model output
                  WoodC[grid_output$i_location[n],grid_output$j_location[n],tmp[t]] = grid_output$wood_gCm2[n,mid_quant,tmp[t]]
             } # loop time steps with obs
             # Wood stock trends
             obs_period_start = tmp[1] ; obs_period_end = tmp[length(tmp)] ; obs_period_years = length(c(obs_period_start:obs_period_end))      
             WoodCobs_trend[n] = (coef(lm(WoodCobs[grid_output$i_location[n],grid_output$j_location[n],obs_period_start:obs_period_end] ~ c(1:obs_period_years)))[2] * 12) # *12 is month to yr adjustment
             WoodCobs_trend_map[grid_output$i_location[n],grid_output$j_location[n]] = WoodCobs_trend[n]
             wood_trend[n] = (coef(lm(grid_output$wood_gCm2[n,mid_quant,obs_period_start:obs_period_end] ~ c(1:obs_period_years)))[2] * 12)
             wood_trend_map[grid_output$i_location[n],grid_output$j_location[n]] = wood_trend[n]
             mean_obs_wood[n] = mean(WoodCobs[grid_output$i_location[n],grid_output$j_location[n],], na.rm=TRUE)
             WoodCobs_mean_CI[n] = mean(drivers$obs[tmp,14], na.rm=TRUE)
             mean_wood[n] = mean(grid_output$mean_wood_gCm2[grid_output$i_location[n],grid_output$j_location[n],mid_quant])
         } # we have more than zero obs
     } # Did this location run
} # Site loop

# Return some information to user
SignalNoise = length(which(abs(WoodCobs_trend*length(run_years)) > as.vector(WoodCobs_mean_CI) & abs(WoodCobs_trend*length(run_years)) > 1)) 
SignalNoise = SignalNoise / length(which(WoodCobs_mean_CI > 0 & abs(WoodCobs_trend*length(run_years)) > 1))
print(paste("Percentage of locations where observed change is greater than CI = ",round(SignalNoise*1e2, digits = 3)," %", sep=""))

###
## Create time series aggregated information for the domain

# Initialise lai for grid annuals
cumarea = 0
lai_grid = array(NA,dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
lai_m2m2 = rep(0,nos_years) ; lai_lower_m2m2 = rep(0,nos_years) ; lai_upper_m2m2 = rep(0,nos_years)
gpp_TgCyr = rep(0,nos_years) ; gpp_lower_TgCyr = rep(0,nos_years) ; gpp_upper_TgCyr = rep(0,nos_years)
rauto_TgCyr = rep(0,nos_years) ; rauto_lower_TgCyr = rep(0,nos_years) ; rauto_upper_TgCyr = rep(0,nos_years)
rhet_TgCyr = rep(0,nos_years) ; rhet_lower_TgCyr = rep(0,nos_years) ; rhet_upper_TgCyr = rep(0,nos_years)
nee_TgCyr = rep(0,nos_years) ; nee_lower_TgCyr = rep(0,nos_years) ; nee_upper_TgCyr = rep(0,nos_years)
nbe_TgCyr = rep(0,nos_years) ; nbe_lower_TgCyr = rep(0,nos_years) ; nbe_upper_TgCyr = rep(0,nos_years)
fire_TgCyr = rep(0,nos_years) ; fire_lower_TgCyr = rep(0,nos_years) ; fire_upper_TgCyr = rep(0,nos_years)
harvest_TgCyr = rep(0,nos_years) ; harvest_lower_TgCyr = rep(0,nos_years) ; harvest_upper_TgCyr = rep(0,nos_years)
# Pool totals
wood_TgC = rep(0,nos_years) ; wood_lower_TgC = rep(0,nos_years) ; wood_upper_TgC = rep(0,nos_years)
lit_TgC = rep(0,nos_years) ; lit_lower_TgC = rep(0,nos_years) ; lit_upper_TgC = rep(0,nos_years)
litwood_TgC = rep(0,nos_years) ; litwood_lower_TgC = rep(0,nos_years) ; litwood_upper_TgC = rep(0,nos_years)
soil_TgC = rep(0,nos_years) ; soil_lower_TgC = rep(0,nos_years) ; soil_upper_TgC = rep(0,nos_years)
# Flux trends
gpp_trend = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
rauto_trend = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
rhet_trend = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
lai_trend = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))

time_vector = seq(0,nos_years, length.out = dim(grid_output$nee_gCm2day)[3])

# Loop through all sites
nos_sites_inc = 0
for (n in seq(1, PROJECT$nosites)) {

     # Extract each sites location within the grid
     i_loc = grid_output$i_location[n] ; j_loc = grid_output$j_location[n]

     # Does an analysis actually exist for this location / model
     if (is.na(i_loc) == FALSE & is.na(j_loc) == FALSE & is.na(landfilter[i_loc,j_loc]) == FALSE) {
         nos_sites_inc = nos_sites_inc + 1
         # Estimate pixel level trends
         gpp_trend[i_loc,j_loc] = coef(lm(grid_output$gpp_gCm2day[n,mid_quant,] ~ time_vector))[2]     # median selected
         rauto_trend[i_loc,j_loc] = coef(lm(grid_output$rauto_gCm2day[n,mid_quant,] ~ time_vector))[2] # median selected
         rhet_trend[i_loc,j_loc] = coef(lm(grid_output$rhet_gCm2day[n,mid_quant,] ~ time_vector))[2]   # median selected
         lai_trend[i_loc,j_loc] = coef(lm(grid_output$lai_m2m2[n,mid_quant,] ~ time_vector))[2]   # median selected
         # Cumulate the total area actually used in the analysis
         cumarea = cumarea + area[i_loc,j_loc]
         lai_grid[i_loc,j_loc,] = rollapply(grid_output$lai_m2m2[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
         lai_m2m2      = lai_m2m2      + lai_grid[i_loc,j_loc,]
         lai_lower_m2m2 = lai_lower_m2m2 + rollapply(grid_output$lai_m2m2[n,low_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
         lai_upper_m2m2 = lai_upper_m2m2 + rollapply(grid_output$lai_m2m2[n,high_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
         # Stocks
         wood_TgC      = wood_TgC      + (rollapply(grid_output$wood_gCm2[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean))
         wood_lower_TgC = wood_lower_TgC + rollapply(grid_output$wood_gCm2[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
         wood_upper_TgC = wood_upper_TgC + rollapply(grid_output$wood_gCm2[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
         lit_TgC       = lit_TgC       + (rollapply(grid_output$lit_gCm2[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean))
         lit_lower_TgC = lit_lower_TgC + rollapply(grid_output$lit_gCm2[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
         lit_upper_TgC = lit_upper_TgC + rollapply(grid_output$lit_gCm2[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
         litwood_TgC   = litwood_TgC   + (rollapply(grid_output$litwood_gCm2[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean))
         litwood_lower_TgC = litwood_lower_TgC + rollapply(grid_output$litwood_gCm2[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
         litwood_upper_TgC = litwood_upper_TgC + rollapply(grid_output$litwood_gCm2[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
         soil_TgC      = soil_TgC      + (rollapply(grid_output$som_gCm2[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean))
         soil_lower_TgC = soil_lower_TgC + rollapply(grid_output$som_gCm2[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
         soil_upper_TgC = soil_upper_TgC + rollapply(grid_output$som_gCm2[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
         # Fluxes
         gpp_TgCyr     = gpp_TgCyr     + (rollapply(grid_output$gpp_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         gpp_lower_TgCyr = gpp_lower_TgCyr + (rollapply(grid_output$gpp_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         gpp_upper_TgCyr = gpp_upper_TgCyr + (rollapply(grid_output$gpp_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         rauto_TgCyr   = rauto_TgCyr   + (rollapply(grid_output$rauto_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         rauto_lower_TgCyr = rauto_lower_TgCyr + (rollapply(grid_output$rauto_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         rauto_upper_TgCyr = rauto_upper_TgCyr + (rollapply(grid_output$rauto_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)         
         rhet_TgCyr    = rhet_TgCyr    + (rollapply(grid_output$rhet_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         rhet_lower_TgCyr = rhet_lower_TgCyr + (rollapply(grid_output$rhet_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         rhet_upper_TgCyr = rhet_upper_TgCyr + (rollapply(grid_output$rhet_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         nee_TgCyr     = nee_TgCyr     + (rollapply(grid_output$nee_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         nee_lower_TgCyr = nee_lower_TgCyr + (rollapply(grid_output$nee_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         nee_upper_TgCyr = nee_upper_TgCyr + (rollapply(grid_output$nee_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         nbe_TgCyr     = nbe_TgCyr     + (rollapply(grid_output$nbe_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         nbe_lower_TgCyr = nbe_lower_TgCyr + (rollapply(grid_output$nbe_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         nbe_upper_TgCyr = nbe_upper_TgCyr + (rollapply(grid_output$nbe_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         fire_TgCyr    = fire_TgCyr    + (rollapply(grid_output$fire_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         fire_lower_TgCyr = fire_lower_TgCyr + (rollapply(grid_output$fire_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         fire_upper_TgCyr = fire_upper_TgCyr + (rollapply(grid_output$fire_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         harvest_TgCyr = harvest_TgCyr + (rollapply(grid_output$harvest_gCm2day[n,mid_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         harvest_lower_TgCyr = harvest_lower_TgCyr + (rollapply(grid_output$harvest_gCm2day[n,low_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         harvest_upper_TgCyr = harvest_upper_TgCyr + (rollapply(grid_output$harvest_gCm2day[n,high_quant,]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
     }
} # loop sites

# LAI averaging
lai_m2m2 = lai_m2m2 / nos_sites_inc
lai_lower_m2m2 = lai_lower_m2m2 / nos_sites_inc
lai_upper_m2m2 = lai_upper_m2m2 / nos_sites_inc
# Now adjust units gC/yr -> TgC/yr
# All AGB
gpp_TgCyr     = gpp_TgCyr * 1e-12
rauto_TgCyr   = rauto_TgCyr * 1e-12
rhet_TgCyr    = rhet_TgCyr * 1e-12
nee_TgCyr     = nee_TgCyr * 1e-12
nbe_TgCyr     = nbe_TgCyr * 1e-12
fire_TgCyr    = fire_TgCyr * 1e-12
harvest_TgCyr = harvest_TgCyr * 1e-12
lit_TgC       = lit_TgC * 1e-12
litwood_TgC   = litwood_TgC * 1e-12
wood_TgC      = wood_TgC * 1e-12
soil_TgC      = soil_TgC * 1e-12
# lower
gpp_lower_TgCyr     = gpp_lower_TgCyr * 1e-12
rauto_lower_TgCyr   = rauto_lower_TgCyr * 1e-12
rhet_lower_TgCyr    = rhet_lower_TgCyr * 1e-12
nee_lower_TgCyr     = nee_lower_TgCyr * 1e-12
nbe_lower_TgCyr     = nbe_lower_TgCyr * 1e-12
fire_lower_TgCyr    = fire_lower_TgCyr * 1e-12
harvest_lower_TgCyr = harvest_lower_TgCyr * 1e-12
lit_lower_TgC       = lit_lower_TgC * 1e-12
litwood_lower_TgC   = litwood_lower_TgC * 1e-12
wood_lower_TgC      = wood_lower_TgC * 1e-12
soil_lower_TgC      = soil_lower_TgC * 1e-12
# upper
gpp_upper_TgCyr     = gpp_upper_TgCyr * 1e-12
rauto_upper_TgCyr   = rauto_upper_TgCyr * 1e-12
rhet_upper_TgCyr    = rhet_upper_TgCyr * 1e-12
nee_upper_TgCyr     = nee_upper_TgCyr * 1e-12
nbe_upper_TgCyr     = nbe_upper_TgCyr * 1e-12
fire_upper_TgCyr    = fire_upper_TgCyr * 1e-12
harvest_upper_TgCyr = harvest_upper_TgCyr * 1e-12
lit_upper_TgC       = lit_upper_TgC * 1e-12
litwood_upper_TgC   = litwood_upper_TgC * 1e-12
wood_upper_TgC      = wood_upper_TgC * 1e-12
soil_upper_TgC      = soil_upper_TgC * 1e-12

###
## C - Budget (TgC/yr)

# Summary C budgets for output to table, NOTE the use of landfilter removes areas outside of the target area
dims = dim(grid_output$mean_gpp_gCm2day)
grid_output$gpp_TgCyr          = apply(grid_output$mean_gpp_gCm2day*array(landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
grid_output$rauto_TgCyr        = apply(grid_output$mean_rauto_gCm2day*array(landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
grid_output$rhet_TgCyr         = apply(grid_output$mean_rhet_gCm2day*array(landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
grid_output$nbe_TgCyr          = apply(grid_output$mean_nbe_gCm2day*array(landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
grid_output$nee_TgCyr          = apply(grid_output$mean_nee_gCm2day*array(landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
grid_output$fire_TgCyr         = apply(grid_output$mean_fire_gCm2day*array(landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
grid_output$harvest_TgCyr      = apply(grid_output$mean_harvest_gCm2day*array(landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
grid_output$labile_TgC         = apply(grid_output$mean_labile_gCm2*array(landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
grid_output$foliage_TgC        = apply(grid_output$mean_foliage_gCm2*array(landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
grid_output$roots_TgC          = apply(grid_output$mean_roots_gCm2*array(landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
grid_output$wood_TgC           = apply(grid_output$mean_wood_gCm2*array(landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
grid_output$lit_TgC            = apply(grid_output$mean_lit_gCm2*array(landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
grid_output$som_TgC            = apply(grid_output$mean_som_gCm2*array(landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
grid_parameters$fNPP           = apply(array(landfilter,dim = dims)*grid_parameters$NPP_foliar_fraction,3,mean, na.rm=TRUE)
grid_parameters$rNPP           = apply(array(landfilter,dim = dims)*grid_parameters$NPP_root_fraction,3,mean, na.rm=TRUE)
grid_parameters$wNPP           = apply(array(landfilter,dim = dims)*grid_parameters$NPP_wood_fraction,3,mean, na.rm=TRUE)
grid_parameters$foliage_MTT_yr = apply(array(landfilter,dim = dims)*grid_parameters$MTT_foliar_years,3,mean, na.rm=TRUE)
grid_parameters$roots_MTT_yr   = apply(array(landfilter,dim = dims)*grid_parameters$MTT_root_years,3,mean, na.rm=TRUE)
grid_parameters$wood_MTT_yr    = apply(array(landfilter,dim = dims)*grid_parameters$MTT_wood_years,3,mean, na.rm=TRUE)
grid_parameters$lit_MTT_yr     = apply(array(landfilter,dim = dims)*grid_parameters$MTT_DeadOrg_years,3,mean, na.rm=TRUE)
grid_parameters$som_MTT_yr     = apply(array(landfilter,dim = dims)*grid_parameters$MTT_som_years,3,mean, na.rm=TRUE)
grid_output$dCfoliage_gCm2yr   = apply(array(landfilter,dim = dims)*(grid_output$final_dCfoliage_gCm2/nos_years),3,mean, na.rm=TRUE)
grid_output$dCroots_gCm2yr     = apply(array(landfilter,dim = dims)*(grid_output$final_dCroots_gCm2/nos_years),3,mean, na.rm=TRUE)
grid_output$dCwood_gCm2yr      = apply(array(landfilter,dim = dims)*(grid_output$final_dCwood_gCm2/nos_years),3,mean, na.rm=TRUE)
grid_output$dClit_gCm2yr       = apply(array(landfilter,dim = dims)*(grid_output$final_dClit_gCm2/nos_years),3,mean, na.rm=TRUE)
grid_output$dCsom_gCm2yr       = apply(array(landfilter,dim = dims)*(grid_output$final_dCsom_gCm2/nos_years),3,mean, na.rm=TRUE)
# Combine output into dataframe
output = data.frame(Quantile = grid_output$num_quantiles, 
                    GPP_TgCyr = grid_output$gpp_TgCyr, Ra_TgCyr = grid_output$rauto_TgCyr, Rhet_TgCyr = grid_output$rhet_TgCyr, 
                    NEE_TgCyr = grid_output$nee_TgCyr, NBE_TgCyr = grid_output$nbe_TgCyr, 
                    Fire_TgCyr = grid_output$fire_TgCyr, Harvest_TgCyr = grid_output$harvest_TgCyr,
                    labile_TgC = grid_output$labile_TgC, foliage_TgC = grid_output$foliage_TgC, 
                    fine_root_TgC = grid_output$roots_TgC, wood_TgC = grid_output$wood_TgC, 
                    litter_TgC = grid_output$lit_TgC, som_TgC = grid_output$som_TgC, 
                    fNPP_fraction = grid_parameters$fNPP, rNPP_fraction = grid_parameters$rNPP, wNPP_fraction = grid_parameters$wNPP,
                    foliage_MTT_yr = grid_parameters$foliage_MTT_yr, fine_roots_MTT_yr = grid_parameters$roots_MTT_yr, 
                    wood_MTT_yr = grid_parameters$wood_MTT_yr,
                    lit_MTT_yr = grid_parameters$lit_MTT_yr, som_MTT_yr = grid_parameters$som_MTT_yr,
                    dCfoliage_gCm2yr = grid_output$dCfoliage_gCm2yr, dCroots_gCm2yr = grid_output$dCroots_gCm2yr,
                    dCwood_gCm2yr = grid_output$dCwood_gCm2yr,
                    dClit_gCm2yr = grid_output$dClit_gCm2yr, dCsom_gCm2yr = grid_output$dCsom_gCm2yr)
# Add any additional variables based on conditional statements
if (length(which(names(grid_output) == "mean_litwood_gCm2")) > 0) {
    # Create variables
    grid_output$litwood_TgC = apply(grid_output$mean_litwood_gCm2*array(landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
    grid_output$dClitwood_gCm2yr = apply(array(landfilter,dim = dims)*(grid_output$final_dClitwood_gCm2/nos_years),3,mean, na.rm=TRUE)
    # Assign to data.frame
    output$wood_litter_TgC = grid_output$litwood_TgC
    output$dClitwood_gCm2yr = grid_output$dClitwood_gCm2yr
} # wood litter is simulated

# Write out C budget
write.table(output, file = paste(out_dir,"/",PROJECT$name,"_C_budget.csv",sep=""), row.names=FALSE, sep=",",append=FALSE)

###
## Determine 1-posterior:prior ratio, i.e. how much have we learned?
###

# Extract parameter prior ranges from source code
prior_ranges = read_src_model_priors(PROJECT)

# Create ratio array
posterior_prior = array(NA, dim=c(dim(grid_parameters$parameters)[1:2],length(prior_ranges$parmin)))
for (n in seq(1, PROJECT$nosites)) {

     # Check that location has run
     if (is.na(grid_output$i_location[n]) == FALSE & is.na(grid_output$j_location[n]) == FALSE) {
         for (p in seq(1,length(prior_ranges$parmin))) {
              tmp = grid_parameters$parameters[grid_output$i_location[n],grid_output$j_location[n],p,high_quant] 
              tmp = tmp - grid_parameters$parameters[grid_output$i_location[n],grid_output$j_location[n],p,low_quant] 
              posterior_prior[grid_output$i_location[n],grid_output$j_location[n],p] = tmp / (prior_ranges$parmax[p]-prior_ranges$parmin[p])
         } # Loop parameters
     } # Does site exist

} # Loop sites

# Generate some summary statistics
#print("====1-(Posterior:Prior====")
#print("===Process parameters===")
#print(summary(apply(1-posterior_prior[,,-c(18:24,28)], 3, mean, na.rm=TRUE)))
#print("===Initial conditions parameters===")
#print(summary(apply(1-posterior_prior[,,c(18:24,28)], 3, mean, na.rm=TRUE)))
print("===All parameters===")
print(summary(apply(1-posterior_prior, 3, mean, na.rm=TRUE)))

# Generate some plots

png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_posterior_prior_reductions.png",sep=""), height = 2000, width = 3000, res = 300)
par(mfrow=c(1,1), mar=c(0.01,1.5,0.3,7),omi=c(0.01,0.1,0.01,0.1))
tmp = area
var1 = apply(1-posterior_prior,c(1,2),mean,na.rm=TRUE)
var1 = raster(vals = t(var1[,dim(var1)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
plot(var1, main="", zlim=c(0,1), col=colour_choices_default, xaxt = "n", yaxt = "n", box = FALSE, bty = "n",
     cex.lab=2, cex.main=2.0, cex.axis = 2, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1))
plot(landmask, add=TRUE)
mtext(expression('Mean posterior reduction (0-1)'), side = 2, cex = 1.6, padj = -0.25, adj = 0.5)
dev.off()

#png(file = "~/WORK/GREENHOUSE/models/CARDAMOM/SECO/figures/posterior_prior_reductions.png", height = (1667/2)*3, width = 4500, res = 300)
#par(mfrow=c(1,2), mar=c(0.5,0.3,2.8,7),omi=c(0.1,0.3,0.1,0.1))
#tmp = area
## Process parameters
#var1 = apply(1-posterior_prior[,,-c(18:24,28)],c(1,2),mean,na.rm=TRUE)
#var1 = raster(vals = t(var1[,dim(var1)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
#plot(var1, main="", zlim=c(0,1), col=colour_choices_default, xaxt = "n", yaxt = "n", box = FALSE, bty = "n",
#     cex.lab=2, cex.main=2.0, cex.axis = 2, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1))
#plot(landmask, add=TRUE)
#mtext(expression('Process parameters'), side = 2, cex = 1.6, padj = -0.15, adj = 0.5)
## Initial conditions
#var1 = apply(1-posterior_prior[,,c(18:24,28)],c(1,2),mean,na.rm=TRUE)
#var1 = raster(vals = t(var1[,dim(var1)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
#plot(var1, main="", zlim=c(0,1), col=colour_choices_default, xaxt = "n", yaxt = "n", box = FALSE, bty = "n",
#     cex.lab=2, cex.main=2.0, cex.axis = 2, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1))
#plot(landmask, add=TRUE)
#mtext(expression('Initial conditions'), side = 2, cex = 1.6, padj = -0.15, adj = 0.5)
#dev.off()

###
## Loading and processing of independent observations
###

# The overall scheme here is to load multiple datasets for each observation type (NBE, GPP, Fire)
# combine them together to provide a mean estimate and an estimate of uncertainty based on the maximum and minimum values for each case

###
## Extract CarbonTracker Europe Inversion (NEE, NBE)

# Extract CarbonTrackerEurope (2000-2017)
CTE = nc_open("/exports/csce/datastore/geos/groups/gcel/AtmosphericInversions/CarbonTrackerEurope/flux1x1_all_years.nc")
cte_nee = ncvar_get(CTE,"bio_flux_opt")  # molC/m2/s (long,lat,ensemble,date) Fire not included
cte_fire = ncvar_get(CTE,"fire_flux_imp") # molC/m2/s 
cte_area = ncvar_get(CTE,"cell_area")    # m2
cte_days_since_2000 = ncvar_get(CTE,"date")
cte_lat = ncvar_get(CTE,"latitude")
cte_long = ncvar_get(CTE,"longitude")
nc_close(CTE)

# Convert cte_nee into nbe
cte_nbe = cte_nee + cte_fire

# Adjust units
cte_nee = cte_nee * 12 * 86400 * 365.25 # gC/m2/yr
cte_fire = cte_fire * 12 * 86400 * 365.25 # gC/m2/yr
cte_nbe = cte_nbe * 12 * 86400 * 365.25 # gC/m2/yr

# Search for africa locations and slot into africa only grid for matching
# Make into CARDAMOM paired masks.
# NOTE: filter 2000-2017 and 0.025,0.25,0.5,0.75,0.975 quantiles
cte_years = c(2000:2017) 
overlap_cte = intersect(cte_years,run_years)
overlap_start = which(cte_years == overlap_cte[1])
overlap_end = which(cte_years == overlap_cte[length(overlap_cte)])
# Create data arrays now
cte_nbe_gCm2yr = array(NA, dim=c(dim(grid_output$mean_lai_m2m2)[1:2],length(overlap_cte)))
cte_nee_gCm2yr = array(NA, dim=c(dim(grid_output$mean_lai_m2m2)[1:2],length(overlap_cte)))
cte_fire_gCm2yr = array(NA, dim=c(dim(grid_output$mean_lai_m2m2)[1:2],length(overlap_cte)))
cte_m2 = array(NA, dim=dim(grid_output$mean_lai_m2m2)[1:2])
for (n in seq(1,PROJECT$nosites)) {
     if (is.na(grid_output$i_location[n]) == FALSE & is.na(grid_output$j_location[n]) == FALSE & is.na(landfilter[grid_output$i_location[n],grid_output$j_location[n]]) == FALSE) {
         output = closest2d(1,cte_lat,cte_long,grid_lat[grid_output$i_location[n],grid_output$j_location[n]],grid_long[grid_output$i_location[n],grid_output$j_location[n]],3)
         i1 = unlist(output)[1] ; j1 = unlist(output)[2]
         cte_nbe_gCm2yr[grid_output$i_location[n],grid_output$j_location[n],] = cte_nbe[i1,j1,overlap_start:overlap_end]
         cte_nee_gCm2yr[grid_output$i_location[n],grid_output$j_location[n],] = cte_nee[i1,j1,overlap_start:overlap_end]
         cte_fire_gCm2yr[grid_output$i_location[n],grid_output$j_location[n],] = cte_fire[i1,j1,overlap_start:overlap_end]
         cte_m2[grid_output$i_location[n],grid_output$j_location[n]] = cte_area[i1,j1] 
     }
}

# Ensure that the timeseries length is consistent between the observed variable and the model analysis
# This assumes that only the timesteps that overlap the model period have been read in the first place,
# so we should only be needing to add extra empty variable space.
tmp = intersect(run_years,cte_years)
if (length(tmp) != length(run_years)) {
    # How many years before the observations need to be added?
    add_beginning = cte_years[1]-run_years[1]
    # How many years after the observations
    add_afterward = run_years[length(run_years)] - cte_years[length(cte_years)]
    if (add_beginning > 0) {
        # Convert these into arrays of the correct shape but empty
        add_beginning = array(NA, dim=c(dim(cte_nbe_gCm2yr)[1:2],add_beginning))
        # Add the extra years 
        cte_nbe_gCm2yr = abind(add_beginning,cte_nbe_gCm2yr, along=3)
    } 
    if (add_afterward > 0) {
        # Convert these into arrays of the correct shape but empty
        add_afterward = array(NA, dim=c(dim(cte_nbe_gCm2yr)[1:2],add_afterward))
        # Add the extra years 
        cte_nbe_gCm2yr = abind(cte_nbe_gCm2yr,add_afterward, along=3)
    }
} # extra years needed

# Calculate domain wide mean net sink for overlap period
cte_domain_nbe_gCm2yr = mean(apply(cte_nbe_gCm2yr,c(1,2),mean, na.rm=TRUE), na.rm=TRUE) # gC/m2/yr
cte_domain_nbe_TgCyr =   sum(apply(cte_nbe_gCm2yr,c(1,2),mean, na.rm=TRUE) * cte_m2 * 1e-12, na.rm=TRUE) # TgC/yr

# Estimate the long term trend
func_lm <-function(var_in) { 
  tmp = length(which(is.na(var_in) == FALSE)) 
  if (tmp > 1) {
      return(coef(lm(var_in ~ c(1:length(var_in))))[2])
  } else {
      return(NA)
  } 
}
cte_domain_nbe_gCm2yr_trend = apply(cte_nbe_gCm2yr,c(1,2),func_lm) # gC/m2/yr

# Tidy up
rm(cte_nee,cte_fire,cte_nbe,cte_lat,cte_long)

###
## Extract CarbonTracker Ensembles
## Two ensembles estimates exist using either flask and oco2 observations

# Make a list of all available biosphere flux estimates
avail_files = list.files("/exports/csce/datastore/geos/groups/gcel/AtmosphericInversions/CarbonTrackerEurope/Gerbrand_ensemble/biofluxopt/", full.names=TRUE)

## First the flask based analysis (2009-2017)

# Select the analyses using flask data
flask_files = avail_files[grepl("flask", avail_files) == TRUE]
# Divide between the biosphere flux and fire flux
# \\ means don't consider . as a wildcard, $ means at the end of the string
flask_files_nee = flask_files[grepl("biofluxopt\\.nc$", flask_files) == TRUE]
flask_files_fire = flask_files[grepl("firefluximp\\.nc$", flask_files) == TRUE]
# Restrict the biosphere fluxes to those which we have the iposed fire emissions, 
# so that we can convert nee into nbe
check_fire_version = unlist(strsplit(x = flask_files_fire, split = "firefluximp.nc"))
tmp = 0
for (i in seq(1,length(check_fire_version))) {
     tmpp = which(grepl(check_fire_version[i],flask_files_nee))
     if (length(tmp) > 0) {
         tmp = append(tmp,tmpp)
     }
}
flask_files_nee = flask_files_nee[tmp]
rm(tmp,tmpp,check_fire_version)
if (length(flask_files_nee) != length(flask_files_fire)) {stop("Oh dear there is a problem with the number of fire and biosphere flux files for flask")}

# Read in the first NEE file to extract spatial and temporal information
flask = nc_open(flask_files_nee[1])
flask_lat = ncvar_get(flask, "latitude")
flask_long = ncvar_get(flask, "longitude")
flask_date = ncvar_get(flask, "date") 
# Read the first values (mol/m2/s)
flask_nee = ncvar_get(flask,"bio_flux_opt")
# Tidy
nc_close(flask)

# Read in the first Fire file to extract spatial information
flask = nc_open(flask_files_fire[1])
# Read the first values (mol/m2/s)
flask_fire = ncvar_get(flask,"fire_flux_imp")
# Tidy
nc_close(flask)

# Estimate step size
flask_step = abs(flask_date[1]-flask_date[2])
# Estimate the number of days in each year since 2000 (the reference point for)
create_years = c(0,2000:2020)
nos_days = 0 ; for (i in seq(2,length(create_years))) { nos_days = append(nos_days,nos_days_in_year(create_years[i])) }
# Convert all into decimal year
for (i in seq(1, length(flask_date))) {
     tmp = which(cumsum(nos_days) > flask_date[i])[1]
     flask_date[i] = create_years[tmp] + ((flask_date[i]-sum(nos_days[1:(tmp-1)]))/nos_days[tmp])
}
# Determine where we will clip the datasets in time
flask_years = floor(flask_date) ; flask_years_keep = 0
for (i in seq(1, length(unique(flask_years)))) {
     if (length(which(flask_years == unique(flask_years)[i])) > 48) {
         flask_years_keep = append(flask_years_keep, which(flask_years == unique(flask_years)[i]))
     } 
}
flask_years_keep = flask_years_keep[-1]
# Select time periods we want only
flask_years = flask_years[flask_years_keep]
overlap_flask = intersect(flask_years,run_years)
flask_nee = flask_nee[,,flask_years_keep]
flask_fire = flask_fire[,,flask_years_keep]

# Restructure to hold all ensemble members
flask_nee = array(flask_nee, dim=c(dim(flask_nee)[1:3],length(flask_files)))
flask_fire = array(flask_fire, dim=c(dim(flask_nee)[1:3],length(flask_files)))
# Loop through all files to get our full ensemble
for (i in seq(2, length(flask_files_nee))) {
     # Open the new file
     flask_bio = nc_open(flask_files_nee[i])
     flask_fir = nc_open(flask_files_fire[i])
     # Read NEE and fire (mol/m2/s)
     tmp_nee = ncvar_get(flask_bio, "bio_flux_opt")
     tmp_fire = ncvar_get(flask_fir, "fire_flux_imp")

     # Close file
     nc_close(flask_bio) ; nc_close(flask_fir)
     # Trim to desired time period
     flask_nee[,,,i] = tmp_nee[,,flask_years_keep]     
     flask_fire[,,,i] = tmp_fire[,,flask_years_keep]     
}

# Now apply units correction (mol/m2/s) to gC/m2/day
flask_nee = flask_nee * 12 * 86400
flask_fire = flask_fire * 12 * 86400
# Create NBE
flask_nbe = flask_nee + flask_fire

# Loop through each year to estimate the annual means
flask_nee_gCm2yr = array(NA, dim=c(dim(flask_nee)[1:2],length(unique(flask_years)),dim(flask_nee)[4]))
flask_nbe_gCm2yr = array(NA, dim=c(dim(flask_nee)[1:2],length(unique(flask_years)),dim(flask_nee)[4]))
flask_fire_gCm2yr = array(NA, dim=c(dim(flask_nee)[1:2],length(unique(flask_years)),dim(flask_nee)[4]))
for (i in seq(1, length(unique(flask_years)))) {
     tmp = which(flask_years == unique(flask_years)[i])
     # Average across each year and scale to annual total
     flask_nee_gCm2yr[,,i,] = apply(flask_nee[,,tmp,],c(1,2,4),mean, na.rm=TRUE) * 365.25
     flask_nbe_gCm2yr[,,i,] = apply(flask_nbe[,,tmp,],c(1,2,4),mean, na.rm=TRUE) * 365.25
     flask_fire_gCm2yr[,,i,] = apply(flask_fire[,,tmp,],c(1,2,4),mean, na.rm=TRUE) * 365.25
}
# Remove the existing output as not needed now
rm(flask_nee,flask_fire,flask_nbe)

# Update flask years to their annuals only
flask_years = unique(flask_years)

# Loop through and extract the correct pixels for the target domain
# At this stage keep the ensemble specific information
# Loop through each year to estimate the annual means
flask_cardamom_nee_gCm2yr = array(NA, dim=c(dim(grid_output$mean_lai_m2m2)[1:2],length(flask_years),dim(flask_fire_gCm2yr)[4]))
flask_cardamom_nbe_gCm2yr = array(NA, dim=c(dim(grid_output$mean_lai_m2m2)[1:2],length(flask_years),dim(flask_fire_gCm2yr)[4]))
flask_cardamom_fire_gCm2yr = array(NA, dim=c(dim(grid_output$mean_lai_m2m2)[1:2],length(flask_years),dim(flask_fire_gCm2yr)[4]))
for (n in seq(1,PROJECT$nosites)) {
     if (is.na(grid_output$i_location[n]) == FALSE & is.na(grid_output$j_location[n]) == FALSE & is.na(landfilter[grid_output$i_location[n],grid_output$j_location[n]]) == FALSE) {
         output = closest2d(1,flask_lat,flask_long,grid_lat[grid_output$i_location[n],grid_output$j_location[n]],grid_long[grid_output$i_location[n],grid_output$j_location[n]],3)
         i1 = unlist(output)[1] ; j1 = unlist(output)[2]
         flask_cardamom_nee_gCm2yr[grid_output$i_location[n],grid_output$j_location[n],,] = flask_nee_gCm2yr[i1,j1,,]
         flask_cardamom_nbe_gCm2yr[grid_output$i_location[n],grid_output$j_location[n],,] = flask_nbe_gCm2yr[i1,j1,,]
         flask_cardamom_fire_gCm2yr[grid_output$i_location[n],grid_output$j_location[n],,] = flask_fire_gCm2yr[i1,j1,,]
     }
}

# Tidy up
rm(flask_nee_gCm2yr,flask_nbe_gCm2yr,flask_fire_gCm2yr)

## Second the OCO-2 based analysis (2015-2017)

# Select the analyses using oco2 data
oco2_files = avail_files[grepl("oco2", avail_files) == TRUE]
# Divide between the biosphere flux and fire flux
# \\ means don't consider . as a wildcard, $ means at the end of the string
oco2_files_nee = oco2_files[grepl("biofluxopt\\.nc$", oco2_files) == TRUE]
oco2_files_fire = oco2_files[grepl("firefluximp\\.nc$", oco2_files) == TRUE]
# Restrict the biosphere fluxes to those which we have the imposed fire emissions, 
# so that we can convert nee into nbe
check_fire_version = unlist(strsplit(x = oco2_files_fire, split = "firefluximp.nc"))
tmp = 0
for (i in seq(1,length(check_fire_version))) {
     tmpp = which(grepl(check_fire_version[i],oco2_files_nee))
     if (length(tmp) > 0) {
         tmp = append(tmp,tmpp)
     }
}
oco2_files_nee = oco2_files_nee[tmp]
rm(tmp,tmpp,check_fire_version)
if (length(oco2_files_nee) != length(oco2_files_fire)) { stop("Oh dear there is a problem with the number of fire and biosphere flux files for oco2") }

# Read in the first file to extract spatial and temporal information
oco2 = nc_open(oco2_files[1])
oco2_lat = ncvar_get(oco2, "latitude")
oco2_long = ncvar_get(oco2, "longitude")
oco2_date = ncvar_get(oco2, "date") 
# Read the first values (mol/m2/s)
oco2_nee = ncvar_get(oco2,"bio_flux_opt")
# Tidy
nc_close(oco2)

# Read in the first Fire file to extract spatial and temporal information
oco2 = nc_open(oco2_files_fire[1])
# Read the first values (mol/m2/s)
oco2_fire = ncvar_get(oco2,"fire_flux_imp")
# Tidy
nc_close(oco2)

# Estimate step size
oco2_step = abs(oco2_date[1]-oco2_date[2])
# Estimate the number of days in each year since 2000 (the reference point for )
create_years = c(0,2000:2020)
nos_days = 0 ; for (i in seq(2,length(create_years))) { nos_days = append(nos_days,nos_days_in_year(create_years[i]))}
# Convert all into decimal year
for (i in seq(1, length(oco2_date))) {
     tmp = which(cumsum(nos_days) > oco2_date[i])[1]
     oco2_date[i] = create_years[tmp] + ((oco2_date[i]-sum(nos_days[1:(tmp-1)]))/nos_days[tmp])
}
# Determine where we will clip the datasets in time
oco2_years = floor(oco2_date) ; oco2_years_keep = 0
for (i in seq(1, length(unique(oco2_years)))) {
     if (length(which(oco2_years == unique(oco2_years)[i])) > 48) {
         oco2_years_keep = append(oco2_years_keep, which(oco2_years == unique(oco2_years)[i]))
     } 
}
oco2_years_keep = oco2_years_keep[-1]
# Select time periods we want only
oco2_years = oco2_years[oco2_years_keep]
oco2_nee = oco2_nee[,,oco2_years_keep]
oco2_fire = oco2_fire[,,oco2_years_keep]
overlap_oco2 = intersect(oco2_years,run_years)

# Restructure to hold all ensemble members
oco2_nee = array(oco2_nee, dim=c(dim(oco2_nee)[1:3],length(oco2_files)))
oco2_fire = array(oco2_fire, dim=c(dim(oco2_nee)[1:3],length(oco2_files)))
# Loop through all files to get our full ensemble
for (i in seq(2, length(oco2_files_nee))) {
     # Open the new file
     oco2_bio = nc_open(oco2_files_nee[i])
     oco2_fir = nc_open(oco2_files_fire[i])
     # Read NEE (mol/m2/s)
     tmp_nee = ncvar_get(oco2_bio, "bio_flux_opt")
     tmp_fire = ncvar_get(oco2_fir, "fire_flux_imp")
     # Close file
     nc_close(oco2_bio) ; nc_close(oco2_fir)
     # Trim to desired time period
     oco2_nee[,,,i] = tmp_nee[,,oco2_years_keep]     
     oco2_fire[,,,i] = tmp_fire[,,oco2_years_keep]     
}

# Now apply units correction (mol/m2/s) to gC/m2/day
oco2_nee = oco2_nee * 12 * 86400
oco2_fire = oco2_fire * 12 * 86400
oco2_nbe = oco2_nee + oco2_fire

# Loop through each year to estimate the annual means
oco2_nee_gCm2yr = array(NA, dim=c(dim(oco2_nee)[1:2],length(unique(oco2_years)),dim(oco2_nee)[4]))
oco2_fire_gCm2yr = array(NA, dim=c(dim(oco2_nee)[1:2],length(unique(oco2_years)),dim(oco2_nee)[4]))
oco2_nbe_gCm2yr = array(NA, dim=c(dim(oco2_nee)[1:2],length(unique(oco2_years)),dim(oco2_nee)[4]))
for (i in seq(1, length(unique(oco2_years)))) {
     tmp = which(oco2_years == unique(oco2_years)[i])
     # Average across each year and scale to annual total
     oco2_nee_gCm2yr[,,i,] = apply(oco2_nee[,,tmp,],c(1,2,4),mean, na.rm=TRUE) * 365.25
     oco2_fire_gCm2yr[,,i,] = apply(oco2_fire[,,tmp,],c(1,2,4),mean, na.rm=TRUE) * 365.25
     oco2_nbe_gCm2yr[,,i,] = apply(oco2_nbe[,,tmp,],c(1,2,4),mean, na.rm=TRUE) * 365.25
}
# Remove the existing output as not needed now
rm(oco2_nee,oco2_fire,oco2_nbe)
# Now update oco2 years to their annuals only
oco2_years = unique(oco2_years)

# Loop through and extract the correct pixels for the target domain
# At this stage keep the ensemble specific information
# Loop through each year to estimate the annual means
oco2_cardamom_nee_gCm2yr = array(NA, dim=c(dim(grid_output$mean_lai_m2m2)[1:2],length(oco2_years),dim(oco2_fire_gCm2yr)[4]))
oco2_cardamom_nbe_gCm2yr = array(NA, dim=c(dim(grid_output$mean_lai_m2m2)[1:2],length(oco2_years),dim(oco2_fire_gCm2yr)[4]))
oco2_cardamom_fire_gCm2yr = array(NA, dim=c(dim(grid_output$mean_lai_m2m2)[1:2],length(oco2_years),dim(oco2_fire_gCm2yr)[4]))
for (n in seq(1,PROJECT$nosites)) {
     if (is.na(grid_output$i_location[n]) == FALSE & is.na(grid_output$j_location[n]) == FALSE & is.na(landfilter[grid_output$i_location[n],grid_output$j_location[n]]) == FALSE) {
         output = closest2d(1,oco2_lat,oco2_long,grid_lat[grid_output$i_location[n],grid_output$j_location[n]],grid_long[grid_output$i_location[n],grid_output$j_location[n]],3)
         i1 = unlist(output)[1] ; j1 = unlist(output)[2]
         oco2_cardamom_nee_gCm2yr[grid_output$i_location[n],grid_output$j_location[n],,] = oco2_nee_gCm2yr[i1,j1,,]
         oco2_cardamom_nbe_gCm2yr[grid_output$i_location[n],grid_output$j_location[n],,] = oco2_nbe_gCm2yr[i1,j1,,]
         oco2_cardamom_fire_gCm2yr[grid_output$i_location[n],grid_output$j_location[n],,] = oco2_fire_gCm2yr[i1,j1,,]
     } # valid value exists
} # loop sites

## Combine the NBE estimates from our datasets

# How many unique years in total
obs_nbe_years = unique(c(flask_years,oco2_years))

# Define the combined timeseries datasets
obs_nbe_gCm2yr = array(NA, dim = c(dim(grid_output$mean_lai_m2m2)[1:2],length(obs_nbe_years),sum(c(dim(oco2_cardamom_nbe_gCm2yr)[4],dim(flask_cardamom_nbe_gCm2yr)[4]))))
obs_nee_gCm2yr = array(NA, dim = c(dim(grid_output$mean_lai_m2m2)[1:2],length(obs_nbe_years),sum(c(dim(oco2_cardamom_nbe_gCm2yr)[4],dim(flask_cardamom_nbe_gCm2yr)[4]))))
obs_fire_gCm2yr = array(NA, dim = c(dim(grid_output$mean_lai_m2m2)[1:2],length(obs_nbe_years),sum(c(dim(oco2_cardamom_nbe_gCm2yr)[4],dim(flask_cardamom_nbe_gCm2yr)[4]))))
for (i in seq(1,length(obs_nbe_years))) {
     # determine whether the flask dataset has any values for this year
     tmp = which(flask_years == obs_nbe_years[i])
     if (length(tmp) > 0) {
         # If there is we shall load this into the output object
         i_s = 1 ; i_e = dim(flask_cardamom_nee_gCm2yr)[4]
         obs_nbe_gCm2yr[,,i,i_s:i_e] = flask_cardamom_nbe_gCm2yr[,,tmp,]
         obs_nee_gCm2yr[,,i,i_s:i_e] = flask_cardamom_nee_gCm2yr[,,tmp,]
         obs_fire_gCm2yr[,,i,i_s:i_e] = flask_cardamom_fire_gCm2yr[,,tmp,]
     }
     # determine whether the oco2 dataset has any values for this year
     tmp = which(oco2_years == obs_nbe_years[i])
     if (length(tmp) > 0) {
         # If there is we shall load this into the output object
         i_s = 1+dim(flask_cardamom_nee_gCm2yr)[4] ; i_e = i_s - 1 + dim(oco2_cardamom_nee_gCm2yr)[4]
         obs_nbe_gCm2yr[,,i,i_s:i_e] = oco2_cardamom_nbe_gCm2yr[,,tmp,]
         obs_nee_gCm2yr[,,i,i_s:i_e] = oco2_cardamom_nee_gCm2yr[,,tmp,]
         obs_fire_gCm2yr[,,i,i_s:i_e] = oco2_cardamom_fire_gCm2yr[,,tmp,]
     }
} # loop years

# Extract the time step mean / min / max for each of these fluxes now
obs_nbe_mean_gCm2yr = apply(obs_nbe_gCm2yr,c(1,2,3),mean, na.rm=TRUE)
obs_nbe_min_gCm2yr = apply(obs_nbe_gCm2yr,c(1,2,3),min, na.rm=TRUE)
obs_nbe_max_gCm2yr = apply(obs_nbe_gCm2yr,c(1,2,3),max, na.rm=TRUE)
obs_nee_mean_gCm2yr = apply(obs_nee_gCm2yr,c(1,2,3),mean, na.rm=TRUE)
obs_nee_min_gCm2yr = apply(obs_nee_gCm2yr,c(1,2,3),min, na.rm=TRUE)
obs_nee_max_gCm2yr = apply(obs_nee_gCm2yr,c(1,2,3),max, na.rm=TRUE)
obs_fire_mean_gCm2yr = apply(obs_fire_gCm2yr,c(1,2,3),mean, na.rm=TRUE)
obs_fire_min_gCm2yr = apply(obs_fire_gCm2yr,c(1,2,3),min, na.rm=TRUE)
obs_fire_max_gCm2yr = apply(obs_fire_gCm2yr,c(1,2,3),max, na.rm=TRUE)
# Filter out the Inf values to NaN
obs_nbe_min_gCm2yr[is.infinite(obs_nbe_min_gCm2yr) == TRUE] = NA
obs_nbe_max_gCm2yr[is.infinite(obs_nbe_max_gCm2yr) == TRUE] = NA
obs_nee_min_gCm2yr[is.infinite(obs_nee_min_gCm2yr) == TRUE] = NA
obs_nee_max_gCm2yr[is.infinite(obs_nee_max_gCm2yr) == TRUE] = NA
obs_fire_min_gCm2yr[is.infinite(obs_fire_min_gCm2yr) == TRUE] = NA
obs_fire_max_gCm2yr[is.infinite(obs_fire_max_gCm2yr) == TRUE] = NA

# Ensure that the timeseries length is consistent between the observed variable and the model analysis
# This assumes that only the timesteps that overlap the model period have been read in the first place,
# so we should only be needing to add extra empty variable space.
tmp = intersect(run_years,obs_nbe_years)
if (length(tmp) != length(run_years)) {
    # How many years before the observations need to be added?
    add_beginning = obs_nbe_years[1]-run_years[1]
    # How many years after the observations
    add_afterward = run_years[length(run_years)] - obs_nbe_years[length(obs_nbe_years)]
    if (add_beginning > 0) {
        # Convert these into arrays of the correct shape but empty
        add_beginning = array(NA, dim=c(dim(obs_nbe_min_gCm2yr)[1:2],add_beginning))
        # Add the extra years 
        obs_nbe_mean_gCm2yr = abind(add_beginning,obs_nbe_mean_gCm2yr, along=3)
        obs_nbe_min_gCm2yr = abind(add_beginning,obs_nbe_min_gCm2yr, along=3)
        obs_nbe_max_gCm2yr = abind(add_beginning,obs_nbe_max_gCm2yr, along=3)
        obs_nee_mean_gCm2yr = abind(add_beginning,obs_nbe_mean_gCm2yr, along=3)
        obs_nee_min_gCm2yr = abind(add_beginning,obs_nbe_min_gCm2yr, along=3)
        obs_nee_max_gCm2yr = abind(add_beginning,obs_nbe_max_gCm2yr, along=3)
        obs_fire_mean_gCm2yr = abind(add_beginning,obs_nbe_mean_gCm2yr, along=3)
        obs_fire_min_gCm2yr = abind(add_beginning,obs_nbe_min_gCm2yr, along=3)
        obs_fire_max_gCm2yr = abind(add_beginning,obs_nbe_max_gCm2yr, along=3)
    } 
    if (add_afterward > 0) {
        # Convert these into arrays of the correct shape but empty
        add_afterward = array(NA, dim=c(dim(obs_nbe_min_gCm2yr)[1:2],add_afterward))
        # Add the extra years 
        obs_nbe_mean_gCm2yr = abind(obs_nbe_mean_gCm2yr,add_afterward, along=3)
        obs_nbe_min_gCm2yr = abind(obs_nbe_min_gCm2yr,add_afterward, along=3)
        obs_nbe_max_gCm2yr = abind(obs_nbe_max_gCm2yr,add_afterward, along=3)
        obs_nee_mean_gCm2yr = abind(obs_nbe_mean_gCm2yr,add_afterward, along=3)
        obs_nee_min_gCm2yr = abind(obs_nbe_min_gCm2yr,add_afterward, along=3)
        obs_nee_max_gCm2yr = abind(obs_nbe_max_gCm2yr, along=3)
        obs_fire_mean_gCm2yr = abind(obs_nbe_mean_gCm2yr,add_afterward, along=3)
        obs_fire_min_gCm2yr = abind(obs_nbe_min_gCm2yr,add_afterward, along=3)
        obs_fire_max_gCm2yr = abind(obs_nbe_max_gCm2yr,add_afterward, along=3)
    }
} # extra years needed

# Generate aggregate values at the domain level
obs_nbe_mean_domain_TgCyr = apply(obs_nbe_mean_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_nbe_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
obs_nbe_min_domain_TgCyr = apply(obs_nbe_min_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_nbe_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
obs_nbe_max_domain_TgCyr = apply(obs_nbe_max_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_nbe_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)

# tidy 
rm(oco2_lat,oco2_long,flask_years,oco2_years,obs_nbe_years,
   oco2_cardamom_nbe_gCm2yr,oco2_cardamom_nee_gCm2yr,oco2_cardamom_fire_gCm2yr,
   flask_cardamom_nbe_gCm2yr,flask_cardamom_nee_gCm2yr,flask_cardamom_fire_gCm2yr)

###
## Extract GPP estimates from Copernicus, FLUXCOM & MODIS

## Extract FLUXCOM GPP CRUJRA (2000-2017)

# Read first file to get additional information
fc_years = c(2000:2017)
fc_years = intersect(fc_years,run_years)

# Read in the first file to extract spatial and temporal information
FC = nc_open(paste("/exports/csce/datastore/geos/groups/gcel/GPP_ESTIMATES/FLUXCOM/CRUJRA/GPP.RS_METEO.FP-ALL.MLM-ALL.METEO-CRUJRA_v1.720_360.monthly.",fc_years[1],".nc",sep=""))
fc_gpp = ncvar_get(FC,"GPP")  # gC/m2/d (long,lat,date) Fire not included
fc_lat = ncvar_get(FC,"lat")
fc_long = ncvar_get(FC,"lon")
nc_close(FC)

# Calculate mean annual from the monthly data
fc_gpp = apply(fc_gpp,c(1,2),mean)

# Update dimensions to accomodate all years about to be read in
fc_gpp = array(fc_gpp, dim=c(dim(fc_gpp)[1:2],length(fc_years)))

# Loop years now reading in each file in turn and adding to the output variable
for (t in seq(2,length(fc_years))) {
     # Read file
     FC = nc_open(paste("/exports/csce/datastore/geos/groups/gcel/GPP_ESTIMATES/FLUXCOM/CRUJRA/GPP.RS_METEO.FP-ALL.MLM-ALL.METEO-CRUJRA_v1.720_360.monthly.",fc_years[t],".nc",sep=""))
     tmp = ncvar_get(FC,"GPP") 
     # Monthly to annual
     fc_gpp[,,t] = apply(tmp,c(1,2),mean)
     # tidy
     nc_close(FC) ; rm(tmp)
}

# Adjust units
fc_gpp = fc_gpp * 365.25 # gC/m2/yr

# Aggregate from 0.5 x 0.5 deg to 1x1
#fc_gpp_tmp = array(NA,dim=c(360,180,length(fc_years)))
#fc_gpp_mad_tmp = array(NA,dim=c(360,180,length(fc_years)))
## Adjust lat / long vectors
#fc_lat = rollapply(fc_lat,by = 2, width = 2, mean)
#fc_long = rollapply(fc_long,by = 2, width = 2, mean)
#x = 1 ; y = 1
#for (i in seq(1,dim(fc_gpp)[1],2)){
#    for (j in seq(1,dim(fc_gpp)[2],2)){
#         fc_gpp_tmp[x,y,] = apply(fc_gpp[i:(i+1),j:(j+1),],3,mean)
#         fc_gpp_mad_tmp[x,y,] = apply(fc_gpp_mad[i:(i+1),j:(j+1),],3,mean)
#         y = y + 1
#    }
#    x = x + 1 ; y = 1
#}
## Once aggregation done replace the original variables and tidy up
#fc_gpp = fc_gpp_tmp ; fc_gpp_mad = fc_gpp_mad_tmp 
#rm(fc_gpp_tmp,fc_gpp_mad_tmp)

# Loop through and extract the correct pixels for the target domain
# At this stage keep the ensemble specific information
# Loop through each year to estimate the annual means
fc_cardamom_gpp_gCm2yr = array(NA, dim=c(dim(grid_output$mean_lai_m2m2)[1:2],length(fc_years)))
for (n in seq(1,PROJECT$nosites)) {
     if (is.na(grid_output$i_location[n]) == FALSE & is.na(grid_output$j_location[n]) == FALSE & is.na(landfilter[grid_output$i_location[n],grid_output$j_location[n]]) == FALSE) {
         output = closest2d(1,fc_lat,fc_long,grid_lat[grid_output$i_location[n],grid_output$j_location[n]],grid_long[grid_output$i_location[n],grid_output$j_location[n]],3)
         i1 = unlist(output)[1] ; j1 = unlist(output)[2]
         fc_cardamom_gpp_gCm2yr[grid_output$i_location[n],grid_output$j_location[n],] = fc_gpp[i1,j1,]
     } # valid value exists
} # loop sites

## Extract Copernicus GPP (2000-2019)

# "PointsOfChange"

# Read first file to get additional information
copernicus_years = c(2000:2019)
copernicus_years = intersect(copernicus_years,run_years)

# Make a list of all available files
copernicus_files = list.files("/exports/csce/datastore/geos/groups/gcel/GPP_ESTIMATES/copernius/global_1deg/", full.names = TRUE)

# Loop years now reading in each file in turn and adding to the output variable
for (y in seq(1,length(copernicus_years))) {
     # Make a list of all files for this year
     infiles = copernicus_files[which(grepl(paste("c_gls_GDMP_",copernicus_years[y],sep=""),copernicus_files) == TRUE)]
     # Loop through all files in the year to store these too
     for (t in seq(1, length(infiles))) {
          # Open file
          copernicus = nc_open(infiles[t])
          # If this is the first year extract some spatial information
          if (y == 1 & t == 1) {
              # Read lat / long
              copernicus_lat = ncvar_get(copernicus,"lat")
              copernicus_long = ncvar_get(copernicus,"lon")         
              # Create output variable we will accumulate into 
              copernicus_gpp = array(NA, dim=c(length(copernicus_long),length(copernicus_lat),length(copernicus_years)))
          } # first year actions only
          # For first step of the year create the new years loading variable
          if (t == 1) { tmp = array(NA, dim=c(length(copernicus_long),length(copernicus_lat),length(infiles))) }
          # Read in variable
          tmp[,,t] = ncvar_get(copernicus, "GDMP")
          # Remove any missing data flags
          tmp[,,t][which(tmp[,,t] == -9999)] = NA
          # Tidy up
          nc_close(copernicus)
     } # loop time step in year
     # Monthly to annual
     copernicus_gpp[,,y] = apply(tmp,c(1,2),mean)
     # Tidy
     rm(tmp)
} # Loop copernicus years

# Adjust units
copernicus_gpp = copernicus_gpp * 365.25 #gC/m2/day -> gC/m2/yr

# Loop through and extract the correct pixels for the target domain
# At this stage keep the ensemble specific information
# Loop through each year to estimate the annual means
copernicus_cardamom_gpp_gCm2yr = array(NA, dim=c(dim(grid_output$mean_lai_m2m2)[1:2],length(copernicus_years)))
for (n in seq(1,PROJECT$nosites)) {
     if (is.na(grid_output$i_location[n]) == FALSE & is.na(grid_output$j_location[n]) == FALSE & is.na(landfilter[grid_output$i_location[n],grid_output$j_location[n]]) == FALSE) {
         output = closest2d(1,copernicus_lat,copernicus_long,grid_lat[grid_output$i_location[n],grid_output$j_location[n]],grid_long[grid_output$i_location[n],grid_output$j_location[n]],3)
         i1 = unlist(output)[1] ; j1 = unlist(output)[2]
         copernicus_cardamom_gpp_gCm2yr[grid_output$i_location[n],grid_output$j_location[n],] = copernicus_gpp[i1,j1,]
     } # valid value exists
} # loop sites



## Extract FluxSat v2 GPP (2000-2019)

# "PointsOfChange"

# Read first file to get additional information
fluxsat_years = c(2000:2019)
fluxsat_years = intersect(fluxsat_years,run_years)

# Make a list of all available files
fluxsat_files = list.files("/exports/csce/datastore/geos/groups/gcel/GPP_ESTIMATES/FluxSat/global_1deg_monthly/", full.names = TRUE)

# Loop years now reading in each file in turn and adding to the output variable
for (y in seq(1,length(fluxsat_years))) {
     # Make a list of all files for this year
     infiles = fluxsat_files[which(grepl(paste("GPP_FluxSat_daily_v2_",fluxsat_years[y],sep=""),fluxsat_files) == TRUE)]
     # Loop through all files in the year to store these too
     for (t in seq(1, length(infiles))) {
          # Open file
          copernicus = nc_open(infiles[t])
          # If this is the first year extract some spatial information
          if (y == 1 & t == 1) {
              # Read lat / long
              fluxsat_lat = ncvar_get(copernicus,"lat")
              fluxsat_long = ncvar_get(copernicus,"lon")         
              # Create output variable we will accumulate into 
              fluxsat_gpp = array(NA, dim=c(length(fluxsat_long),length(fluxsat_lat),length(fluxsat_years)))
          } # first year actions only
          # For first step of the year create the new years loading variable
          if (t == 1) { tmp = array(NA, dim=c(length(fluxsat_long),length(fluxsat_lat),length(infiles))) }
          # Read in variable
          tmp[,,t] = ncvar_get(copernicus, "GPP")
          # Remove any missing data flags
          tmp[,,t][which(tmp[,,t] == -9999)] = NA
          # Tidy up
          nc_close(copernicus)
     } # loop time step in year
     # Monthly to annual
     fluxsat_gpp[,,y] = apply(tmp,c(1,2),mean)
     # Tidy
     rm(tmp)
} # Loop copernicus years

# Adjust units
fluxsat_gpp = fluxsat_gpp * 365.25 #gC/m2/day -> gC/m2/yr

# Loop through and extract the correct pixels for the target domain
# At this stage keep the ensemble specific information
# Loop through each year to estimate the annual means
fluxsat_cardamom_gpp_gCm2yr = array(NA, dim=c(dim(grid_output$mean_lai_m2m2)[1:2],length(fluxsat_years)))
fluxsat_cardamom_gpp_gCm2yr_trend = array(NA, dim=c(dim(grid_output$mean_lai_m2m2)[1:2]))
for (n in seq(1,PROJECT$nosites)) {
     if (is.na(grid_output$i_location[n]) == FALSE & is.na(grid_output$j_location[n]) == FALSE & is.na(landfilter[grid_output$i_location[n],grid_output$j_location[n]]) == FALSE) {
         output = closest2d(1,fluxsat_lat,fluxsat_long,grid_lat[grid_output$i_location[n],grid_output$j_location[n]],grid_long[grid_output$i_location[n],grid_output$j_location[n]],3)
         i1 = unlist(output)[1] ; j1 = unlist(output)[2]
         fluxsat_cardamom_gpp_gCm2yr[grid_output$i_location[n],grid_output$j_location[n],] = fluxsat_gpp[i1,j1,]
         # Estimate the GPP trend at this time too
         if (length(which(is.na(fluxsat_gpp[i1,j1,]) == FALSE)) > 1) {
             fluxsat_cardamom_gpp_gCm2yr_trend[grid_output$i_location[n],grid_output$j_location[n]] = coef(lm(fluxsat_gpp[i1,j1,]~fluxsat_years))[2]
         } 
     } # valid value exists
} # loop sites

## Combine the GPP estimates from our datasets

# How many unique years in total
obs_gpp_years = unique(c(fc_years,copernicus_years,fluxsat_years))

# Define the combined timeseries datasets
nos_gpp_databases = 3
obs_gpp_gCm2yr = array(NA, dim = c(dim(grid_output$mean_lai_m2m2)[1:2],length(obs_gpp_years),nos_gpp_databases))
for (i in seq(1,length(obs_gpp_years))) {
     # determine whether the fluxcom dataset has any values for this year
     tmp = which(fc_years == obs_gpp_years[i])
     if (length(tmp) > 0) {
         # If there is we shall load this into the output object
         i_s = 1 ; i_e = 1
         obs_gpp_gCm2yr[,,i,i_s:i_e] = fc_cardamom_gpp_gCm2yr[,,tmp]
     }
     # determine whether the copernicus dataset has any values for this year
     tmp = which(copernicus_years == obs_gpp_years[i])
     if (length(tmp) > 0) {
         # If there is we shall load this into the output object
         i_s = 1+1 ; i_e = i_s - 1 + 1
         obs_gpp_gCm2yr[,,i,i_s:i_e] = copernicus_cardamom_gpp_gCm2yr[,,tmp]
     }
     # determine whether the fluxsatv2 dataset has any values for this year
     tmp = which(fluxsat_years == obs_gpp_years[i])
     if (length(tmp) > 0) {
         # If there is we shall load this into the output object
         i_s = 1+1+1 ; i_e = i_s - 1 + 1
         obs_gpp_gCm2yr[,,i,i_s:i_e] = fluxsat_cardamom_gpp_gCm2yr[,,tmp]
     }
} # loop years

# Extract the time step mean / min / max for each of these fluxes now
obs_gpp_mean_gCm2yr = apply(obs_gpp_gCm2yr,c(1,2,3),mean, na.rm=TRUE)
obs_gpp_min_gCm2yr = apply(obs_gpp_gCm2yr,c(1,2,3),min, na.rm=TRUE)
obs_gpp_max_gCm2yr = apply(obs_gpp_gCm2yr,c(1,2,3),max, na.rm=TRUE)
# Filter out the Inf values to NaN
obs_gpp_min_gCm2yr[is.infinite(obs_gpp_min_gCm2yr) == TRUE] = NA
obs_gpp_max_gCm2yr[is.infinite(obs_gpp_max_gCm2yr) == TRUE] = NA

# Ensure that the timeseries length is consistent between the observed variable and the model analysis
# This assumes that only the timesteps that overlap the model period have been read in the first place,
# so we should only be needing to add extra empty variable space.
tmp = intersect(run_years,obs_gpp_years)
if (length(tmp) != length(run_years)) {
    # How many years before the observations need to be added?
    add_beginning = obs_gpp_years[1]-run_years[1]
    # How many years after the observations
    add_afterward = run_years[length(run_years)] - obs_gpp_years[length(obs_gpp_years)]
    if (add_beginning > 0) {
        # Convert these into arrays of the correct shape but empty
        add_beginning = array(NA, dim=c(dim(obs_gpp_min_gCm2yr)[1:2],add_beginning))
        # Add the extra years 
        obs_gpp_mean_gCm2yr = abind(add_beginning,obs_gpp_mean_gCm2yr, along=3)
        obs_gpp_min_gCm2yr = abind(add_beginning,obs_gpp_min_gCm2yr, along=3)
        obs_gpp_max_gCm2yr = abind(add_beginning,obs_gpp_max_gCm2yr, along=3)
    } 
    if (add_afterward > 0) {
        # Convert these into arrays of the correct shape but empty
        add_afterward = array(NA, dim=c(dim(obs_gpp_min_gCm2yr)[1:2],add_afterward))
        # Add the extra years 
        obs_gpp_mean_gCm2yr = abind(obs_gpp_mean_gCm2yr,add_afterward, along=3)
        obs_gpp_min_gCm2yr = abind(obs_gpp_min_gCm2yr,add_afterward, along=3)
        obs_gpp_max_gCm2yr = abind(obs_gpp_max_gCm2yr,add_afterward, along=3)
    }
} # extra years needed

apply(copernicus_cardamom_gpp_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(copernicus_cardamom_gpp_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)

# Generate aggregate values at the domain level
obs_gpp_mean_domain_TgCyr = apply(obs_gpp_mean_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_gpp_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
obs_gpp_min_domain_TgCyr = apply(obs_gpp_min_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_gpp_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
obs_gpp_max_domain_TgCyr = apply(obs_gpp_max_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_gpp_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)

###
## Independent fire emissions estimate
## Two estimates available, gfed and gfas
###

## Extract GFEDv4.1s (2001-2016)

# Read first file to get additional information
gfed_years = c(2001:2016)
gfed_years = intersect(gfed_years,run_years)

gfed = nc_open(paste("/exports/csce/datastore/geos/groups/gcel/GFED4/1.0deg/monthly_emissions/GFED4.1s_",gfed_years[1],".nc",sep=""))
gfed_fire = ncvar_get(gfed,"emissions")  # gC/m2/month (long,lat,date) 
gfed_lat = ncvar_get(gfed,"lat")
gfed_long = ncvar_get(gfed,"lon")
nc_close(gfed)

# Aggregate to annual flux
gfed_fire = apply(gfed_fire,c(1,2),sum,na.rm=TRUE)
# correct dimension
gfed_fire = array(gfed_fire, dim=c(dim(gfed_fire)[1:2],length(gfed_years)))
for (t in seq(2,length(gfed_years))) {
     gfed = nc_open(paste("/exports/csce/datastore/geos/groups/gcel/GFED4/1.0deg/monthly_emissions/GFED4.1s_",gfed_years[t],".nc",sep=""))
     tmp = ncvar_get(gfed,"emissions")  # gC/m2/month (long,lat,date) 
     tmp = apply(tmp,c(1,2),sum,na.rm=TRUE)
     gfed_fire[,,t] = tmp
     nc_close(gfed)
}
# Search for africa locations and slot into africa only grid for matching
# Make into CARDAMOM paired masks.
gfed_cardamom_fire_gCm2yr = array(NA, dim=c(dim(grid_output$mean_lai_m2m2)[1:2],length(gfed_years)))
gfed_m2 = array(NA, dim=dim(grid_output$mean_lai_m2m2)[1:2])
for (n in seq(1,PROJECT$nosites)) {
     if (is.na(grid_output$i_location[n]) == FALSE & is.na(grid_output$j_location[n]) == FALSE & is.na(landfilter[grid_output$i_location[n],grid_output$j_location[n]]) == FALSE) {
         output = closest2d(1,gfed_lat,gfed_long,grid_lat[grid_output$i_location[n],grid_output$j_location[n]],grid_long[grid_output$i_location[n],grid_output$j_location[n]],3)
         i1 = unlist(output)[1] ; j1 = unlist(output)[2]
         gfed_cardamom_fire_gCm2yr[grid_output$i_location[n],grid_output$j_location[n],] = gfed_fire[i1,j1,]
     }
}

## Annual Emissions: (units: Teragrams (Tg) carbon)
# 2004-2017. GFAS
# Based on MODIS radiative energy emissions and converted to biomass based on PFT conversion factors

# Time period
gfas_years = c(2004:2017)
gfas_years = intersect(gfas_years,run_years)
gfas = nc_open(paste("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/cssp_brazil/GFAS_annual/cams_gfas_co2fire_",gfas_years[1],".nc",sep=""))
gfas_fire = ncvar_get(gfas, "AnnualFire")
gfas_lat = ncvar_get(gfas,"latitude")
gfas_long = ncvar_get(gfas,"longitude")

# Make space for all years
gfas_fire = array(gfas_fire, dim=c(dim(gfas_fire)[1:2],length(gfas_years)))
# Read in all the files and collect into single array
for (t in seq(2, length(gfas_years))) {
     gfas = nc_open(paste("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/cssp_brazil/GFAS_annual/cams_gfas_co2fire_",gfas_years[t],".nc",sep=""))
     tmp = ncvar_get(gfas, "AnnualFire")
     gfas_fire[,,t] = tmp
}
# Tidy
nc_close(gfas)

# Search for africa locations and slot into africa only grid for matching onto CARDAMOM
gfas_cardamom_fire_gCm2yr = array(NA, dim=c(dim(grid_output$mean_lai_m2m2)[1:2],length(gfas_years)))
for (n in seq(1,PROJECT$nosites)) {
     if (is.na(grid_output$i_location[n]) == FALSE & is.na(grid_output$j_location[n]) == FALSE & is.na(landfilter[grid_output$i_location[n],grid_output$j_location[n]]) == FALSE) {
         output = closest2d(1,gfas_lat,gfas_long,grid_lat[grid_output$i_location[n],grid_output$j_location[n]],grid_long[grid_output$i_location[n],grid_output$j_location[n]],3)
         i1 = unlist(output)[1] ; j1 = unlist(output)[2]
         gfas_cardamom_fire_gCm2yr[grid_output$i_location[n],grid_output$j_location[n],] = gfas_fire[i1,j1,]
     }
}

## Combine the Fire estimates from our datasets

# How many unique years in total
obs_fire_years = unique(c(gfed_years,gfas_years))

# Define the combined timeseries datasets
nos_fire_databases = 2
obs_fire_gCm2yr = array(NA, dim = c(dim(grid_output$mean_lai_m2m2)[1:2],length(obs_fire_years),nos_fire_databases))
for (i in seq(1,length(obs_fire_years))) {
     # determine whether the flask dataset has any values for this year
     tmp = which(gfed_years == obs_fire_years[i])
     if (length(tmp) > 0) {
         # If there is we shall load this into the output object
         i_s = 1 ; i_e = 1
         obs_fire_gCm2yr[,,i,i_s:i_e] = gfed_cardamom_fire_gCm2yr[,,tmp]
     }
     # determine whether the oco2 dataset has any values for this year
     tmp = which(gfas_years == obs_fire_years[i])
     if (length(tmp) > 0) {
         # If there is we shall load this into the output object
         i_s = 1+1 ; i_e = i_s - 1 + 1
         obs_fire_gCm2yr[,,i,i_s:i_e] = gfas_cardamom_fire_gCm2yr[,,tmp]
     }
} # loop years

# Extract the time step mean / min / max for each of these fluxes now
obs_fire_mean_gCm2yr = apply(obs_fire_gCm2yr,c(1,2,3),mean, na.rm=TRUE)
obs_fire_min_gCm2yr = apply(obs_fire_gCm2yr,c(1,2,3),min, na.rm=TRUE)
obs_fire_max_gCm2yr = apply(obs_fire_gCm2yr,c(1,2,3),max, na.rm=TRUE)
# Filter out the Inf values to NaN
obs_fire_min_gCm2yr[is.infinite(obs_fire_min_gCm2yr) == TRUE] = NA
obs_fire_max_gCm2yr[is.infinite(obs_fire_max_gCm2yr) == TRUE] = NA

# Ensure that the timeseries length is consistent between the observed variable and the model analysis
# This assumes that only the timesteps that overlap the model period have been read in the first place,
# so we should only be needing to add extra empty variable space.
tmp = intersect(run_years,obs_fire_years)
if (length(tmp) != length(run_years)) {
    # How many years before the observations need to be added?
    add_beginning = obs_fire_years[1]-run_years[1]
    # How many years after the observations
    add_afterward = run_years[length(run_years)] - obs_fire_years[length(obs_fire_years)]
    if (add_beginning > 0) {
        # Convert these into arrays of the correct shape but empty
        add_beginning = array(NA, dim=c(dim(obs_fire_min_gCm2yr)[1:2],add_beginning))
        # Add the extra years 
        obs_fire_mean_gCm2yr = abind(add_beginning,obs_fire_mean_gCm2yr, along=3)
        obs_fire_min_gCm2yr = abind(add_beginning,obs_fire_min_gCm2yr, along=3)
        obs_fire_max_gCm2yr = abind(add_beginning,obs_fire_max_gCm2yr, along=3)
    } 
    if (add_afterward > 0) {
        # Convert these into arrays of the correct shape but empty
        add_afterward = array(NA, dim=c(dim(obs_fire_min_gCm2yr)[1:2],add_afterward))
        # Add the extra years 
        obs_fire_mean_gCm2yr = abind(obs_fire_mean_gCm2yr,add_afterward, along=3)
        obs_fire_min_gCm2yr = abind(obs_fire_min_gCm2yr,add_afterward, along=3)
        obs_fire_max_gCm2yr = abind(obs_fire_max_gCm2yr,add_afterward, along=3)
    }
} # extra years needed

# Generate aggregate values at the domain level
obs_fire_mean_domain_TgCyr = apply(obs_fire_mean_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_fire_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
obs_fire_min_domain_TgCyr = apply(obs_fire_min_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_fire_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
obs_fire_max_domain_TgCyr = apply(obs_fire_max_gCm2yr*array(area, dim=c(dim(area)[1:2],dim(obs_fire_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)

# Tidy
rm(gfas_cardamom_fire_gCm2yr,gfas_long,gfas_lat,gfas_years,gfas_fire,
   gfed_cardamom_fire_gCm2yr,gfed_long,gfed_lat,gfed_years,gfed_fire)

###
## Estimate spatial consistency of DALEC analyses and independent estimates

# Determine locations where we are consistent with independent evaluation
# DALEC consistency with range described by independent estimates
nbe_sig_latitude = 0   ; nbe_sig_longitude = 0
gpp_sig_latitude = 0   ; gpp_sig_longitude = 0
fire_sig_latitude = 0  ; fire_sig_longitude = 0
# Create objects for mean CARDAMOM observational overlap 
grid_output$gpp_obs_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
grid_output$nbe_obs_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
grid_output$fire_obs_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
# Currently coded to find consistency at the 95 % CI level
nos_site_inc = 0
for (i in seq(1, PROJECT$long_dim)) {
     for (j in seq(1,PROJECT$lat_dim)) {
          if (is.na(grid_output$mean_nbe_gCm2day[i,j,high_quant]) == FALSE & is.na(landfilter[i,j]) == FALSE) {
              nos_site_inc = nos_site_inc + 1
              # Determine correct site location for extracting time period specific information
              n = which(grid_output$i_location == i & grid_output$j_location == j)
              ## Does DALEC ensemble and independent estimates overlap?

              ## Where are we consistent with GPP observations
              tmp = array(NA, dim=c(length(wanted_quant),nos_years))
              # Determine mean annual flux per quantile
              for (q in seq(1, length(wanted_quant))) {
                   tmp[q,] = rollapply(grid_output$gpp_gCm2day[n,wanted_quant[q],], by = steps_per_year, width = steps_per_year, FUN = mean)
              }
              # scale to annual value
              tmp = tmp * 365.25 ; nobs = 0 ; npdf = 0 ; grid_output$gpp_obs_overlap_fraction[i,j] = 0
              # Loop through time to assess model overlap with observations
              nobs = 0 ; grid_output$gpp_obs_overlap_fraction[i,j] = 0
              for (t in seq(1, nos_years)) {
                   if (is.na(obs_gpp_mean_gCm2yr[i,j,t]) == FALSE) {
                       if ((obs_gpp_min_gCm2yr[i,j,t] - obs_gpp_max_gCm2yr[i,j,t]) != 0 ) {
                           # Create list object containing each observations distributions
                           hist_list = list(o = c(obs_gpp_min_gCm2yr[i,j,t],obs_gpp_max_gCm2yr[i,j,t]),
                                            m = tmp[,t])
                           # Estimate average model ensemble within observated range
                           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
                           if (tmp2 > 0) {npdf = npdf + 1}
                           grid_output$gpp_obs_overlap_fraction[i,j] = grid_output$gpp_obs_overlap_fraction[i,j] + tmp2
                           nobs = nobs + 1
                       } 
                   }
              } # looping for cal period
              if (nobs > 0) {
                  grid_output$gpp_obs_overlap_fraction[i,j] = grid_output$gpp_obs_overlap_fraction[i,j] / nobs
                  npdf = npdf / nobs
              } else {
                  grid_output$gpp_obs_overlap_fraction[i,j] = 0
              }
              # Where are we consistent with Fluxcom GPP and its uncertainty
              if (npdf > 0.90) {
                  gpp_sig_latitude = append(gpp_sig_latitude,j-0.5)
                  gpp_sig_longitude = append(gpp_sig_longitude,i-0.5)
              } # GPP

              ## Where are we consistent with GPP observations
              tmp = array(NA, dim=c(length(wanted_quant),nos_years))
              # Determine mean annual flux per quantile
              for (q in seq(1, length(wanted_quant))) {
                   tmp[q,] = rollapply(grid_output$nbe_gCm2day[n,wanted_quant[q],], by = steps_per_year, width = steps_per_year, FUN = mean)
              }
              tmp = tmp * 365.25 ; nobs = 0 ; npdf = 0 ; grid_output$nbe_obs_overlap_fraction[i,j] = 0
              for (t in seq(1, nos_years)) {
                   if (is.na(obs_nbe_mean_gCm2yr[i,j,t]) == FALSE){
                       if ((obs_nbe_min_gCm2yr[i,j,t] - obs_nbe_max_gCm2yr[i,j,t]) != 0 ) {
                           # Create list object containing each observations distributions
                           hist_list = list(o = c(obs_nbe_min_gCm2yr[i,j,t],obs_nbe_max_gCm2yr[i,j,t]),
                                            m = tmp[,t])
                           # Estimate average model ensemble within observated range
                           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
                           if (tmp2 > 0) {npdf = npdf + 1}
                           grid_output$nbe_obs_overlap_fraction[i,j] = grid_output$nbe_obs_overlap_fraction[i,j] + tmp2
                           nobs = nobs + 1
                       } 
                   }
              } # looping for cal period
              if (nobs > 0) {
                  grid_output$nbe_obs_overlap_fraction[i,j] = grid_output$nbe_obs_overlap_fraction[i,j] / nobs
                  npdf = npdf / nobs
              } else {
                  grid_output$nbe_obs_overlap_fraction[i,j] = 0
              }
              # Where are we consistent with CTE NBE ensemble
              if (npdf > 0.90) {
                  nbe_sig_latitude = append(nbe_sig_latitude,j-0.5)
                  nbe_sig_longitude = append(nbe_sig_longitude,i-0.5)
              } # NEE

              ## Where are we consistent with Fire observations
              tmp = array(NA, dim=c(length(wanted_quant),nos_years))
              # Determine mean annual flux per quantile
              for (q in seq(1, length(wanted_quant))) {
                   tmp[q,] = rollapply(grid_output$fire_gCm2day[n,wanted_quant[q],], by = steps_per_year, width = steps_per_year, FUN = mean)
              }
              tmp = tmp * 365.25 ; nobs = 0 ; npdf = 0 ; grid_output$fire_obs_overlap_fraction[i,j] = 0
              for (t in seq(1, nos_years)) {
                   if (is.na(obs_fire_mean_gCm2yr[i,j,t]) == FALSE) {
                       if (obs_fire_mean_gCm2yr[i,j,t] > 0.01 & (obs_fire_min_gCm2yr[i,j,t] - obs_fire_max_gCm2yr[i,j,t]) != 0 ) {
                           # Create list object containing each observations distributions
                           hist_list = list(o = c(obs_fire_min_gCm2yr[i,j,t],obs_fire_max_gCm2yr[i,j,t]),
                                            m = tmp[,t])
                           # Estimate average model ensemble within observated range
                           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
                           if (tmp2 > 0) {npdf = npdf + 1}
                           grid_output$fire_obs_overlap_fraction[i,j] = grid_output$fire_obs_overlap_fraction[i,j] + tmp2
                           nobs = nobs + 1
                       } else if (obs_fire_mean_gCm2yr[i,j,t] <= 0.1 & tmp[,t] <= 0.1) {
                           nobs = nobs + 1 ; npdf = npdf + 1
                       } 
                   }
              } # looping for cal period
              if (nobs > 0) {
                  grid_output$fire_obs_overlap_fraction[i,j] = grid_output$fire_obs_overlap_fraction[i,j] / nobs
                  npdf = npdf / nobs
              } else {
                  grid_output$fire_obs_overlap_fraction[i,j] = 0
              }
              # Where are we consistent for fire (2003-2016)
              if (npdf > 0.90) {
                  fire_sig_latitude = append(fire_sig_latitude,j-0.5)
                  fire_sig_longitude = append(fire_sig_longitude,i-0.5)
              } # fire obs

          } # NA
     } # j 
} # i
## Remove initial value
# Uncertainy / range of independent estimate overlaps
nbe_sig_latitude    = nbe_sig_latitude[-1]
nbe_sig_longitude   = nbe_sig_longitude[-1]
gpp_sig_latitude    = gpp_sig_latitude[-1]
gpp_sig_longitude   = gpp_sig_longitude[-1]
fire_sig_latitude   = fire_sig_latitude[-1]
fire_sig_longitude  = fire_sig_longitude[-1]

# Print % of pixels that are consistent
print(paste("NBE pixel consistent =",round(100*(length(nbe_sig_latitude) / nos_site_inc), digits=3)," %",sep=" "))
print(paste("GPP pixel consistent =",round(100*(length(gpp_sig_latitude) / nos_site_inc), digits=3)," %",sep=" "))
print(paste("Fire pixel consistent =",round(100*(length(fire_sig_latitude) / nos_site_inc), digits=3)," %",sep=" "))

# Determine locations where we have confidence of net source / sink of C
nbp_sig_latitude = 0    ; nbp_sig_longitude = 0
dCwood_sig_latitude = 0 ; dCwood_sig_longitude = 0
dCsom_sig_latitude = 0  ; dCsom_sig_longitude = 0
# Loop through locations
nos_site_inc = 0
for (i in seq(1, PROJECT$long_dim)) {
     for (j in seq(1,PROJECT$lat_dim)) {
          if (is.na(grid_output$mean_nbe_gCm2day[i,j,high_quant]) == FALSE & is.na(landfilter[i,j]) == FALSE) {
              nos_site_inc = nos_site_inc + 1
              # Is NBE confidently a source or sink?
              if ((grid_output$mean_nbe_gCm2day[i,j,high_quant] > 0 & grid_output$mean_nbe_gCm2day[i,j,low_quant] > 0) | 
                  (grid_output$mean_nbe_gCm2day[i,j,high_quant] < 0 & grid_output$mean_nbe_gCm2day[i,j,low_quant] < 0)) {
                  nbp_sig_latitude = append(nbp_sig_latitude,j-0.5)
                  nbp_sig_longitude = append(nbp_sig_longitude,i-0.5)
              }
              # Is wood chance confidently a source or sink?
              if ((grid_output$mean_dCwood_gCm2[i,j,high_quant] > 0 & grid_output$mean_dCwood_gCm2[i,j,low_quant] > 0) | 
                  (grid_output$mean_dCwood_gCm2[i,j,high_quant] < 0 & grid_output$mean_dCwood_gCm2[i,j,low_quant] < 0)) {
                  dCwood_sig_latitude = append(dCwood_sig_latitude,j-0.5)
                  dCwood_sig_longitude = append(dCwood_sig_longitude,i-0.5)
              }
              # Is soil confidently a source or sink?
              if ((grid_output$mean_dCsom_gCm2[i,j,high_quant] > 0 & grid_output$mean_dCsom_gCm2[i,j,low_quant] > 0) | 
                  (grid_output$mean_dCsom_gCm2[i,j,high_quant] < 0 & grid_output$mean_dCsom_gCm2[i,j,low_quant] < 0)) {
                  dCsom_sig_latitude = append(dCsom_sig_latitude,j-0.5)
                  dCsom_sig_longitude = append(dCsom_sig_longitude,i-0.5)
              }
          }
     } # j
} # i
# remove the initial value
nbp_sig_latitude = nbp_sig_latitude[-1]
nbp_sig_longitude = nbp_sig_longitude[-1]
dCwood_sig_latitude = dCwood_sig_latitude[-1]
dCwood_sig_longitude = dCwood_sig_longitude[-1]
dCsom_sig_latitude = dCsom_sig_latitude[-1]
dCsom_sig_longitude = dCsom_sig_longitude[-1]

# Fraction of locations with significant change
# NBE
print(paste("NBE sig  = ",round((length(nbp_sig_latitude) / nos_site_inc)*100, digits=3)," %",sep=""))
# dCwood
print(paste("dCwood sig  = ",round((length(dCwood_sig_latitude) / nos_site_inc)*100, digits=3)," %",sep=""))
# dCsom
print(paste("dCsom sig  = ",round((length(dCsom_sig_latitude) / nos_site_inc)*100, digits=3)," %",sep=""))

# Statisical correlation between NBP and wood change
print(paste("NBP ~ dCwood R2 = ",round(summary(lm(as.vector(-grid_output$mean_nbe_gCm2[,,mid_quant]) ~ as.vector(grid_output$mean_dCwood_gCm2[,,mid_quant])))$adj.r.squared,digits=3),sep=""))
# Statisical correlation between NBE and som change
print(paste("NBP ~ dCsom R2  = ",round(summary(lm(as.vector(-grid_output$mean_nbe_gCm2[,,mid_quant]) ~ as.vector(grid_output$mean_dCsom_gCm2[,,mid_quant])))$adj.r.squared,digits=3),sep=""))

###
## Plot Observations

# Compare analyses against the observational constraints (LAI, Soil C prior, Cwood stock, potAGB)
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_compare_observation.png",sep=""), height = 4000, width = 4500, res = 300)
par(mfrow=c(2,2), mar=c(3,4.2,3,2), omi = c(0.35,0.4,0.1,0.1))
# Plot LAI mean annual
var1 = as.vector(LAIobs) # as.vector(LAIobs*array(landfilter,dim=dim(LAIobs))) ; var1 = var1[which(is.na(var1) == FALSE)]
var2 = as.vector(lai_grid) #; var2 = var2[which(is.na(var2) == FALSE)]
plot(var2 , var1, col=model_colours[1],
     pch=1, cex = 1.6, cex.lab=2.4, cex.axis = 2.4, cex.main=2.0, ylab="", xlab="", main="")
mtext(expression(paste('CARDAMOM',sep="")), side = 1, cex = 2.4, padj = 1.85)
mtext(expression(paste('Annual LAI (',m^2,'/',m^2,')',sep="")), side = 2, cex = 2.4, padj = -1.05)
abline(0,1, col="grey", lwd=3)
# Plot wood
plot(as.vector(1e-2*WoodCobs) ~ as.vector(1e-2*WoodC), pch=1, cex = 1.6, cex.lab=2.0, cex.axis = 2.3, cex.main=2.0, ylab="", xlab="", main="", col=model_colours[1])
mtext(expression(paste('CARDAMOM',sep="")), side = 1, cex = 2.4, padj = 1.85)
mtext(expression(paste('Wood stocks (MgC/ha)',sep="")), side = 2, cex = 2.4, padj = -1.00)
abline(0,1, col="grey", lwd=3)
# Now plot LAI time series
var3  = apply(LAIobs*array(landfilter,dim=dim(LAIobs)),3,mean,na.rm=TRUE)
var4  = lai_m2m2
zrange = range(c(var3,var4), na.rm=TRUE) * c(0.8,1.2)
plot(var3~run_years, main="", cex.lab=2.4, cex.main=2, cex.axis=2.4, ylim=zrange,
      col="black", type="l", lwd=4, ylab="", xlab="")
lines(var4~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var4~run_years, col=model_colours[1], pch=16)
legend("topleft", legend = c("Copernicus","CARDAMOM"), col = c("black",model_colours[1]), lty = c(1,2), pch=c(NA,NA), horiz = FALSE, bty = "n", cex=2.1, lwd=3, ncol = 2)
mtext(expression(paste('Year',sep="")), side = 1, cex = 2.4, padj = 1.85)
mtext(expression(paste('Analysis-wide LAI (',m^2,'/',m^2,')',sep="")), side = 2, cex = 2.4, padj = -1.05)
abline(0,1, col="grey", lwd=3)
# Now plot initial soil
plot(as.vector(1e-2*grid_parameters$parameters[,,23,mid_quant]) ~ as.vector(1e-2*SoilCPrior), pch=1, cex = 1.6, cex.lab=2.4, cex.axis = 2.4, cex.main=2.0, ylab="", xlab="", main="", col=model_colours[1])
mtext(expression(paste('CARDAMOM',sep="")), side = 1, cex = 2.4, padj = 1.85)
mtext(expression(paste('Initial soil C (MgC/ha)',sep="")), side = 2, cex = 2.4, padj = -1.00)
abline(0,1, col="grey", lwd=3)
dev.off()

# Determine whether we have any observed wood trend information, i.e. do we have more than 1 wood stock
if (length(which(is.na(WoodCobs_trend) == FALSE)) > 0) {

    png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_observed_wood_timeseries.png",sep=""), height = 2000, width = 3500, res = 300)
    par(mfrow=c(1,1), mar=c(4.2,4.7,2.8,2),omi=c(0.01,0.01,0.01,0.01))
    # Now plot LAI time series
    var2 = rollapply(apply(WoodCobs_CI*1e-2*array(landfilter,dim=dim(WoodCobs))**2,3,mean,na.rm=TRUE), FUN = mean, by = 12, width = 12, na.rm=TRUE)
    var2 = sqrt(var2)
    var3 = rollapply(apply(WoodCobs*1e-2*array(landfilter,dim=dim(WoodCobs)),3,mean,na.rm=TRUE), FUN = mean, by = 12, width = 12, na.rm=TRUE)
    var4 = rollapply(apply(WoodC*1e-2*array(landfilter,dim=dim(WoodC)),3,mean,na.rm=TRUE), FUN = mean, by = 12, width = 12, na.rm=TRUE)
    zrange = range(c(var3,var4), na.rm=TRUE) * c(0.85,1.15)
    plotCI(x = run_years, y = var3, uiw = var2, main="", cex.lab=2.4, cex.main=2, cex.axis=2.4, ylim=zrange,
          col="black", lwd=4, ylab="", xlab="")
    lines(var4~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var4~run_years, col=model_colours[1], pch=16)
    legend("topleft", legend = c("Obs","CARDAMOM"), col = c("black",model_colours[1]), lty = c(1,2), pch=c(NA,NA), 
           horiz = FALSE, bty = "n", cex=2.1, lwd=3, ncol = 2)
    mtext(expression(paste('Year',sep="")), side = 1, cex = 2.4, padj = 1.85)
    mtext(expression(paste('Wood stocks (MgC/ha)',sep="")), side = 2, cex = 2.4, padj = -1.3)
    dev.off()

    # restricted axis version
    png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_wood_trend_CI_comparison_restricted_axes_heatmap.png",sep=""), 
        height = 1200, width = 4000, res = 300)
    par(mfrow=c(1,3), mar=c(4.2,5.4,2.8,2),omi=c(0.01,0.01,0.01,0.01))
    # X~Y scatter
    yrange = c(-1,1) * quantile(abs(c(WoodCobs_trend,wood_trend)), prob=c(0.999), na.rm=TRUE) * 1e-2
    smoothScatter((WoodCobs_trend*1e-2) ~ as.vector(1e-2*wood_trend), xlim=yrange, ylim=yrange, 
         ylab = expression(paste("Obs wood trend (MgC h",a^-1,"",y^-1,")",sep="")), 
         main = " ", xlab = expression(paste("Model wood trend (MgC h",a^-1,"",y^-1,")",sep="")), 
         pch=16, cex = 1.5, cex.main=2, cex.axis = 2.2, cex.lab=2.2,
         transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), colramp=smoothScatter_colours, nrpoints = 0,
         nbin = 1500)
    abline(0,1,col="red", lwd=3) ; abline(0,0,col="grey", lwd=2) ; abline(v = 0,col="grey", lwd=2)
    # Observed wood change vs stock
    yrange = c(-1,1) * quantile((WoodCobs_trend*length(run_years)), prob=c(0.999), na.rm=TRUE) * 1e-2
    plot((1e-2*WoodCobs_trend*length(run_years)) ~ as.vector(1e-2*mean_obs_wood), ylab=expression(paste("Obs total AGB change (MgC h",a^-1,")",sep="",)), 
         xlab = expression(paste("Mean obs wood stock (MgC h",a^-1,")",sep="")), 
         ylim=yrange, xlim=c(0,max(mean_obs_wood*1e-2,na.rm=TRUE)), pch = 16, cex = 1.5, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
    abline(0,0,col="grey", lwd=2)
    # Observed wood change vs CI
    plot((1e-2*WoodCobs_trend*length(run_years))  ~ as.vector(1e-2*WoodCobs_mean_CI), 
         ylab=expression(paste("Obs total AGB change (MgC h",a^-1,")",sep="",)), 
         xlab = expression(paste("Obs mean CI (MgC h",a^-1,")",sep="")), 
         ylim=yrange, xlim=c(0,max(1e-2*WoodCobs_mean_CI, na.rm=TRUE)), pch = 16, cex = 1.5, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
    lines(c(0:max(1e-2*WoodCobs_mean_CI, na.rm=TRUE))~c(0:max(1e-2*WoodCobs_mean_CI, na.rm=TRUE)), col="red", lwd=2)
    lines(c(0:-max(1e-2*WoodCobs_mean_CI, na.rm=TRUE))~c(0:max(1e-2*WoodCobs_mean_CI, na.rm=TRUE)), col="red", lwd=2)
    abline(0,0,col="grey", lwd=2) 
    dev.off()

    # restricted axis version
    png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_wood_trend_CI_comparison_restricted_axes_heatmap_plus_maps.png",sep=""), 
        height = 2400, width = 4000, res = 300)
    par(mfrow=c(2,3), mar=c(4.2,5.4,2.8,2),omi=c(0.01,0.01,0.01,0.01))
    # X~Y scatter
    yrange = c(-1,1) * quantile(abs(c(WoodCobs_trend,wood_trend)), prob=c(0.999), na.rm=TRUE) * 1e-2
    smoothScatter((WoodCobs_trend*1e-2) ~ as.vector(1e-2*wood_trend), xlim=yrange, ylim=yrange, 
         ylab = expression(paste("Obs wood trend (MgC h",a^-1,"",y^-1,")",sep="")), 
         main = " ", xlab = expression(paste("Model wood trend (MgC h",a^-1,"",y^-1,")",sep="")), 
         pch=16, cex = 1.5, cex.main=2, cex.axis = 2.2, cex.lab=2.2,
         transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), colramp=smoothScatter_colours, nrpoints = 0,
         nbin = 1500)
    abline(0,1,col="red", lwd=3) ; abline(0,0,col="grey", lwd=2) ; abline(v = 0,col="grey", lwd=2)
    # Observed wood change vs stock
    yrange = c(-1,1) * quantile((WoodCobs_trend*length(run_years)), prob=c(0.999), na.rm=TRUE) * 1e-2
    plot((1e-2*WoodCobs_trend*length(run_years)) ~ as.vector(1e-2*mean_obs_wood), ylab=expression(paste("Obs total AGB change (MgC h",a^-1,")",sep="",)), 
         xlab = expression(paste("Mean obs wood stock (MgC h",a^-1,")",sep="")), 
         ylim=yrange, xlim=c(0,max(mean_obs_wood*1e-2,na.rm=TRUE)), pch = 16, cex = 1.5, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
    abline(0,0,col="grey", lwd=2)
    # Observed wood change vs CI
    plot((1e-2*WoodCobs_trend*length(run_years))  ~ as.vector(1e-2*WoodCobs_mean_CI), 
         ylab=expression(paste("Obs total AGB change (MgC h",a^-1,")",sep="",)), 
         xlab = expression(paste("Obs mean CI (MgC h",a^-1,")",sep="")), 
         ylim=yrange, xlim=c(0,max(1e-2*WoodCobs_mean_CI, na.rm=TRUE)), pch = 16, cex = 1.5, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
    lines(c(0:max(1e-2*WoodCobs_mean_CI, na.rm=TRUE))~c(0:max(1e-2*WoodCobs_mean_CI, na.rm=TRUE)), col="red", lwd=2)
    lines(c(0:-max(1e-2*WoodCobs_mean_CI, na.rm=TRUE))~c(0:max(1e-2*WoodCobs_mean_CI, na.rm=TRUE)), col="red", lwd=2)
    abline(0,0,col="grey", lwd=2) 
    # Calculate variables needed
    var2 = WoodCobs_trend_map*1e-2 # gCm2yr -> MgChayr total
    filter = quantile(var2, prob=c(0.025, 0.975), na.rm=TRUE) 
    var2[var2 < filter[1]] = filter[1] ; var2[var2 > filter[2]] = filter[2]
    var2 = raster(vals = t(landfilter[,dim(area)[2]:1]*var2[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    var3 = abs((WoodCobs_trend_map*length(run_years)) / apply(WoodCobs_CI,c(1,2),mean,na.rm=TRUE)) # signal:uncertainty
    filter = quantile(var3, prob=c(0.025, 0.975), na.rm=TRUE) 
    var3[var3 < filter[1]] = filter[1] ; var3[var3 > filter[2]] = filter[2]
    var4 = var3 ; var4[var4 > 1] = 1 ; var4[var4 < 1] = 0
    var3 = raster(vals = t(landfilter[,dim(area)[2]:1]*var3[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    var4 = raster(vals = t(landfilter[,dim(area)[2]:1]*var4[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    # Correct spatial area and mask
    var2 = crop(var2, landmask) ; var3 = crop(var3, landmask) ; var4 = crop(var4, landmask)
    var2 = mask(var2, landmask) ; var3 = mask(var3, landmask) ; var4 = mask(var4, landmask)
    # create axis
    zrange2 = c(-1,1) * max(abs(range(values(var2),na.rm=TRUE)), na.rm=TRUE)
    zrange3 = c(0,1) * max(abs(range(values(var3),na.rm=TRUE)), na.rm=TRUE)
    plot(var2, main="",col = colour_choices_sign, zlim=zrange2, xaxt = "n", yaxt = "n", box = FALSE, bty = "n",
         cex.lab=2.6, cex.main=2.6, cex.axis = 2, legend.width = 2.3, axes = FALSE, axis.args=list(cex.axis=2.6,hadj=0.1))
    mtext(expression(paste("Obs ",Delta,"wood (MgC h",a^-1,"y",r^-1,")",sep="")), side = 3, cex = 1.6, padj = +0.3, adj = 0.5)
    plot(landmask, add=TRUE)
    plot(var3, main="",col = colour_choices_CI, zlim=zrange3, xaxt = "n", yaxt = "n", box = FALSE, bty = "n",
         cex.lab=2.6, cex.main=2.6, cex.axis = 2, legend.width = 2.3, axes = FALSE, axis.args=list(cex.axis=2.6,hadj=0.1))
    mtext(expression(paste("Obs ",Delta,"wood:CI",sep="")), side = 3, cex = 1.6, padj = +0.3, adj = 0.5)
    plot(landmask, add=TRUE)
    plot(var4, main="",col = colour_choices_default, zlim=c(0,1), xaxt = "n", yaxt = "n", box = FALSE, bty = "n",
         cex.lab=2.6, cex.main=2.6, cex.axis = 2, legend.width = 2.3, axes = FALSE, axis.args=list(cex.axis=2.6,hadj=0.1))
    mtext(expression(paste("Signal:Noise Catagorical",sep="")), side = 3, cex = 1.6, padj = +0.3, adj = 0.5)
    plot(landmask, add=TRUE)
    dev.off()

    # restricted axis version
    png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_wood_trend_CI_comparison_heatmap.png",sep=""), height = 1200, width = 4000, res = 300)
    par(mfrow=c(1,3), mar=c(4.2,5.4,2.8,2),omi=c(0.01,0.01,0.01,0.01))
    # X~Y scatter
    yrange = c(-1,1) * quantile(abs(c(WoodCobs_trend,wood_trend)), prob=c(1), na.rm=TRUE) * 1e-2
    smoothScatter((WoodCobs_trend*1e-2) ~ as.vector(1e-2*wood_trend), xlim=yrange, ylim=yrange, 
         ylab = expression(paste("Obs wood trend (MgC h",a^-1,"",y^-1,")",sep="")), 
         main = " ", xlab = expression(paste("Model wood trend (MgC h",a^-1,"",y^-1,")",sep="")), 
         pch=16, cex = 1.5, cex.main=2, cex.axis = 2.2, cex.lab=2.2,
         transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), colramp=smoothScatter_colours, nrpoints = 0,
         nbin = 1500)
    abline(0,1,col="red", lwd=3) ; abline(0,0,col="grey", lwd=2) ; abline(v = 0,col="grey", lwd=2)
    # Observed wood change vs stock
    yrange = c(-1,1) * quantile((WoodCobs_trend*length(run_years)), prob=c(0.999), na.rm=TRUE) * 1e-2
    plot((1e-2*WoodCobs_trend*length(run_years)) ~ as.vector(1e-2*mean_obs_wood), ylab=expression(paste("Obs total AGB change (MgC h",a^-1,")",sep="",)), 
         xlab = expression(paste("Mean obs wood stock (MgC h",a^-1,")",sep="")), 
         ylim=yrange, xlim=c(0,max(mean_obs_wood*1e-2,na.rm=TRUE)), pch = 16, cex = 1.5, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
    abline(0,0,col="grey", lwd=2)
    # Observed wood change vs CI
    plot((1e-2*WoodCobs_trend*length(run_years))  ~ as.vector(1e-2*WoodCobs_mean_CI), 
         ylab=expression(paste("Obs total AGB change (MgC h",a^-1,")",sep="",)), 
         xlab = expression(paste("Obs mean CI (MgC h",a^-1,")",sep="")), 
         ylim=yrange, xlim=c(0,max(1e-2*WoodCobs_mean_CI, na.rm=TRUE)), pch = 16, cex = 1.5, cex.main=2, cex.axis = 2.2, cex.lab=2.2)
    lines(c(0:max(1e-2*WoodCobs_mean_CI, na.rm=TRUE))~c(0:max(1e-2*WoodCobs_mean_CI, na.rm=TRUE)), col="red", lwd=2)
    lines(c(0:-max(1e-2*WoodCobs_mean_CI, na.rm=TRUE))~c(0:max(1e-2*WoodCobs_mean_CI, na.rm=TRUE)), col="red", lwd=2)
    abline(0,0,col="grey", lwd=2) 
    dev.off()

} # multiple wood stocks available for trend analysis?

###
## Independent evaluation plots

# Are CARDAMOM models consistent with the range described by CTE NBE ensemble, FC GPP ensemble and GFED / GFAS Fire products
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NBE_GPP_FIRE_evaluation_stippling.png",sep=""), height = 1000, width = 4000, res = 300)
# Plot differences
par(mfrow=c(1,3), mar=c(0.05,0.9,0.05,7.2), omi = c(0.01,0.2,0.3,0.1))
var1 = raster(vals = t(landfilter[,dim(area)[2]:1]*365.25*1e-2*grid_output$mean_nbe_gCm2day[,dim(area)[2]:1,mid_quant]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = raster(vals = t(landfilter[,dim(area)[2]:1]*365.25*1e-2*grid_output$mean_gpp_gCm2day[,dim(area)[2]:1,mid_quant]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t(landfilter[,dim(area)[2]:1]*365.25*1e-2*grid_output$mean_fire_gCm2day[,dim(area)[2]:1,mid_quant]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Correct spatial area and mask
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask)
# create axis
zrange1 = c(-1,1) * max(abs(range(values(var1),na.rm=TRUE)), na.rm=TRUE)
zrange2 = c(0,max(values(var2), na.rm=TRUE))
zrange3 = c(0,max(values(var3), na.rm=TRUE))
plot(var1, main="",col = rev(colour_choices_default), zlim=zrange1, xaxt = "n", yaxt = "n", box = FALSE, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, legend.width = 2.3, axes = FALSE, axis.args=list(cex.axis=2.6,hadj=0.1))
mtext(expression(paste("NBE (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = +0.3, adj = 0.5)
points(grid_long[nbe_sig_longitude+0.5,1],grid_lat[1,nbe_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE)
plot(var2, main="",col = colour_choices_gain, zlim=zrange2, xaxt = "n", yaxt = "n",  box = FALSE, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, legend.width = 2.3, axes = FALSE, axis.args=list(cex.axis=2.6,hadj=0.1))
mtext(expression(paste("GPP (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = +0.3, adj = 0.5)
points(grid_long[gpp_sig_longitude+0.5,1],grid_lat[1,gpp_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE)
plot(var3, main="",col = (colour_choices_loss), zlim=zrange3, xaxt = "n", yaxt = "n",  box = FALSE, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, legend.width = 2.3, axes = FALSE, axis.args=list(cex.axis=2.6,hadj=0.1))
mtext(expression(paste("Fire (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = +0.3, adj = 0.5)
points(grid_long[fire_sig_longitude+0.5,1],grid_lat[1,fire_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE)
dev.off()

# Are CARDAMOM models consistent with the range described by CTE NBE ensemble, FC GPP ensemble and GFED / GFAS Fire products
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NBE_GPP_FIRE_evaluation_bias_stippling.png",sep=""), height = 1000, width = 4000, res = 300)
# Plot differences
par(mfrow=c(1,3), mar=c(0.05,0.9,0.05,7.2), omi = c(0.01,0.2,0.3,0.1))
var1 = (365.25*grid_output$mean_nbe_gCm2day[,,mid_quant]) - apply(obs_nbe_mean_gCm2yr, c(1,2), mean, na.rm=TRUE)
var1 = raster(vals = t(landfilter[,dim(area)[2]:1]*1e-2*var1[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = (365.25*grid_output$mean_gpp_gCm2day[,,mid_quant]) - apply(obs_gpp_mean_gCm2yr, c(1,2), mean, na.rm=TRUE)
var2 = raster(vals = t(landfilter[,dim(area)[2]:1]*1e-2*var2[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = (365.25*grid_output$mean_fire_gCm2day[,,mid_quant]) - apply(obs_fire_mean_gCm2yr, c(1,2), mean, na.rm=TRUE)
var3 = raster(vals = t(landfilter[,dim(area)[2]:1]*1e-2*var3[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Correct spatial area and mask
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask)
# create axis
zrange1 = c(-1,1) * max(abs(range(values(var1),na.rm=TRUE)), na.rm=TRUE)
zrange2 = c(-1,1) * max(abs(range(values(var2),na.rm=TRUE)), na.rm=TRUE)
zrange3 = c(-1,1) * max(abs(range(values(var3),na.rm=TRUE)), na.rm=TRUE)
plot(var1, main="",col = colour_choices_sign, zlim=zrange1, xaxt = "n", yaxt = "n", box = FALSE, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, legend.width = 2.3, axes = FALSE, axis.args=list(cex.axis=2.6,hadj=0.1))
mtext(expression(paste("NBE (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = +0.3, adj = 0.5)
points(grid_long[nbe_sig_longitude+0.5,1],grid_lat[1,nbe_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE)
plot(var2, main="",col = colour_choices_sign, zlim=zrange2, xaxt = "n", yaxt = "n",  box = FALSE, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, legend.width = 2.3, axes = FALSE, axis.args=list(cex.axis=2.6,hadj=0.1))
mtext(expression(paste("GPP (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = +0.3, adj = 0.5)
points(grid_long[gpp_sig_longitude+0.5,1],grid_lat[1,gpp_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE)
plot(var3, main="",col = (colour_choices_sign), zlim=zrange3, xaxt = "n", yaxt = "n",  box = FALSE, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, legend.width = 2.3, axes = FALSE, axis.args=list(cex.axis=2.6,hadj=0.1))
mtext(expression(paste("Fire (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = +0.3, adj = 0.5)
points(grid_long[fire_sig_longitude+0.5,1],grid_lat[1,fire_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE)
dev.off()

# Are CARDAMOM models consistent with the range described by CTE NBE ensemble, FC GPP ensemble and GFED / GFAS Fire products
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NBE_GPP_FIRE_evaluation_fraction_overlap.png",sep=""), height = 1000, width = 4000, res = 300)
# Plot differences
par(mfrow=c(1,3), mar=c(0.05,0.9,0.05,7.2), omi = c(0.01,0.2,0.3,0.1))
# C1
var1 = raster(vals = t(landfilter[,dim(area)[2]:1]*grid_output$nbe_obs_overlap_fraction[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = raster(vals = t(landfilter[,dim(area)[2]:1]*grid_output$gpp_obs_overlap_fraction[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t(landfilter[,dim(area)[2]:1]*grid_output$fire_obs_overlap_fraction[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Correct spatial area and mask
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask)
# create axis
zrange1 = c(0,1)
zrange2 = c(0,1)
zrange3 = c(0,1)
plot(var1, main="",col = colour_choices_gain, zlim=zrange1, xaxt = "n", yaxt = "n", box = FALSE, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, legend.width = 2.3, axes = FALSE, axis.args=list(cex.axis=2.6,hadj=0.1))
mtext(expression(paste("NBE overlap fraction",sep="")), side = 3, cex = 1.8, padj = +0.3, adj = 0.5)
plot(landmask, add=TRUE)
plot(var2, main="",col = colour_choices_gain, zlim=zrange2, xaxt = "n", yaxt = "n",  box = FALSE, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, legend.width = 2.3, axes = FALSE, axis.args=list(cex.axis=2.6,hadj=0.1))
mtext(expression(paste("GPP overlap fraction",sep="")), side = 3, cex = 1.8, padj = +0.3, adj = 0.5)
plot(landmask, add=TRUE)
plot(var3, main="",col = (colour_choices_gain), zlim=zrange3, xaxt = "n", yaxt = "n",  box = FALSE, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, legend.width = 2.3, axes = FALSE, axis.args=list(cex.axis=2.6,hadj=0.1))
mtext(expression(paste("Fire overlap fraction",sep="")), side = 3, cex = 1.8, padj = +0.3, adj = 0.5)
plot(landmask, add=TRUE)
dev.off()

# Domain wide NBE (yaxis) model (xaxis), include independent estimates
model_flags=c("CARDAMOM")
obs_flags=c("CTE","FC/Copernicus/FluxSatv2","GFEDv4.1s / GFAS")
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NBE_GPP_Fire_timeseries_comparison.png",sep=""), height=3800, width=2500, res=300)
par(mfrow=c(3,1),mai=c(0.3,0.65,0.3,0.2),omi=c(0.2,0.2,0.3,0.005))
# Now plot NBE, annual time series TgC/yr
dims = dim(cte_nbe_gCm2yr)
var1  = c(apply(cte_nbe_gCm2yr * array(cte_m2, dim=dims),c(3),sum, na.rm=TRUE) * 1e-12)
var2  = cbind(cbind(c(obs_nbe_mean_domain_TgCyr),c(obs_nbe_min_domain_TgCyr)),c(obs_nbe_max_domain_TgCyr))
var3  = nbe_TgCyr
zrange = range(c(var1,var2,var3), na.rm=TRUE)
zrange[2] = zrange[2] + 500
plot(var1~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=obs_colours[1], type="l", lwd=4, ylab="", xlab="")
plotconfidence(var2,run_years,2,obs_colours[1])
lines(var1~run_years, col=obs_colours[1], lwd=4)
lines(var3~run_years, col=model_colours[1], lwd=3, lty = 1) ; points(var3~run_years, col=model_colours[1], pch=16)
abline(0,0,col="grey", lwd=2)
legend("topleft", legend = c(obs_flags,model_flags), col = c(obs_colours[1:3],model_colours), 
       lty = c(rep(1,length(obs_flags)),rep(1,length(model_flags))), pch=rep(NA,length(c(obs_flags,model_flags))), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)
mtext(expression(paste("Net Biome Exchange (TgC y",r^-1,")",sep="")), side=2, padj=-2.65,cex=1.5)
#mtext("Year", side=1, padj=2.0,cex=1.6)

# Now plot GPP
var3  = cbind(cbind(c(obs_gpp_mean_domain_TgCyr),c(obs_gpp_min_domain_TgCyr)),c(obs_gpp_max_domain_TgCyr))
var4  = gpp_TgCyr   
zrange = range(c(var3,var4), na.rm=TRUE)*c(0.9,1.0)
plot(var4~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd = 4, ylab="", xlab="", lty = 2)
plotconfidence(var3,run_years,2,obs_colours[2])
lines(var4~run_years, col=model_colours[1], lwd = 4, lty = 1) ; points(var4~run_years, col=model_colours[1], pch=16)
#legend("bottomright", legend = c(obs_flags[-5],model_flags), col = c(obs_colours[1:4],model_colours), 
#       lty = c(rep(1,length(obs_flags[-5])),rep(2,length(model_flags))), pch=rep(NA,length(c(obs_flags[-5],model_flags))), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)
#mtext("Year", side=1, padj=2.0,cex=1.6)
mtext(expression(paste("Gross Primary Productivity (TgC y",r^-1,")",sep="")), side=2, padj=-2.65, cex=1.5)

# Now plot fire
var3  = cbind(cbind(c(obs_fire_mean_domain_TgCyr),c(obs_fire_min_domain_TgCyr)),c(obs_fire_max_domain_TgCyr))
var4  = fire_TgCyr 
zrange = range(c(var3,var4), na.rm=TRUE)*c(0.9,1.1)
plot(var4~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, lty=2, ylab="", xlab="")
plotconfidence(var3,run_years,2,obs_colours[3])
lines(var4~run_years, col=model_colours[1], lwd=4, lty = 1) ; points(var4~run_years, col=model_colours[1], pch=16)
mtext("Year", side=1, padj=2.0,cex=1.6)
mtext(expression(paste("Fire Emissions (TgC y",r^-1,")",sep="")), side=2, padj=-2.65,cex=1.5)
dev.off()

# Domain wide NBE (yaxis) model (xaxis), include independent estimates
model_flags=c("CARDAMOM")
obs_flags=c("CTE","FC/Copernicus/FluxSatv2","GFEDv4.1s / GFAS")
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NBE_GPP_Fire_timeseries_comparison_plusCI.png",sep=""), height=3800, width=2500, res=300)
par(mfrow=c(3,1),mai=c(0.3,0.65,0.3,0.2),omi=c(0.2,0.2,0.3,0.005))
# Now plot NBE, annual time series TgC/yr
dims = dim(cte_nbe_gCm2yr)
var1  = c(apply(cte_nbe_gCm2yr * array(cte_m2, dim=dims),c(3),sum, na.rm=TRUE) * 1e-12)
var2  = cbind(cbind(c(obs_nbe_mean_domain_TgCyr),c(obs_nbe_min_domain_TgCyr)),c(obs_nbe_max_domain_TgCyr))
var3  = nbe_TgCyr ; var4  = nbe_lower_TgCyr ; var5  = nbe_upper_TgCyr
zrange = range(c(var1,var2,var3,var4,var5), na.rm=TRUE)
zrange[2] = zrange[2] + 500
plot(var3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, ylab="", xlab="", lty=1)
plotconfidence(var2,run_years,2,obs_colours[1])
lines(var3~run_years, col=model_colours[1], lwd=3, lty = 1) ; points(var3~run_years, col=model_colours[1], pch=16)
lines(var4~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var4~run_years, col=model_colours[1], pch=16)
lines(var5~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var5~run_years, col=model_colours[1], pch=16)
abline(0,0,col="grey", lwd=2)
legend("topleft", legend = c(obs_flags,model_flags), col = c(obs_colours[1:3],model_colours), 
       lty = c(rep(1,length(obs_flags)),rep(1,length(model_flags))), pch=rep(NA,length(c(obs_flags,model_flags))), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)
mtext(expression(paste("Net Biome Exchange (TgC y",r^-1,")",sep="")), side=2, padj=-2.65,cex=1.5)
#mtext("Year", side=1, padj=2.0,cex=1.6)

# Now plot GPP
var3  = cbind(cbind(c(obs_gpp_mean_domain_TgCyr),c(obs_gpp_min_domain_TgCyr)),c(obs_gpp_max_domain_TgCyr))
var4  = gpp_TgCyr ; var5  = gpp_lower_TgCyr ; var6  = gpp_upper_TgCyr   
zrange = range(c(var3,var4,var5,var6), na.rm=TRUE)*c(0.9,1.0)
plot(var4~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd = 4, ylab="", xlab="", lty = 2)
plotconfidence(var3,run_years,2,obs_colours[2])
lines(var4~run_years, col=model_colours[1], lwd = 4, lty = 1) ; points(var4~run_years, col=model_colours[1], pch=16)
lines(var5~run_years, col=model_colours[1], lwd = 4, lty = 2) ; points(var5~run_years, col=model_colours[1], pch=16)
lines(var6~run_years, col=model_colours[1], lwd = 4, lty = 2) ; points(var6~run_years, col=model_colours[1], pch=16)
#legend("bottomright", legend = c(obs_flags[-5],model_flags), col = c(obs_colours[1:4],model_colours), 
#       lty = c(rep(1,length(obs_flags[-5])),rep(2,length(model_flags))), pch=rep(NA,length(c(obs_flags[-5],model_flags))), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)
#mtext("Year", side=1, padj=2.0,cex=1.6)
mtext(expression(paste("Gross Primary Productivity (TgC y",r^-1,")",sep="")), side=2, padj=-2.65, cex=1.5)

# Now plot fire
var3  = cbind(cbind(c(obs_fire_mean_domain_TgCyr),c(obs_fire_min_domain_TgCyr)),c(obs_fire_max_domain_TgCyr))
var4  = fire_TgCyr  ; var5  = fire_lower_TgCyr ; var6  = fire_upper_TgCyr
zrange = range(c(var3,var4,var5,var6), na.rm=TRUE)*c(0.9,1.1)
plot(var4~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, lty=2, ylab="", xlab="")
plotconfidence(var3,run_years,2,obs_colours[3])
lines(var4~run_years, col=model_colours[1], lwd=4, lty = 1) ; points(var4~run_years, col=model_colours[1], pch=16)
lines(var5~run_years, col=model_colours[1], lwd=4, lty = 2) ; points(var5~run_years, col=model_colours[1], pch=16)
lines(var6~run_years, col=model_colours[1], lwd=4, lty = 2) ; points(var5~run_years, col=model_colours[1], pch=16)
mtext("Year", side=1, padj=2.0,cex=1.6)
mtext(expression(paste("Fire Emissions (TgC y",r^-1,")",sep="")), side=2, padj=-2.65,cex=1.5)
dev.off()

###
## Statistical significance / trend maps for C-budget terms

# Comparison of NBE, wood and soil stock change over the analysis period by model
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NBP_dCwood_dCsom.png",sep=""), height = 1000, width = 3500, res = 300)
# Plot differences
par(mfrow=c(1,3), mar=c(0.05,1,0.05,7.0), omi = c(0.01,0.4,0.3,0.05))
# Create raster
var1 = raster(vals = t(365.25*-grid_output$mean_nbe_gCm2day[,dim(area)[2]:1,mid_quant]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = raster(vals = t((1/nos_years)*grid_output$final_dCwood_gCm2[,dim(area)[2]:1,mid_quant]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((1/nos_years)*grid_output$final_dCsom_gCm2[,dim(area)[2]:1,mid_quant]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Crop to size
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask)
## Restrict parameter range to +/- 800 gC/m2/yr
#var1[var1 > 800] = 800 ; var1[var1 < -800] = -800
#var2[var2 > 800] = 800 ; var2[var2 < -800] = -800
#var3[var3 > 800] = 800 ; var3[var3 < -800] = -800
# Convert Units gC/m2/yr -> MgC/ha/yr
var1 = var1 * 1e-2   ; var2 = var2 * 1e-2 ; var3 = var3 * 1e-2 
tmp = c(maxValue(var1),minValue(var1),maxValue(var2),minValue(var2))
tmp1 = c(maxValue(var3),minValue(var3))
zrange = max(abs(range(tmp, na.rm=TRUE))) * c(-1,1)
zrange1 = max(abs(c(range(tmp1, na.rm=TRUE)))) * c(-1,1)
# C1 Mean annual NBP, dCwood, dCsom
plot(var1, ylab="", xlab="", main="",  box = FALSE, bty = "n",
     xaxt = "n", yaxt = "n", zlim=zrange,
     col=colour_choices_sign, cex.lab=2, cex.main=2.2, cex.axis = 2, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1))
plot(landmask, add=TRUE)
points(grid_long[nbp_sig_longitude+0.5,1],grid_lat[1,nbp_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
mtext(expression(paste("NBP (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.7, padj = +0.5, adj = 0.5)
plot(var2, ylab="", xlab="", main="",  box = FALSE, bty = "n",
     xaxt = "n", yaxt = "n", zlim=zrange,
     col=colour_choices_sign, cex.lab=2, cex.main=2.2, cex.axis = 2, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1))
plot(landmask, add=TRUE)
points(grid_long[dCwood_sig_longitude+0.5,1],grid_lat[1,dCwood_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
mtext(expression(paste("Wood Change (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.7, padj = +0.5, adj = 0.5)
plot(var3, ylab="", xlab="", main="",  box = FALSE, bty = "n",
     xaxt = "n", yaxt = "n", zlim=zrange1,
     col=colour_choices_sign, cex.lab=2, cex.main=2.2, cex.axis = 2, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1))
plot(landmask, add=TRUE)
points(grid_long[dCsom_sig_longitude+0.5,1],grid_lat[1,dCsom_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
mtext(expression(paste("Soil Change (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.7, padj = +0.5, adj = 0.5)
dev.off()

# GPP, Rauto, Rhet trend
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_GPP_Rauto_Rhet_LAI_trend.png",sep=""), height = 2800, width = 3300, res = 300)
# Plot differences
par(mfrow=c(2,2), mar=c(0.1,0.1,1.0,8),omi=c(0.05,0.1,0.2,0.15))
# C1 - GPP. Rauto, Rhet and LAI trends
var1 = raster(vals = t(1e-2*365.25*gpp_trend[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = raster(vals = t(1e-2*365.25*rauto_trend[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t(1e-2*365.25*rhet_trend[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var16 = raster(vals = t(lai_trend[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Crop to size
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask)
var16 = crop(var16, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask)
var16 = mask(var16, landmask) 
# Determine gradients
zrange = max(abs(quantile(c(values(var1),values(var2),values(var3)), prob = c(0.001,0.999), na.rm=TRUE))) * c(-1,1)
zrange1 = max(abs(quantile(values(var16), prob = c(0.001,0.999), na.rm=TRUE))) * c(-1,1)
# Restrict to sensible bounds
#var1[var1 < zrange[1]] = zrange[1] ; var1[var1 > zrange[2]] = zrange[2]
#var2[var2 < zrange[1]] = zrange[1] ; var2[var2 > zrange[2]] = zrange[2]
#var3[var3 < zrange[1]] = zrange[1] ; var3[var3 > zrange[2]] = zrange[2]
# C1 GPP, Rauto, Rhet trend (MgC/ha/yr2)
plot(var1, ylab="", xlab="", main="", box = FALSE, bty = "n",
     xaxt = "n", yaxt = "n", zlim=zrange,
     col=colour_choices_sign, cex.lab=2, cex.main=2.2, cex.axis = 2, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1))
plot(landmask, add=TRUE)
mtext(expression(paste('GPP Trend (MgC h',a^-1,'y',r^-2,')',sep="")), cex=1.8, padj = +0.25)
plot(var2, ylab="", xlab="", main="",  box = FALSE, bty = "n",
     xaxt = "n", yaxt = "n", zlim=zrange,
     col=colour_choices_sign, cex.lab=2, cex.main=2.2, cex.axis = 2, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1))
plot(landmask, add=TRUE)
mtext(expression(paste(R[auto],' Trend (MgC h',a^-1,'y',r^-2,')',sep="")), cex=1.8, padj = +0.25)
plot(var3, ylab="", xlab="", main="",  box = FALSE, bty = "n",
           xaxt = "n", yaxt = "n", zlim=zrange,
           col=colour_choices_sign, cex.lab=2, cex.main=2.2, cex.axis = 2, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1))
plot(landmask, add=TRUE)
mtext(expression(paste(R[het],' Trend (MgC h',a^-1,'y',r^-2,')',sep="")), cex=1.8, padj = +0.25)
plot(var16, ylab="", xlab="", main="",  box = FALSE, bty = "n",
     xaxt = "n", yaxt = "n", zlim=zrange1,
     col=colour_choices_sign, cex.lab=2, cex.main=2.2, cex.axis = 2, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1))
plot(landmask, add=TRUE)
mtext(expression(paste('LAI Trend (',m^2,'',m^-2,' y',r^-2,')',sep="")), cex=1.8, padj = +0.25)
dev.off()

###
## Plot carbon fluxes / uncertainty

# C fluxes
# Assign variables
var1 = grid_output$mean_nbe_gCm2day[,,mid_quant]*1e-2*365.25
var2 = grid_output$mean_gpp_gCm2day[,,mid_quant]*1e-2*365.25
var3 = grid_output$mean_reco_gCm2day[,,mid_quant]*1e-2*365.25 
var4 = grid_output$mean_fire_gCm2day[,,mid_quant]*1e-2*365.25 
var5 = (grid_output$mean_nbe_gCm2day[,,high_quant]-grid_output$mean_nbe_gCm2day[,,low_quant])*1e-2*365.25
var6 = (grid_output$mean_gpp_gCm2day[,,high_quant]-grid_output$mean_gpp_gCm2day[,,low_quant])*1e-2*365.25
var7 = (grid_output$mean_reco_gCm2day[,,high_quant]-grid_output$mean_reco_gCm2day[,,low_quant])*1e-2*365.25
var8 = (grid_output$mean_fire_gCm2day[,,high_quant]-grid_output$mean_fire_gCm2day[,,low_quant])*1e-2*365.25
# Apply filter
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var4[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var5[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var6[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var7[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var8[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Determine ranges
zrange1 = c(-1,1)*max(abs(range(values(var1),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(values(var2),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(values(var3),na.rm=TRUE)))
zrange4 = c(0,1)*max(abs(range(values(var4),na.rm=TRUE)))
zrange5 = c(0,1)*max(abs(range(values(var5),na.rm=TRUE)))
zrange6 = c(0,1)*max(abs(range(values(var6),na.rm=TRUE)))
zrange7 = c(0,1)*max(abs(range(values(var7),na.rm=TRUE)))
zrange8 = c(0,1)*max(abs(range(values(var8),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_C_fluxes_median_CI.png",sep=""), height = 2100, width = 5000, res = 300)
par(mfrow=c(2,4), mar=c(0.5,0.5,2.8,7),omi=c(0.1,0.4,0.2,0.2))
# Mean annual median estimates
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste('NBE (MgC h',a^-1,' y',r^-1,')',sep="")), col=rev(colour_choices_default))
plot(landmask, add=TRUE)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste('GPP (MgC h',a^-1,' y',r^-1,')',sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste('Reco (MgC h',a^-1,' y',r^-1,')',sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var4, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste('Fire (MgC h',a^-1,' y',r^-1,')',sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
# Mean annual estimates uncertainty
plot(var5, zlim=zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste('NBE CI (MgC h',a^-1,' y',r^-1,')',sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste('GPP CI (MgC h',a^-1,' y',r^-1,')',sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var7, zlim=zrange7, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste('Reco CI (MgC h',a^-1,' y',r^-1,')',sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var8, zlim=zrange8, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste('Fire CI (MgC h',a^-1,' y',r^-1,')',sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
dev.off()

# C fluxes
# Assign variables
var1 = grid_output$mean_nbe_gCm2day[,,mid_quant]*1e-2*365.25
var2 = grid_output$mean_gpp_gCm2day[,,mid_quant]*1e-2*365.25
var3 = grid_output$mean_reco_gCm2day[,,mid_quant]*1e-2*365.25 
var4 = grid_output$mean_fire_gCm2day[,,mid_quant]*1e-2*365.25 
var5 = (grid_output$mean_nbe_gCm2day[,,high_quant]-grid_output$mean_nbe_gCm2day[,,low_quant])*1e-2*365.25
var6 = (grid_output$mean_gpp_gCm2day[,,high_quant]-grid_output$mean_gpp_gCm2day[,,low_quant])*1e-2*365.25
var7 = (grid_output$mean_reco_gCm2day[,,high_quant]-grid_output$mean_reco_gCm2day[,,low_quant])*1e-2*365.25
var8 = (grid_output$mean_fire_gCm2day[,,high_quant]-grid_output$mean_fire_gCm2day[,,low_quant])*1e-2*365.25
# Apply filter
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var4[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var5[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var6[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var7[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var8[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = raster(vals = t((var7)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = raster(vals = t((var8)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Determine ranges
zrange1 = c(-1,1)*max(abs(range(values(var1),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(values(var2),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(values(var3),na.rm=TRUE)))
zrange4 = c(0,1)*max(abs(range(values(var4),na.rm=TRUE)))
zrange5 = c(0,1)*max(abs(range(c(values(var5),values(var6),values(var7),values(var8)),na.rm=TRUE)))
zrange6 = zrange5
zrange7 = zrange5
zrange8 = zrange5
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_C_fluxes_median_CI_axis_matched.png",sep=""), height = 2100, width = 5000, res = 300)
par(mfrow=c(2,4), mar=c(0.5,0.5,2.8,7),omi=c(0.1,0.4,0.2,0.2))
# Mean annual median estimates
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste('NBE (MgC h',a^-1,' y',r^-1,')',sep="")), col=rev(colour_choices_default))
plot(landmask, add=TRUE)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste('GPP (MgC h',a^-1,' y',r^-1,')',sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste('Reco (MgC h',a^-1,' y',r^-1,')',sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var4, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste('Fire (MgC h',a^-1,' y',r^-1,')',sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
# Mean annual estimates uncertainty
plot(var5, zlim=zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste('NBE CI (MgC h',a^-1,' y',r^-1,')',sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste('GPP CI (MgC h',a^-1,' y',r^-1,')',sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var7, zlim=zrange7, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste('Reco CI (MgC h',a^-1,' y',r^-1,')',sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var8, zlim=zrange8, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste('Fire CI (MgC h',a^-1,' y',r^-1,')',sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
dev.off()

###
## Plot the final C stocks, change and uncertainty

# Final stocks
# Assign variables
var1 = grid_output$final_totalC_gCm2[,,mid_quant]*1e-2
var2 = grid_output$final_biomass_gCm2[,,mid_quant]*1e-2 
var3 = grid_output$final_dom_gCm2[,,mid_quant]*1e-2 
var4 = (grid_output$final_totalC_gCm2[,,high_quant]-grid_output$final_totalC_gCm2[,,low_quant])*1e-2 
var5 = (grid_output$final_biomass_gCm2[,,high_quant]-grid_output$final_biomass_gCm2[,,low_quant])*1e-2  
var6 = (grid_output$final_dom_gCm2[,,high_quant]-grid_output$final_dom_gCm2[,,low_quant])*1e-2  
# Apply filter
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var4[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var5[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var6[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(values(var1),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(values(var2),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(values(var3),na.rm=TRUE)))
zrange4 = c(0,1)*max(abs(range(c(values(var4),values(var5),values(var6)),na.rm=TRUE)))
zrange5 = zrange4
zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_Final_stock_median_CI.png",sep=""), height = 2700, width = 4900, res = 300)
par(mfrow=c(2,3), mar=c(0.6,0.4,2.9,7),omi=c(0.1,0.4,0.18,0.2))
# Final C stocks, median estimate
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Total (MgC h",a^-1,")",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Biomass (MgC h",a^-1,")",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("DOM (MgC h",a^-1,")",sep="")), col=colour_choices_gain)
plot(landmask, add=TRUE)
# Final C stocks, confidence interval
plot(var4, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Total CI (MgC h",a^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var5, zlim=zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Biomass CI (MgC h",a^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("DOM CI (MgC h",a^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
dev.off()

# Change stocks
# Assign variables
var1 = grid_output$final_dCtotalC_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var2 = grid_output$final_dCbio_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var3 = grid_output$final_dCdom_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var4 = (grid_output$final_dCtotalC_gCm2[,,high_quant]-grid_output$final_dCtotalC_gCm2[,,low_quant])*1e-2*(1/nos_years)
var5 = (grid_output$final_dCbio_gCm2[,,high_quant]-grid_output$final_dCbio_gCm2[,,low_quant])*1e-2*(1/nos_years)
var6 = (grid_output$final_dCdom_gCm2[,,high_quant]-grid_output$final_dCdom_gCm2[,,low_quant])*1e-2*(1/nos_years)
# Apply filter
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var4[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var5[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var6[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(-1,1)*max(abs(range(values(var1),na.rm=TRUE)))
zrange2 = c(-1,1)*max(abs(range(values(var2),na.rm=TRUE)))
zrange3 = c(-1,1)*max(abs(range(values(var3),na.rm=TRUE)))
zrange4 = c(0,1)*max(abs(range(c(values(var4),values(var5),values(var6)),na.rm=TRUE)))
zrange5 = zrange4
zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_Final_stock_change_median_CI.png",sep=""), height = 2700, width = 4900, res = 300)
par(mfrow=c(2,3), mar=c(0.5,0.4,2.8,7),omi=c(0.1,0.4,0.2,0.2))
# Final stock changes, median estimates
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste(Delta,"Total (MgC h",a^-1,"y",r^-1,")",sep="")), col=(colour_choices_default))
plot(landmask, add=TRUE)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste(Delta,"Biomass (MgC h",a^-1,"y",r^-1,")",sep="")), col=(colour_choices_default))
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste(Delta,"DOM (MgC h",a^-1,"y",r^-1,")",sep="")), col=(colour_choices_default))
plot(landmask, add=TRUE)
# Final stock changes, confidence interval
plot(var4, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste(Delta,"Total CI (MgC h",a^-1,"y",r^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var5, zlim=zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste(Delta,"Biomass CI (MgC h",a^-1,"y",r^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste(Delta,"DOM CI (MgC h",a^-1,"y",r^-1,")",sep="")), col=colour_choices_CI)
plot(landmask, add=TRUE)
dev.off()

###
## Plot the MRTwood and NPPwood

# Traits
# Assign variables
wood_mrt_limit = 60
print(paste("Wood MRT plotting range has been limited to ",wood_mrt_limit," years",sep=""))
var1 = pmin(wood_mrt_limit,grid_parameters$MTT_wood_years[,,mid_quant])
var1 = array(var1, dim = dim(grid_parameters$MTT_wood_years)[1:2])
var2 = grid_parameters$NPP_wood_fraction[,,mid_quant]
var3 = grid_parameters$SS_wood_gCm2[,,mid_quant]*1e-2
var4 = pmin(wood_mrt_limit*2,grid_parameters$MTT_wood_years[,,high_quant]) - pmin(wood_mrt_limit*2,grid_parameters$MTT_wood_years[,,low_quant])
var4 = array(var4, dim = dim(grid_parameters$MTT_wood_years)[1:2])
var5 = grid_parameters$NPP_wood_fraction[,,high_quant] - grid_parameters$NPP_wood_fraction[,,low_quant]
var6 = (grid_parameters$SS_wood_gCm2[,,high_quant]*1e-2) - (grid_parameters$SS_wood_gCm2[,,low_quant]*1e-2)
# Apply filter
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var4[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var5[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var6[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# ranges
zrange1 = c(0,1)*max(abs(range(values(var1),na.rm=TRUE)))
zrange2 = c(0,1)
zrange3 = c(0,1)*max(abs(range(values(var3),na.rm=TRUE)))
zrange4 = c(0,1)*max(abs(range(values(var4),na.rm=TRUE)))
zrange5 = c(0,1)
zrange6 = c(0,1)*max(abs(range(values(var6),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NPP_MRT_SS_median_CI.png",sep=""), height = 2700, width = 4900, res = 300)
par(mfrow=c(2,3), mar=c(0.5,0.3,2.8,8),omi=c(0.1,0.3,0.2,0.2))
# Ecosystem traits, median estimates
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "MRT Wood (years)", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var2, zlim=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "NPP wood (0-1)", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "SS wood (MgC/ha)", col=colour_choices_gain)
plot(landmask, add=TRUE)
# Ecosystem traits, confidence intervals
plot(var4, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "MRT Wood CI (years)", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var5, zlim=zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "NPP wood CI (0-1)", col=colour_choices_gain)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "SS wood CI (MgC/ha)", col=colour_choices_gain)
plot(landmask, add=TRUE)
dev.off()

###
## Partition the importance of disturbance on residence times

# Estimate the proportion of turnover determined by natural
wood_turn = (grid_parameters$MTT_wood_years[,,mid_quant]**-1)
var1 = grid_parameters$MTTnatural_wood_years[,,mid_quant]**-1 
var2 = grid_parameters$MTTfire_wood_years[,,mid_quant]**-1
var3 = grid_parameters$MTTharvest_wood_years[,,mid_quant]**-1
# Now make proportional
var1 = var1 / wood_turn ; var2 = var2 / wood_turn ; var3 = var3 / wood_turn
# Filter for the miombo AGB map locations
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Filter for the miombo AGB map locations
var1[which(landfilter == 1 & is.na(var1) == TRUE)] = 0
var2[which(landfilter == 1 & is.na(var2) == TRUE)] = 0
var3[which(landfilter == 1 & is.na(var3) == TRUE)] = 0
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# specify ranges
zrange1 = c(0,1)
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_Wood_turnover_contribution.png",sep=""), height = 1300, width = 4900, res = 300)
par(mfrow=c(1,3), mar=c(0.5,0.4,3.0,7),omi=c(0.1,0.3,0.1,0.2))
# Partitioning of wood turnover, median estimate
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "Natural MRT comp (0-1)", col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var2, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "Fire MRT comp (0-1)", col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = "Biomass removal MRT comp (0-1)", col=colour_choices_loss)
plot(landmask, add=TRUE)
dev.off()

# Estimate the proportion of turnover determined by natural
wood_turn = (grid_parameters$MTT_wood_years[,,mid_quant]**-1)
var1 = grid_parameters$MTTnatural_wood_years[,,mid_quant]**-1 
var2 = grid_parameters$MTTfire_wood_years[,,mid_quant]**-1
var3 = grid_parameters$MTTharvest_wood_years[,,mid_quant]**-1
var4 = grid_parameters$MTTnatural_wood_years[,,mid_quant]**-1
var5 = grid_parameters$MTTfire_wood_years[,,mid_quant]**-1
var6 = grid_parameters$MTTharvest_wood_years[,,mid_quant]**-1
# Now make proportional
var1 = var1 / wood_turn ; var2 = var2 / wood_turn ; var3 = var3 / wood_turn
# Filter for the miombo AGB map locations
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var4[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var5[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var6[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Filter for the miombo AGB map locations
var1[which(landfilter == 1 & is.na(var1) == TRUE)] = 0
var2[which(landfilter == 1 & is.na(var2) == TRUE)] = 0
var3[which(landfilter == 1 & is.na(var3) == TRUE)] = 0
var4[which(landfilter == 1 & is.na(var4) == TRUE)] = 0
var5[which(landfilter == 1 & is.na(var5) == TRUE)] = 0
var6[which(landfilter == 1 & is.na(var6) == TRUE)] = 0
# Convert to raster
var1 = raster(vals = t((var1)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = raster(vals = t((var2)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = raster(vals = t((var3)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = raster(vals = t((var4)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = raster(vals = t((var5)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = raster(vals = t((var6)[,dim(area)[2]:1]), ext = extent(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# specify ranges
zrange1 = c(0,1)
zrange4 = c(0,1) * max(abs(range(c(values(var4),values(var5),values(var6)), na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_Wood_turnover_contribution_turnover.png",sep=""), height = 2500, width = 4900, res = 300)
par(mfrow=c(2,3), mar=c(0.5,0.4,3.0,7),omi=c(0.1,0.3,0.1,0.2))
# Partitioning of wood turnover, median estimate
plot(var1, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Natural MRT comp (0-1)",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var2, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Fire MRT comp (0-1)",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var3, zlim=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Biomass removal MRT comp (0-1)",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var4, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Natural turnover (",y^-1,")",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var5, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Fire turnover (",y^-1,")",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
plot(var6, zlim=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, box = FALSE, bty = "n",
     cex.axis = 2.5, legend.width = 2.2, axes = FALSE, axis.args=list(cex.axis=2.0,hadj=0.1),
     main = expression(paste("Biomass removal turnover (",y^-1,")",sep="")), col=colour_choices_loss)
plot(landmask, add=TRUE)
dev.off()

# Plot Foliage, fine root, wood, litter(foliar+fine root+wood?), soil mean residence times against main meteorology
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_MRT_meteorology_association.png",sep=""), height = 2200, width = 4500, res = 300)
par(mfrow=c(3,5), mar=c(4,2,1.4,1), omi = c(0.1,0.2,0.1,0.1))
# Temperature
plot(grid_parameters$MTT_foliar_years[,,mid_quant]~(grid_parameters$mean_temperature_C), main="Foliar MRT (yrs)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
#mtext(expression('C1'), side = 2, cex = 1.6, padj = -2.5, adj = 0.5)
plot(grid_parameters$MTT_root_years[,,mid_quant]~(grid_parameters$mean_temperature_C), main="Root MRT (yrs)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$MTT_wood_years[,,mid_quant]~(grid_parameters$mean_temperature_C), main="Wood MRT (yrs)", ylab="", xlab="Mean Temperature (C)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
plot(grid_parameters$MTT_DeadOrg_years[,,mid_quant]~(grid_parameters$mean_temperature_C), main="DeadOrg MRT (yrs)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$MTT_som_years[,,mid_quant]~(grid_parameters$mean_temperature_C), main="Soil MRT (yrs)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
# Precipitation
plot(grid_parameters$MTT_foliar_years[,,mid_quant]~(grid_parameters$mean_precipitation_kgm2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$MTT_root_years[,,mid_quant]~(grid_parameters$mean_precipitation_kgm2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$MTT_wood_years[,,mid_quant]~(grid_parameters$mean_precipitation_kgm2yr), main="", ylab="", xlab="Mean precipitation (mm/yr)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
plot(grid_parameters$MTT_DeadOrg_years[,,mid_quant]~(grid_parameters$mean_precipitation_kgm2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$MTT_som_years[,,mid_quant]~(grid_parameters$mean_precipitation_kgm2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
# Vapour pressure deficit
plot(grid_parameters$MTT_foliar_years[,,mid_quant]~(grid_parameters$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$MTT_root_years[,,mid_quant]~(grid_parameters$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$MTT_wood_years[,,mid_quant]~(grid_parameters$mean_vpd_Pa), main="", ylab="", xlab="Mean VPD (Pa)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
plot(grid_parameters$MTT_DeadOrg_years[,,mid_quant]~(grid_parameters$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$MTT_som_years[,,mid_quant]~(grid_parameters$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
dev.off()

# Plot Foliage, fine root, wood, litter(foliar+fine root+wood?), soil mean residence times against main disturbance
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_MRT_disturbance_association.png",sep=""), height = 2200, width = 4500, res = 300)
par(mfrow=c(3,5), mar=c(4,2,1.4,1), omi = c(0.1,0.2,0.1,0.1))
# Mean annual number of fires
plot(grid_parameters$MTT_foliar_years[,,mid_quant]~(FireFreq), main="Foliar MRT (yrs)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
#mtext(expression('C1'), side = 2, cex = 1.6, padj = -2.5, adj = 0.5)
plot(grid_parameters$MTT_root_years[,,mid_quant]~(FireFreq), main="Root MRT (yrs)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$MTT_wood_years[,,mid_quant]~(FireFreq), main="Wood MRT (yrs)", ylab="", xlab="No. annual fires", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
plot(grid_parameters$MTT_DeadOrg_years[,,mid_quant]~(FireFreq), main="DeadOrg MRT (yrs)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$MTT_som_years[,,mid_quant]~(FireFreq), main="Soil MRT (yrs)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
# Mean annual burned fraction
plot(grid_parameters$MTT_foliar_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$MTT_root_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$MTT_wood_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab="Annual burned fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
plot(grid_parameters$MTT_DeadOrg_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$MTT_som_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
# Mean annual forest harvest fraction
plot(grid_parameters$MTT_foliar_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$MTT_root_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$MTT_wood_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab="Annual harvested fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
plot(grid_parameters$MTT_DeadOrg_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$MTT_som_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
dev.off()

# Plot Foliage, fine root, wood NPP allocation fractions main meteorology
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NPP_meteorology_association.png",sep=""), height = 2200, width = 2800, res = 300)
par(mfrow=c(3,3), mar=c(4,2,1.4,1), omi = c(0.1,0.2,0.1,0.1))
# Temperature
plot(grid_parameters$NPP_foliar_fraction[,,mid_quant]~(grid_parameters$mean_temperature_C), main="Foliar NPP (0-1)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
#mtext(expression('C1'), side = 2, cex = 1.6, padj = -2.5, adj = 0.5)
plot(grid_parameters$NPP_root_fraction[,,mid_quant]~(grid_parameters$mean_temperature_C), main="Root NPP (0-1)", ylab="", xlab="Mean Temperature (C)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$NPP_wood_fraction[,,mid_quant]~(grid_parameters$mean_temperature_C), main="Wood NPP (0-1)", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
# Precipitation
plot(grid_parameters$NPP_foliar_fraction[,,mid_quant]~(grid_parameters$mean_precipitation_kgm2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$NPP_root_fraction[,,mid_quant]~(grid_parameters$mean_precipitation_kgm2yr), main="", ylab="", xlab="Mean precipitation (mm/yr)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$NPP_wood_fraction[,,mid_quant]~(grid_parameters$mean_precipitation_kgm2yr), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
# Vapour pressure deficit
plot(grid_parameters$NPP_foliar_fraction[,,mid_quant]~(grid_parameters$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$NPP_root_fraction[,,mid_quant]~(grid_parameters$mean_vpd_Pa), main="", ylab="", xlab="Mean VPD (Pa)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$NPP_wood_fraction[,,mid_quant]~(grid_parameters$mean_vpd_Pa), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
dev.off()

# Plot Foliage, fine root, wood NPP allocation fractions against main disturbance
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NPP_disturbance_association.png",sep=""), height = 2200, width = 2800, res = 300)
par(mfrow=c(3,3), mar=c(4,2,1.4,1), omi = c(0.1,0.2,0.1,0.1))
# Temperature
plot(grid_parameters$NPP_foliar_fraction[,,mid_quant]~(FireFreq), main="Foliar NPP (0-1)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
#mtext(expression('C1'), side = 2, cex = 1.6, padj = -2.5, adj = 0.5)
plot(grid_parameters$NPP_root_fraction[,,mid_quant]~(FireFreq), main="Root NPP (0-1)", ylab="", xlab="No. annual fires", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$NPP_wood_fraction[,,mid_quant]~(FireFreq), main="Wood NPP (0-1)", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
# Precipitation
plot(grid_parameters$NPP_foliar_fraction[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$NPP_root_fraction[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab="Annual burned fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$NPP_wood_fraction[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
# Vapour pressure deficit
plot(grid_parameters$NPP_foliar_fraction[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$NPP_root_fraction[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab="Annual harvested fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_parameters$NPP_wood_fraction[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
dev.off()

# Plot Foliage, fine root, wood, litter(foliar+fine root+wood?), soil mean residence times against main meteorology
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_MRT_meteorology_association_heatmap.png",sep=""), height = 2200, width = 4500, res = 300)
par(mfrow=c(3,5), mar=c(4,2,1.4,3.8), omi = c(0.1,0.2,0.1,0.1))
fudgeit.leg.lab=""
# Temperature
smoothScatter(grid_parameters$MTT_foliar_years[,,mid_quant]~(grid_parameters$mean_temperature_C), main="Foliar MRT (yrs)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
#mtext(expression('C1'), side = 2, cex = 1.6, padj = -2.5, adj = 0.5)
smoothScatter(grid_parameters$MTT_root_years[,,mid_quant]~(grid_parameters$mean_temperature_C), main="Root MRT (yrs)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
smoothScatter(grid_parameters$MTT_wood_years[,,mid_quant]~(grid_parameters$mean_temperature_C), main="Wood MRT (yrs)", ylab="", 
     xlab="Mean Temperature (C)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours,   
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
smoothScatter(grid_parameters$MTT_DeadOrg_years[,,mid_quant]~(grid_parameters$mean_temperature_C), main="DeadOrg MRT (yrs)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_parameters$MTT_som_years[,,mid_quant]~(grid_parameters$mean_temperature_C), main="Soil MRT (yrs)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
# Precipitation
fudgeit.leg.lab=""
smoothScatter(grid_parameters$MTT_foliar_years[,,mid_quant]~(grid_parameters$mean_precipitation_kgm2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
smoothScatter(grid_parameters$MTT_root_years[,,mid_quant]~(grid_parameters$mean_precipitation_kgm2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
smoothScatter(grid_parameters$MTT_wood_years[,,mid_quant]~(grid_parameters$mean_precipitation_kgm2yr), main="", ylab="", 
     xlab="Mean precipitation (mm/yr)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
smoothScatter(grid_parameters$MTT_DeadOrg_years[,,mid_quant]~(grid_parameters$mean_precipitation_kgm2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_parameters$MTT_som_years[,,mid_quant]~(grid_parameters$mean_precipitation_kgm2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
# Vapour pressure deficit
fudgeit.leg.lab=""
smoothScatter(grid_parameters$MTT_foliar_years[,,mid_quant]~(grid_parameters$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
smoothScatter(grid_parameters$MTT_root_years[,,mid_quant]~(grid_parameters$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
smoothScatter(grid_parameters$MTT_wood_years[,,mid_quant]~(grid_parameters$mean_vpd_Pa), main="", ylab="", xlab="Mean VPD (Pa)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
smoothScatter(grid_parameters$MTT_DeadOrg_years[,,mid_quant]~(grid_parameters$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_parameters$MTT_som_years[,,mid_quant]~(grid_parameters$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
dev.off()

# Plot Foliage, fine root, wood, litter(foliar+fine root+wood?), soil mean residence times against main disturbance
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_MRT_disturbance_association_heatmap.png",sep=""), height = 2200, width = 4500, res = 300)
par(mfrow=c(3,5), mar=c(4,2,1.4,3.8), omi = c(0.1,0.2,0.1,0.1))
# Mean annual number of fires
fudgeit.leg.lab=""
smoothScatter(grid_parameters$MTT_foliar_years[,,mid_quant]~(FireFreq), main="Foliar MRT (yrs)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
#mtext(expression('C1'), side = 2, cex = 1.6, padj = -2.5, adj = 0.5)
smoothScatter(grid_parameters$MTT_root_years[,,mid_quant]~(FireFreq), main="Root MRT (yrs)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
smoothScatter(grid_parameters$MTT_wood_years[,,mid_quant]~(FireFreq), main="Wood MRT (yrs)", ylab="", xlab="No. annual fires", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
smoothScatter(grid_parameters$MTT_DeadOrg_years[,,mid_quant]~(FireFreq), main="DeadOrg MRT (yrs)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_parameters$MTT_som_years[,,mid_quant]~(FireFreq), main="Soil MRT (yrs)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
# Mean annual burned fraction
fudgeit.leg.lab=""
smoothScatter(grid_parameters$MTT_foliar_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
smoothScatter(grid_parameters$MTT_root_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
smoothScatter(grid_parameters$MTT_wood_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab="Annual burned fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
smoothScatter(grid_parameters$MTT_DeadOrg_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_parameters$MTT_som_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
# Mean annual forest harvest fraction
fudgeit.leg.lab=""
smoothScatter(grid_parameters$MTT_foliar_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
smoothScatter(grid_parameters$MTT_root_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
smoothScatter(grid_parameters$MTT_wood_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab="Annual harvested fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
smoothScatter(grid_parameters$MTT_DeadOrg_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_parameters$MTT_som_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
dev.off()

# Plot Foliage, fine root, wood NPP allocation fractions main meteorology
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NPP_meteorology_association_heatmap.png",sep=""), height = 2200, width = 2800, res = 300)
par(mfrow=c(3,3), mar=c(4,2,1.4,3.8), omi = c(0.1,0.2,0.1,0.1))
fudgeit.leg.lab=""
# Temperature
smoothScatter(grid_parameters$NPP_foliar_fraction[,,mid_quant]~(grid_parameters$mean_temperature_C), main="Foliar NPP (0-1)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_parameters$mean_temperature_C, na.rm=TRUE)))
smoothScatter(grid_parameters$NPP_root_fraction[,,mid_quant]~(grid_parameters$mean_temperature_C), main="Root NPP (0-1)", 
     ylab="", xlab="Mean Temperature (C)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_parameters$mean_temperature_C, na.rm=TRUE)))
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_parameters$NPP_wood_fraction[,,mid_quant]~(grid_parameters$mean_temperature_C), main="Wood NPP (0-1)", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_parameters$mean_temperature_C, na.rm=TRUE)))
# Precipitation
fudgeit.leg.lab=""
smoothScatter(grid_parameters$NPP_foliar_fraction[,,mid_quant]~(grid_parameters$mean_precipitation_kgm2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_parameters$mean_precipitation_kgm2yr, na.rm=TRUE)))
smoothScatter(grid_parameters$NPP_root_fraction[,,mid_quant]~(grid_parameters$mean_precipitation_kgm2yr), main="", 
     ylab="", xlab="Mean precipitation (mm/yr)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_parameters$mean_precipitation_kgm2yr, na.rm=TRUE)))
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_parameters$NPP_wood_fraction[,,mid_quant]~(grid_parameters$mean_precipitation_kgm2yr), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_parameters$mean_precipitation_kgm2yr, na.rm=TRUE)))
# Vapour pressure deficit
fudgeit.leg.lab=""
smoothScatter(grid_parameters$NPP_foliar_fraction[,,mid_quant]~(grid_parameters$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_parameters$mean_vpd_Pa, na.rm=TRUE)))
smoothScatter(grid_parameters$NPP_root_fraction[,,mid_quant]~(grid_parameters$mean_vpd_Pa), main="", ylab="", xlab="Mean VPD (Pa)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_parameters$mean_vpd_Pa, na.rm=TRUE)))
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_parameters$NPP_wood_fraction[,,mid_quant]~(grid_parameters$mean_vpd_Pa), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_parameters$mean_vpd_Pa, na.rm=TRUE)))
dev.off()

# Plot Foliage, fine root, wood NPP allocation fractions against main disturbance
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NPP_disturbance_association_heatmap.png",sep=""), height = 2200, width = 2800, res = 300)
par(mfrow=c(3,3), mar=c(4,2,1.4,3.8), omi = c(0.1,0.2,0.1,0.1))
fudgeit.leg.lab=""
# Temperature
smoothScatter(grid_parameters$NPP_foliar_fraction[,,mid_quant]~(FireFreq), main="Foliar NPP (0-1)", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
smoothScatter(grid_parameters$NPP_root_fraction[,,mid_quant]~(FireFreq), main="Root NPP (0-1)", ylab="", xlab="No. annual fires", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_parameters$NPP_wood_fraction[,,mid_quant]~(FireFreq), main="Wood NPP (0-1)", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
# Precipitation
fudgeit.leg.lab=""
smoothScatter(grid_parameters$NPP_foliar_fraction[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
smoothScatter(grid_parameters$NPP_root_fraction[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab="Annual burned fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_parameters$NPP_wood_fraction[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
# Vapour pressure deficit
fudgeit.leg.lab=""
smoothScatter(grid_parameters$NPP_foliar_fraction[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
smoothScatter(grid_parameters$NPP_root_fraction[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab="Annual harvested fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_parameters$NPP_wood_fraction[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
dev.off()
