
###
## Process CARDAMOM into Global Carbon Project C-budget domains.
## This file will create an ascii file containing the C-budget terms
## for the globe, tropics, north (>N30 degree) and south (>S30 degree).
## Terms are net biome productivity, including all simulated terms.
## Units PgC/yr, annual estimates. 
## The median and propogated uncertainty information are provided.
## Created 29/07/2022 by T. Luke Smallman
###

###
## Job specific information

# set working directory
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")

# set input and output directories
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca_gpp_fire/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca_CsomPriorNCSCD/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca_nbe_CsomPriorNCSDC3m/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Trendyv9_historical/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/reccap2_permafrost_1deg_C7_isimip3a_agb_lca_gpp_fire_nbe/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Mexico_1deg_C7_agb_lca_gpp_fire_nbe/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/Miombo_0.5deg_allWood/"
input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/global_1deg_C7_GCP_LCA_AGB/"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC_CDEA_ACM2_BUCKET_MHMCMC/global_1deg_C7_GCP_LCA_AGB_GPP_FIRE/"

# Specify any extra information for the filename
output_prefix = "CARDAMOM_S3_" # follow with "_"
#output_prefix = "CARDAMOM_S2_" # follow with "_"
output_suffix = "" # begin with "_"

###
## Load libraries, functions and data needed

# load needed libraries
library(ncdf4)
library(raster)
library(compiler)

# Load needed functions 
source("~/WORK/GREENHOUSE/models/CARDAMOM/R_functions/load_all_cardamom_functions.r")

# load the CARDAMOM files
load(paste(input_dir,"/infofile.RData",sep=""))
load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))

# Create output directory to aid storage management
out_dir = paste(PROJECT$results_processedpath,"/trendy_output",sep="")
if (dir.exists(out_dir) == FALSE) {
    dir.create(out_dir)
}

###
## Begin creating information for processing and subsequent saving to files

# Time information
years = c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year))
nos_years = length(years)
steps_per_year = dim(grid_output$lai_m2m2)[3] / nos_years

# create lat / long axes, assumes regular WGS-84 grid
output = determine_lat_long_needed(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution,PROJECT$grid_type,PROJECT$waterpixels)
# NOTE: rev due to CARDAMOM grid being inverse of what comes out of raster function. Should consider changing this at some point.
longitude = output$obs_long_grid[,1] ; latitude = rev(output$obs_lat_grid[1,]) 
# Tidy up
rm(output) ; gc(reset=TRUE,verbose=FALSE)

# Extract the available quantiles 
quantiles_wanted = grid_output$num_quantiles
nos_quantiles = length(quantiles_wanted)
# Check that the quantiles we want to use are available
if (length(which(quantiles_wanted == 0.05)) == 1) {
    low_quant = which(quantiles_wanted == 0.05)
} else {
    stop("Desired low quantile cannot be found")
}
if (length(which(quantiles_wanted == 0.5)) == 1) {
    mid_quant = which(quantiles_wanted == 0.5)
} else {
    stop("Median quantile cannot be found")
}
if (length(which(quantiles_wanted == 0.95)) == 1) {
    high_quant = which(quantiles_wanted == 0.95)
} else {
    stop("Desired high quantile cannot be found")
}

###
## Loop through grid to aggregate the net biome productivity (NBP = GPP - Ra - Rh - Fire - Harvest)

# Median estimate
NBP_global_PgCyr = rep(0, nos_years)
NBP_north_PgCyr = rep(0, nos_years)
NBP_tropics_PgCyr = rep(0, nos_years)
NBP_south_PgCyr = rep(0, nos_years)
# Lower CI
NBP_global_PgCyr_lowCI = rep(0, nos_years)
NBP_north_PgCyr_lowCI = rep(0, nos_years)
NBP_tropics_PgCyr_lowCI = rep(0, nos_years)
NBP_south_PgCyr_lowCI = rep(0, nos_years)
# Upper CI
NBP_global_PgCyr_highCI = rep(0, nos_years)
NBP_north_PgCyr_highCI = rep(0, nos_years)
NBP_tropics_PgCyr_highCI = rep(0, nos_years)
NBP_south_PgCyr_highCI = rep(0, nos_years)

# Fill the output arrays
for (n in seq(1, length(PROJECT$sites))) {

     # Ensure the site has been processed
     if (is.na(grid_output$i_location[n]) == FALSE) {

         # Extract grid position
         i = grid_output$i_location[n]
         j = grid_output$j_location[n]

         # Read in site specific drivers
         drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep=""))
         
         # Determine nos days per time step         
         deltat = drivers$met[,1] # run_day
         deltat[2:length(deltat)] = diff(deltat)

         # Accumulate global
         NBP_global_PgCyr = NBP_global_PgCyr + (rollapply(grid_output$nbp_gCm2day[n,mid_quant,]*deltat, width = steps_per_year, by = steps_per_year, FUN = sum, na.rm=TRUE)*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NBP_global_PgCyr_lowCI = NBP_global_PgCyr_lowCI + (rollapply(grid_output$nbp_gCm2day[n,low_quant,]*deltat, width = steps_per_year, by = steps_per_year, FUN = sum, na.rm=TRUE)*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         NBP_global_PgCyr_highCI = NBP_global_PgCyr_highCI + (rollapply(grid_output$nbp_gCm2day[n,high_quant,]*deltat, width = steps_per_year, by = steps_per_year, FUN = sum, na.rm=TRUE)*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
         
         # Where appropriate accumulate tropical
         if (grid_output$lat[i,j] <= 30 & grid_output$lat[i,j] >= -30) {
             NBP_tropics_PgCyr = NBP_tropics_PgCyr + (rollapply(grid_output$nbp_gCm2day[n,mid_quant,]*deltat, width = steps_per_year, by = steps_per_year, FUN = sum, na.rm=TRUE)*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_tropics_PgCyr_lowCI = NBP_tropics_PgCyr_lowCI + (rollapply(grid_output$nbp_gCm2day[n,low_quant,]*deltat, width = steps_per_year, by = steps_per_year, FUN = sum, na.rm=TRUE)*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_tropics_PgCyr_highCI = NBP_tropics_PgCyr_highCI + (rollapply(grid_output$nbp_gCm2day[n,high_quant,]*deltat, width = steps_per_year, by = steps_per_year, FUN = sum, na.rm=TRUE)*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)         
         }

         # Where appropriate accumulate north
         if (grid_output$lat[i,j] > 30) {
             NBP_north_PgCyr = NBP_north_PgCyr + (rollapply(grid_output$nbp_gCm2day[n,mid_quant,]*deltat, width = steps_per_year, by = steps_per_year, FUN = sum, na.rm=TRUE)*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_north_PgCyr_lowCI = NBP_north_PgCyr_lowCI + (rollapply(grid_output$nbp_gCm2day[n,low_quant,]*deltat, width = steps_per_year, by = steps_per_year, FUN = sum, na.rm=TRUE)*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_north_PgCyr_highCI = NBP_north_PgCyr_highCI + (rollapply(grid_output$nbp_gCm2day[n,high_quant,]*deltat, width = steps_per_year, by = steps_per_year, FUN = sum, na.rm=TRUE)*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)         
         }

         # Where appropriate accumulate south
         if (grid_output$lat[i,j] < -30) {
             NBP_south_PgCyr = NBP_south_PgCyr + (rollapply(grid_output$nbp_gCm2day[n,mid_quant,]*deltat, width = steps_per_year, by = steps_per_year, FUN = sum, na.rm=TRUE)*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_south_PgCyr_lowCI = NBP_south_PgCyr_lowCI + (rollapply(grid_output$nbp_gCm2day[n,low_quant,]*deltat, width = steps_per_year, by = steps_per_year, FUN = sum, na.rm=TRUE)*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)
             NBP_south_PgCyr_highCI = NBP_south_PgCyr_highCI + (rollapply(grid_output$nbp_gCm2day[n,high_quant,]*deltat, width = steps_per_year, by = steps_per_year, FUN = sum, na.rm=TRUE)*grid_output$area_m2[i,j]*grid_output$land_fraction[i,j]*1e-15)         
         }

     } # Does the file exist / has it been processed

} # site loop

###
## Combined together the data into an output variable

output = data.frame(Year = years, Global = NBP_global_PgCyr, North = NBP_north_PgCyr, Tropics = NBP_tropics_PgCyr, South = NBP_south_PgCyr,
                                  Global_5pc = NBP_global_PgCyr_lowCI, North_5pc = NBP_north_PgCyr_lowCI, Tropics_5pc = NBP_tropics_PgCyr_lowCI, South_5pc = NBP_south_PgCyr_lowCI,
                                  Global_95pc = NBP_global_PgCyr_highCI, North_95pc = NBP_north_PgCyr_highCI, Tropics_95pc = NBP_tropics_PgCyr_highCI, South_95pc = NBP_south_PgCyr_highCI)
###
## Begin writing out the file

write.table(output,file = paste(out_dir,"/",output_prefix,"NBP",output_suffix,".txt",sep=""), sep=" ", row.names = FALSE)
