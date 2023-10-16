
###
## Process CARDAMOM-DALEC output files into NetCDF files
## consistent with the TRENDYv11 / GCP model intercomparison structure.
## NOTE: unlike the sibling script this script does not put each variable into a single file
## and is thus not consistent with the latest guidance. But adapted to be more useful
## at convening CARDAMOM outputs
###

##TODO:
# COnvert back to standard names and units for CARDAMOM

###
## Job specific information

# set working directory
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")

# set input and output directories
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/reccap2_permafrost_1deg_dalec2_isimip3a_agb_lca_nbe_gpp_CsomPriorNCSDC3m"
input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/Miombo_0.5deg_allWood"
#input_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_1deg_dalec4_GCP_LCA_AGB_GPP"

# Specify any extra information for the filename
output_prefix = "" # follow with "_"
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
nos_years = length(c(as.numeric(PROJECT$start_year):as.numeric(PROJECT$end_year)))
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

###
## Begin defining variables
###

###
## Variables we assume always exist

## At model time step variables
# DRIVERS
AIRT_MIN = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
AIRT_MAX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
SWRAD = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
CO2 = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
DOY = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
PRECIP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
FLOSS_FRAC = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days))) 
BURNT_FRAC = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
WINDSPD = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
VPD = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
# STATE OBSERVATIONS (gC/m2, except LAI = m2/m2)
LAI_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
LAI_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
WOOD_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
WOOD_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
SOIL_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
SOIL_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
# FLUX OBSERVATIONS (gC/m2/day)
GPP_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
GPP_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
NEE_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
NEE_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
NBE_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
NBE_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
FIRE_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))
FIRE_UNC_OBS = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,length(PROJECT$model$timestep_days)))

###
## Variables based on their presence

# C STATES
if (exists(x = "lai_m2m2", where = grid_output)) {LAI = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "Ctotal_gCm2", where = grid_output)) {TOT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "labile_gCm2", where = grid_output)) {LAB = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "foliage_gCm2", where = grid_output)) {FOL = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "roots_gCm2", where = grid_output)) {ROOT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "wood_gCm2", where = grid_output)) {WOOD = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "litter_gCm2", where = grid_output)) {LIT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "som_gCm2", where = grid_output)) {SOIL = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "woodlitter_gCm2", where = grid_output)) {WLIT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "dom_gCm2", where = grid_output)) {DOM = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "biomass_gCm2", where = grid_output)) {BIO = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# C STATE CHANGE ESTIMATES
if (exists(x = "dCbiomass_gCm2", where = grid_output)) {dBIO = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "dCdom_gCm2", where = grid_output)) {dDOM = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# C STATES - mean
if (exists(x = "mean_lai_m2m2", where = grid_output)) {MLAI = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_Ctotal_gCm2", where = grid_output)) {MTOT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_labile_gCm2", where = grid_output)) {MLAB = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_foliage_gCm2", where = grid_output)) {MFOL = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_roots_gCm2", where = grid_output)) {MROOT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_wood_gCm2", where = grid_output)) {MWOOD = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "litter_gCm2", where = grid_output)) {MLIT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "som_gCm2", where = grid_output)) {MSOIL = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_woodlitter_gCm2", where = grid_output)) {MWLIT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_dom_gCm2", where = grid_output)) {MDOM = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_biomass_gCm2", where = grid_output)) {MBIO = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
# C FLUXES
if (exists(x = "gpp_gCm2day", where = grid_output)) {GPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "rauto_gCm2day", where = grid_output)) {RAU = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "rhet_gCm2day", where = grid_output)) {RHE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "npp_gCm2day", where = grid_output)) {NPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "fire_gCm2day", where = grid_output)) {FIR = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "harvest_gCm2day", where = grid_output)) {HARV = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "reco_gCm2day", where = grid_output)) {RECO = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "nee_gCm2day", where = grid_output)) {NEE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "nbe_gCm2day", where = grid_output)) {NBE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "nbp_gCm2day", where = grid_output)) {NBP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# C FLUXES - ANNUAL
if (exists(x = "mean_annual_gpp_gCm2day", where = grid_output)) {AGPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_rauto_gCm2day", where = grid_output)) {ARAU = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_rhet_gCm2day", where = grid_output)) {ARHE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_npp_gCm2day", where = grid_output)) {ANPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_fire_gCm2day", where = grid_output)) {AFIR = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_harvest_gCm2day", where = grid_output)) {AHARV = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_reco_gCm2day", where = grid_output)) {ARECO = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_nee_gCm2day", where = grid_output)) {ANEE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_nbe_gCm2day", where = grid_output)) {ANBE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
if (exists(x = "mean_annual_nbp_gCm2day", where = grid_output)) {ANBP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,nos_years))}
# C FLUXES - means
if (exists(x = "mean_gpp_gCm2day", where = grid_output)) {MGPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_rauto_gCm2day", where = grid_output)) {MRAU = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_rhet_gCm2day", where = grid_output)) {MRHE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_npp_gCm2day", where = grid_output)) {MNPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_fire_gCm2day", where = grid_output)) {MFIR = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_harvest_gCm2day", where = grid_output)) {MHARV = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_reco_gCm2day", where = grid_output)) {MRECO = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_nee_gCm2day", where = grid_output)) {MNEE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_nbe_gCm2day", where = grid_output)) {MNBE = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
if (exists(x = "mean_nbp_gCm2day", where = grid_output)) {MNBP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles))}
# H2O FLUXES
if (exists(x = "ET_kgH2Om2day", where = grid_output)) {ET = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "Etrans_kgH2Om2day", where = grid_output)) {Etrans = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "Esoil_kgH2Om2day", where = grid_output)) {Esoil = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "Ewetcanopy_kgH2Om2day", where = grid_output)) {Ewetcanopy = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "total_drainage_kgH2Om2day", where = grid_output)) {total_drainage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "runoff_kgH2Om2day", where = grid_output)) {runoff = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "underflow_kgH2Om2day", where = grid_output)) {underflow = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# H2O STATES
if (exists(x = "SurfWater_kgH2Om2", where = grid_output)) {SurfWater = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "snow_kgH2Om2", where = grid_output)) {SNOW = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# BIOPHYSICAL
if (exists(x = "CiCa", where = grid_output)) {CiCa = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "wue_eco_gCkgH2O", where = grid_output)) {wue_eco = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "wue_plant_gCkgH2O", where = grid_output)) {wue_plant = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "wSWP_MPa", where = grid_output)) {wSWP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "APAR_MJm2day", where = grid_output)) {APAR = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "gs_demand_supply_ratio", where = grid_output)) {gs_demand_supply_ratio = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "gs_mmolH2Om2day", where = grid_output)) {gs = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "gb_mmolH2Om2day", where = grid_output)) {gb = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Direct allocation of NPP (labile, foliar, fine root, wood, gC/m2/day)
if (exists(x = "alloc_labile_gCm2day", where = grid_output)) {NPP_labile_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "alloc_foliage_gCm2day", where = grid_output)) {NPP_foliage_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "alloc_wood_gCm2day", where = grid_output)) {NPP_wood_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "alloc_roots_gCm2day", where = grid_output)) {NPP_root_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Combined (direct / indirect) allocation of C to foliage (gC/m2/day)
if (exists(x = "combined_alloc_foliage_gCm2day", where = grid_output)) {NPP_combinedfoliage_FLX = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Allocation fraction of expressed NPP (foliar, fine root, wood)
if (exists(x = "NPP_foliage_fraction", where = grid_output)) {fNPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "NPP_wood_fraction", where = grid_output)) {wNPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "NPP_roots_fraction", where = grid_output)) {rNPP = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# MRT (labile, foliar, wood, fine root, litter, soil, biomass, dead organic matter; years)
if (exists(x = "MTT_annual_labile_years", where = grid_output)) {labMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "MTT_annual_foliage_years", where = grid_output)) {folMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "MTT_annual_wood_years", where = grid_output)) {wooMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "MTT_annual_roots_years", where = grid_output)) {rooMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "MTT_annual_litter_years", where = grid_output)) {litMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "MTT_annual_woodlitter_years", where = grid_output)) {wlitMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "MTT_annual_som_years", where = grid_output)) {somMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "MTT_annual_biomass_years", where = grid_output)) {bioMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "MTT_annual_dom_years", where = grid_output)) {domMRT = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Total outflux from these pools (labile, foliar, wood, fine root, litter, soil, biomass, dead organic matter; gC/m2/day)
if (exists(x = "outflux_labile_gCm2day", where = grid_output)) {outflux_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "outflux_foliage_gCm2day", where = grid_output)) {outflux_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "outflux_wood_gCm2day", where = grid_output)) {outflux_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "outflux_roots_gCm2day", where = grid_output)) {outflux_root = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "outflux_litter_gCm2day", where = grid_output)) {outflux_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "outflux_woodlitter_gCm2day", where = grid_output)) {outflux_woodlitter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "outflux_som_gCm2day", where = grid_output)) {outflux_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "outflux_biomass_gCm2day", where = grid_output)) {outflux_bio = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "outflux_dom_gCm2day", where = grid_output)) {outflux_dom = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Fire mortality outflux (labile, foliar, wood, fine root, litter, biomass, dead organic matter; gC/m2/day)
if (exists(x = "FIRElitter_labile_gCm2day", where = grid_output)) {FIRElitter_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIRElitter_foliage_gCm2day", where = grid_output)) {FIRElitter_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIRElitter_wood_gCm2day", where = grid_output)) {FIRElitter_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIRElitter_roots_gCm2day", where = grid_output)) {FIRElitter_root = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIRElitter_litter_gCm2day", where = grid_output)) {FIRElitter_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIRElitter_woodlitter_gCm2day", where = grid_output)) {FIRElitter_woodlitter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIRElitter_biomass_gCm2day", where = grid_output)) {FIRElitter_bio = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Fire combustion outflux (labile, foliar, wood, fine root, litter, soil, biomass, dead organic matter; gC/m2/day)
if (exists(x = "FIREemiss_labile_gCm2day", where = grid_output)) {FIREemiss_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_foliage_gCm2day", where = grid_output)) {FIREemiss_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_wood_gCm2day", where = grid_output)) {FIREemiss_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_roots_gCm2day", where = grid_output)) {FIREemiss_root = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_litter_gCm2day", where = grid_output)) {FIREemiss_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_woodlitter_gCm2day", where = grid_output)) {FIREemiss_woodlitter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_som_gCm2day", where = grid_output)) {FIREemiss_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_biomass_gCm2day", where = grid_output)) {FIREemiss_bio = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "FIREemiss_dom_gCm2day", where = grid_output)) {FIREemiss_dom = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Harvest mortality outflux (labile, foliar, wood, fine root, litter, soil, biomass, dead organic matter; gC/m2/day)
if (exists(x = "HARVESTlitter_labile_gCm2day", where = grid_output)) {HARVESTlitter_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTlitter_foliage_gCm2day", where = grid_output)) {HARVESTlitter_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTlitter_wood_gCm2day", where = grid_output)) {HARVESTlitter_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTlitter_roots_gCm2day", where = grid_output)) {HARVESTlitter_root = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTlitter_biomass_gCm2day", where = grid_output)) {HARVESTlitter_bio = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
# Harvest extracted outflux (labile, foliar, wood, fine root, litter, soil, biomass, dead organic matter; gC/m2/day)
if (exists(x = "HARVESTextracted_labile_gCm2day", where = grid_output)) {HARVESTextracted_labile = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTextracted_foliage_gCm2day", where = grid_output)) {HARVESTextracted_foliage = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTextracted_wood_gCm2day", where = grid_output)) {HARVESTextracted_wood = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTextracted_roots_gCm2day", where = grid_output)) {HARVESTextracted_root = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTextracted_litter_gCm2day", where = grid_output)) {HARVESTextracted_litter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTextracted_woodlitter_gCm2day", where = grid_output)) {HARVESTextracted_woodlitter = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTextracted_som_gCm2day", where = grid_output)) {HARVESTextracted_som = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTextracted_biomass_gCm2day", where = grid_output)) {HARVESTextracted_bio = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}
if (exists(x = "HARVESTextracted_dom_gCm2day", where = grid_output)) {HARVESTextracted_dom = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim,nos_quantiles,length(PROJECT$model$timestep_days)))}

# Fill the output arrays
for (n in seq(1, length(PROJECT$sites))) {

     # Ensure the site has been processed
     if (is.na(grid_output$i_location[n]) == FALSE) {

         # Extract grid position
         i = grid_output$i_location[n]
         j = grid_output$j_location[n]

         # Read in site specific drivers
         drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep=""))

         # DRIVERS
         AIRT_MIN[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,2]+273.15     # mint C -> K
         AIRT_MAX[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,3]+273.15     # maxt C -> K
         SWRAD[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,4]*1e6*(1/86400) # SWRAD MJ/m2/day -> W/m2
         CO2[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,5]                 # CO2 ppm
         DOY[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,6]                 # Julian day of year
         PRECIP[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,7]              # Precipitation kgH2O/m2/s
         FLOSS_FRAC[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,8]          # Forest loss fraction
         BURNT_FRAC[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,9]          # Burned fraction
         WINDSPD[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,15]            # Wind speed m/s
         VPD[grid_output$i_location[n],grid_output$j_location[n],] = drivers$met[,16]                # Vapour pressure deficit Pa

         # STATE OBSERVATIONS
         LAI_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,3]     # LAI m2/m2
         LAI_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,4] # LAI UNC m2/m2
         # ...first assign priors to the 1st time step to simplify our storage
         WOOD_OBS[grid_output$i_location[n],grid_output$j_location[n],1] = drivers$parpriors[21] # Wood gC/m2
         WOOD_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],1] = drivers$parpriorunc[21] # Wood UNC gC/m2
         SOIL_OBS[grid_output$i_location[n],grid_output$j_location[n],1] = drivers$parpriors[23] # SOM gC/m2
         SOIL_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],1] = drivers$parpriorunc[23] # SOM UNC gC/m2
         # ...second assign time series information if any exists
         WOOD_OBS[grid_output$i_location[n],grid_output$j_location[n],2:length(PROJECT$model$timestep_days)] = drivers$obs[2:length(PROJECT$model$timestep_days),13]
         WOOD_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],2:length(PROJECT$model$timestep_days)] = drivers$obs[2:length(PROJECT$model$timestep_days),14]
         SOIL_OBS[grid_output$i_location[n],grid_output$j_location[n],2:length(PROJECT$model$timestep_days)] = drivers$obs[2:length(PROJECT$model$timestep_days),19]
         SOIL_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],2:length(PROJECT$model$timestep_days)] = drivers$obs[2:length(PROJECT$model$timestep_days),20]
         # FLUX OBSERVATIONS (gC/m2/day)
         GPP_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,1]
         GPP_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,2]
         NEE_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,5]
         NEE_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,6]
         NBE_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,35]
         NBE_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,36]
         FIRE_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,7]
         FIRE_UNC_OBS[grid_output$i_location[n],grid_output$j_location[n],] = drivers$obs[,8]

         ## At model time step
         # STATES (NOTE: unit conversions gC/m2 -> kgC/m2, except LAI)
         if (exists("LAI")) {LAI[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$lai_m2m2[n,,]}
         if (exists("TOT")) {TOT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$Ctotal_gCm2[n,,]*1e-3}
         if (exists("LAB")) {LAB[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$labile_gCm2[n,,]*1e-3}
         if (exists("FOL")) {FOL[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$foliage_gCm2[n,,]*1e-3}
         if (exists("ROOT")) {ROOT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$roots_gCm2[n,,]*1e-3}
         if (exists("WOOD")) {WOOD[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$wood_gCm2[n,,]*1e-3}
         if (exists("LIT")) {LIT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$litter_gCm2[n,,]*1e-3}
         if (exists("SOIL")) {SOIL[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$som_gCm2[n,,]*1e-3}
         if (exists("WLIT")) {WLIT[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$woodlitter_gCm2[n,,]*1e-3}
         if (exists("DOM")) {DOM[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$dom_gCm2[n,,]*1e-3}
         if (exists("BIO")) {BIO[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$biomass_gCm2[n,,]*1e-3}
         # Change in stocks
         if (exists("dBIO")) {dBIO[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$dCbiomass_gCm2[n,,]*1e-3}
         if (exists("dDOM")) {dDOM[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$dCdom_gCm2[n,,]*1e-3}
         # STATES - means (NOTE: unit conversions gC/m2 -> kgC/m2, except LAI)
         if (exists("MLAI")) {MLAI[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_lai_m2m2[grid_output$i_location[n],grid_output$j_location[n],]}
         if (exists("MTOT")) {MTOT[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_Ctotal_gCm2[grid_output$i_location[n],grid_output$j_location[n],]*1e-3}
         if (exists("MLAB")) {MLAB[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_labile_gCm2[grid_output$i_location[n],grid_output$j_location[n],]*1e-3}
         if (exists("MFOL")) {MFOL[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_foliage_gCm2[grid_output$i_location[n],grid_output$j_location[n],]*1e-3}
         if (exists("MROOT")) {MROOT[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_roots_gCm2[grid_output$i_location[n],grid_output$j_location[n],]*1e-3}
         if (exists("MWOOD")) {MWOOD[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_wood_gCm2[grid_output$i_location[n],grid_output$j_location[n],]*1e-3}
         if (exists("MLIT")) {MLIT[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_litter_gCm2[grid_output$i_location[n],grid_output$j_location[n],]*1e-3}
         if (exists("MSOIL")) {MSOIL[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_som_gCm2[grid_output$i_location[n],grid_output$j_location[n],]*1e-3}
         if (exists("MWLIT")) {MWLIT[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_woodlitter_gCm2[grid_output$i_location[n],grid_output$j_location[n],]*1e-3}
         if (exists("MDOM")) {MDOM[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_dom_gCm2[grid_output$i_location[n],grid_output$j_location[n],]*1e-3}
         if (exists("MBIO")) {MBIO[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_biomass_gCm2[grid_output$i_location[n],grid_output$j_location[n],]*1e-3}
         # FLUXES (NOTE; unit conversion gC/m2/day -> kgC/m2/s)
         if (exists("GPP")) {GPP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$gpp_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("RAU")) {RAU[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$rauto_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("RHE")) {RHE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$rhet_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NPP")) {NPP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$npp_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("FIR")) {FIR[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$fire_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("HARV")) {HARV[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$harvest_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("RECO")) {RECO[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$reco_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NEE")) {NEE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$nee_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NBE")) {NBE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$nbe_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NBP")) {NBP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$nbp_gCm2day[n,,]* 1e-3 * (1/86400)}
         # FLUXES - ANNUAL (NOTE; unit conversion gC/m2/day -> kgC/m2/s)
         if (exists("AGPP")) {AGPP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_gpp_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("ARAU")) {ARAU[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_rauto_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("ARHE")) {ARHE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_rhet_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("ANPP")) {ANPP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_npp_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("AFIR")) {AFIR[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_fire_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("AHARV")) {AHARV[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_harvest_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("ARECO")) {ARECO[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_reco_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("ANEE")) {ANEE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_nee_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("ANBE")) {ANBE[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_nbe_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("ANBP")) {ANBP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$mean_annual_nbp_gCm2day[n,,]* 1e-3 * (1/86400)}
         # C FLUXES - means
         if (exists("MGPP")) {MGPP[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_gpp_gCm2day[grid_output$i_location[n],grid_output$j_location[n],]* 1e-3 * (1/86400)}
         if (exists("MRAU")) {MRAU[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_rauto_gCm2day[grid_output$i_location[n],grid_output$j_location[n],]* 1e-3 * (1/86400)}
         if (exists("MRHE")) {MRHE[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_rhet_gCm2day[grid_output$i_location[n],grid_output$j_location[n],]* 1e-3 * (1/86400)}
         if (exists("MNPP")) {MNPP[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_npp_gCm2day[grid_output$i_location[n],grid_output$j_location[n],]* 1e-3 * (1/86400)}
         if (exists("MFIR")) {MFIR[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_fire_gCm2day[grid_output$i_location[n],grid_output$j_location[n],]* 1e-3 * (1/86400)}
         if (exists("MHARV")) {MHARV[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_harvest_gCm2day[grid_output$i_location[n],grid_output$j_location[n],]* 1e-3 * (1/86400)}
         if (exists("MRECO")) {MRECO[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_reco_gCm2day[grid_output$i_location[n],grid_output$j_location[n],]* 1e-3 * (1/86400)}
         if (exists("MNEE")) {MNEE[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_nee_gCm2day[grid_output$i_location[n],grid_output$j_location[n],]* 1e-3 * (1/86400)}
         if (exists("MNBE")) {MNBE[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_nbe_gCm2day[grid_output$i_location[n],grid_output$j_location[n],]* 1e-3 * (1/86400)}
         if (exists("MNBP")) {MNBP[grid_output$i_location[n],grid_output$j_location[n],] = grid_output$mean_nbp_gCm2day[grid_output$i_location[n],grid_output$j_location[n],]* 1e-3 * (1/86400)}
         # H2O FLUXES
         if (exists("ET")) {ET[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$ET_kgH2Om2day[n,,] * (1/86400)}
         if (exists("Etrans")) {Etrans[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$Etrans_kgH2Om2day[n,,] * (1/86400)}
         if (exists("Esoil")) {Esoil[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$Esoil_kgH2Om2day[n,,] * (1/86400)}
         if (exists("Ewetcanopy")) {Ewetcanopy[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$Ewetcanopy_kgH2Om2day[n,,] * (1/86400)}
         if (exists("runoff")) {runoff[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$runoff_kgH2Om2day[n,,] * (1/86400)}
         if (exists("underflow")) {underflow[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$underflow_kgH2Om2day[n,,] * (1/86400)}
         if (exists("total_drainage")) {total_drainage[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$total_drainage_kgH2Om2day[n,,] * (1/86400)}
         # H2O STATES
         if (exists("SurfWater")) {SurfWater[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$SurfWater_kgH2Om2[n,,]}
         if (exists("SNOW")) {SNOW[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$snow_kgH2Om2[n,,]}
         # BIOPHYSICAL
         if (exists("CiCa")) {CiCa[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$CiCa[n,,]}
         if (exists("wue_eco")) {wue_eco[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$wue_eco_gCkgH2O[n,,]}
         if (exists("wue_plant")) {wue_plant[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$wue_plant_gCkgH2O[n,,]}
         if (exists("wSWP")) {wSWP[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$wSWP_MPa[n,,]}
         if (exists("APAR")) {APAR[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$APAR_MJm2day[n,,]}
         if (exists("gs_demand_supply_ratio")) {gs_demand_supply_ratio[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$gs_demand_supply_ratio[n,,]}
         if (exists("gs")) {gs[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$gs_mmolH2Om2day[n,,]}
         if (exists("gb")) {gb[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$gb_mmolH2Om2day[n,,]}
         # NPP (foliar, root, wood; gC/m2/day -> kgC/m2/s)
         if (exists("NPP_labile_FLX")) {NPP_labile_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$alloc_labile_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NPP_root_FLX")) {NPP_root_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$alloc_roots_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NPP_wood_FLX")) {NPP_wood_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$alloc_wood_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NPP_combinedfoliage_FLX")) {NPP_combinedfoliage_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$combined_alloc_foliage_gCm2day[n,,]* 1e-3 * (1/86400)}
         if (exists("NPP_foliage_FLX")) {NPP_foliage_FLX[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$alloc_foliage_gCm2day[n,,]* 1e-3 * (1/86400)}

         # NPP (fraction) and MRT years are requested to have same number of time steps as stocks and fluxes
         # This is awkward as no easy way to repeat specific elements without loop for variables which have no meaningful value at sub-annual timescales
         # (and are therefore calculated as annuals)
         for (q in seq(1, nos_quantiles)) {
              # MRT
              if (exists("labMRT")) {labMRT[grid_output$i_location[n],grid_output$j_location[n],q,]  = rep(grid_output$MTT_annual_labile_years[n,q,], each = steps_per_year)}
              if (exists("folMRT")) {folMRT[grid_output$i_location[n],grid_output$j_location[n],q,]  = rep(grid_output$MTT_annual_foliage_years[n,q,], each = steps_per_year)}
              if (exists("rooMRT")) {rooMRT[grid_output$i_location[n],grid_output$j_location[n],q,]  = rep(grid_output$MTT_annual_roots_years[n,q,], each = steps_per_year)}
              if (exists("wooMRT")) {wooMRT[grid_output$i_location[n],grid_output$j_location[n],q,]  = rep(grid_output$MTT_annual_wood_years[n,q,], each = steps_per_year)}
              if (exists("litMRT")) {litMRT[grid_output$i_location[n],grid_output$j_location[n],q,]  = rep(grid_output$MTT_annual_litter_years[n,q,], each = steps_per_year)}
              if (exists("wlitMRT")) {wlitMRT[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_output$MTT_annual_woodlitter_years[n,q,], each = steps_per_year)}
              if (exists("somMRT")) {somMRT[grid_output$i_location[n],grid_output$j_location[n],q,]  = rep(grid_output$MTT_annual_som_years[n,q,], each = steps_per_year)}
              if (exists("bioMRT")) {bioMRT[grid_output$i_location[n],grid_output$j_location[n],q,]  = rep(grid_output$MTT_annual_biomass_years[n,q,], each = steps_per_year)}
              if (exists("domMRT")) {domMRT[grid_output$i_location[n],grid_output$j_location[n],q,]  = rep(grid_output$MTT_annual_dom_years[n,q,], each = steps_per_year)}
              # NPP fractional allocation
              if (exists("fNPP")) {fNPP[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_output$NPP_foliage_fraction[grid_output$i_location[n],grid_output$j_location[n],q], each = length(PROJECT$model$timestep_days))}
              if (exists("rNPP")) {rNPP[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_output$NPP_roots_fraction[grid_output$i_location[n],grid_output$j_location[n],q], each = length(PROJECT$model$timestep_days))}
              if (exists("wNPP")) {wNPP[grid_output$i_location[n],grid_output$j_location[n],q,] = rep(grid_output$NPP_wood_fraction[grid_output$i_location[n],grid_output$j_location[n],q], each = length(PROJECT$model$timestep_days))}
         } # loop quantiles

         # Total outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day -> kgC/m2/s)
         if (exists("outflux_labile")) {outflux_labile[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_labile_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("outflux_foliage")) {outflux_foliage[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_foliage_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("outflux_wood")) {outflux_wood[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_wood_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("outflux_root")) {outflux_root[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_roots_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("outflux_litter")) {outflux_litter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_litter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("outflux_woodlitter")) {outflux_woodlitter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_woodlitter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("outflux_som")) {outflux_som[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_som_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("outflux_bio")) {outflux_bio[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_biomass_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("outflux_dom")) {outflux_dom[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$outflux_dom_gCm2day[n,,] * 1e-3 * (1/86400)}

         # Fire mortality outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day -> kgC/m2/s)
         if (exists("FIRElitter_labile")) {FIRElitter_labile[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_labile_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_foliage")) {FIRElitter_foliage[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_foliage_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_wood")) {FIRElitter_wood[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_wood_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_root")) {FIRElitter_root[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_roots_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_litter")) {FIRElitter_litter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_litter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_woodlitter")) {FIRElitter_woodlitter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_woodlitter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_bio")) {FIRElitter_bio[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIRElitter_biomass_gCm2day[n,,] * 1e-3 * (1/86400)}

         # Fire combustion outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day -> kgC/m2/s)
         if (exists("FIREemiss_labile")) {FIREemiss_labile[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_labile_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_foliage")) {FIREemiss_foliage[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_foliage_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_wood")) {FIREemiss_wood[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_wood_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_root")) {FIREemiss_root[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_roots_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_litter")) {FIREemiss_litter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_litter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_woodlitter")) {FIREemiss_woodlitter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_woodlitter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_som")) {FIREemiss_som[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_som_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_bio")) {FIREemiss_bio[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_biomass_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIREemiss_dom")) {FIREemiss_dom[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$FIREemiss_dom_gCm2day[n,,] * 1e-3 * (1/86400)}

         # HARVEST mortality outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day -> kgC/m2/s)
         if (exists("FIRElitter_labile")) {HARVESTlitter_labile[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTlitter_labile_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_foliage")) {HARVESTlitter_foliage[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTlitter_foliage_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_wood")) {HARVESTlitter_wood[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTlitter_wood_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_root")) {HARVESTlitter_root[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTlitter_roots_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("FIRElitter_bio")) {HARVESTlitter_bio[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTlitter_biomass_gCm2day[n,,] * 1e-3 * (1/86400)}

         # HARVEST extraction outflux (labile, foliar, wood, fine root, litter, soil; gC/m2/day -> kgC/m2/s)
         if (exists("HARVESTextracted_labile")) {HARVESTextracted_labile[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_labile_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("HARVESTextracted_foliage")) {HARVESTextracted_foliage[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_foliage_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("HARVESTextracted_wood")) {HARVESTextracted_wood[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_wood_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("HARVESTextracted_root")) {HARVESTextracted_root[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_roots_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("HARVESTextracted_litter")) {HARVESTextracted_litter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_litter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("HARVESTextracted_woodlitter")) {HARVESTextracted_woodlitter[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_woodlitter_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("HARVESTextracted_som")) {HARVESTextracted_som[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_som_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("HARVESTextracted_bio")) {HARVESTextracted_bio[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_biomass_gCm2day[n,,] * 1e-3 * (1/86400)}
         if (exists("HARVESTextracted_dom")) {HARVESTextracted_dom[grid_output$i_location[n],grid_output$j_location[n],,] = grid_output$HARVESTextracted_dom_gCm2day[n,,] * 1e-3 * (1/86400)}

     } # Does the file exist / has it been processed

} # site loop

###
## Define dimensions which will be used across all files

## define dimension
lat_dimen <- ncdim_def( "lat", units="degree north (-90->90)", latitude )
long_dimen <- ncdim_def( "lon", units="degree east (-180->180)", longitude )
time_dimen <- ncdim_def( "time", units="", 1:length(PROJECT$model$timestep_days))
quantile_dimen <- ncdim_def( "quantile", units="-", quantiles_wanted)
year_dimen <- ncdim_def( "year", units="", 1:nos_years)
npar_dimen <- ncdim_def( "nos_parameters", units="", 1:(max(PROJECT$model$nopars)+1)) # NOTE: +1 is to account for the log-likelihood

###
## Create C STATES file with timing information, landsea fraction and grid area
###

## define output variable
var0 = ncvar_def("Time", units = "d", longname = paste("Monthly time step given in days since 01/01/",PROJECT$start_year,sep=""),
                 dim=list(time_dimen), missval = -99999, prec="double", compression = 9)
var1 = ncvar_def("grid_area", units = "m2", longname = paste("Pixel area",sep=""),
                 dim=list(long_dimen,lat_dimen), missval = -99999, prec="double", compression = 9)
var2 = ncvar_def("land_fraction", units = "1", longname = paste("Fraction of pixel which is land",sep=""),
                 dim=list(long_dimen,lat_dimen), missval = -99999, prec="double", compression = 9)
# Define the output file name
output_name = paste(PROJECT$results_processedpath,output_prefix,"CSTOCK_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
# Delete if the file currently exists
if (file.exists(output_name)) {file.remove(output_name)}
# Create the empty file space
new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
# Load first variable into the file
# TIMING
ncvar_put(new_file, var0, drivers$met[,1])
# Grid area
ncvar_put(new_file, var1, grid_output$area_m2)
# Land fraction
ncvar_put(new_file, var2, grid_output$land_fraction)

# Close the existing file to ensure its written to file
nc_close(new_file)

###
## Re-open the file so that we can add to it a variable at a time
###

new_file <- nc_open( output_name, write=TRUE )

###
## ADD STATE VARIABLES
###

# LAI
if(exists("LAI")) {
   # Median
   var_new  = ncvar_def("lai", unit="m2.m-2", longname = "Leaf Area Index - Median ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  LAI)
}

# Labile
if(exists("LAB")) {
   var_new  = ncvar_def("cLabile_ensemble", unit="kg.m-2", longname = "Carbon in labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  LAB)
}

# Foliar
if(exists("FOL")) {
   var_new = ncvar_def("cLeaf_ensemble", unit="kg.m-2", longname = "Carbon in leaves - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  FOL)
}

# Fine root
if(exists("ROOT")) {
   var_new = ncvar_def("cFineRoot_ensemble", unit="kg.m-2", longname = "Carbon in fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  ROOT)
}

# Wood
if(exists("WOOD")) {
   var_new = ncvar_def("cWoodTotal_ensemble", unit="kg.m-2", longname = "Carbon in (AGB + BGB) wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  WOOD)
}

# Foliar + fine root litter
if(exists("LIT")) {
   var_new = ncvar_def("cLeafFineRootLitter_ensemble", unit="kg.m-2", longname = "Carbon in (Foliar + fine root) litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  LIT)
}

# Wood litter
if(exists("WLIT")) {
   var_new = ncvar_def("cWoodLitter_ensemble", unit="kg.m-2", longname = "Carbon in (wood) litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  WLIT)
}

# Soil organic matter
if(exists("SOIL")) {
   var_new = ncvar_def("cSOM_ensemble", unit="kg.m-2", longname = "Carbon in soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  SOIL)
}

# Dead organic matter
if(exists("DOM")) {
   var_new = ncvar_def("cDOM_ensemble", unit="kg.m-2", longname = "Carbon in leaf, fine root, wood litter, and soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  DOM)
}

# Biomass
if(exists("BIO")) {
   var_new = ncvar_def("cVeg_ensemble", unit="kg.m-2", longname = "Carbon in live biomass - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, BIO)
}

# Change in Biomass
if(exists("dBIO")) {
   var_new = ncvar_def("dcVeg_ensemble", unit="kg.m-2", longname = "Change in Carbon in live biomass since t=1 - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, dBIO)
}

# Change in Biomass
if(exists("dDOM")) {
   var_new = ncvar_def("dcDOM_ensemble", unit="kg.m-2", longname = "Change in Carbon in leaf, fine root, wood litter, and soil organic matter since t=1 - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, dDOM)
}

# TotalC
if(exists("TOT")) {
   var_new = ncvar_def("cTotal_ensemble", unit="kg.m-2", longname = "Carbon in live and dead organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, TOT)
}

###
## ADD STATE VARIABLES - mean
###

# Mean LAI
if(exists("MLAI")) {
   # Median
   var_new  = ncvar_def("mean_lai", unit="m2.m-2", longname = "Mean Leaf Area Index - Median ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )	# NOTE this returns a modified netcdf file handle
   ncvar_put(new_file, var_new,  MLAI)
}

# Mean Labile
if(exists("MLAB")) {
   var_new  = ncvar_def("mean_cLabile_ensemble", unit="kg.m-2", longname = "Mean carbon in labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  MLAB)
}

# Mean Foliar
if(exists("MFOL")) {
   var_new = ncvar_def("mean_cLeaf_ensemble", unit="kg.m-2", longname = "Mean carbon in leaves - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  MFOL)
}

# Mean Fine root
if(exists("MROOT")) {
   var_new = ncvar_def("mean_cFineRoot_ensemble", unit="kg.m-2", longname = "Mean carbon in fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  MROOT)
}

# Wood
if(exists("WOOD")) {
   var_new = ncvar_def("mean_cWoodTotal_ensemble", unit="kg.m-2", longname = "Mean carbon in (AGB + BGB) wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  MWOOD)
}

# Mean Foliar + fine root litter
if(exists("MLIT")) {
   var_new = ncvar_def("mean_cLeafFineRootLitter_ensemble", unit="kg.m-2", longname = "Mean carbon in (Foliar + fine root) litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  MLIT)
}

# Mean Wood litter
if(exists("MWLIT")) {
   var_new = ncvar_def("mean_cWoodLitter_ensemble", unit="kg.m-2", longname = "Mean carbon in (wood) litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  MWLIT)
}

# Mean Soil organic matter
if(exists("MSOIL")) {
   var_new = ncvar_def("mean_cSOM_ensemble", unit="kg.m-2", longname = "Mean carbon in soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  MSOIL)
}

# Mean Dead organic matter
if(exists("MDOM")) {
   var_new = ncvar_def("mean_cDOM_ensemble", unit="kg.m-2", longname = "Mean carbon in leaf, fine root, wood litter, and soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new,  MDOM)
}

# Mean Biomass
if(exists("MBIO")) {
   var_new = ncvar_def("mean_cVeg_ensemble", unit="kg.m-2", longname = "Mean carbon in live biomass - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, MBIO)
}

# Mean TotalC
if(exists("MTOT")) {
   var_new = ncvar_def("mean_cTotal_ensemble", unit="kg.m-2", longname = "Mean carbon in live and dead organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, MTOT)
}

###
## close the file to write to disk
###

nc_close(new_file)

###
## Create file for C FLUXES with timing information
###

# Define the output file name
output_name = paste(PROJECT$results_processedpath,output_prefix,"CFLUX_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
# Delete if the file currently exists
if (file.exists(output_name)) {file.remove(output_name)}
# Create the empty file space
new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
# Load first variable into the file
# TIMING
ncvar_put(new_file, var0, drivers$met[,1])
# Grid area
ncvar_put(new_file, var1, grid_output$area_m2)
# Land fraction
ncvar_put(new_file, var2, grid_output$land_fraction)

# Close the existing file to ensure its written to file
nc_close(new_file)

###
## Re-open the file so that we can add to it a variable at a time
###

new_file <- nc_open( output_name, write=TRUE )

###
## FLUXES - at time step
###

# GPP
if(exists("GPP")) {
   var_new  = ncvar_def("gpp_ensemble", unit="kg.m-2.s-1", longname = "Gross Primary Productivity - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, GPP)
}

# Autotrophic respiration
if(exists("RAU")) {
   var_new  = ncvar_def("ra_ensemble", unit="kg.m-2.s-1", longname = "Autotrophic (Plant) Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, RAU)
}

# Heterotrophic respiration
if(exists("RHE")) {
   var_new  = ncvar_def("rh_ensemble", unit="kg.m-2.s-1", longname = "Heterotrophic Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, RHE)
}

# Ecosystem respiration
if(exists("RECO")) {
   var_new  = ncvar_def("reco_ensemble", unit="kg.m-2.s-1", longname = "Ecosystem (Ra + Rh) Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, RECO)
}

# Net Primary Productivity
if(exists("NPP")) {
   var_new  = ncvar_def("npp_ensemble", unit="kg.m-2.s-1", longname = "Net Primary Productivity - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, NPP)
}

# Net Ecosystem Exchange
if(exists("NEE")) {
   var_new  = ncvar_def("nee_ensemble", unit="kg.m-2.s-1", longname = "Net Ecosystem Exchange - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, NEE)
}

# Net Biome Exchange
if(exists("NBE")) {
   var_new  = ncvar_def("nbe_ensemble", unit="kg.m-2.s-1", longname = "Net Biome Exchange (NEE + Fire) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, NBE)
}

# Net Biome Productivity
if(exists("NBP")) {
   var_new  = ncvar_def("nbp_ensemble", unit="kg.m-2.s-1", longname = "Net Biome Productivity (-NEE - Fire - fLuc) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, NBP)
}

# Fire emissions
if(exists("FIR")) {
   var_new  = ncvar_def("fFire_ensemble", unit="kg.m-2.s-1", longname = "Fire - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIR)
}

# Flux from forest loss
if(exists("HARV")) {
   var_new = ncvar_def("fLuc_ensemble", unit="kg.m-2.s-1", longname = "Forest harvest - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, HARV)
}

# Total outflux from labile
if(exists("outflux_labile")) {
   var_new  = ncvar_def("outflux_cLabile_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, outflux_labile)
}

# Total outflux from foliage
if(exists("outflux_foliage")) {
   var_new  = ncvar_def("outflux_cLeaf_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, outflux_foliage)
}

# Total outflux from root
if(exists("outflux_root")) {
   var_new  = ncvar_def("outflux_cFineRoot_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, outflux_root)
}

# Total outflux from wood
if(exists("outflux_wood")) {
   var_new  = ncvar_def("outflux_cWoodTotal_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, outflux_wood)
}

# Total outflux from litter
if(exists("outflux_litter")) {
   var_new  = ncvar_def("outflux_cLeafFineRootlitter_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from foliar and fine root litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, outflux_litter)
}

# Total outflux from wood litter
if(exists("outflux_woodlitter")) {
   var_new  = ncvar_def("outflux_cWoodlitter_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from wood litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, outflux_woodlitter)
}

# Total outflux from som
if(exists("outflux_som")) {
   var_new  = ncvar_def("outflux_cSOM_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, outflux_som)
}

# Biomass
if(exists("outflux_bio")) {
   var_new  = ncvar_def("outflux_cVeg_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from vegetation - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, outflux_bio)
}

# Dead Organic Matter
if(exists("outflux_dom")) {
   var_new  = ncvar_def("outflux_cDOM_ensemble", unit="kg.m-2.s-1", longname = "Total C output flux from dead organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, outflux_dom)
}

# Fire mortality outflux from labile
if(exists("FIRElitter_labile")) {
   var_new  = ncvar_def("FIRElitter_cLabile_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIRElitter_labile)
}

# Fire mortality outflux from foliage
if(exists("FIRElitter_foliage")) {
   var_new  = ncvar_def("FIRElitter_cLeaf_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIRElitter_foliage)
}

# Fire mortality outflux from root
if(exists("FIRElitter_root")) {
   var_new  = ncvar_def("FIRElitter_cFineRoot_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIRElitter_root)
}

# Fire mortality outflux from wood
if(exists("FIRElitter_wood")) {
   var_new  = ncvar_def("FIRElitter_cWoodTotal_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIRElitter_wood)
}

# Fire mortality outflux from litter
if(exists("FIRElitter_litter")) {
   var_new  = ncvar_def("FIRElitter_cLeafFineRootlitter_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from foliar and fine root litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIRElitter_litter)
}

# Fire mortality outflux from litter
if(exists("FIRElitter_woodlitter")) {
   var_new  = ncvar_def("FIRElitter_cWoodlitter_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from wood litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIRElitter_woodlitter)
}

# Fire mortality outflux from som
if(exists("FIRElitter_bio")) {
   var_new  = ncvar_def("FIRElitter_cVeg_ensemble", unit="kg.m-2.s-1", longname = "Fire mortality C output flux from Vegetation - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIRElitter_bio)
}

# Fire mortality outflux from labile
if(exists("FIREemiss_labile")) {
   var_new  = ncvar_def("FIREemiss_cLabile_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIREemiss_labile)
}

# Fire mortality outflux from foliage
if(exists("FIREemiss_foliage")) {
   var_new  = ncvar_def("FIREemiss_cLeaf_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIREemiss_foliage)
}

# Fire mortality outflux from root
if(exists("FIREemiss_root")) {
   var_new  = ncvar_def("FIREemiss_cFineRoot_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIREemiss_root)
}

# Fire mortality outflux from wood
if(exists("FIREemiss_wood")) {
   var_new  = ncvar_def("FIREemiss_cWoodTotal_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIREemiss_wood)
}

# Fire mortality outflux from litter
if(exists("FIREemiss_litterLAI")) {
   var_new  = ncvar_def("FIREemiss_cLeafFineRootlitter_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from foliar and fine root litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIREemiss_litter)
}

# Fire mortality outflux from wood litter
if(exists("FIREemiss_woodlitter")) {
   var_new  = ncvar_def("FIREemiss_cWoodlitter_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from wood litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIREemiss_woodlitter)
}

# Fire mortality outflux from som
if(exists("FIREemiss_som")) {
   var_new  = ncvar_def("FIREemiss_cSOM_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIREemiss_som)
}

# Fire mortality outflux from vegetation
if(exists("FIREemiss_bio")) {
   var_new  = ncvar_def("FIREemiss_cVeg_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from vegetation - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIREemiss_bio)
}

# Fire mortality outflux from dom
if(exists("FIREemiss_dom")) {
   var_new  = ncvar_def("FIREemiss_cDOM_ensemble", unit="kg.m-2.s-1", longname = "Fire combusted C output flux from dead organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, FIREemiss_dom)
}

# Harvest litter outflux from labile
if(exists("HARVESTlitter_labile")) {
   var_new  = ncvar_def("HARVESTlitter_cLabile_ensemble", unit="kg.m-2.s-1", longname = "Harvest litter C output flux from labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, HARVESTlitter_labile)
}

# Harvest litter outflux from foliage
if(exists("HARVESTlitter_foliage")) {
   var_new  = ncvar_def("HARVESTlitter_cLeaf_ensemble", unit="kg.m-2.s-1", longname = "Harvest litter C output flux from foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, HARVESTlitter_foliage)
}

# Harvest litter outflux from root
if(exists("HARVESTlitter_root")) {
   var_new  = ncvar_def("HARVESTlitter_cFineRoot_ensemble", unit="kg.m-2.s-1", longname = "Harvest litter C output flux from fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, HARVESTlitter_root)
}

# Harvest litter outflux from wood
if(exists("HARVESTlitter_wood")) {
   var_new  = ncvar_def("HARVESTlitter_cWoodTotal_ensemble", unit="kg.m-2.s-1", longname = "Harvest litter C output flux from wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, HARVESTlitter_wood)
}

# Harvest litter outflux from vegetation
if(exists("HARVESTlitter_bio")) {
   var_new  = ncvar_def("HARVESTlitter_cVeg_ensemble", unit="kg.m-2.s-1", longname = "Harvest litter C output flux from Vegetation - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, HARVESTlitter_bio)
}

# Harvest extracted outflux from labile
if(exists("HARVESTextracted_labile")) {
   var_new  = ncvar_def("HARVESTextracted_cLabile_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, HARVESTextracted_labile)
}

# Harvest extracted outflux from foliage
if(exists("HARVESTextracted_foliage")) {
   var_new  = ncvar_def("HARVESTextracted_cLeaf_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, HARVESTextracted_foliage)
}

# Harvest extracted outflux from root
if(exists("HARVESTextracted_root")) {
   var_new  = ncvar_def("HARVESTextracted_cFineRoot_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, HARVESTextracted_root)
}

# Harvest extracted outflux from wood
if(exists("HARVESTextracted_wood")) {
   var_new  = ncvar_def("HARVESTextracted_cWoodTotal_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, HARVESTextracted_wood)
}

# Harvest extracted outflux from litter
if(exists("HARVESTextracted_litter")) {
   var_new  = ncvar_def("HARVESTextracted_cLeafFineRootlitter_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from foliar and fine root litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, HARVESTextracted_litter)
}

# Harvest extracted outflux from wood litter
if(exists("HARVESTextracted_woodlitter")) {
   var_new  = ncvar_def("HARVESTextracted_cWoodlitter_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from wood litter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, HARVESTextracted_woodlitter)
}

# Harvest extracted outflux from som
if(exists("HARVESTextracted_som")) {
   var_new  = ncvar_def("HARVESTextracted_cSOM_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from soil organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, HARVESTextracted_som)
}

# Harvest extracted outflux from vegetation
if(exists("HARVESTextracted_bio")) {
   var_new  = ncvar_def("HARVESTextracted_cVeg_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from vegetation - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, HARVESTextracted_bio)
}

# Harvest extracted outflux from dom
if(exists("HARVESTextracted_dom")) {
   var_new  = ncvar_def("HARVESTextracted_cDOM_ensemble", unit="kg.m-2.s-1", longname = "Harvest extracted C output flux from dead organic matter - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, HARVESTextracted_dom)
}

###
## FLUXES - at annual time step
###

# Annual GPP
if(exists("AGPP")) {
   var_new  = ncvar_def("annual_gpp_ensemble", unit="kg.m-2.s-1", longname = "Mean Annual Gross Primary Productivity - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, AGPP)
}

# Annual Autotrophic respiration
if(exists("ARAU")) {
   var_new  = ncvar_def("annual_ra_ensemble", unit="kg.m-2.s-1", longname = "Mean Annual Autotrophic (Plant) Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, ARAU)
}

# Annual Heterotrophic respiration
if(exists("ARHE")) {
   var_new  = ncvar_def("annual_rh_ensemble", unit="kg.m-2.s-1", longname = "Mean Annual Heterotrophic Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, ARHE)
}

# Annual Ecosystem respiration
if(exists("ARECO")) {
   var_new  = ncvar_def("annual_reco_ensemble", unit="kg.m-2.s-1", longname = "Mean Annual Ecosystem (Ra + Rh) Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, ARECO)
}

# Annual Net Primary Productivity
if(exists("ANPP")) {
   var_new  = ncvar_def("annual_npp_ensemble", unit="kg.m-2.s-1", longname = "Mean Annual Net Primary Productivity - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, ANPP)
}

# Annual Net Ecosystem Exchange
if(exists("ANEE")) {
   var_new  = ncvar_def("annual_nee_ensemble", unit="kg.m-2.s-1", longname = "Mean Annual Net Ecosystem Exchange - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, ANEE)
}

# Annual Net Biome Exchange
if(exists("ANBE")) {
   var_new  = ncvar_def("annual_nbe_ensemble", unit="kg.m-2.s-1", longname = "Mean Annual Net Biome Exchange (NEE + Fire) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, ANBE)
}

# Annual Net Biome Productivity
if(exists("ANBP")) {
   var_new  = ncvar_def("annual_nbp_ensemble", unit="kg.m-2.s-1", longname = "Mean Annual Net Biome Productivity (-NEE - Fire - fLuc) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, ANBP)
}

# Annual Fire emissions
if(exists("AFIR")) {
   var_new  = ncvar_def("annual_fFire_ensemble", unit="kg.m-2.s-1", longname = "Mean Annual Fire - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, AFIR)
}

# Annual Flux from forest loss
if(exists("AHARV")) {
   var_new = ncvar_def("annual_fLuc_ensemble", unit="kg.m-2.s-1", longname = "Mean Annual Forest harvest - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,year_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, AHARV)
}

###
## FLUXES - timeseries mean
###

# Mean GPP
if(exists("MGPP")) {
   var_new  = ncvar_def("mean_gpp_ensemble", unit="kg.m-2.s-1", longname = "Mean Gross Primary Productivity - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, MGPP)
}

# Mean Autotrophic respiration
if(exists("MRAU")) {
   var_new  = ncvar_def("mean_ra_ensemble", unit="kg.m-2.s-1", longname = "Mean Autotrophic (Plant) Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, MRAU)
}

# Mean Heterotrophic respiration
if(exists("MRHE")) {
   var_new  = ncvar_def("mean_rh_ensemble", unit="kg.m-2.s-1", longname = "Mean Heterotrophic Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, MRHE)
}

# Mean Ecosystem respiration
if(exists("MRECO")) {
   var_new  = ncvar_def("mean_reco_ensemble", unit="kg.m-2.s-1", longname = "Mean Ecosystem (Ra + Rh) Respiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, MRECO)
}

# Mean Net Primary Productivity
if(exists("MNPP")) {
   var_new  = ncvar_def("mean_npp_ensemble", unit="kg.m-2.s-1", longname = "Mean Net Primary Productivity - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, MNPP)
}

# Mean Net Ecosystem Exchange
if(exists("MNEE")) {
   var_new  = ncvar_def("mean_nee_ensemble", unit="kg.m-2.s-1", longname = "Mean Net Ecosystem Exchange - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, MNEE)
}

# Mean Net Biome Exchange
if(exists("MNBE")) {
   var_new  = ncvar_def("mean_nbe_ensemble", unit="kg.m-2.s-1", longname = "Mean Net Biome Exchange (NEE + Fire) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, MNBE)
}

# Mean Net Biome Productivity
if(exists("MNBP")) {
   var_new  = ncvar_def("mean_nbp_ensemble", unit="kg.m-2.s-1", longname = "Mean Net Biome Productivity (-NEE - Fire - fLuc) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, MNBP)
}

# Mean Fire emissions
if(exists("MFIR")) {
   var_new  = ncvar_def("mean_fFire_ensemble", unit="kg.m-2.s-1", longname = "Mean Fire - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, MFIR)
}

# Mean Flux from forest loss
if(exists("MHARV")) {
   var_new = ncvar_def("mean_fLuc_ensemble", unit="kg.m-2.s-1", longname = "Mean Forest harvest - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, MHARV)
}


###
## close the file to write to disk
###

nc_close(new_file)

###
## Create file for NPP and MRT with timing information
###

# Define the output file name
output_name = paste(PROJECT$results_processedpath,output_prefix,"NPP_MRT_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
# Delete if the file currently exists
if (file.exists(output_name)) {file.remove(output_name)}
# Create the empty file space
new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
# Load first variable into the file
# TIMING
ncvar_put(new_file, var0, drivers$met[,1])
# Grid area
ncvar_put(new_file, var1, grid_output$area_m2)
# Land fraction
ncvar_put(new_file, var2, grid_output$land_fraction)

# Close the existing file to ensure its written to file
nc_close(new_file)

###
## Re-open the file so that we can add to it a variable at a time
###

new_file <- nc_open( output_name, write=TRUE )

###
## Mean residence times and NPP allocation / fluxes
###

## Mean Residence Times
# Labile
if(exists("labMRT")) {
   var_new = ncvar_def("MTT_lab_ensemble", unit="year", longname = "Mean Labile Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, labMRT)
}

# Foliar
if(exists("folMRT")) {
   var_new = ncvar_def("MTT_fol_ensemble", unit="year", longname = "Mean Foliar Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, folMRT)
}

# Fine root
if(exists("rooMRT")) {
   var_new = ncvar_def("MTT_root_ensemble", unit="year", longname = "Mean fine root Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, rooMRT)
}

# Wood
if(exists("wooMRT")) {
   var_new = ncvar_def("MTT_wood_ensemble", unit="year", longname = "Mean wood Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, wooMRT)
}

# Fine litter (fol + fine root)
if(exists("litMRT")) {
   var_new = ncvar_def("MTT_lit_ensemble", unit="year", longname = "Mean foliar + fine root litter Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, litMRT)
}

# Wood litter
if(exists("wlitMRT")) {
   var_new = ncvar_def("MTT_wlit_ensemble", unit="year", longname = "Mean wood litter Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, wlitMRT)
}

# Soil
if(exists("somMRT")) {
   var_new = ncvar_def("MTT_som_ensemble", unit="year", longname = "Mean Soil Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, somMRT)
}

# Dead Organic Matter
if(exists("domMRT")) {
   var_new = ncvar_def("MTT_dom_ensemble", unit="year", longname = "Mean Soil + litter Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, domMRT)
}

# Biomass / vegetation
if(exists("bioMRT")) {
   var_new = ncvar_def("MTT_veg_ensemble", unit="year", longname = "Mean biomass / vegetation Transit Time - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, bioMRT)
}

## NPP allocation fractions
# Foliar
if(exists("fNPP")) {
   var_new = ncvar_def("NPP_fol_ensemble", unit="1", longname = "Fraction of Net Primary Productivity to foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, fNPP)
}

# Fine root
if(exists("rNPP")) {
   var_new = ncvar_def("NPP_root_ensemble", unit="1", longname = "Fraction of Net Primary Productivity to fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, rNPP)
}

# Wood
if(exists("wNPP")) {
   var_new = ncvar_def("NPP_wood_ensemble", unit="1", longname = "Fraction of Net Primary Productivity to wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, wNPP)
}


## NPP allocation fluxes
# Labile
if(exists("NPP_labile_FLX")) {
   var_new= ncvar_def("NPP_labile_flx_ensemble", unit="kg.m-2.s-1", longname = "Net Primary Productivity to labile - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, NPP_labile_FLX)
}

# Direct allocation to Foliar
if(exists("NPP_foliage_FLX")) {
   var_new = ncvar_def("NPP_direct_fol_flx_ensemble", unit="kg.m-2.s-1", longname = "Net Primary Productivity direct to foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, NPP_foliage_FLX)
}

# Combined direct and via labile allocation of NPP to foliage
if(exists("NPP_combinedfoliage_FLX")) {
   var_new = ncvar_def("NPP_fol_flx_ensemble", unit="kg.m-2.s-1", longname = "Both direct and via labile Net Primary Productivity to foliage - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, NPP_combinedfoliage_FLX)
}

# Fine root
if(exists("NPP_root_FLX")) {
   var_new= ncvar_def("NPP_root_flx_ensemble", unit="kg.m-2.s-1", longname = "Net Primary Productivity to fine root - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, NPP_root_FLX)
}

# Wood
if(exists("NPP_wood_FLX")) {
   var_new = ncvar_def("NPP_wood_flx_ensemble", unit="kg.m-2.s-1", longname = "Net Primary Productivity to wood - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, NPP_wood_FLX)
}
###
## close the file to write to disk
###

nc_close(new_file)

###
## Create file for parameters and clusters
###

# Define the output file name
output_name = paste(PROJECT$results_processedpath,output_prefix,"PARS_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
# Delete if the file currently exists
if (file.exists(output_name)) {file.remove(output_name)}
# Create the parameter variable definition
var_pars = ncvar_def("parameters", unit="-", longname = "Parameter array + log-likelihood, see specific model version for desccription and units - Ensemble", dim=list(long_dimen,lat_dimen,npar_dimen,quantile_dimen), missval = -99999, prec="double",compression = 9)
# Create the empty file space, using the already defined timing variable
new_file=nc_create(filename=output_name, vars=list(var_pars), force_v4 = TRUE)
# Load first variable into the file
# PARAMETERS
ncvar_put(new_file, var_pars, grid_output$parameters)

# Close the existing file to ensure its written to file
nc_close(new_file)

###
## Re-open the file so that we can add to it a variable at a time
###

new_file <- nc_open( output_name, write=TRUE )

###
## Other parameter or cluster related information
###

## Cluster by affinity propogation process
# parameters only, i.e. initial conditions not included
if(exists(x = "pars_clusters", where = grid_output)) {
   var_new = ncvar_def("pars_clusters", unit="-", longname = "Affinity propogation clustering based on process parameters (i.e. not including initial conditions)", dim=list(long_dimen,lat_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, grid_output$pars_clusters)
}

## Cluster by affinity propogation process
# parameters only, i.e. initial conditions not included
if(exists(x = "clusters", where = grid_output)) {
   var_new = ncvar_def("clusters", unit="-", longname = "Affinity propogation clustering based on process parameters and initial conditions", dim=list(long_dimen,lat_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, grid_output$clusters)
}

###
## close the file to write to disk
###

nc_close(new_file)

###
## Create file for DRIVERS and OBSERVATIONS with timing information
###

# Define the output file name
output_name = paste(PROJECT$results_processedpath,output_prefix,"DRIVERS_OBS_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
# Delete if the file currently exists
if (file.exists(output_name)) {file.remove(output_name)}
# Create the empty file space
new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
# Load first variable into the file
# TIMING
ncvar_put(new_file, var0, drivers$met[,1])
# Grid area
ncvar_put(new_file, var1, grid_output$area_m2)
# Land fraction
ncvar_put(new_file, var2, grid_output$land_fraction)

# Close the existing file to ensure its written to file
nc_close(new_file)

###
## Re-open the file so that we can add to it a variable at a time
###

new_file <- nc_open( output_name, write=TRUE )

## DRIVERS
# Minimum air temperature
var_new = ncvar_def("tas_min", unit="K", longname = "Mean daily minimum near surface air temperature", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, AIRT_MIN)

# Maximum air temperature
var_new = ncvar_def("tas_max", unit="K", longname = "Mean daily maximum near surface air temperature", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, AIRT_MAX)

# Short Wave radiation
var_new = ncvar_def("rsds", unit="W.m-2", longname = "Mean downwelling short wave radiation", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, SWRAD)

# Atmospheric CO2 concentration
var_new = ncvar_def("co2", unit="ppm", longname = "Mean atmospheric CO2 concentration", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, CO2)

# Precipitation
var_new = ncvar_def("pr", unit="kg.m-2.s-1", longname = "Mean precipitation - combined liquid and solid phase", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, PRECIP)

# Forest loss fraction
var_new = ncvar_def("forest_loss_fraction", unit="1", longname = "Forest loss fraction", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, FLOSS_FRAC)

# Burnt fraction
var_new = ncvar_def("burntArea", unit="1", longname = "Burnt fraction", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, BURNT_FRAC)

# Wind Speed
var_new = ncvar_def("wsp", unit="m.s-1", longname = "Mean wind speed", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, WINDSPD)

# Vapour pressure deficit
var_new = ncvar_def("vpd", unit="Pa", longname = "Mean vapour pressure deficit", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, VPD)

## STATE OBSERVATIONS
# Leaf area index
var_new = ncvar_def("LAI_OBS", unit="m-2.m-2", longname = "Observed Leaf area index", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, LAI_OBS)

# Leaf area index Uncertainty
var_new = ncvar_def("LAI_UNC_OBS", unit="m-2.m-2", longname = "Uncertainty on observed Leaf area index", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, LAI_UNC_OBS)

# Wood stock
var_new = ncvar_def("WOOD_OBS", unit="g.m-2", longname = "Observed wood stock C (above + below + coarse root)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, WOOD_OBS)

# Wood stock uncertainty
var_new = ncvar_def("WOOD_UNC_OBS", unit="g.m-2", longname = "Uncertainty on observed wood stock C (above + below + coarse root)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, WOOD_UNC_OBS)

# Soil stock
var_new = ncvar_def("SOIL_OBS", unit="g.m-2", longname = "Observed soil stock C (assumed to include wood litter)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, SOIL_OBS)

# Soil stock uncertainty
var_new = ncvar_def("SOIL_UNC_OBS", unit="g.m-2", longname = "Uncertainty on observed wood stock C (assumed to include wood litter)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, SOIL_UNC_OBS)

## FLUX OBSERVATIONS
# GPP
var_new = ncvar_def("GPP_OBS", unit="g.m-2.d-1", longname = "Observed gross primary productivity of C (GPP)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, GPP_OBS)

# GPP uncertainty
var_new = ncvar_def("GPP_UNC_OBS", unit="g.m-2.d-1", longname = "Uncertainty on observed gross primary productivity of C (GPP)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, GPP_UNC_OBS)

# NEE
var_new = ncvar_def("NEE_OBS", unit="g.m-2.d-1", longname = "Observed net ecosystem exchange of C (NEE = Reco - GPP)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, NEE_OBS)

# NEE uncertainty
var_new = ncvar_def("NEE_UNC_OBS", unit="g.m-2.d-1", longname = "Uncertainty on observed net ecosystem exchange of C (NEE = Reco - GPP)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, NEE_UNC_OBS)

# NBE
var_new = ncvar_def("NBE_OBS", unit="g.m-2.d-1", longname = "Observed net biome exchange of C (NEE = Reco - GPP + FIRE)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, NBE_OBS)

# NBE uncertainty
var_new = ncvar_def("NBE_UNC_OBS", unit="g.m-2.d-1", longname = "Uncertainty on observed net biome exchange of C (NEE = Reco - GPP + FIRE)", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, NBE_UNC_OBS)

# Fire
var_new = ncvar_def("FIRE_OBS", unit="g.m-2.d-1", longname = "Observed C emission due to fire combustion", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, FIRE_OBS)

# Fire uncertainty
var_new = ncvar_def("FIRE_UNC_OBS", unit="g.m-2.d-1", longname = "Uncertainty on observed C emission due to fire combustion", dim=list(long_dimen,lat_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
new_file <- ncvar_add( new_file, var_new )
ncvar_put(new_file, var_new, FIRE_UNC_OBS)

## close the file to write to disk
nc_close(new_file)

###
## Create file for H2O FLUXES with timing information
###

# Define the output file name
output_name = paste(PROJECT$results_processedpath,output_prefix,"H2OFLUX_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
# Delete if the file currently exists
if (file.exists(output_name)) {file.remove(output_name)}
# Create the empty file space
new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
# Load first variable into the file
# TIMING
ncvar_put(new_file, var0, drivers$met[,1])
# Grid area
ncvar_put(new_file, var1, grid_output$area_m2)
# Land fraction
ncvar_put(new_file, var2, grid_output$land_fraction)

# Close the existing file to ensure its written to file
nc_close(new_file)

###
## Re-open the file so that we can add to it a variable at a time
###

new_file <- nc_open( output_name, write=TRUE )

###
## H2O FLUXES
###

# Evapotranspiration
if(exists("ET")) {
   var_new  = ncvar_def("evapotrans_ensemble", unit="kg.m-2.s-1", longname = "Evapotranspiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, ET)
}

# Transpiration
if(exists("Etrans")) {
   var_new  = ncvar_def("trans_ensemble", unit="kg.m-2.s-1", longname = "Transpiration - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, Etrans)
}

# Soil evaporation
if(exists("Esoil")) {
   var_new  = ncvar_def("evapsoil_ensemble", unit="kg.m-2.s-1", longname = "Soil evaporation - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, Esoil)
}

# Canopy intercepted rainfall evaporation
if(exists("Ewetevap")) {
   var_new  = ncvar_def("evapwetcanopy_ensemble", unit="kg.m-2.s-1", longname = "Canopy intercepted rainfall evaporation - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, Ewetcanopy)
}

# Surface runoff
if(exists("runoff")) {
   var_new  = ncvar_def("runoff_ensemble", unit="kg.m-2.s-1", longname = "Soil surface runoff - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, runoff)
}

# Drainage from bottom of soil column
if(exists("underflow")) {
   var_new  = ncvar_def("underflow_ensemble", unit="kg.m-2.s-1", longname = "Drainage from bottom of soil column - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, underflow)
}

# Total drainage from soil surface and bottom of soil column
if(exists("total_drainage")) {
   var_new  = ncvar_def("mrro_ensemble", unit="kg.m-2.s-1", longname = "Total drainage from soil surface and bottom of soil column - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, total_drainage)
}

###
## close the file to write to disk
###

nc_close(new_file)

###
## Create file for H2O STATES with timing information
###

# Define the output file name
output_name = paste(PROJECT$results_processedpath,output_prefix,"H2OSTATES_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
# Delete if the file currently exists
if (file.exists(output_name)) {file.remove(output_name)}
# Create the empty file space
new_file=nc_create(filename=output_name, vars=list(var0,var1,var2), force_v4 = TRUE)
# Load first variable into the file
# TIMING
ncvar_put(new_file, var0, drivers$met[,1])
# Grid area
ncvar_put(new_file, var1, grid_output$area_m2)
# Land fraction
ncvar_put(new_file, var2, grid_output$land_fraction)

# Close the existing file to ensure its written to file
nc_close(new_file)

###
## Re-open the file so that we can add to it a variable at a time
###

new_file <- nc_open( output_name, write=TRUE )

###
## H2O STATES
###

# Soil surface water content
if(exists("SurfWater")) {
   var_new  = ncvar_def("SurfWater_ensemble", unit="kg.m-2", longname = "Soil water content (0-30 cm) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, SurfWater)
}

# Soil surface snow cover
if(exists("SNOW")) {
   var_new  = ncvar_def("SNOW_ensemble", unit="kg.m-2", longname = "Soil surface snow cover - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, SNOW)
}

###
## close the file to write to disk
###

nc_close(new_file)

###
## Create file for Biophysical diagnostics with timing information
###

# Define the output file name
output_name = paste(PROJECT$results_processedpath,output_prefix,"BIOPHYSDIAG_",PROJECT$start_year,"_",PROJECT$end_year,output_suffix,".nc",sep="")
# Delete if the file currently exists
if (file.exists(output_name)) {file.remove(output_name)}
# Create the empty file space, using the already defined timing variable
new_file=nc_create(filename=output_name, vars=list(var0), force_v4 = TRUE)
# Load first variable into the file
# TIMING
ncvar_put(new_file, var0, drivers$met[,1])

# Close the existing file to ensure its written to file
nc_close(new_file)

###
## Re-open the file so that we can add to it a variable at a time
###

new_file <- nc_open( output_name, write=TRUE )

###
## Biophysical diagnostics
###

# CiCa
if(exists("CiCa")) {
   var_new = ncvar_def("CiCa_ensemble", unit="1", longname = "Internal:Ambiant CO2 ratio - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, CiCa)
}

# Ecosystem WUE (GPP/ET)
if(exists("wue_eco")) {
   var_new = ncvar_def("wue_eco_ensemble", unit="g.kg", longname = "Ecosystem water use efficiency (GPP/ET) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, wue_eco)
}

# Plant WUE (GPP/Etrans)
if(exists("wue_eco")) {
   var_new = ncvar_def("wue_plant_ensemble", unit="g.kg", longname = "Plant water use efficiency (GPP/Etrans) - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, wue_eco)
}

# Soil water potential weighted by root access
if(exists("wSWP")) {
   var_new  = ncvar_def("wSWP_ensemble", unit="MPa", longname = "Soil water potential weighted by root access - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, wSWP)
}

# Canopy absorbed photosynthetically active radiation
if(exists("APAR")) {
   var_new  = ncvar_def("APAR_ensemble", unit="MJ.m2.d-1", longname = "Canopy absorbed photosynthatically active radiation - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, APAR)
}

# Ratio of actual stomatal conductance to that allowed by water supply to canopy
if(exists("gs_demand_supply_ratio")) {
   var_new  = ncvar_def("gs_demand_supply_ratio_ensemble", unit="1", longname = "Ratio of actual stomatal conductance to potential allowed by water supply to canopy - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, gs_demand_supply_ratio)
}

# Stomatal conductance for water
if(exists("gs")) {
   var_new  = ncvar_def("gs_ensemble", unit="mmol.m-2.d-1", longname = "Stomatal conductance for water - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, gs)
}

# Boundary layer conductance for water
if(exists("gb")) {
   var_new  = ncvar_def("gb_ensemble", unit="mmol.m-2.d-1", longname = "Boundary conductance for water - Ensemble", dim=list(long_dimen,lat_dimen,quantile_dimen,time_dimen), missval = -99999, prec="double",compression = 9)
   new_file <- ncvar_add( new_file, var_new )
   ncvar_put(new_file, var_new, gb)
}

###
## close the file to write to disk
###

nc_close(new_file)
