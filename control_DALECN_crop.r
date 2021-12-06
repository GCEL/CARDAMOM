
###
## Script to control the creation of files and submission to eddie for 
## ACM-TESSEL-DALEC
###

local({
## Prepare
## Compute
setwd("/exports/csce/datastore/geos/groups/ATEC/INTERIM_WORKINGS/gcel_github/")
## Print result
})

###
## Options

## Load needed libraries and internal functions
source("./R_functions/load_all_cardamom_functions.r")

## use parallel functions?
use_parallel=TRUE
numWorkers=5 # number of cores to assign to parallel job
## about you
username="arevill"
home_computer="racadal.geos.ed.ac.uk"
## projname
# Give a runid
projname="DALECN_Crop_single_point" 

## Language
# i.e. "Fortran", "C"
language="Fortran"

## Compiler options (Fortan only)
compiler="gfortran" # "ifort" or "gfortran"
timing=FALSE
debug=FALSE

## Model    
# i.e. "DALEC_CDEA", "DALEC_CDEA_FR", "DALECcrop", "AT_DALEC", "DALEC_GSI_FR"
model="DALECN_CROP"
pft_specific_parameters=FALSE

## MDF method
# i.e. MHMCMC or other
method="MHMCMC"

## Land cover map
# which land cover map to use
use_lcm="ECMWF" # choices are "CORINE2006", "LCM2007", "CORINE2006_1km", "ECMWF"
pft_wanted=FALSE

## Met paths
#path_to_met_source="/disk/scratch/local.2/lsmallma/ECMWF/ERA-Interim/0.25deg_global/"
path_to_met_source="./input_data/"
#path_to_met_source="/home/lsmallma/gcel/CHESS_v1.0/met_data/"
#path_to_met_source="/home/lsmallma/gcel/princeton/1.0deg/3hourly/"
#path_to_lai="/home/lsmallma/gcel/MODIS_LAI/global_2001_2016_0.25deg/"
#path_to_lai="/disk/scratch/local.2/copernicus/LAI_0.25x0.25/"
#path_to_lai="/disk/scratch/local.2/copernicus/LAI_1km_linked/"
#path_to_crop_management=""
#path_to_sand_clay="/exports/csce/datastore/geos/groups/gcel/SoilGrids/processed/global_5km/"
#path_to_Csom="/exports/csce/datastore/geos/groups/gcel/SoilGrids/processed/global_5km/"
#path_to_Csom="/home/lsmallma/gcel/HWSD/processed_file/global_0.25deg/"
#path_to_sand_clay="/home/lsmallma/gcel/HWSD/processed_file/global_0.25deg/"
#path_to_forestry="/home/lsmallma/gcel/GlobalForestWatch/global_0.125_degree/"
#path_to_Cwood=" "
#path_to_Cwood_initial="/home/lsmallma/gcel/lsmallma_BIOMASS_maps/Avitabile_AGB_0.25/"
#path_to_Cwood_potential=" "
#path_to_biomass="/home/lsmallma/gcel/lsmallma_BIOMASS_maps/Avitabile_AGB_0.25/"
#path_to_gleam="/home/lsmallma/gcel/GLEAM/v3.1a_soilmoisture_prior_TLS/"
#path_to_burnt_area="/exports/csce/datastore/geos/groups/gcel/BurnedArea/MCD64A1/global_0.125_degree/"
#path_to_nbe = " "
#path_to_lca = "/exports/csce/datastore/geos/groups/gcel/TraitMaps/Butler/LCA/global_1deg/"
path_to_site_obs="./input_data/"
met_interp=FALSE

## Data streams
met_source="site_specific" # "PRINCETON" or "ECMWF"
lai_source="site_specific" # "MODIS" or "site_specific"
Csom_source=" " #"HWSD" # "HWSD" or "site_specific"
sand_clay_source="SoilGrids" # HWSD or "site_specific
soilwater_initial_source = " "
nbe_source = " "
Evap_source=" "
woodinc_source=" " 	# " " or "site_specific"
GPP_source=" " 	# " " or "site_specific"
Reco_source=" " 	# " " or "site_specific"
NEE_source=" " # "site_specific" 	# " " or "site_specific"
# i.e. single value valid for beginning of simulation
Cfol_initial_source=" " #"site_specific" 	# " " or "site_specific"
Cwood_initial_source=" " #"site_specific" 	# " " or "site_specific"
Croots_initial_source=" " #"site_specific" 	# " " or "site_specific"
Clit_initial_source=" " #"site_specific"  	# " " or "site_specific"
# i.e. time series of stock estimates
Cfol_stock_source=" " 	# " " or "site_specific"
Cfolmax_stock_source=" " 	# " " or "site_specific"
Cwood_stock_source="site_specific" 	# " " or "site_specific" or "mpi_biomass"
Cstem_stock_source=" "      # " " or "site_specific"
Cbranch_stock_source=" "      # " " or "site_specific"
Cagb_stock_source=" " 	# " " or "site_specific"
Ccoarseroot_stock_source=" " 	# " " or "site_specific"
Croots_stock_source=" " 	# " " or "site_specific"
Clit_stock_source=" "  	# " " or "site_specific"
Csom_stock_source="site_specific"  	# " " or "site_specific"
lca_source = " " # "Butler" or " " or "site_specific"
# Steady state attractor
Cwood_potential_source = " " # "site_specific" or ""
# Management drivers
burnt_area_source="site_specific"# " " or "GFED4" or "site_specific"
deforestation_source=" " # " ", "site_specific" or "maryland_forestry_commission"
crop_management_source=" " # "_" or "site_specific" or "sacks_crop_calendar"
snow_source=" "

## sites for analysis
# start year
years_to_do=as.character(c(2020:2021)) # c("2000","2001","2002","2003","2004","2005","2006","2007","2008","2009")
# is this run "site" level or over a "grid"?
cardamom_type="site"
# if type = "grid" then what resolution in (m for UK), (degrees for global)?
cardamom_resolution=0.25
# which grid are we on? "UK" or "wgs84"
cardamom_grid_type="wgs84"

# site names if specific locations e.g. "UKGri"
sites_cardamom="west_fortune"
# lat/long of sites, if type = "grid"then these these are bottom left and top right corners
sites_cardamom_lat=56.0044 # c(52,58) , 
sites_cardamom_long=-2.7460 # c(-3,0.5)
# timestep mode, currently "daily" or "monthly" or "weekly"
timestep_type="daily"

## Stage 5: Driver modifications
# Define a proportional change to weather drivers.
airt_factor  = 1  # NOTE: this is also impact air temperature used in GSI model
swrad_factor = 1
co2_factor   = 1
rainfall_factor = 1
wind_spd_factor = 1
vpd_factor      = 1
# Define a proportional change on disturbance drivers.
# NOTE: at this point the only the intensity is impacted rather than frequency of events
deforestation_factor = 1
burnt_area_factor    = 1

## Define the project setup
# NOTE: if these are not set CARDAMOM will ask you for them
request_nos_chains = 3        # Number of chains CARDAMOM should run for each location
request_nos_samples = 100000   # Total number of parameter samples / iterations to be explored
request_nos_subsamples = 1e3  # Number of parameter sets to be sub-sampled from the chain
request_use_server = FALSE     # Use remote server? Currently coded for UoE Eddie.
request_runtime = 12          # How many hours of compute to request per job. Only applied for running on remote server
request_compile_server = FALSE # Copy and compile current source code on remote server
request_use_EDCs = TRUE       # Use EDCs
   
## Stage
# stage -1 : Fix or create project first time (load source to eddie)
# stage  1 : Create met / obs containing files for the specifc project
# stage  2 : Submit the project to eddie
# stage  3 : Copy back results and process vectors
# stage  4 ; Do some standard result checking
# stage  5 : Prepare parameters for DTESSEL???
stage=-1
repair=1 # to force (=1) re-run processed results or driver files if they already exist
grid_override=FALSE # force site specific files to be saved and figures to be generated when in "grid" operation

##
# Call CARDAMOM with specific stages
cardamom(projname,model,method,stage)

rm(list=ls()) ; gc() ; gc()

#find -name "*PARS" -size -160k -delete
#qstat | grep "UK_public" | wc -l
