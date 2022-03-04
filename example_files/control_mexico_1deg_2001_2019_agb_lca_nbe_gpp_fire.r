
###
## Script to control the creation of files and submission to eddie for
## CARDAMOM-DALEC
###

local({
## Prepare
## Compute
setwd("<your_cardamom_directory_here>")
## Print result
})

###
## Options

## Load needed libraries and internal functions
source("./R_functions/load_all_cardamom_functions.r")

## use parallel functions?
use_parallel=TRUE
numWorkers=12 # number of cores to assign to parallel job
## about you
username="<uun>"
home_computer="ssh.geos.ed.ac.uk"
## projname
# Give a runid
projname="Mexico_1deg_agb_lca_gpp_fire_nbe_example"

## Language
# i.e. "Fortran", "C"
language="Fortran"

## Compiler options (Fortan only)
compiler="ifort" # "ifort" or "gfortran"
timing=FALSE
debug=FALSE

## Model
# i.e. "DALEC_CDEA", "DALEC_CDEA_FR", "DALECcrop", "AT_DALEC", "DALEC_GSI_FR"
model="DALEC_CDEA_ACM2_BUCKET"
pft_specific_parameters=FALSE

## MDF method
# i.e. MHMCMC or other
method="MHMCMC"

## Land cover map
# which land cover map to use
use_lcm="ECMWF" # choices are "CORINE2006", "LCM2007", "CORINE2006_1km", "ECMWF"
pft_wanted=FALSE

## Met paths
#path_to_met_source="/exports/csce/datastore/geos/groups/gcel/ECMWF/ERA5/0.125deg_global/"
#path_to_met_source="/exports/csce/datastore/geos/groups/gcel/ISI-MIP/isimip3a/historical/meteorology/global_0.5deg/"
path_to_met_source="/exports/csce/datastore/geos/groups/gcel/Trendy_v9_met/monthly/"
path_to_lai="/disk/scratch/local.2/copernicus/LAI_0.125deg/" # located on racadal
path_to_crop_management="/exports/csce/datastore/geos/groups/gcel/Crop_calendar_dataset_Sacks/netCDF_5_min/Wheat/"
path_to_sand_clay="/exports/csce/datastore/geos/groups/gcel/SoilGrids/processed/global_5km/"
path_to_Csom="/exports/csce/datastore/geos/groups/gcel/SoilGrids/processed/global_5km/"
path_to_Cwood_inc = "/exports/csce/datastore/geos/groups/gcel/cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg/wood_productivity/"
path_to_Cwood_mortality = "/exports/csce/datastore/geos/groups/gcel/cssp_rainfor_amazon_brazil/rainfor_leeds_data/modified_for_CARDAMOM/1deg/wood_mortality/"
path_to_Cwood="/exports/csce/datastore/geos/groups/gcel/AGB/ESA_CCI_BIOMASS/ESA_CCI_AGB_0.125deg/"
path_to_Cwood_initial="/exports/csce/datastore/geos/groups/gcel/AGB/ESA_CCI_BIOMASS/ESA_CCI_AGB_5km/"
path_to_Cwood_potential=" "
path_to_gleam="/exports/csce/datastore/geos/groups/gcel/GLEAM/v3.1a_soilmoisture_prior_TLS/"
path_to_nbe = "/exports/csce/datastore/geos/groups/gcel/AtmosphericInversions/combined_nbe/global_1deg_monthly/"
path_to_gpp = "/exports/csce/datastore/geos/groups/gcel/GPP_ESTIMATES/combined_gpp/global_1deg_monthly/"
path_to_fire = "/exports/csce/datastore/geos/groups/gcel/FIRE_ESTIMATES/combined_fire/global_1deg_monthly/"
path_to_forestry="/exports/csce/datastore/geos/groups/gcel/GlobalForestWatch/global_0.125deg/"
path_to_burnt_area="/exports/csce/datastore/geos/groups/gcel/BurnedArea/MCD64A1/global_0.125deg/"
path_to_lca = "/exports/csce/datastore/geos/groups/gcel/TraitMaps/Butler/LCA/global_1deg/"
path_to_landsea = "default" # "default" or "filepath/to/raster/covering/the/target/area.tif"
path_to_site_obs=" "
met_interp=FALSE # linear interpolation if needed

## Data streams
met_source="trendy_v9" # "trendy_v9" or "ERA" or "isimip3a"
lai_source="COPERNICUS" # "MODIS" or "site_specific"
Csom_source="SoilGrids" # "HWSD" or "SoilGrids" or "site_specific"
soilwater_initial_source = " " # initial soil water fraction (m3/m3)
sand_clay_source="SoilGrids" # "HWSD" or "site_specific" or "SoilGrids"
Evap_source=" "
Cwood_inc_source = " " # " " or "site_specific" or "Rainfor"
Cwood_mortality_source = " " # "site_specific" or " " or "Rainfor"
GPP_source="Global_Combined" 	# " " or "site_specific"
fire_source="Global_Combined" # " " or "site_specific" or "Global_Combined"
Reco_source=" " 	# " " or "site_specific"
NEE_source=" " # "site_specific" 	# " " or "site_specific"
nbe_source = "Global_Combined" # "site_specific" or "GEOSCHEM" or "Global_Combined" or " "
# i.e. single value valid for beginning of simulation
Cfol_initial_source=" " #"site_specific" 	# " " or "site_specific"
Cwood_initial_source=" " #"site_specific" 	# " " or "site_specific"
Croots_initial_source=" " #"site_specific" 	# " " or "site_specific"
Clit_initial_source=" " #"site_specific"  	# " " or "site_specific"
# i.e. time series of stock estimates
Cfol_stock_source=" " 	# " " or "site_specific"
Cfolmax_stock_source=" " 	# " " or "site_specific"
Cwood_stock_source="ESA_CCI_Biomass" 	# " " or "site_specific" or "ESA_CCI_Biomass"
Cstem_stock_source=" "      # " " or "site_specific"
Cbranch_stock_source=" "      # " " or "site_specific"
Cagb_stock_source=" " 	# " " or "site_specific"
Ccoarseroot_stock_source=" " 	# " " or "site_specific"
Croots_stock_source=" " 	# " " or "site_specific"
Clit_stock_source=" "  	# " " or "site_specific"
Csom_stock_source=" "  	# " " or "site_specific"
lca_source = "Butler" # "Butler" or " " or "site_specific"
# Steady state attractor
Cwood_potential_source = " " # "site_specific" or ""
# Management drivers
burnt_area_source="MCD64A1" # " " or "MCD64A1" or "GFED4" or "site_specific"
deforestation_source="GFW" # " ", "site_specific" or "GFW"
crop_management_source=" " # "_" or "site_specific" or "sacks_crop_calendar"
snow_source=" "

## sites for analysis
# start year
years_to_do=as.character(c(2001:2019)) 
# is this run "site" level or over a "grid"?
cardamom_type="grid"
# if type = "grid" then what resolution in (m for UK), (degrees for global)
cardamom_resolution=1
# which grid are we on? "UK" or "wgs84"
cardamom_grid_type="wgs84"

# site names if specific locations e.g. "UKGri"
sites_cardamom="Mexico"
# lat/long of sites, if type = "grid"then these these are bottom left and top right corners
sites_cardamom_lat=c(12,34) # c(52,58)
sites_cardamom_long=c(-120,-84) # c(-3,0.5)
# timestep mode, currently "daily" or "monthly" or "weekly"
timestep_type="monthly"
select_country = TRUE # use function available_functions() for valid counties in sites_cardamom

## Define the project setup
# NOTE: if these are not set CARDAMOM will ask you for them
request_nos_chains = 3        # Number of chains CARDAMOM should run for each location
request_nos_samples = 100e6   # Total number of parameter samples / iterations to be explored
request_nos_subsamples = 1e3  # Number of parameter sets to be sub-sampled from the chain
request_use_server = TRUE     # Use remote server? Currently coded for UoE Eddie.
request_runtime = 48          # How many hours of compute to request per job. Only applied for running on remote server
request_compile_server = TRUE # Copy and compile current source code on remote server
request_use_EDCs = TRUE       # Use EDCs

## Stage
# stage -1 : Fix or create project first time (load source to eddie)
# stage  1 : Create met / obs containing files for the specifc project
# stage  2 : Submit the project to eddie
# stage  3 : Copy back results and process vectors
# stage  4 : Do some standard result checking
stage=-1
repair=0 # to force (=1) re-run processed results or driver files if they already exist
grid_override=FALSE # force site specific files to be saved and figures to be generated when in "grid" operation

##
# Call CARDAMOM with specific stages
cardamom(projname,model,method,stage)

rm(list=ls()) ; gc() ; gc()

#find -name "*PARS" -size -160k -delete
#qstat | grep "UK_public" | wc -l
