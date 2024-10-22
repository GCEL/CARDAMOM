
###
## Example script to control CARDAMOM for a single site run
###

# Set working directory in which the CARDAMOM code base can be found
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")

###
## Options

## Load needed libraries and internal functions
source("./R_functions/load_all_cardamom_functions.r")

## projname
# Give a runid
projname="FI-Hyy_example"
#projname="FI-Hyy_example_O3"
#projname="FI-Hyy_example_Ofast"

## Language
# i.e. "Fortran", "C"
language="Fortran"
## Compiler options (Fortan only)
compiler="ifort" #"ifort", "gfortran"
compiler_optimisation = "-Ofast"
timing=FALSE
debug=FALSE

## about you (only valid if working on UoE remote server)
username="lsmallma" # put your Edinburgh uun here
home_computer="ssh.geos.ed.ac.uk"
sshpass_key_home = "~/.ssh/id_rsa_geos.pub" # location of passkey on remote server

## use parallel functions?
use_parallel = TRUE
numWorkers = 3 # number of cores to assign to parallel job

## Model - which DALEC 
# see "MODEL_DESCRIPTIONS.md" for available models
model="DALEC.4."
pft_specific_parameters=FALSE # impacts crop model only

## MDF method
# i.e. MHMCMC or other. Note that while MHMCMC remains the only available choice this has been updated to an AP-MCMC
method="MHMCMC"

## Land cover map
# which land cover map to use. Used to determine the land / sea mask
use_lcm="ECMWF" # coded choices exist for other maps however only "ECMWF" map is provided with source code
pft_wanted=FALSE # Impacts crop model only
path_to_landsea="default" # If gridded analysis, any raster layer with >0 values will be taken as the mask area. To ignore = "default"

## Met paths
path_to_met_source=" "
#path_to_met_source="/exports/csce/datastore/geos/groups/gcel/Trendy_v11_met/monthly/"
#path_to_met_source="/exports/csce/datastore/geos/groups/gcel/ECMWF/ERA5/0.125deg_global/"
path_to_lai=" " #"/exports/csce/datastore/geos/groups/gcel/LAI_ESTIMATES/MCD15A2H.061/global_0.0625deg/"
path_to_fapar = " "
path_to_crop_management=" "
path_to_sand_clay=" " #"/exports/csce/datastore/geos/groups/gcel/SoilGrids/version2/processed/global_5km/"
path_to_Csom=" " #"/exports/csce/datastore/geos/groups/gcel/SoilGrids/version2/processed/global_5km/"
path_to_Cwood_inc = ""
path_to_Cwood_mortality = ""
path_to_Cwood=" " #"/exports/csce/datastore/geos/groups/gcel/AGB/ESA_CCI_BIOMASS/ESA_CCI_AGB_0.125deg/"
path_to_Cwood_initial=" "
path_to_Cwood_potential=" "
path_to_gleam=" "
path_to_nbe = " "
path_to_gpp = " "
path_to_fire = " "
path_to_forestry=" " #"/exports/csce/datastore/geos/groups/gcel/GlobalForestWatch/global_0.0625deg/"
path_to_burnt_area=" " #"/exports/csce/datastore/geos/groups/gcel/BurnedArea/MCD64A1/global_0.0625deg/"
path_to_lca = "/exports/csce/datastore/geos/groups/gcel/TraitMaps/Butler/LCA/global_1deg/"
path_to_co2 = "/exports/csce/datastore/geos/groups/gcel/Trendy_v11_met/global_CO2/"
path_to_site_obs="./example_files/inputs/"
path_to_landsea = "default"
met_interp=TRUE

## Data streams - The currently coded data streams which can be used to drive or constrain the models
met_source="site_specific" # "trendy_v9" or "trendy_v11" or "ERA" or "isimip3a" or "site_specific"
lai_source="site_specific" # "COPERNICUS" or "MODIS" or "site_specific"
fapar_source=" " # "COPERNICUS" or "MODIS" or "site_specific"
Csom_source="site_specific" # "SoilGrids" or "SoilGrids_v2" or "HWSD" or "site_specific
sand_clay_source="site_specific" # "SoilGrids" or "SoilGrids_v2" or "HWSD" or "site_specific
soilwater_initial_source = " " # initial soil water fraction (m3/m3)
Evap_source="site_specific"        # " " or "site_specific"
Cwood_inc_source = " " # "site_specific" or " " or "Rainfor"
Cwood_mortality_source = " " # "site_specific" or " " or "Rainfor"
fire_source=" " # " " or "site_specific" or "Global_Combined"
GPP_source=" " 	# " " or "site_specific" or "Global_Combined"
Reco_source=" " 	# " " or "site_specific"
NEE_source="site_specific" # " " or "site_specific"
nbe_source = " " # " " or "site_specific" or "Global_Combined" or "GEOSCHEM" or "OCO2MIP"
harvest_source = ""
foliage_to_litter_source = " " # " " or "site_specific"
# i.e. single value valid for beginning of simulation
Cfol_initial_source=" " #"site_specific" 	# " " or "site_specific"
Cwood_initial_source=" " #"site_specific" 	# " " or "site_specific"
Croots_initial_source=" " #"site_specific" 	# " " or "site_specific"
Clit_initial_source=" " #"site_specific"  	# " " or "site_specific"
# i.e. time series of stock estimates
Cfol_stock_source=" " 	# " " or "site_specific"
Cfolmax_stock_source=" " 	# " " or "site_specific"
Cwood_stock_source="site_specific" 	# " " or "site_specific" or "McNicol" or "Saatchi_2021" or "ESA_CCI_Biomass"
Cstem_stock_source=" "      # " " or "site_specific"
Cagb_stock_source=" " 	# " " or "site_specific"
Ccoarseroot_stock_source=" " 	# " " or "site_specific"
Croots_stock_source=" " 	# " " or "site_specific"
Clit_stock_source=" "  	# " " or "site_specific"
Csom_stock_source=" "  	# " " or "site_specific"
# Parameter priors
lca_source = " " # "Butler" or " " or "site_specific"
frac_Cwood_coarse_root_source = "" # " " or "site_specific"
minLWP_source = "" # " " or "site_specific"
# Steady state attractor
Cwood_potential_source = " " # "site_specific" or ""
# Management drivers
burnt_area_source=" " # " " or "MCD64A1" or "GFED4" or "site_specific"
deforestation_source=" " # " ", "site_specific" or "GFW"
crop_management_source=" " # "_" or "site_specific" or "sacks_crop_calendar"
snow_source=" "

## sites for analysis
# start year
years_to_do=as.character(c(1999:2014)) 
# is this run "site" level or over a "grid"?
cardamom_type="site"
cardamom_grid_type=" " # "UK" or "wgs84", no value needed if site run
# if type = "grid" then at what spatial resolution (UK = m, wgs84 = degree)?
cardamom_resolution=1e5

# site names if specific locations e.g. "UKGri"
sites_cardamom=c("FI-Hyy")
# lat/long of sites, if type = "grid" then these these are bottom left and top right corners
sites_cardamom_lat=61.84741
sites_cardamom_long=24.29477
timestep_type="monthly"
select_country = FALSE # If gridded run and path_to_landsea = "default", 
                       # select country based on site_cardamom. Use function
                       # available_countries() for compatible country names.

## Define the project setup
# NOTE: if these are not set CARDAMOM will ask you for them
request_nos_chains = 3        # Number of chains CARDAMOM should run for each location
request_nos_samples = 10e6   # Total number of parameter samples / iterations to be explored
request_nos_subsamples = 1e3  # Number of parameter sets to be sub-sampled from the chain
request_use_server = FALSE    # Use remote server? Currently coded for UoE Eddie.
request_use_local_slurm = FALSE# Only applies if request_use_server == FALSE
request_runtime = 48          # How many hours of compute to request per job. Only applied for running on remote server
request_compile_server = FALSE# Copy and compile current source code on remote server
request_compile_local = TRUE  # Compile local copy of the source code 
request_use_EDCs = TRUE       # Use EDCs
request_extended_mcmc = FALSE # Extend the current MCMC by adding a further request_nos_extended_samples + request_nos_samples
request_nos_extended_samples = 90e6 # If request_extened_mcmc == TRUE then this is the number of additional proposals to be made
request_cost_function_scaling = 2 # 0 = Default, no normaliation of the likelihood score
                                  # 1 = Normalisation of the likelihood score by sample size
                                  # 2 = Normalisation of the likelihood score by sqrt(sample size)
                                  # 3 = Normalisation of the likelihood score by log(sample size) 

## Stage
# stage -1 : Create project first time (load source to eddie)
# stage  0 : Re-compile source code, if needed but without re-creating the PROJECT related infofile.RData
# stage  1 : Create met / obs containing files for the specifc project
# stage  2 : Submit the project to eddie
# stage  3 : Copy back results and process vectors
# stage  4 : Do some standard figure creation (and further processing for gridded analysis)
# stage  5 : Generic dump of RESULTS_PROCESSED files to netcdf
stage=4
repair=1 # to force (=1) re-run processed results or driver files if they already exist
grid_override=FALSE # force site specific files to be saved and figures to be generated when in "grid" operation

##
# Call CARDAMOM with specific stages
cardamom(projname,model,method,stage)


