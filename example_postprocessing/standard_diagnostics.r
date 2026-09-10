
###
## Script to control various functions defined in "functions_for_standard_diagnostics.r"
## Author(s): T. Luke Smallman (t.l.smallman@ed.ac.uk) 
## Created: 29/05/2025
## Modification History:
## 1) ... 
###

# Assumptions:
# 1) That all put file are in formats which could be read into CARDAMOM.
# 2) Consistent with 1) functions used to read files are drawn from CARDAMOM.

###
## Define various file paths

# Set working directory for CARDAMOM codebase
cardamom_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/"
# Set directory for the CARDAMOM project to be analysed
project_dir = "~/gcel_ceph/cardamom_analyses/lsmallma/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.004_MHMCMC/global_0.5deg_dalec4_trendyv15_LCA_TWB_GPP_fAPAR_hashimoto_SGDB/"
#project_dir = "~/gcel_ceph/cardamom_analyses/lsmallma/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.004_MHMCMC/global_0.5deg_dalec4_trendyv14_LCA_TWB_GPP_fAPAR/"
# Set the current working directory for this script
script_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/example_postprocessing/"
# Set output directory, if left empty then the FIGURES directory of the specific project will be assumed to be the location
output_dir = ""

###
## Settings for the analysis

# Median estimate array index
mid_quant = 5
# Lower and lower quantiles array index
low_quant = 1  ; high_quant = 9 # 95 %
#low_quant = 2 ; high_quant = 8 # 68, i.e. 1-SD %

###
## Load functions in memory

# Load pre-defined CARDAMOM functions
setwd(cardamom_dir)
source(paste(cardamom_dir,"R_functions/load_all_cardamom_functions.r",sep=""))
# Load pre-defined diagnositc functions
source(paste(script_dir,"functions_for_standard_diagnostics.r",sep=""))

###
## Load CARDAMOM project into memory

# Load project information
load(paste(project_dir,"infofile.RData",sep=""))
# Load project gridded output file
load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))

###
## Calculate repeat use information

# Useful values from the PROJECT object
calc_useful()

###
## Load forcing information

# Load additional variables not calculated based on the forcings
# and observations going into the analysis
load_forcings_to_grid()

###
## Begin diagnostics

## Set sink for printed information
sink(file = paste(output_dir,"standard_diagnostics_output_",PROJECT$name,"_",Sys.Date(),".txt",sep=""))

## Do analyses

# Generate global and zonal carbon budgets, consistent with TRENDY model intercomparison
global_zonal_budget()

# Conduct K-means cluster analysis, save output and plot.
# Useful for exploring what 'PFTs' you could have with constrained analyses
k_means_clustering(nos_clusters=9) # default is number of JULES PFTs
# Plot summary information based on clusters
plot_clustering()

# Conduct between pixel correlations for parameters and then parameters~key variables
tmp = grid_output$landmask ; tmp[tmp == 0] = NA
key_variables_parameters_spatial_correlation(tmp,"Global")
                  
# Generate time series plots of the aggregated fluxes and anomalies - modelled
# These should be made to include each of the reccap regions in a single panel plot
var_and_units    = c("wSWP_MPa","SurfWater_kgH2Om2","CiCa","lai_m2m2","nbp_PgCyr", "nbe_PgCyr", "nee_PgCyr","npp_PgCyr", "gpp_PgCyr", "reco_PgCyr","rhet_PgCyr", "rhet_litter_PgCyr", "rhet_som_PgCyr", 
                    "rauto_PgCyr", "fire_PgCyr", "harvest_PgCyr","grazing_PgCyr","combined_alloc_foliage_PgCyr","alloc_roots_PgCyr","dnbp_PgCyr", "dnbe_PgCyr", 
                    "dnee_PgCyr","dnpp_PgCyr", "dgpp_PgCyr","dreco_PgCyr","drhet_PgCyr", "drhet_litter_PgCyr","drhet_som_PgCyr", 
                    "drauto_PgCyr", "dfire_PgCyr", "dharvest_PgCyr","dgrazing_PgCyr","biomass_PgC", "dom_PgC", "labile_PgC","foliage_PgC", "roots_PgC", 
                    "litter_PgC", "som_PgC","dCbiomass_PgC","dCdom_PgC", "dClabile_PgC","dCfoliage_PgC","dCroots_PgC", "dClitter_PgC","dCsom_PgC",
                    "dCiCa","dwSWP_MPa","dSurfWater_kgH2Om2")
outfile_var_name    = c("wSWP","SurfWater","CiCa","LAI","NBP", "NBE", "NEE", "NPP", "GPP", "Reco","Rhet", "Rhet_litter", "Rhet_som", "Rauto", "Fire", "Harvest",
                        "Grazing","NPPflux_foliage","NPPflux_roots","NBP_anomaly", "NBE_anomaly", "NEE_anomaly", "NPP_anomaly", 
                        "GPP_anomaly", "Reco_anomaly","Rhet_anomaly", "Rhet_litter_anomaly","Rhet_som_anomaly", "Rauto_anomaly", "Fire_anomaly", "Harvest_anomaly",
                        "Grazing_anomaly","Biomass", "DOM", "Labile", "Foliage", "FineRoots", "Litter", "SOM","Biomass_anomaly", "DOM_anomaly", "Labile_anomaly",
                        "Foliage_anomaly", "FineRoots_anomaly", "Litter_anomaly", "SOM_anomaly","CiCa_anomaly","wSWP_anomaly","SurfWater_anomaly")
outfile_var_units    = c(expression('(MPa)'),expression(paste('kgH2Om'^-2,sep="")),expression('(0-1)'),expression(paste('m'^2,'m'^-2,sep="")),expression('(PgC/yr)'),
                         expression('(PgC/yr)'),expression('(PgC/yr)'),expression('(PgC/yr)'),expression('(PgC/yr)'),expression('(PgC/yr)'),
                         expression('(PgC/yr)'),expression('(PgC/yr)'),expression('(PgC/yr)'),expression('(PgC/yr)'),expression('(PgC/yr)'),
                         expression('(PgC/yr)'),expression('(PgC/yr)'),expression('(PgC/yr)'),expression('(PgC/yr)'),
                         expression('(PgC/yr)'),expression('(PgC/yr)'),expression('(PgC/yr)'),expression('(PgC/yr)'),
                         expression('(PgC/yr)'),expression('(PgC/yr)'),expression('(PgC/yr)'),expression('(PgC/yr)'),
                         expression('(PgC/yr)'),expression('(PgC/yr)'),expression('(PgC/yr)'),expression('(PgC/yr)'),
                         expression('(PgC/yr)'),expression('(PgC)'),expression('(PgC)'),expression('(PgC)'),
                         expression('(PgC)'),expression('(PgC)'),expression('(PgC)'),
                         expression('(PgC)'),expression('(PgC)'),expression('(PgC)'),expression('(PgC)'),expression('(PgC)'),
                         expression('(PgC)'),expression('(PgC)'),expression('(PgC)'),expression('(0-1)'),expression('(MPa)'),
                         expression(paste('kgH2Om'^-2,sep="")))
for (i in seq(1, length(var_and_units))) { 
     create_spatially_aggregate_mean_annual_timeseries_and_anomaly(do_global = TRUE, global_lab="Global", do_obs = TRUE, outfile_prefix = "zonal",
                                                                   region_list[-3], region_list[-3], 
                                                                   var_and_units[i], outfile_var_name[i], outfile_var_units[i]) 
     create_spatially_aggregate_mean_annual_timeseries_and_anomaly_per_unit_area(do_global = TRUE, global_lab="Global", do_obs = TRUE, outfile_prefix = "zonal",
                                                                                 region_list[-3], region_list[-3], 
                                                                                 var_and_units[i], outfile_var_name[i], outfile_var_units[i]) 

}

# Do seasonal cycles for model outputs
# These should be made to include each of the reccap regions in a single panel plot
var_and_units    = c("wSWP_MPa","SurfWater_kgH2Om2","CiCa","lai_m2m2","nbp_PgCday", "nbe_PgCday", "nee_PgCday","npp_PgCday", "gpp_PgCday", 
                     "reco_PgCday","rhet_PgCday", "rhet_litter_PgCday", "rhet_som_PgCday", "rauto_PgCday", "fire_PgCday", "harvest_PgCday",
                     "grazing_PgCday","combined_alloc_foliage_PgCday","alloc_roots_PgCday")
outfile_var_name    = c("wSWP","SurfWater","CiCa","LAI","NBP", "NBE", "NEE", "NPP", "GPP", "Reco","Rhet", "Rhet_litter", "Rhet_som", "Rauto", "Fire", "Harvest",
                        "Grazing","NPPflux_foliage","NPPflux_roots")
outfile_var_units    = c(expression('(MPa)'),expression(paste('kgH2Om'^-2,sep="")),expression('(0-1)'),expression(paste('m'^2,'m'^-2,sep="")),expression('(PgC/day)'),
                         expression('(PgC/day)'),expression('(PgC/day)'),expression('(PgC/day)'),expression('(PgC/day)'),expression('(PgC/day)'),
                         expression('(PgC/day)'),expression('(PgC/day)'),expression('(PgC/day)'),expression('(PgC/day)'),
                         expression('(PgC/day)'),expression('(PgC/day)'),expression('(PgC/day)'),expression('(PgC/day)'),
                         expression('(PgC/day)'))
for (i in seq(1, length(var_and_units))) { 
     create_spatially_aggregate_seasonal_cycles(do_global = TRUE, global_lab="Global", do_obs = TRUE, outfile_prefix = "zonal",
                                                region_list[-3], region_list[-3], 
                                                var_and_units[i], outfile_var_name[i], outfile_var_units[i]) 
     create_spatially_aggregate_seasonal_cycles_per_unit_area(do_global = TRUE, global_lab="Global", do_obs = TRUE, outfile_prefix = "zonal",
                                                              region_list[-3], region_list[-3], 
                                                              var_and_units[i], outfile_var_name[i], outfile_var_units[i]) 
}

# Do timeseries and anomaly plots - forcings
var_and_units     = c("daily_min_temperature_C","daily_max_temperature_C","mean_temperature_C","sw_radiation_MJm2day","precipitation_kgH2Om2s","mean_vpd_Pa","biomass_removal_fraction",
                     "burned_fraction","mean_wind_speed_ms")
outfile_var_name  = c("Min_Air_Temperature","Max_Air_Temperature","Air_Temperature","SW_radiation","Precipitation","VPD","LUC","BA","Wind_Speed")
outfile_var_units = c(expression('(Celsius)'),expression('(Celsius)'),expression('(Celsius)'),expression(paste('MJm'^-2,'d'^-1,sep="")),
                      expression(paste('kgH2Om'^-2,'s'^-1,sep="")),expression('(Pa)'),expression('(0-1)'),expression('(0-1)'),expression(paste('ms'^-1,sep="")))
for (i in seq(1, length(var_and_units))) { 
     create_spatially_aggregate_mean_annual_timeseries_and_anomaly_forcings(do_global = FALSE, do_obs = TRUE, outfile_prefix = "zonal",
                                                                            zonal_names[c(-3,-6)], outfile_zonal_names[c(-3,-6)], 
                                                                            var_and_units[i], outfile_var_name[i], outfile_var_units[i]) 
}

# NPP
npp_plots()
# MRT
mrt_plots()
# Disturbance impacts on mrt
disturbance_and_mrt()
# General summary plots
summary_plots()

###
## Load independent observations

# The intention in this section is to load independent datasets which we can use later for evaluation.
# The code assumes that these have been converted into the CARDAMOM consistent formats including units and naming conventions,
# allowing use of existing CARDAMOM R code functions.

## Time varying with uncertainty 

# Net Biome Exchange (gC/m2/day)
nbe_name = "OCO2MIPv10"
nbe_source = "Gridded_nc" # "Gridded_nc" or "Gridded_tif" or " "
path_to_nbe = "/exports/geos.ed.ac.uk/gcel/spatial_datasets/net_biome_exchange/OCO2MIP/v10/global_1deg/" # Where are the files located
# Gross Primary Production (gC/m2/day)
gpp_name = "InvTOMCAT" #"MOD17A3HGF"
gpp_source = "Gridded_tif" # "Gridded_tif" # "Gridded_nc" or "Gridded_tif" or " "
path_to_gpp = "/home/lsmallma/gcel_ceph/spatial_datasets/gross_primary_production/InvTOMCAT/global_2.8deg/" # "/exports/geos.ed.ac.uk/gcel/spatial_datasets/gross_primary_production/MOD17A3HGF/collection_6.1/global_0.5deg/" # Where are the files located
# Evapotranspiration (kgH2O/m2/day)
et_name = " "
et_source = " " # "Gridded_nc" or "Gridded_tif" or " "
path_to_et = "" # Where are the files located
# Fire carbon emissions (gC/m2/day)
fire_name = "GFEDv5_GFASv1.2"
fire_source = "Gridded_nc" # "Gridded_nc" or "Gridded_tif" or " "
path_to_fire = "/exports/geos.ed.ac.uk/gcel/spatial_datasets/fire_C_emissions/merged_GFED_GFAS/global_0.5deg/" # Where are the files located
# Wood stock (gC/m2)
Cwood_name = " " # "Xu2021"
Cwood_stock_source = " " # "Gridded_tif" # " " or "Gridded_nc" or "Gridded_tif" or "site_specific"
path_to_Cwood      = " " # "/exports/geos.ed.ac.uk/gcel/spatial_datasets/total_woody_biomass/Xu2021/all/" # Where are the files located

## Static with uncertainty

# Leaf Carbon per unit leaf Area (LCA, gC/m2)
## NOT IN USE
lca_name = "Butler2016"
lca_source         = "Gridded_tif" # " " or "Gridded_nc" or "Gridded_tif" or "site_specific"
path_to_lca        = "/exports/geos.ed.ac.uk/gcel/spatial_datasets/LCA/Butler/global_0.5deg/" # Where are the files located

# Call the actual function to load independent evaluation datasets
eval_data = load_evaluation_datasets()
# Independent evaluaton
do_independent_evaluation(eval_data)

###
## Save work done for later

save(grid_output,eval_data, file = paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""), 
     compress = "gzip", compression_level = 9)

# End sink
sink()




