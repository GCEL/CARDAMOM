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
# Function to check the for the presence of options from control files,
# where appropriate set default values and provide warnings for obvious errors
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

check_control_file_defaults<-function(paths) {

  ## Set defaults incase missing, NOTE: <<- to assign global

  
  # Analysis options
  if (exists("select_country") == FALSE)                {select_country <<- FALSE}
  if (exists("met_interp") == FALSE)                    {met_interp <<- FALSE}
  if (exists("pft_specific_parameters") == FALSE)       {pft_specific_parameters <<- FALSE}
  if (exists("cardamom_grid_type") == FALSE)            {cardamom_grid_type <<- "epsg:4326"}
  if (cardamom_grid_type == " " | cardamom_grid_type == "") {cardamom_grid_type <<- "epsg:4326"}
  # Site combined forcings and assimilated data path
  if (exists("path_to_site_obs") == FALSE)              {path_to_site_obs <<- " "}
  # Forcings data paths 
  if (exists("path_to_landsea") == FALSE)               {path_to_landsea <<- "default"}
  if (exists("path_to_met_source") == FALSE)            {path_to_met_source <<- " "}
  if (exists("path_to_sand_clay") == FALSE)             {path_to_sand_clay <<- " "}
  if (exists("path_to_forestry") == FALSE)              {path_to_forestry <<- " "}
  if (exists("path_to_burnt_area") == FALSE)            {path_to_burnt_area <<- " "}
  if (exists("path_to_co2") == FALSE)                   {path_to_co2 <<- "./R_functions/"}
  if (exists("path_to_lai_change") == FALSE)            {path_to_lai_change <<- " "}
  # Assimilated data paths 
  if (exists("path_to_lai") == FALSE)                   {path_to_lai <<- " "}
  if (exists("path_to_fapar") == FALSE)                 {path_to_fapar <<- " "}
  if (exists("path_to_crop_management") == FALSE)       {path_to_crop_management <<- " "}
  if (exists("path_to_RhetQ10") == FALSE)               {path_to_RhetQ10 <<- " "}  
  if (exists("path_to_MTTsom") == FALSE)                {path_to_MTTsom <<- " "}  
  if (exists("path_to_Csom") == FALSE)                  {path_to_Csom <<- " "}
  if (exists("path_to_Cwood_change") == FALSE)          {path_to_Cwood_change <<- " "}
  if (exists("path_to_Cwood_growth") == FALSE)          {path_to_Cwood_growth <<- " "}
  if (exists("path_to_Cwood_mortality") == FALSE)       {path_to_Cwood_mortality <<- " "}
  if (exists("path_to_Cwood") == FALSE)                 {path_to_Cwood <<- " "}
  if (exists("path_to_Cagb_change") == FALSE)           {path_to_Cagb_change <<- " "}
  if (exists("path_to_Cagb_growth") == FALSE)           {path_to_Cagb_growth <<- " "}
  if (exists("path_to_Cagb_mortality") == FALSE)        {path_to_Cagb_mortality <<- " "}
  if (exists("path_to_Cagb") == FALSE)                  {path_to_Cagb <<- " "}
  if (exists("path_to_Cwood_initial") == FALSE)         {path_to_Cwood_initial <<- " "}
  if (exists("path_to_Cwood_potential") == FALSE)       {path_to_Cwood_potential <<- " "}
  if (exists("path_to_soilwater") == FALSE)             {path_to_soilwater <<- " "}
  if (exists("path_to_nbe") == FALSE)                   {path_to_nbe <<- " "}
  if (exists("path_to_gpp") == FALSE)                   {path_to_gpp <<- " "}
  if (exists("path_to_et") == FALSE)                    {path_to_et <<- " "}
  if (exists("path_to_fire") == FALSE)                  {path_to_fire <<- " "}
  if (exists("path_to_lca") == FALSE)                   {path_to_lca <<- " "}
  if (exists("path_to_leaflifespan") == FALSE)          {path_to_leaflifespan <<- " "}
  if (exists("path_to_labile_release_period") == FALSE) {path_to_labile_release_period <<- " "}
  if (exists("path_to_labile_release_timing") == FALSE) {path_to_labile_release_timing <<- " "}
  if (exists("path_to_MaxRootDepth") == FALSE)          {path_to_MaxRootDepth <<- " "}  
  if (exists("path_to_frac_Cwood_coarse_root") == FALSE) {path_to_frac_Cwood_coarse_root <<- " "}  
  # Forcings data options
  if (exists("met_source") == FALSE)                    {met_source <<- " "}
  if (exists("burnt_area_source") == FALSE)             {burnt_area_source <<- " "}
  if (exists("deforestation_source") == FALSE)          {deforestation_source <<- " "}
  if (exists("sand_clay_source") == FALSE)              {sand_clay_source <<- " "}
  if (exists("lai_change_source") == FALSE)             {lai_change_source <<- " "}
  # Assimilated data options 
  if (exists("lai_source") == FALSE)                    {lai_source <<- " "}
  if (exists("fapar_source") == FALSE)                  {fapar_source <<- " "}
  if (exists("RhetQ10_source") == FALSE)                {RhetQ10_source <<- " "}
  if (exists("MTTsom_source") == FALSE)                 {MTTsom_source <<- " "}  
  if (exists("Csom_source") == FALSE)                   {Csom_source <<- " "} 
  if (exists("soilwater_source") == FALSE)              {soilwater_source <<- " "}
  if (exists("et_source") == FALSE)                     {et_source <<- " "}
  if (exists("Cwood_change_source") == FALSE)           {Cwood_change_source <<- " "}
  if (exists("Cwood_growth_source") == FALSE)           {Cwood_growth_source <<- " "}
  if (exists("Cwood_mortality_source") == FALSE)        {Cwood_mortality_source <<- " "}
  if (exists("Cagb_change_source") == FALSE)            {Cagb_change_source <<- " "}
  if (exists("Cagb_growth_source") == FALSE)            {Cagbgrowth_source <<- " "}
  if (exists("Cagb_mortality_source") == FALSE)         {Cagb_mortality_source <<- " "}
  if (exists("gpp_source") == FALSE)                    {gpp_source <<- " "}
  if (exists("fire_source") == FALSE)                   {fire_source <<- " "}
  if (exists("Reco_source") == FALSE)                   {Reco_source <<- " "}
  if (exists("nee_source") == FALSE)                    {nee_source <<- " "}
  if (exists("nbe_source") == FALSE)                    {nbe_source <<- " "}
  if (exists("harvest_source") == FALSE)                {harvest_source <<- " "}
  if (exists("foliage_to_litter_source") == FALSE)      {foliage_to_litter_source <<- " "}
  if (exists("Cfol_initial_source") == FALSE)           {Cfol_initial_source <<- " "}
  if (exists("Cwood_initial_source") == FALSE)          {Cwood_initial_source <<- " "}
  if (exists("Croots_initial_source") == FALSE)         {Croots_initial_source <<- " "}
  if (exists("Clit_initial_source") == FALSE)           {Clit_initial_source <<- " "}
  if (exists("Cfol_stock_source") == FALSE)             {Cfol_stock_source <<- " "}
  if (exists("Cfolmax_stock_source") == FALSE)          {Cfolmax_stock_source <<- " "}
  if (exists("Cwood_stock_source") == FALSE)            {Cwood_stock_source <<- " "}
  if (exists("Cstem_stock_source") == FALSE)            {Cstem_stock_source <<- " "}
  if (exists("Cbranch_stock_source") == FALSE)          {Cbranch_stock_source <<- " "}
  if (exists("Cagb_stock_source") == FALSE)             {Cagb_stock_source <<- " "}
  if (exists("Ccoarseroot_stock_source") == FALSE)      {Ccoarseroot_stock_source <<- " "}
  if (exists("Croots_stock_source") == FALSE)           {Croots_stock_source <<- " "}
  if (exists("Clit_stock_source") == FALSE)             {Clit_stock_source <<- " "}
  if (exists("lca_source") == FALSE)                    {lca_source <<- " "}
  if (exists("leaflifespan_source") == FALSE)           {leaflifespan_source <<- " "}
  if (exists("labile_release_timing_source") == FALSE)  {labile_release_timing_source <<- " "}
  if (exists("labile_release_period_source") == FALSE)  {labile_release_period_source <<- " "}
  if (exists("frac_Cwood_coarse_root_source") == FALSE) {frac_Cwood_coarse_root_source <<- " "}
  if (exists("minLWP_source") == FALSE)                 {minLWP_source <<- " "}
  if (exists("MaxRootDepth_source") == FALSE)           {MaxRootDepth_source <<- " "}
  if (exists("Cwood_potential_source") == FALSE)        {Cwood_potential_source <<- " "}
  if (exists("crop_management_source") == FALSE)        {crop_management_source <<- " "}
  if (exists("snow_source") == FALSE)                   {snow_source <<- " "}
  # Algorithm options
  if (exists("request_nos_chains") == FALSE)            {request_nos_chains <<- 3}
  if (exists("request_nos_samples") == FALSE)           {request_nos_samples <<- 10e6}
  if (exists("request_nos_subsamples") == FALSE)        {request_nos_subsamples <<- 1e3}
  if (exists("request_use_EDCs") == FALSE)              {request_use_EDCs <<- TRUE}
  if (exists("request_extended_mcmc") == FALSE)         {request_extended_mcmc <<- 10e6}
  if (exists("request_cost_function_scaling") == FALSE) {request_cost_function_scaling <<- 0}
  if (exists("method") == FALSE)                        {method <<- "APMCMC"} # "APMCMC" or "DEMCz" 
  # Computer running / compiler locations
  if (exists("request_use_server") == FALSE)            {request_use_server <<- FALSE}
  if (exists("request_use_local_slurm") == FALSE)       {request_use_local_slurm <<- FALSE} 
  if (exists("request_runtime") == FALSE)               {request_runtime <<- 48}
  if (exists("request_compile_server") == FALSE)        {request_compile_server <<- FALSE}
  if (exists("request_compile_local") == FALSE)         {request_compile_local <<- TRUE}
  # Computer defaults
  #if (request_use_server & exists("sshpass_key_home") == FALSE) {print("CARDAMOM R code base assumes that access to the remote service is managed using ssh passkey. 
  #                                                                     The sshpass_key_home variable is missing. 
  #                                                                     This please create a passkey for both your home and server machines and specify using 
  #                                                                     sshpass_key_home and sshpass_key_server as appropriate")}
  #if (request_use_server & exists("sshpass_key_server") == FALSE) {print("CARDAMOM R code base assumes that access to the remote service is managed using ssh passkey. 
  #                                                                        The sshpass_key_server variable is missing. 
  #                                                                        This please create a passkey for both your home and server machines and specify using 
  #                                                                        sshpass_key_home and sshpass_key_server as appropriate")}
  if (request_use_server & exists("home_computer") == FALSE) {print("A request to run CARDAMOM on a remote server has been made, but not home_computer has been set.
                                                                     Please specify the address of the home computer used for moving files to and from the remote server")}
  # Compiler options
  if (exists("language") == FALSE) {language <<- "Fortran"} # Assume that the language is Fortran, currently all that actually works
  if (exists("compiler") == FALSE) {compiler <<- "gfortran"} # "gfortran" (GNU) or "ifx" (Intel); passed to CMake as CMAKE_Fortran_COMPILER
  if (exists("compiler_optimisation") == FALSE) {compiler_optimisation <<- "-O2"} # tested options are current models are -O2, -O3, -Ofast
  # Slurm options
  if (exists("slurm_concurrent_cpus") == FALSE) {slurm_concurrent_cpus = 60} # number of concurrent cpus allowed by slurm
  if (exists("slurm_max_run_time") == FALSE) {slurm_max_run_time = 12} # Number of hours which will be asked for task in stage 3 
  if (request_use_local_slurm & exists("slurm_account") == FALSE) {print("The slurm local cluster has been requested for use but no slurm account has been specified (e.g. slurm_account = 'geos_research')")} 
  # Any R library options which need setting
  if (exists("cardamom_temporary_directory") == FALSE) { cardamom_temporary_directory <<- paste(paths$cardamom_outputs,"/temporary_files",sep="") }
  if (dir.exists(cardamom_temporary_directory) == FALSE) {dir.create(cardamom_temporary_directory)}
  ## Check for obvious combination errors
  # Forcings datasets
  if (met_source != "site_specific" & path_to_met_source == " ")                                       {stop(paste("specified 'met_source' and 'path_to_met_source' incompatible"))}
  if (burnt_area_source != "site_specific" & burnt_area_source != " " & path_to_burnt_area == " ")     {stop(paste("specified 'burnt_area_source' and 'path_to_burnt_area' incompatible"))}
  if (deforestation_source != "site_specific" & deforestation_source != " " & path_to_forestry == " ") {stop(paste("specified 'deforestation_source' and 'path_to_forestry' incompatible"))}
  if (sand_clay_source != "site_specific" & sand_clay_source != " " & path_to_sand_clay == " ")        {stop(paste("specified 'sand_clay_source' and 'path_to_sand_clay' incompatible"))}
  if (lai_change_source != "site_specific" & lai_change_source != " " & path_to_lai_change == " ")     {stop(paste("specified 'lai_change_source' and 'path_to_lai_change' incompatible"))}
  # Assimilated datasets, note these only consider those attached to gridded datasets
  if (lai_source != "site_specific" & lai_source != " " & path_to_lai == " ")                                    {stop(paste("specified 'lai_source' and 'path_to_lai' incompatible"))}
  if (fapar_source != "site_specific" & fapar_source != " " & path_to_fapar == " ")                              {stop(paste("specified 'fapar_source' and 'path_to_fapar' incompatible"))}
  if (crop_management_source != "site_specific" & crop_management_source != " " & path_to_crop_management == " "){stop(paste("specified 'crop_management_source' and 'path_to_crop_management' incompatible"))}
  if (Csom_source != "site_specific" & Csom_source != " " & path_to_Csom == " ")                                 {stop(paste("specified 'Csom_source' and 'path_to_Csom' incompatible"))}
  if (MTTsom_source != "site_specific" & MTTsom_source != " " & path_to_MTTsom == " ")                           {stop(paste("specified 'MTTsom_source' and 'path_to_MTTsom' incompatible"))}  
  if (RhetQ10_source != "site_specific" & RhetQ10_source != " " & path_to_RhetQ10 == " ")                        {stop(paste("specified 'RhetQ10_source' and 'path_to_RhetQ10' incompatible"))}  
  if (Cwood_change_source != "site_specific" & Cwood_change_source != " " & path_to_Cwood_change == " ")         {stop(paste("specified 'Cwood_change_source' and 'path_to_Cwood_change' incompatible"))}
  if (Cwood_growth_source != "site_specific" & Cwood_growth_source != " " & path_to_Cwood_growth == " ")         {stop(paste("specified 'Cwood_growth_source' and 'path_to_Cwood_growth' incompatible"))}
  if (Cwood_mortality_source != "site_specific" & Cwood_mortality_source != " " & path_to_Cwood_mortality == " "){stop(paste("specified 'Cwood_mortality_source' and 'path_to_Cwood_mortality' incompatible"))}
  if (Cwood_stock_source != "site_specific" & Cwood_stock_source != " " & path_to_Cwood == " ")                  {stop(paste("specified 'Cwood_stock_source' and 'path_to_Cwood' incompatible"))}
  if (Cagb_change_source != "site_specific" & Cagb_change_source != " " & path_to_Cagb_change == " ")            {stop(paste("specified 'Cagb_change_source' and 'path_to_Cagb_change' incompatible"))}
  if (Cagb_growth_source != "site_specific" & Cagb_growth_source != " " & path_to_Cagb_growth == " ")            {stop(paste("specified 'Cagb_growth_source' and 'path_to_Cagb_growth' incompatible"))}
  if (Cagb_mortality_source != "site_specific" & Cagb_mortality_source != " " & path_to_Cagb_mortality == " ")   {stop(paste("specified 'Cagb_mortality_source' and 'path_to_Cagb_mortality' incompatible"))}
  if (Cagb_stock_source != "site_specific" & Cagb_stock_source != " " & path_to_Cagb == " ")                     {stop(paste("specified 'Cagb_stock_source' and 'path_to_Cagb' incompatible"))}
  if (Cwood_initial_source != "site_specific" & Cwood_initial_source != " " & path_to_Cwood_initial == " ")      {stop(paste("specified 'Cwood_initial_source' and 'path_to_Cwood_initial' incompatible"))}
  if (Cwood_potential_source != "site_specific" & Cwood_potential_source != " " & path_to_Cwood_potential == " "){stop(paste("specified 'Cwood_potential_source' and 'path_to_Cwood_potential' incompatible"))}
  if (soilwater_source != "site_specific" & soilwater_source != " " & path_to_soilwater == " ")                  {stop(paste("specified 'soilwater_source' and 'path_to_soilwater' incompatible"))}
  if (nbe_source != "site_specific" & nbe_source != " " & path_to_nbe == " ")                                    {stop(paste("specified 'nbe_source' and 'path_to_nbe' incompatible"))}
  if (gpp_source != "site_specific" & gpp_source != " " & path_to_gpp == " ")                                    {stop(paste("specified 'gpp_source' and 'path_to_gpp' incompatible"))}
  if (fire_source != "site_specific" & fire_source != " " & path_to_fire == " ")                                 {stop(paste("specified 'fire_source' and 'path_to_fire' incompatible"))}
  if (lca_source != "site_specific" & lca_source != " " & path_to_lca == " ")                                    {stop(paste("specified 'lca_source' and 'path_to_lca' incompatible"))}
  if (leaflifespan_source != "site_specific" & leaflifespan_source != " " & path_to_leaflifespan == " ")         {stop(paste("specified 'leaflifespan_source' and 'path_to_leaflifespan' incompatible"))}
  if (labile_release_timing_source != "site_specific" & labile_release_timing_source != " " & path_to_labile_release_timing == " ") {stop(paste("specified 'labile_release_timing_source' and 'path_to_labile_release_timing' incompatible"))}
  if (labile_release_period_source != "site_specific" & labile_release_period_source != " " & path_to_labile_release_period == " ") {stop(paste("specified 'labile_release_period_source' and 'path_to_labile_release_period' incompatible"))}
  if (MaxRootDepth_source != "site_specific" & MaxRootDepth_source != " " & path_to_MaxRootDepth == " ")         {stop(paste("specified 'MaxRootDepth_source' and 'path_to_MaxRootDepth' incompatible"))}  
  if (frac_Cwood_coarse_root_source != "site_specific" & frac_Cwood_coarse_root_source != " " & path_to_frac_Cwood_coarse_root == " ")         {stop(paste("specified 'frac_Cwood_coarse_root_source' and 'path_to_frac_Cwood_coarse_root' incompatible"))}  

  # NOTE: current assimilated variables which have currently only been done on site scale, without gridded datasets
#  if (exists("Evap_source") == FALSE)                   {Evap_source <<- " "}
#  if (exists("fire_source") == FALSE)                   {fire_source <<- " "}
#  if (exists("Reco_source") == FALSE)                   {Reco_source <<- " "}
#  if (exists("nee_source") == FALSE)                    {nee_source <<- " "}
#  if (exists("harvest_source") == FALSE)                {harvest_source <<- " "}
#  if (exists("foliage_to_litter_source") == FALSE)      {foliage_to_litter_source <<- " "}
#  if (exists("Cfol_initial_source") == FALSE)           {Cfol_initial_source <<- " "}
#  if (exists("Croots_initial_source") == FALSE)         {Croots_initial_source <<- " "}
#  if (exists("Clit_initial_source") == FALSE)           {Clit_initial_source <<- " "}
#  if (exists("Cfol_stock_source") == FALSE)             {Cfol_stock_source <<- " "}
#  if (exists("Cfolmax_stock_source") == FALSE)          {Cfolmax_stock_source <<- " "}
#  if (exists("Cstem_stock_source") == FALSE)            {Cstem_stock_source <<- " "}
#  if (exists("Cbranch_stock_source") == FALSE)          {Cbranch_stock_source <<- " "}
#  if (exists("Cagb_stock_source") == FALSE)             {Cagb_stock_source <<- " "}
#  if (exists("Ccoarseroot_stock_source") == FALSE)      {Ccoarseroot_stock_source <<- " "}
#  if (exists("Croots_stock_source") == FALSE)           {Croots_stock_source <<- " "}
#  if (exists("Clit_stock_source") == FALSE)             {Clit_stock_source <<- " "}
#  if (exists("frac_Cwood_coarse_root_source") == FALSE) {frac_Cwood_coarse_root_source <<- " "}
#  if (exists("minLWP_source") == FALSE)                 {minLWP_source <<- " "}
#  if (exists("crop_management_source") == FALSE)        {crop_management_source <<- " "}
#  if (exists("snow_source") == FALSE)                   {snow_source <<- " "}

  return(print("Standard checks on control file completed"))

} # end function check_control_file_defaults

## Use byte compile
check_control_file_defaults<-cmpfun(check_control_file_defaults)


