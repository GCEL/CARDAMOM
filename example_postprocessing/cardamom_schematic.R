# =============================================================================
# cardamom_schematic.R
#
# Reproduces the CARDAMOM terrestrial ecosystem carbon-budget schematic.
#
# All numeric values (medians, confidence interval bounds) are declared as
# named variables in Section 1 below.  To generate a new schematic from a
# different model run, only Section 1 needs to be updated — the drawing
# code in Section 3 onwards reads exclusively from those variables.
#
# All cex (character expansion / text size) parameters are centralised in
# Section 2 as named scalars.  Adjusting a single value there propagates
# to every label that uses it throughout the figure.
#
# Typical workflow for multiple runs
# ------------------------------------
#   source("cardamom_schematic.R")  # draws and saves the PNG
#
# The script is structured as a callable function `cardamom_schematic()`
# so it can be sourced and called programmatically (see bottom of file).
# =============================================================================

# Function calculates specific values needed for the schematic
prepare_cardamom_values<-function(PROJECT, site_nos) {

   # Create empty output list object
   output = list()

   # Set number decimal places
   dp = 1

   # Do something different depending on whether it is a grid or site
   if (PROJECT$spatial_type == "grid") {
   
       # A gridded analysis
    
       # Specify quantiles we want to extract, should only be 3 (low/median/high).
       # From gridded analysis this must be from those available in the file. 
       # Check grid_output$num_quantiles for available
       quantiles_wanted = c(1,5,8)

       # Load the processed site file
       load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))

       # Calculate the pixel area weightings for the spatial average
       area_weighting = grid_output$land_fraction
       area_weighting[is.na(grid_output$mean_gpp_gCm2day[,,1])] = 0 # filters for missing pixels
       area_weighting = area_weighting*PROJECT$area_m2 # multiple by the pixel areas
       area_weighting = area_weighting / sum(area_weighting) # complete scaling to make this a product scaler for aggregating with sum
       area_weighting = array(area_weighting, dim = c(dim(area_weighting),length(quantiles_wanted)))

       # Extract or calculate required derived values
       # NOTE: unit conversion from gC/m2/day to MgC/ha/yr

       # NATURAL FLUXES
       output$gpp_gCm2yr = format(round(apply(grid_output$mean_gpp_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$rauto_gCm2yr = format(round(apply(grid_output$mean_rauto_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$rhet_litter_gCm2yr = format(round(apply(grid_output$mean_rhet_litter_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$rhet_som_gCm2yr = format(round(apply(grid_output$mean_rhet_som_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$rhet_gCm2yr = format(round(apply(grid_output$mean_rhet_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$npp_gCm2yr = format(round(apply(grid_output$mean_npp_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$npp_labile_gCm2yr = format(round(apply(grid_output$mean_alloc_labile_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$npp_foliage_gCm2yr = format(round(apply(grid_output$mean_alloc_foliage_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$labile_to_foliage_gCm2yr = format(round(apply(grid_output$mean_labile_to_foliage_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$npp_roots_gCm2yr = format(round(apply(grid_output$mean_alloc_roots_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$npp_wood_gCm2yr = format(round(apply(grid_output$mean_alloc_wood_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$foliage_to_litter_gCm2yr = format(round(apply(grid_output$mean_foliage_to_litter_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$roots_to_litter_gCm2yr = format(round(apply(grid_output$mean_roots_to_litter_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$wood_to_litter_gCm2yr = format(round(apply(grid_output$mean_wood_to_litter_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$litter_to_som_gCm2yr = format(round(apply(grid_output$mean_litter_to_som_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       # FIRE FLUXES
       output$FIRElitter_labile_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_labile_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIRElitter_foliage_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_foliage_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIRElitter_roots_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_roots_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIRElitter_wood_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_wood_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIRElitter_litter_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_litter_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIREemiss_labile_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_labile_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIREemiss_foliage_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_foliage_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIREemiss_roots_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_roots_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIREemiss_wood_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_wood_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIREemiss_litter_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_litter_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIREemiss_som_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_som_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$fire_gCm2yr = format(round(apply(grid_output$mean_fire_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       # Net fluxes
       output$nbe_gCm2yr = format(round(apply(grid_output$mean_nbe_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$nee_gCm2yr = format(round(apply(grid_output$mean_nee_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$nbp_gCm2yr = format(round(apply(grid_output$mean_nbp_gCm2day[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       # STOCKS
       output$labile_gCm2 = format(round(apply(grid_output$mean_labile_gCm2[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 1e-2, digits = dp), nsmall = dp)
       output$foliage_gCm2 = format(round(apply(grid_output$mean_foliage_gCm2[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 1e-2, digits = dp), nsmall = dp)
       output$roots_gCm2 = format(round(apply(grid_output$mean_roots_gCm2[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 1e-2, digits = dp), nsmall = dp)
       output$wood_gCm2 = format(round(apply(grid_output$mean_wood_gCm2[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 1e-2, digits = dp), nsmall = dp)
       output$litter_gCm2 = format(round(apply(grid_output$mean_litter_gCm2[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 1e-2, digits = dp), nsmall = dp)
       output$som_gCm2 = format(round(apply(grid_output$mean_som_gCm2[,,quantiles_wanted]*area_weighting,3,sum, na.rm=TRUE) * 1e-2, digits = dp), nsmall = dp)
  
   } else if (PROJECT$spatial_type == "site") {
   
       # A site analysis
        
       # Load the processed site file
       load(paste(PROJECT$results_processedpath,PROJECT$sites[site_nos],".RData",sep=""))

       # Specify quantiles we want to extract, should only be 3 (low/median/high)
       quantiles_wanted = c(0.025,0.50,0.975)    
       #quantiles_wanted = c(0.05,0.50,0.95)    
       # Restrict the time bounds?
       s = 1 ; f = dim(states_all$lai_m2m2)[2]
#       s = f-(17*12) ; f = dim(states_all$lai_m2m2)[2]
    
       # Extract or calculate required derived values
       # NOTE: unit conversion from gC/m2/day to MgC/ha/yr

       # NATURAL FLUXES
       output$gpp_gCm2yr = format(round(quantile(apply(states_all$gpp_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$rauto_gCm2yr = format(round(quantile(apply(states_all$rauto_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$rhet_litter_gCm2yr = format(round(quantile(apply(states_all$rhet_litter_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$rhet_som_gCm2yr = format(round(quantile(apply(states_all$rhet_som_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$rhet_gCm2yr = format(round(quantile(apply(states_all$rhet_litter_gCm2day[,s:f]+states_all$rhet_som_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$npp_gCm2yr = format(round(quantile(apply(states_all$gpp_gCm2day[,s:f]-states_all$rauto_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$npp_labile_gCm2yr = format(round(quantile(apply(states_all$alloc_labile_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$npp_foliage_gCm2yr = format(round(quantile(apply(states_all$alloc_foliage_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$labile_to_foliage_gCm2yr = format(round(quantile(apply(states_all$labile_to_foliage_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$npp_roots_gCm2yr = format(round(quantile(apply(states_all$alloc_roots_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$npp_wood_gCm2yr = format(round(quantile(apply(states_all$alloc_wood_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$foliage_to_litter_gCm2yr = format(round(quantile(apply(states_all$foliage_to_litter_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$roots_to_litter_gCm2yr = format(round(quantile(apply(states_all$roots_to_litter_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$wood_to_litter_gCm2yr = format(round(quantile(apply(states_all$wood_to_litter_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$litter_to_som_gCm2yr = format(round(quantile(apply(states_all$litter_to_som_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       # FIRE FLUXES
       output$FIRElitter_labile_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_labile_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIRElitter_foliage_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_foliage_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIRElitter_roots_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_roots_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIRElitter_wood_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_wood_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIRElitter_litter_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_litter_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIREemiss_labile_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_labile_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIREemiss_foliage_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_foliage_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIREemiss_roots_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_roots_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIREemiss_wood_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_wood_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIREemiss_litter_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_litter_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$FIREemiss_som_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_som_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$fire_gCm2yr = format(round(quantile(apply(states_all$fire_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       # Net fluxes
       output$nbe_gCm2yr = format(round(quantile(apply((states_all$fire_gCm2day+states_all$rhet_litter_gCm2day+states_all$rhet_som_gCm2day+states_all$rauto_gCm2day)[,s:f]-states_all$gpp_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$nee_gCm2yr = format(round(quantile(apply((states_all$rhet_litter_gCm2day+states_all$rhet_som_gCm2day+states_all$rauto_gCm2day)[,s:f]-states_all$gpp_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       output$nbp_gCm2yr = format(round(quantile(apply((states_all$nbp_gCm2day)[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
       # STOCKS
       output$labile_gCm2 = format(round(quantile(apply(states_all$labile_gCm2[,s:f],1,mean) * 1e-2, prob=quantiles_wanted), digits = dp), nsmall = dp)
       output$foliage_gCm2 = format(round(quantile(apply(states_all$foliage_gCm2[,s:f],1,mean) * 1e-2, prob=quantiles_wanted), digits = dp), nsmall = dp)
       output$roots_gCm2 = format(round(quantile(apply(states_all$roots_gCm2[,s:f],1,mean), prob=quantiles_wanted) * 1e-2, digits = dp), nsmall = dp)
       output$wood_gCm2 = format(round(quantile(apply(states_all$wood_gCm2[,s:f],1,mean), prob=quantiles_wanted) * 1e-2, digits = dp), nsmall = dp)
       output$litter_gCm2 = format(round(quantile(apply(states_all$litter_gCm2[,s:f],1,mean), prob=quantiles_wanted) * 1e-2, digits = dp), nsmall = dp)
       output$som_gCm2 = format(round(quantile(apply(states_all$som_gCm2[,s:f],1,mean), prob=quantiles_wanted) * 1e-2, digits = dp), nsmall = dp)
    
   } else {
       # We have a compatibility problem
       stop("PROJECT$spatial_type does not have a compatible value, i.e. grid or site")
   } # end if is_grid

   # Return 
   return(output)

} # end function prepare_cardamom_schematic

# A function to draw a schematic of CARDAMOM, but in R language
cardamom_schematic <- function(list_in, out_file = "cardamom_schematic.png") {

   # ===========================================================================
   # SECTION 1 — All numeric values (medians and 95 % CI bounds)
   #
   # Each variable name follows the pattern:
   #   <component>_med   : median estimate
   #   <component>_lo    : 2.5th percentile (lower CI bound)
   #   <component>_hi    : 97.5th percentile (upper CI bound)
   #
   # Units: MgC ha-1 yr-1 for fluxes; MgC ha-1 for stocks (pools).
   # Replace these default values with outputs read from your model files.
   # ===========================================================================

   # ── Gross primary production ─────────────────────────────────────────────────
   GPP_med = list_in$gpp_gCm2yr[2]
   GPP_lo  = list_in$gpp_gCm2yr[1]
   GPP_hi  = list_in$gpp_gCm2yr[3]

   # ── Autotrophic respiration ───────────────────────────────────────────────────
   Ra_med  = list_in$rauto_gCm2yr[2]
   Ra_lo   = list_in$rauto_gCm2yr[1]
   Ra_hi   = list_in$rauto_gCm2yr[3]

   # ── Net primary production (total) ───────────────────────────────────────────
   NPP_med = list_in$npp_gCm2yr[2]
   NPP_lo  = list_in$npp_gCm2yr[1]
   NPP_hi  = list_in$npp_gCm2yr[3]

   # ── NPP allocation fluxes ─────────────────────────────────────────────────
   NPP_lab_med  = list_in$npp_labile_gCm2yr[2] ; NPP_lab_lo  = list_in$npp_labile_gCm2yr[1] ; NPP_lab_hi  = list_in$npp_labile_gCm2yr[3]
   NPP_fol_med  = list_in$npp_foliage_gCm2yr[2]; NPP_fol_lo  = list_in$npp_foliage_gCm2yr[1]; NPP_fol_hi  = list_in$npp_foliage_gCm2yr[3]
   NPP_root_med = list_in$npp_roots_gCm2yr[2]  ; NPP_root_lo = list_in$npp_roots_gCm2yr[1]  ; NPP_root_hi = list_in$npp_roots_gCm2yr[3]
   NPP_wood_med = list_in$npp_wood_gCm2yr[2]   ; NPP_wood_lo = list_in$npp_wood_gCm2yr[1]   ; NPP_wood_hi = list_in$npp_wood_gCm2yr[3]

   # ── Labile-to-foliage transfer ───────────────────────────────────────────────
   lab2fol_med  = list_in$labile_to_foliage_gCm2yr[2] ; lab2fol_lo  = list_in$labile_to_foliage_gCm2yr[1] ; lab2fol_hi  = list_in$labile_to_foliage_gCm2yr[3]
   
   # ── Carbon pools (stocks) ─────────────────────────────────────────────────────
   Clab_med    = list_in$labile_gCm2[2]  ; Clab_lo    = list_in$labile_gCm2[1]  ; Clab_hi    = list_in$labile_gCm2[3]
   Cfol_med    = list_in$foliage_gCm2[2] ; Cfol_lo    = list_in$foliage_gCm2[1] ; Cfol_hi    = list_in$foliage_gCm2[3]
   Croot_med   = list_in$roots_gCm2[2]   ; Croot_lo   = list_in$roots_gCm2[1]   ; Croot_hi   = list_in$roots_gCm2[3]
   Cwood_med   = list_in$wood_gCm2[2]    ; Cwood_lo   = list_in$wood_gCm2[1]    ; Cwood_hi   = list_in$wood_gCm2[3]
   Clitter_med = list_in$litter_gCm2[2]  ; Clitter_lo = list_in$litter_gCm2[1]  ; Clitter_hi = list_in$litter_gCm2[3]
   Csom_med    = list_in$som_gCm2[2]     ; Csom_lo    = list_in$som_gCm2[1]     ; Csom_hi    = list_in$som_gCm2[3]

   # ── Mortality fluxes — biogenic (black) ──────────────────────────────────────
   #Mort_lab_med  = 0.00 ; Mort_lab_lo   = 0.00 ; Mort_lab_hi   = 0.00
   Mort_fol_med  = list_in$foliage_to_litter_gCm2yr[2] ; Mort_fol_lo  = list_in$foliage_to_litter_gCm2yr[1] ; Mort_fol_hi  = list_in$foliage_to_litter_gCm2yr[3]
   Mort_root_med = list_in$roots_to_litter_gCm2yr[2]   ; Mort_root_lo = list_in$roots_to_litter_gCm2yr[1]   ; Mort_root_hi = list_in$roots_to_litter_gCm2yr[3]
   Mort_wood_med = list_in$wood_to_litter_gCm2yr[2]    ; Mort_wood_lo = list_in$wood_to_litter_gCm2yr[1]    ; Mort_wood_hi = list_in$wood_to_litter_gCm2yr[3]

   # ── Mortality fluxes — fire-driven (red) ─────────────────────────────────────
   Mort_lab_fire_med  = list_in$FIRElitter_labile_gCm2yr[2]  ; Mort_lab_fire_lo  = list_in$FIRElitter_labile_gCm2yr[1]  ; Mort_lab_fire_hi  = list_in$FIRElitter_labile_gCm2yr[3]
   Mort_fol_fire_med  = list_in$FIRElitter_foliage_gCm2yr[2] ; Mort_fol_fire_lo  = list_in$FIRElitter_foliage_gCm2yr[1] ; Mort_fol_fire_hi  = list_in$FIRElitter_foliage_gCm2yr[3]
   Mort_root_fire_med = list_in$FIRElitter_roots_gCm2yr[2]   ; Mort_root_fire_lo = list_in$FIRElitter_roots_gCm2yr[1]   ; Mort_root_fire_hi = list_in$FIRElitter_roots_gCm2yr[3]
   Mort_wood_fire_med = list_in$FIRElitter_wood_gCm2yr[2]   ; Mort_wood_fire_lo = list_in$FIRElitter_wood_gCm2yr[1]    ; Mort_wood_fire_hi = list_in$FIRElitter_wood_gCm2yr[3]

   # ── Decomposition: litter → SOM ──────────────────────────────────────────────
   decomp_med      = list_in$litter_to_som_gCm2yr[2]     ; decomp_lo      = list_in$litter_to_som_gCm2yr[1]     ; decomp_hi      = list_in$litter_to_som_gCm2yr[3]
   decomp_fire_med = list_in$FIRElitter_litter_gCm2yr[2] ; decomp_fire_lo = list_in$FIRElitter_litter_gCm2yr[1] ; decomp_fire_hi = list_in$FIRElitter_litter_gCm2yr[3]

   # ── Heterotrophic respiration ─────────────────────────────────────────────────
   Rh_litter_med = list_in$rhet_litter_gCm2yr[2] ; Rh_litter_lo = list_in$rhet_litter_gCm2yr[1] ; Rh_litter_hi = list_in$rhet_litter_gCm2yr[3]
   Rh_som_med    = list_in$rhet_som_gCm2yr[2]    ; Rh_som_lo    = list_in$rhet_som_gCm2yr[1]    ; Rh_som_hi    = list_in$rhet_som_gCm2yr[3]
   Rh_med        = list_in$rhet_gCm2yr[2]        ; Rh_lo        = list_in$rhet_gCm2yr[1]        ; Rh_hi        = list_in$rhet_gCm2yr[3]

   # ── Fire emissions — per pool ─────────────────────────────────────────────────
   E_lab_med    = list_in$FIREemiss_labile_gCm2yr[2] ; E_lab_lo    = list_in$FIREemiss_labile_gCm2yr[1] ; E_lab_hi    = list_in$FIREemiss_labile_gCm2yr[3]
   E_fol_med    = list_in$FIREemiss_foliage_gCm2yr[2]; E_fol_lo    = list_in$FIREemiss_foliage_gCm2yr[1]; E_fol_hi    = list_in$FIREemiss_foliage_gCm2yr[3]
   E_root_med   = list_in$FIREemiss_roots_gCm2yr[2]  ; E_root_lo   = list_in$FIREemiss_roots_gCm2yr[1]  ; E_root_hi   = list_in$FIREemiss_roots_gCm2yr[3]
   E_wood_med   = list_in$FIREemiss_wood_gCm2yr[2]   ; E_wood_lo   = list_in$FIREemiss_wood_gCm2yr[1]   ; E_wood_hi   = list_in$FIREemiss_wood_gCm2yr[3]
   E_litter_med = list_in$FIREemiss_litter_gCm2yr[2] ; E_litter_lo = list_in$FIREemiss_litter_gCm2yr[1] ; E_litter_hi = list_in$FIREemiss_litter_gCm2yr[3]
   E_som_med    = list_in$FIREemiss_som_gCm2yr[2]    ; E_som_lo    = list_in$FIREemiss_som_gCm2yr[1]    ; E_som_hi    = list_in$FIREemiss_som_gCm2yr[3]
   E_total_med  = list_in$fire_gCm2yr[2]             ; E_total_lo  = list_in$fire_gCm2yr[1]             ; E_total_hi  = list_in$fire_gCm2yr[3]

   # ── Net ecosystem / biome productivity ───────────────────────────────────────
   NEE_med = list_in$nee_gCm2yr[2] ; NEE_lo = list_in$nee_gCm2yr[1] ; NEE_hi = list_in$nee_gCm2yr[3]
   NBE_med = list_in$nbe_gCm2yr[2] ; NBE_lo = list_in$nbe_gCm2yr[1] ; NBE_hi = list_in$nbe_gCm2yr[3]
   NBP_med = list_in$nbp_gCm2yr[2] ; NBP_lo = list_in$nbp_gCm2yr[1] ; NBP_hi = list_in$nbp_gCm2yr[3]   

   # ── Caption period string (update when running a different time window) ───────
   #caption_period = "2003\u20132024"

   # ===========================================================================
   # SECTION 2 — Typography / size parameters (cex = character expansion)
   #
   # Changing any value here propagates to EVERY label that uses that role.
   # Increase all values proportionally to scale up text for a larger figure.
   # ===========================================================================

   cex_title   = 0.72+0.6   # section headings, pool/flux names, NBP/NEE labels
   cex_median  = 0.62+0.3   # median value labels alongside arrows and in boxes
   cex_ci      = 0.50+0.3   # confidence-interval bound labels  (lo / hi)
   cex_fire    = 0.52+0.3   # fire-emission median labels (slightly smaller)
   cex_fire_ci = 0.44+0.3   # fire-emission CI labels
   cex_caption = 0.44+0.3   # bottom caption text
   cex_alloc   = 0.52+0.3   # rotated NPP-allocation box label

   # ===========================================================================
   # SECTION 3 — Internal formatting helpers
   #
   # These functions are defined inside the main function so they close over
   # the cex_* parameters automatically — no need to pass them explicitly.
   # ===========================================================================

   # fire / emission colour (matches LaTeX \color{red})
   red <- "#FF8C00"  # "#CC0000" = red, #FF8C00 = dark orange

   # Format a single number to 2 decimal places for display.
   fmt <- function(x) { formatC(x, format = "f", digits = 2)}

   # Format a CI pair as "lo / hi".
   fmt_ci <- function(lo, hi) { paste0("(",fmt(lo), " / ", fmt(hi),")") }

   # Draw a solid rectangle.
   draw_box <- function(x0, y0, w, h, border = "black", lwd = 1.2, lty = 1, bg = "white") {
      rect(x0, y0, x0 + w, y0 + h, col = bg, border = border, lwd = lwd, lty = lty)
   } # end function draw_box

   # Draw a dashed rectangle (LaTeX \dashbox).
   dash_box <- function(x0, y0, w, h, border = "black", lwd = 1.0) {
     rect(x0, y0, x0 + w, y0 + h, col = NA, border = border, lwd = lwd, lty = 2)
   }

   # Arrowhead line.
   arrow <- function(x0, y0, x1, y1, col = "black", lwd = 1.4, len = 0.10) {
      arrows(x0, y0, x1, y1, length = len, angle = 20, code = 2, col = col, lwd = lwd)
   }

   # Plain line segment.
   seg <- function(x0, y0, x1, y1, col = "black", lwd = 1.2, lty = 1) {
      segments(x0, y0, x1, y1, col = col, lwd = lwd, lty = lty)
   }

   # Text label (bottom-left anchor matches LaTeX \put default).
   lbl <- function(x, y, txt, col = "black", cex = cex_title, adj = c(0, 0), font = 1, srt = 0) {
      text(x, y, txt, col = col, cex = cex, adj = adj, font = font, srt = srt)
   }

   # Composite: flux label = median on one row, CI on next row.
   # Uses cex_median and cex_ci by default.
   flux_lbl <- function(x_med, y_med, med, x_ci, y_ci, lo, hi, col = "black") {
     lbl(x_med, y_med, fmt(med), col = col, cex = cex_median)
     lbl(x_ci, y_ci, fmt_ci(lo, hi), col = col, cex = cex_ci)
   }

   # Fire-emission flux label (uses cex_fire / cex_fire_ci).
   fire_lbl <- function(x_med, y_med, med, x_ci,  y_ci,  lo, hi) {
     lbl(x_med, y_med, fmt(med),       col = red, cex = cex_fire)
     lbl(x_ci,  y_ci,  fmt_ci(lo, hi), col = red, cex = cex_fire_ci)
   }

   # Draw a pool box with name, median and CI pre-positioned inside.
   pool_box <- function(x0, y0, w = 2.5, h = 2.05, pool_expr, med, lo, hi, bg = "white") {
     draw_box(x0, y0, w, h, bg = bg, lwd = 1.4)
     lbl(x0 + 0.45, y0 + h - 0.60, pool_expr, cex = cex_title, font = 2)
     lbl(x0 + 0.55, y0 + h - 1.20, fmt(med), cex = cex_median)
     lbl(x0 + 0.55, y0 + h - 1.80, fmt_ci(lo, hi), cex = cex_ci)
   }

   # ===========================================================================
   # SECTION 4 — Open the graphics device and initialise the canvas
   # ===========================================================================

   png(out_file, width = 2800, height = 1400, res = 200)
   on.exit(dev.off())   # guarantee device is closed even if an error occurs

   # Prepare the plotting area
   par(mar = c(0.2, 0.2, 0.2, 0.2), bg = "white")
   plot.new()
   plot.window(xlim = c(0, 28), ylim = c(0, 14), asp = 1)

   # Outer figure border (LaTeX \fbox).
   rect(0.05, 0.05, 28.20, 13.95, col = NA, border = "black", lwd = 1.5)

   # ===========================================================================
   # SECTION 5 — Static structural elements (labels, boxes, routing lines)
   #             that carry no numeric data values.
   # ===========================================================================

   # ── Input / Output section headings ──────────────────────────────────────────
   lbl(0.1, 13.5, "Input",   font = 3, cex = cex_title * 0.95)
   lbl(0.1, 13.0, "carbon",  font = 3, cex = cex_title * 0.95)
   lbl(0.1, 12.5, "rates",   font = 3, cex = cex_title * 0.95)
   lbl(26.7, 13.5, "Output", font = 3, cex = cex_title * 0.95)
   lbl(26.7, 13.0, "carbon", font = 3, cex = cex_title * 0.95)
   lbl(26.7, 12.5, "rates",  font = 3, cex = cex_title * 0.95)

   # ── Internal carbon rates dashed boundary box ─────────────────────────────────
   dash_box(2.95, 0.15, 20.55, 11.2, lwd = 1.2)
   lbl(3.0, 10.70, "Internal carbon rates", font = 3, cex = cex_title)

   # ── CUE dashed box + routing lines ───────────────────────────────────────────
   dash_box(3.20, 3.70, 1.75, 1.65)
   lbl(3.57, 4.35, "CUE", font = 2, cex = cex_title * 0.95)
   seg(4.075, 3.70, 4.075, 0.65)       # vertical drop from CUE box
   arrow(4.075, 0.65, 26.0, 0.65)      # horizontal run to Ra

   # ── NPP allocation dashed box ─────────────────────────────────────────────────
   dash_box(7.95, 1.0, 1.2, 7.0, lwd = 1.0)
   lbl(8.55, 4.50, "NPP allocation", font = 3, cex = cex_alloc,
       srt = 90, adj = c(0.5, 0.5))  

   # ── Mortality routing lines (biogenic — go to litter or SOM) ─────────────────
   # Labile → litter
   seg(15.25, 7.0, 19.70, 7.0)
   arrow(19.70, 7.0, 20.30, 6.40)
   # Foliage → litter
   seg(15.25, 9.75, 19.70, 9.75)
   arrow(19.70, 9.75, 20.95, 7.20)
   # Fine root → litter
   seg(15.25, 4.50, 19.70, 4.50)
   arrow(19.70, 4.50, 20.30, 5.20)
   # Wood → SOM
   seg(15.25, 2.0, 19.70, 2.0)
   arrow(19.70, 2.0, 20.30, 1.80)

   # ── Decomposition arrow (litter → SOM) ───────────────────────────────────────
   arrow(21.50, 5.20, 21.50, 3.75)

   # ── Heterotrophic respiration routing lines ───────────────────────────────────
   # Rhet - litter
   seg(22.75, 6.0,  24.85, 6.0)
   arrow(24.85, 6.0, 25.70, 4.60)
   # Rhet - som
   seg(22.75, 2.75, 24.85, 2.75)
   arrow(24.85, 2.75, 25.85, 3.60)

   # ── Fire emission vertical lines (red) ───────────────────────────────────────
   # E_lab: horizontal tap from labile pool side
   seg(12.5, 7.5, 12.7, 7.5, col = red, lwd = 1.0)
   seg(12.5, 7.5, 12.5, 13.5, col = red, lwd = 1.0)
   # E_fol: from foliage pool top
   seg(14.1, 10.80, 14.1, 13.5, col = red, lwd = 1.0)
   # E_root: horizontal tap then vertical
   seg(11.85, 5.0,  12.70, 5.0,  col = red, lwd = 1.0)
   seg(11.85, 5.0,  11.85, 13.5, col = red, lwd = 1.0)
   # E_wood: horizontal tap then vertical
   seg(15.20, 2.70, 15.65, 2.70, col = red, lwd = 1.0)
   seg(15.65, 2.70, 15.65, 13.5, col = red, lwd = 1.0)
   # E_litter: from litter pool top
   seg(20.95, 7.20, 20.95, 13.5, col = red, lwd = 1.0)
   # E_som: horizontal tap then vertical
   seg(22.75, 3.40, 23.15, 3.40, col = red, lwd = 1.0)
   seg(23.15, 3.40, 23.15, 13.5, col = red, lwd = 1.0)
   # E_total: horizontal arrow collecting all emissions
   arrow(11.85, 13.5, 24.10, 13.5, col = red, lwd = 1.6, len = 0.12)

   # ── NBP and NEE dashed summary boxes ─────────────────────────────────────────
   dash_box(2.95, 11.5, 2.2, 2.2, lwd = 1.0)
   dash_box(5.95, 11.5, 2.2, 2.2, lwd = 1.0)

   # ===========================================================================
   # SECTION 6 — All text labels that carry numeric data values
   #             Every hard-coded number has been replaced by a variable.
   # ===========================================================================

   # ── GPP ──────────────────────────────────────────────────────────────────────
   arrow(0.05, 4.50, 3.25, 4.50)
   lbl(0.76, 4.70, "GPP", font = 2, cex = cex_title)
   flux_lbl(0.15, 4.0, GPP_med, 0.15+0.9, 4.0+0.024, GPP_lo, GPP_hi) # Suggested gap in x = 0.9, y = 0.024

   # ── Ra ───────────────────────────────────────────────────────────────────────
   lbl(26.0, 0.78, expression(R[a]), font = 2, cex = cex_title)
   flux_lbl(25.2, 0.20, Ra_med, 25.2+0.9, 0.20+0.024, Ra_lo, Ra_hi) 

   # ── NPP ──────────────────────────────────────────────────────────────────────
   arrow(4.95, 4.50, 7.95, 4.50)
   lbl(5.70, 4.70, "NPP", font = 2, cex = cex_title)
   flux_lbl(5.05, 4.0,  NPP_med, 5.05+0.9, 4.0+0.024, NPP_lo, NPP_hi)

   # ── NPP to labile ─────────────────────────────────────────────────────────────
   arrow(9.15, 6.70, 12.75, 6.70)
   lbl(9.35, 6.90, expression(NPP[lab]), cex = cex_median)
   flux_lbl(9.3, 6.2, NPP_lab_med, 9.3+0.9, 6.2+0.024, NPP_lab_lo, NPP_lab_hi)

   # ── NPP to foliage (angled ~30°) ─────────────────────────────────────────────
   arrow(9.15, 7.8, 12.70, 10.0, lwd = 1.4)
   lbl(9.25, 8.0, expression(NPP[fol]), cex = cex_median, srt = 31)
   lbl(9.40, 7.5, fmt(NPP_fol_med),     cex = cex_median, srt = 31)
   lbl(9.40+0.9, 7.5+0.574, fmt_ci(NPP_fol_lo, NPP_fol_hi), cex = cex_ci, srt = 31)

   # ── NPP to fine root ──────────────────────────────────────────────────────────
   arrow(9.15, 4.50, 12.70, 4.50)
   lbl(9.35, 4.70, expression(NPP[root]), cex = cex_median)
   flux_lbl(9.3, 4.0, NPP_root_med, 9.3+0.9, 4.0+0.024, NPP_root_lo, NPP_root_hi)

   # ── NPP to wood ───────────────────────────────────────────────────────────────
   arrow(9.15, 2.0, 12.70, 2.0)
   lbl(9.35, 2.20, expression(NPP[wood]), cex = cex_median)
   flux_lbl(9.3, 1.5, NPP_wood_med, 9.3+0.9, 1.5+0.024, NPP_wood_lo, NPP_wood_hi)

   # ── Labile → foliage transfer ─────────────────────────────────────────────────
   arrow(13.8, 8.05, 13.8, 8.80)
   flux_lbl(13.00, 8.25, lab2fol_med, 13.0+0.9, 8.25+0.024, lab2fol_lo, lab2fol_hi)

   # ── Pool boxes ────────────────────────────────────────────────────────────────
   # Colours: #EFF6FF = blue, #ECFDF5 = green, #FFF7ED = orange, #FDF4FF = purple, #FFFBEB = yellow
   #          #CAFF70 = light olive green, #A2CD5A = darker olive green 
   pool_box(12.70, 6.0, pool_expr = expression(C[lab]), med = Clab_med, lo = Clab_lo, hi = Clab_hi, bg = "#CAFF70")
   pool_box(12.70, 8.75,pool_expr = expression(C[fol]), med = Cfol_med, lo = Cfol_lo, hi = Cfol_hi, bg = "#CAFF70")
   pool_box(12.70, 3.5, pool_expr = expression(C[root]), med = Croot_med, lo = Croot_lo, hi = Croot_hi, bg = "#CAFF70")
   pool_box(12.70, 1.0, pool_expr = expression(C[wood]), med = Cwood_med, lo = Cwood_lo, hi = Cwood_hi, bg = "#CAFF70")
   pool_box(20.25, 5.2, w = 2.5, h = 2.0, pool_expr = expression(C[litter]), med = Clitter_med, lo = Clitter_lo, hi = Clitter_hi, bg = "#A2CD5A")

   # SOM box is taller (h = 2.75) so draw manually to keep correct offsets.
   draw_box(20.25, 1.0, 2.5, 2.75, bg = "#A2CD5A", lwd = 1.4)
   lbl(20.55, 2.90, expression(C[som]), font = 2, cex = cex_title)
   lbl(20.55, 2.10, fmt(Csom_med), cex = cex_median)
   lbl(20.55, 1.45, fmt_ci(Csom_lo, Csom_hi), cex = cex_ci)

   # ── Mortality: labile ─────────────────────────────────────────────────────────
   lbl(15.70, 7.20, expression(Mort[lab]), cex = cex_median)
   # biogenic natural (near-zero for labile; still show for completeness)
   #flux_lbl(15.70, 6.50, Mort_lab_med, 15.70+0.9, 6.50+0.024, Mort_lab_lo, Mort_lab_hi)
   # fire
   fire_lbl(17.3, 7.2, Mort_lab_fire_med, 17.3+0.9, 7.20+0.024, Mort_lab_fire_lo, Mort_lab_fire_hi)

   # ── Mortality: foliage ────────────────────────────────────────────────────────
   lbl(15.70, 9.95, expression(Mort[fol]), cex = cex_median)
   flux_lbl(15.7, 9.25, Mort_fol_med, 15.7+0.9, 9.25+0.024, Mort_fol_lo, Mort_fol_hi)
   fire_lbl(17.3, 9.95, Mort_fol_fire_med, 17.3+0.9, 9.95+0.024, Mort_fol_fire_lo, Mort_fol_fire_hi)

   # ── Mortality: fine root ──────────────────────────────────────────────────────
   lbl(15.70, 4.65, expression(Mort[root]), cex = cex_median)
   flux_lbl(15.70, 3.95, Mort_root_med, 15.7+0.9, 3.95+0.024, Mort_root_lo, Mort_root_hi)
   fire_lbl(17.3, 4.65, Mort_root_fire_med, 17.3+0.90, 4.65+0.024, Mort_root_fire_lo, Mort_root_fire_hi)

   # ── Mortality: wood ───────────────────────────────────────────────────────────
   lbl(15.70, 2.20, expression(Mort[wood]), cex = cex_median)
   flux_lbl(15.70, 1.5, Mort_wood_med, 15.70+0.9, 1.5+0.024, Mort_wood_lo, Mort_wood_hi)
   fire_lbl(17.3, 2.2, Mort_wood_fire_med, 17.3+0.9, 2.2+0.024, Mort_wood_fire_lo, Mort_wood_fire_hi)

   # ── Decomposition: litter → SOM ──────────────────────────────────────────────
   # biogenic
   flux_lbl(21.55, 4.65, decomp_med, 21.55, 4.10, decomp_lo, decomp_hi)
   # fire-driven combusted litter to SOM
   fire_lbl(20.60, 4.65, decomp_fire_med, 19.95, 4.10, decomp_fire_lo, decomp_fire_hi)

   # ── Heterotrophic respiration (Rlit, Rsom, Rh)────────────────────────────────
   lbl(23.65, 6.20, expression(R["h-litter"]), cex = cex_title)
   flux_lbl(24.1, 5.55, Rh_litter_med, 24.1+1.1, 5.55+0.024, Rh_litter_lo, Rh_litter_hi)
   lbl(23.65, 2.9, expression(R["h-som"]), cex = cex_title)
   flux_lbl(24.1, 2.30, Rh_som_med, 24.1+1.1, 2.30+0.024, Rh_som_lo, Rh_som_hi)
   lbl(26.0, 4.60, expression(R[h]), font = 2, cex = cex_title)
   flux_lbl(25.2, 4.00, Rh_med, 25.2+0.9, 4.00+0.024, Rh_lo, Rh_hi)

   # ── Fire emission labels ──────────────────────────────────────────────────────
   # E_lab
   lbl(12.55, 12.70, expression(E[lab]),  col = red, cex = cex_fire, srt = 0)
   lbl(12.55, 12.40, fmt(E_lab_med),     col = red, cex = cex_fire)
   lbl(12.55, 12.10, fmt_ci(E_lab_lo, E_lab_hi), col = red, cex = cex_fire_ci)
   # E_fol
   lbl(14.15, 12.70, expression(E[fol]), col = red, cex = cex_fire, srt = 0)
   lbl(14.15, 12.40, fmt(E_fol_med),     col = red, cex = cex_fire)
   lbl(14.15, 12.10, fmt_ci(E_fol_lo, E_fol_hi), col = red, cex = cex_fire_ci)
   # E_root
   lbl(11.15, 12.70, expression(E[root]),col = red, cex = cex_fire, srt = 0)
   lbl(11.15, 12.40, fmt(E_root_med),    col = red, cex = cex_fire)
   lbl(10.30, 12.10, fmt_ci(E_root_lo, E_root_hi), col = red, cex = cex_fire_ci)
   # E_wood
   lbl(15.70, 12.70, expression(E[wood]),col = red, cex = cex_fire, srt = 0)
   lbl(15.70, 12.40, fmt(E_wood_med),    col = red, cex = cex_fire)
   lbl(15.70, 12.10, fmt_ci(E_wood_lo, E_wood_hi), col = red, cex = cex_fire_ci)
   # E_litter
   lbl(21.05, 12.70, expression(E[litter]), col = red, cex = cex_fire, srt = 0)
   lbl(21.05, 12.40, fmt(E_litter_med),  col = red, cex = cex_fire)
   lbl(21.05, 12.10, fmt_ci(E_litter_lo, E_litter_hi), col = red, cex = cex_fire_ci)
   # E_som
   lbl(23.20, 12.70, expression(E[som]), col = red, cex = cex_fire, srt = 0)
   lbl(23.20, 12.40, fmt(E_som_med),     col = red, cex = cex_fire)
   lbl(23.20, 12.10, fmt_ci(E_som_lo, E_som_hi), col = red, cex = cex_fire_ci)
   # E_total
   lbl(24.15, 13.35, expression(E[total]), col = red, font = 2, cex = cex_title)
   lbl(24.00, 13.00, fmt(E_total_med),   col = red, cex = cex_median)
   lbl(24.90, 13.00+0.024, fmt_ci(E_total_lo, E_total_hi), col = red, cex = cex_ci)

   # ── NBP ───────────────────────────────────────────────────────────────────────
   lbl(3.50, 13.05, "NBP", font = 2, cex = cex_title)
   flux_lbl(3.65, 12.40, NBP_med, 3.25, 11.80, NBP_lo, NBP_hi)

   # ── NEE ───────────────────────────────────────────────────────────────────────
   lbl(6.50, 13.05, "NEE", font = 2, cex = cex_title)
   flux_lbl(6.65, 12.40, NEE_med, 6.25, 11.80, NEE_lo, NEE_hi)

   # ── Caption ───────────────────────────────────────────────────────────────────
#   caption <- paste0("CARDAMOM global terrestrial C-budget (", caption_period, "). ",
#                     "Stocks in MgC ha\u207B\u00B9; fluxes in MgC ha\u207B\u00B9 y\u207B\u00B9. ",
#                     "Black = biogenic; red = fire-driven." )
#   lbl(14.0, 0.12, caption, cex = cex_caption, adj = c(0.5, 0), col = "grey30")

   # on.exit(dev.off()) handles device closure
   invisible(out_file)

} # end function cardamom_schematic

# =============================================================================

# Load project file
load("~/gcel_ceph/cardamom_analyses/lsmallma/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.004_MHMCMC/global_0.5deg_dalec4_trendyv15_LCA_TWB_GPP_fAPAR_hashimoto_SGDB/infofile.RData")

# If a site run with multiple sites, which site number to use?
# Load from the output object and prepare the wanted values
output = prepare_cardamom_values(PROJECT, site_nos = 1)

# Now generate the schematic
cardamom_schematic(output, out_file = "cardamom_schematic.png")
cat("Saved: cardamom_schematic.png\n")

