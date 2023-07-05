
###
## Script for the creation of a latex script to generate C-budgets for DALEC
## Author: T L Smallman (t.l.smallman@ed.ac.uk)
## Created: 12/07/2022
## Last updated: 14/07/2022
## Contributing authors: 
## NOTES: 
## 1) This version presents a distinct labile and foliage pools. 
## 2) Assumes a CDEA style phenology with both direct GPP to foliage and via labile pathways
###

## TO DO
#1) Add distinct labile / foliage pools

###
## Prepare the work space

# Set working directory
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/")

# Define the output file name for the created script
outfilename_prefix = "miombo_latex_C_budget"
outfilename = paste(outfilename_prefix,".tex",sep="")

###
## Load files from which C-budget is extracted, prepare C-budget values

# Load information file
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/reccap2_permafrost_1deg_dalec4_isimip3a_agb_lca_nbe_CsomPriorNCSDC3m/infofile.RData")
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/Miombo_0.5deg_allWood/infofile.RData")

# If this is a site analysis 
# which site are we using?
site_nos = 1

###
## Define figure labels

# The caption must be written in correct latex syntax
# NOTE: if the latex language requires use of "\" ensure it is a "\\".
figure_caption = "CARDAMOM C-cycle analysis of Miombo region, monthly time step between 2006-2017. Numbers show median estimate of fluxes (alongside arrows) and of stocks (in boxes). Units are MgC ha$^{-1}$ for stocks and MgC ha$^{-1}$ y$^{-1}$ for fluxes. 95\\% confidence intervals are shown in a fractional form with 2.5 and 97.5 percentiles as numerator and denominator. Black fluxes are biogenic, including net primary production ($NPP$), mortality ($Mort$), autotrophic respiration ($Ra$) and heterotrophic respiration ($Rh$). $NEE = Ra + Rh - GPP$. $NBP = -NEE - E_{total} - Forest_loss$ (not shown). Red fluxes are fire-driven emissions ($E$)."
# The label will be used for referencing the figure in the latex document
figure_label = "SIFig:global_budget"
# Desired precision, i.e. decimal places
dp = 1

if (PROJECT$spatial_type == "grid") {
    # A gridded analysis
    
    # Specify quantiles we want to extract, should only be 3 (low/median/high).
    # From gridded analysis this must be from those available in the file. 
    # Check grid_output$num_quantiles for available
    quantiles_wanted = c(1,4,7)

    # Load the processed site file
    load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))

    # Extract or calculate required derived values
    # NOTE: unit conversion from gC/m2/day to MgC/ha/yr

    # NATURAL FLUXES
    gpp_gCm2yr = format(round(apply(grid_output$mean_gpp_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    rauto_gCm2yr = format(round(apply(grid_output$mean_rauto_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    rhet_litter_gCm2yr = format(round(apply(grid_output$mean_rhet_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    rhet_som_gCm2yr = format(round(apply(grid_output$mean_rhet_som_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    rhet_gCm2yr = format(round(apply(grid_output$mean_rhet_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    npp_gCm2yr = format(round(apply(grid_output$mean_npp_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    npp_labile_gCm2yr = format(round(apply(grid_output$mean_alloc_labile_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    npp_foliage_gCm2yr = format(round(apply(grid_output$mean_alloc_foliage_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    labile_to_foliage_gCm2yr = format(round(apply(grid_output$mean_labile_to_foliage_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    npp_roots_gCm2yr = format(round(apply(grid_output$mean_alloc_roots_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    npp_wood_gCm2yr = format(round(apply(grid_output$mean_alloc_wood_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    foliage_to_litter_gCm2yr = format(round(apply(grid_output$mean_foliage_to_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    roots_to_litter_gCm2yr = format(round(apply(grid_output$mean_roots_to_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    wood_to_litter_gCm2yr = format(round(apply(grid_output$mean_wood_to_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    litter_to_som_gCm2yr = format(round(apply(grid_output$mean_litter_to_som_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    # FIRE FLUXES
    FIRElitter_labile_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_labile_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIRElitter_foliage_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_foliage_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIRElitter_roots_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_roots_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIRElitter_wood_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_wood_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIRElitter_litter_gCm2yr = format(round(apply(grid_output$mean_FIRElitter_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIREemiss_labile_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_labile_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIREemiss_foliage_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_foliage_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIREemiss_roots_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_roots_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIREemiss_wood_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_wood_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIREemiss_litter_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_litter_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIREemiss_som_gCm2yr = format(round(apply(grid_output$mean_FIREemiss_som_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    fire_gCm2yr = format(round(apply(grid_output$mean_fire_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    # Net fluxes
    nbe_gCm2yr = format(round(apply(grid_output$mean_nbe_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    nee_gCm2yr = format(round(apply(grid_output$mean_nee_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    nbp_gCm2yr = format(round(apply(grid_output$mean_nbp_gCm2day[,,quantiles_wanted],3,mean, na.rm=TRUE) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    # STOCKS
    labile_gCm2 = format(round(apply(grid_output$mean_labile_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE) * 1e-2, digits = dp), nsmall = dp)
    foliage_gCm2 = format(round(apply(grid_output$mean_foliage_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE) * 1e-2, digits = dp), nsmall = dp)
    roots_gCm2 = format(round(apply(grid_output$mean_roots_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE) * 1e-2, digits = dp), nsmall = dp)
    wood_gCm2 = format(round(apply(grid_output$mean_wood_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE) * 1e-2, digits = dp), nsmall = dp)
    litter_gCm2 = format(round(apply(grid_output$mean_litter_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE) * 1e-2, digits = dp), nsmall = dp)
    som_gCm2 = format(round(apply(grid_output$mean_som_gCm2[,,quantiles_wanted],3,mean, na.rm=TRUE) * 1e-2, digits = dp), nsmall = dp)
  
} else if (PROJECT$spatial_type == "site") {
    # A site analysis
        
    # Load the processed site file
    load(paste(PROJECT$results_processedpath,PROJECT$sites[site_nos],".RData",sep=""))

    # Specify quantiles we want to extract, should only be 3 (low/median/high)
    quantiles_wanted = c(0.025,0.50,0.975)    
    # Restrict the time bounds?
    s = 1 ; f = dim(states_all$lai_m2m2)[2]
#    s = f-(17*12) ; f = dim(states_all$lai_m2m2)[2]
    
    # Extract or calculate required derived values
    # NOTE: unit conversion from gC/m2/day to MgC/ha/yr

    # NATURAL FLUXES
    gpp_gCm2yr = format(round(quantile(apply(states_all$gpp_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    rauto_gCm2yr = format(round(quantile(apply(states_all$rauto_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    rhet_litter_gCm2yr = format(round(quantile(apply(states_all$rhet_litter_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    rhet_som_gCm2yr = format(round(quantile(apply(states_all$rhet_som_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    rhet_gCm2yr = format(round(quantile(apply(states_all$rhet_litter_gCm2day[,s:f]+states_all$rhet_som_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    npp_gCm2yr = format(round(quantile(apply(states_all$gpp_gCm2day[,s:f]-states_all$rauto_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    npp_labile_gCm2yr = format(round(quantile(apply(states_all$alloc_labile_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    npp_foliage_gCm2yr = format(round(quantile(apply(states_all$alloc_foliage_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    labile_to_foliage_gCm2yr = format(round(quantile(apply(states_all$labile_to_foliage_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    npp_roots_gCm2yr = format(round(quantile(apply(states_all$alloc_roots_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    npp_wood_gCm2yr = format(round(quantile(apply(states_all$alloc_wood_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    foliage_to_litter_gCm2yr = format(round(quantile(apply(states_all$foliage_to_litter_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    roots_to_litter_gCm2yr = format(round(quantile(apply(states_all$roots_to_litter_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    wood_to_litter_gCm2yr = format(round(quantile(apply(states_all$wood_to_litter_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    litter_to_som_gCm2yr = format(round(quantile(apply(states_all$litter_to_som_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    # FIRE FLUXES
    FIRElitter_labile_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_labile_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIRElitter_foliage_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_foliage_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIRElitter_roots_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_roots_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIRElitter_wood_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_wood_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIRElitter_litter_gCm2yr = format(round(quantile(apply(states_all$FIRElitter_litter_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIREemiss_labile_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_labile_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIREemiss_foliage_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_foliage_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIREemiss_roots_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_roots_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIREemiss_wood_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_wood_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIREemiss_litter_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_litter_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    FIREemiss_som_gCm2yr = format(round(quantile(apply(states_all$FIREemiss_som_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    fire_gCm2yr = format(round(quantile(apply(states_all$fire_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    # Net fluxes
    nbe_gCm2yr = format(round(quantile(apply((states_all$fire_gCm2day+states_all$rhet_litter_gCm2day+states_all$rhet_som_gCm2day+states_all$rauto_gCm2day)[,s:f]-states_all$gpp_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    nee_gCm2yr = format(round(quantile(apply((states_all$rhet_litter_gCm2day+states_all$rhet_som_gCm2day+states_all$rauto_gCm2day)[,s:f]-states_all$gpp_gCm2day[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    nbp_gCm2yr = format(round(quantile(apply((states_all$nbp_gCm2day)[,s:f],1,mean), prob=quantiles_wanted) * 365.25 * 1e-2, digits = dp), nsmall = dp)
    # STOCKS
    labile_gCm2 = format(round(quantile(apply(states_all$labile_gCm2[,s:f],1,mean) * 1e-2, prob=quantiles_wanted), digits = dp), nsmall = dp)
    foliage_gCm2 = format(round(quantile(apply(states_all$labile_gCm2[,s:f] + states_all$foliage_gCm2[,s:f],1,mean) * 1e-2, prob=quantiles_wanted), digits = dp), nsmall = dp)
    roots_gCm2 = format(round(quantile(apply(states_all$roots_gCm2[,s:f],1,mean), prob=quantiles_wanted) * 1e-2, digits = dp), nsmall = dp)
    wood_gCm2 = format(round(quantile(apply(states_all$wood_gCm2[,s:f],1,mean), prob=quantiles_wanted) * 1e-2, digits = dp), nsmall = dp)
    litter_gCm2 = format(round(quantile(apply(states_all$litter_gCm2[,s:f],1,mean), prob=quantiles_wanted) * 1e-2, digits = dp), nsmall = dp)
    som_gCm2 = format(round(quantile(apply(states_all$som_gCm2[,s:f],1,mean), prob=quantiles_wanted) * 1e-2, digits = dp), nsmall = dp)
    
} else {
    # We have a compatibility problem
    stop("PROJECT$spatial_type does not have a compatible value, i.e. grid or site")
} # end if is_grid

###
## Begin writing the latex code

col_sep = ""
nos_cols = 20

write(    c("\\documentclass{article}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = FALSE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\title{C budget template figure}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage[T1]{fontenc}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{verbatim}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{color}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{hyperref}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{cleveref}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{fixmath}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{ulem}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{lscape}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{subfigure}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{array,multirow,graphicx}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{chngcntr}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage[landscape, margin=2cm]{geometry}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\usepackage{helvet}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\renewcommand{\\familydefault}{\\sfdefault}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\begin{document}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("% Biogenic flux and emissions figure with budgets"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\begin{figure}[]"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("% NOTE: \\put(x coord,y coord){ ... } where to put ..."), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("%       \\vector(x slope,y slope){length} used for arrows"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("%       \\line(x slope,y slope){length} used for lines"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("%       \\framebox(x,y){...} puts ... at box centre"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("   \\centering"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("   \\fbox{"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\setlength{\\unitlength}{0.60cm}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\begin{picture}(28,14)"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % GPP"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(0.0,4.50){\\vector(1,0){3.25}} % arrow for GPP"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(0.76,4.7){$GPP$}               % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(0.0,3.9){\\small ",gpp_gCm2yr[2],"}                 % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(1.66,3.925){\\scriptsize $\\frac{",gpp_gCm2yr[1],"}{",gpp_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Ra"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(26.0,0.78){$Ra$}              % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(25.2,0.1){\\small ",rauto_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(26.9,0.1){\\scriptsize $\\frac{",rauto_gCm2yr[1],"}{",rauto_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % CUE and associated arrows         "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(3.37,4.35){$CUE$} % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(3.20,3.70){\\dashbox{0.2}(1.75,1.65)} % Box around CUE"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(4.075,3.70){\\line(0,-1){3.05}} % Vertical line down from box"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(4.075,0.65){\\vector(1,0){21.925}} % horizontal line to Ra"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)         
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Add labels and box for the internal C-cycle dynamics"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(3.0,10.70){\\textit{Internal carbon rates}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(2.95,0){\\dashbox{0.2}(20.55,11.2)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Add labels for input / output"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(0.1,13.5){\\textit{Input}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(0.1,13.0){\\textit{carbon}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(0.1,12.5){\\textit{rates}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(26.5,13.5){\\textit{Output}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(26.5,13.0){\\textit{carbon}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(26.5,12.5){\\textit{rates}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Total NPP"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(4.95,4.5){\\vector(1,0){3.0}} % arrow for NPP"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(5.7,4.7){$NPP$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(5.05,3.9){\\small ",npp_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(6.6,3.925){\\scriptsize $\\frac{",npp_gCm2yr[1],"}{",npp_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Add partitioning point in graphic"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(8.40,4.75){\\small \\rotatebox[origin=c]{90}{NPP allocation}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(7.95,1.0){\\dashbox{0.2}(1.2,7)} % Add box around label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % NPP to labile"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.15,6.3){\\vector(1,0){3.6}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.35,6.5){\\small $NPP_{lab}$}    % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(9.55,5.70){\\small ",npp_labile_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(10.72,5.75){\\scriptsize $\\frac{",npp_labile_gCm2yr[1],"}{",npp_labile_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % NPP to foliage"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(8.92,8.6){\\rotatebox[origin=c]{30}{\\vector(1,0){4.0}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.15,8.95){\\small \\rotatebox[origin=c]{30}{$NPP_{fol}$}}    % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(9.30,8.2){\\small \\rotatebox[origin=c]{30}{",npp_foliage_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(10.60,9.0){\\scriptsize \\rotatebox[origin=c]{30}{$\\frac{",npp_foliage_gCm2yr[1],"}{",npp_foliage_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % labile to foliage"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(13.95,8.05){\\vector(0,1){0.75}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(12.65,8.20){\\small ",labile_to_foliage_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(14.1,8.30){\\scriptsize $\\frac{",labile_to_foliage_gCm2yr[1],"}{",labile_to_foliage_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % NPP to fine root"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.15,4.5){\\vector(1,0){3.6}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.35,4.7){\\small $NPP_{root}$} % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(9.55,3.9){\\small ",npp_roots_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(10.72,3.95){\\scriptsize $\\frac{",npp_roots_gCm2yr[1],"}{",npp_roots_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % NPP to wood"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.15,2){\\vector(1,0){3.6}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(9.35,2.2){\\small $NPP_{wood}$} % Label"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(9.55,1.45){\\small ",npp_wood_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(10.72,1.50){\\scriptsize $\\frac{",npp_wood_gCm2yr[1],"}{",npp_wood_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Labile C pools"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.70,6.0){\\framebox(2.5,2.05)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(13.20,7.45){$C_{lab}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)              
write(paste("         \\put(13.30,6.80){\\small ",labile_gCm2[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(13.40,6.20){\\scriptsize $\\frac{",labile_gCm2[1],"}{",labile_gCm2[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Foliage C pools"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.70,8.75){\\framebox(2.5,2.05)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(13.20,10.2){$C_{fol}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)              
write(paste("         \\put(13.30,9.55){\\small ",foliage_gCm2[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(13.40,8.95){\\scriptsize $\\frac{",foliage_gCm2[1],"}{",foliage_gCm2[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Fine root C pools"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.70,3.5){\\framebox(2.5,2.05)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(13.20,4.95){$C_{root}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)              
write(paste("         \\put(13.30,4.35){\\small ",roots_gCm2[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(13.40,3.75){\\scriptsize $\\frac{",roots_gCm2[1],"}{",roots_gCm2[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Wood C pools"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.70,1){\\framebox(2.5,2.05)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(13.1,2.45){$C_{wood}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)              
write(paste("         \\put(13.10,1.85){\\small ",wood_gCm2[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(13.20,1.25){\\scriptsize $\\frac{",wood_gCm2[1],"}{",wood_gCm2[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % litter/mortality, natural and fire driven (red)"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Labile - fire only, labile to foliage carried out above in allocation"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.7,7.2){\\small $Mort_{lab}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(18.0,6.4){\\color{red}{\\small ",FIRElitter_labile_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(18.95,6.45){\\color{red}{\\scriptsize $\\frac{",FIRElitter_labile_gCm2yr[1],"}{",FIRElitter_labile_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Assign arrow to Clitter"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.25,7.0){\\line(1,0){4.45}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(19.70,7.0){\\vector(1,-1){0.6}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Foliage"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.7,9.95){\\small $Mort_{fol}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(15.65,9.15){\\small ",foliage_to_litter_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(16.85,9.20){\\scriptsize $\\frac{",foliage_to_litter_gCm2yr[1],"}{",foliage_to_litter_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(18.0,9.15){\\color{red}{\\small ",FIRElitter_foliage_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(18.95,9.20){\\color{red}{\\scriptsize $\\frac{",FIRElitter_foliage_gCm2yr[1],"}{",FIRElitter_foliage_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Assign arrow to Clitter"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.25,9.75){\\line(1,0){4.45}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(19.70,9.75){\\vector(1,-2){1.25}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Fine root"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.7,4.65){\\small $Mort_{root}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(15.65,3.9){\\small ",roots_to_litter_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(16.85,3.95){\\scriptsize $\\frac{",roots_to_litter_gCm2yr[1],"}{",roots_to_litter_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(18.0,3.9){\\color{red}{\\small ",FIRElitter_roots_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(18.95,3.95){\\color{red}{\\scriptsize $\\frac{",FIRElitter_roots_gCm2yr[1],"}{",FIRElitter_roots_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Assign arrow to Clitter"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.25,4.5){\\line(1,0){4.45}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(19.70,4.5){\\vector(1,2){0.6}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Wood"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.70,2.2){\\small $Mort_{wood}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(15.65,1.45){\\small ",wood_to_litter_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(16.85,1.50){\\scriptsize $\\frac{",wood_to_litter_gCm2yr[1],"}{",wood_to_litter_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(18.0,1.45){\\color{red}{\\small ",FIRElitter_wood_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(18.95,1.50){\\color{red}{\\scriptsize $\\frac{",FIRElitter_wood_gCm2yr[1],"}{",FIRElitter_wood_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Assign arrow to Csom"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.25,2.0){\\line(1,0){4.45}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(19.70,2.0){\\vector(1,1){0.6}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Fire emission fluxes"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Labile"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.0,12.4){\\small \\rotatebox[origin=c]{90}{$E_{lab}$}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.5,7.5){\\color{red}{\\line(0,1){6.0}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(12.5,7.5){\\color{red}{\\line(1,0){0.2}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(12.7,12.7){\\color{red}{\\small ",FIREemiss_labile_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(12.7,12.1){\\color{red}{\\scriptsize $\\frac{",FIREemiss_labile_gCm2yr[1],"}{",FIREemiss_labile_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Foliage"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(13.57,12.4){\\small \\rotatebox[origin=c]{90}{$E_{fol}$}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(14.2,10.8){\\color{red}{\\line(0,1){2.7}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(14.3,12.7){\\color{red}{\\small ",FIREemiss_foliage_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(14.3,12.1){\\color{red}{\\scriptsize $\\frac{",FIREemiss_foliage_gCm2yr[1],"}{",FIREemiss_foliage_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Fine roots"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(10.45,12.4){\\small \\rotatebox[origin=c]{90}{$E_{root}$}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(11.85,5.0){\\color{red}{\\line(1,0){0.85}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(11.85,5.0){\\color{red}{\\line(0,1){8.5}}} "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(11.00,12.7){\\color{red}{\\small ",FIREemiss_roots_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(10.95,12.1){\\color{red}{\\scriptsize $\\frac{",FIREemiss_roots_gCm2yr[1],"}{",FIREemiss_roots_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Wood"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.05,12.4){\\small \\rotatebox[origin=c]{90}{$E_{wood}$}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.2,2.7){\\color{red}{\\line(1,0){0.35}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(15.55,2.7){\\color{red}{\\line(0,1){10.8}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(15.65,12.7){\\color{red}{\\small ",FIREemiss_wood_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(15.70,12.1){\\color{red}{\\scriptsize $\\frac{",FIREemiss_wood_gCm2yr[1],"}{",FIREemiss_wood_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Litter"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20.45,12.4){\\small \\rotatebox[origin=c]{90}{$E_{litter}$}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20.95,7.2){\\color{red}{\\line(0,1){6.3}}} "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(21.05,12.7){\\color{red}{\\small ",FIREemiss_litter_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(21.1,12.1){\\color{red}{\\scriptsize $\\frac{",FIREemiss_litter_gCm2yr[1],"}{",FIREemiss_litter_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % SOM"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(22.5,12.4){\\small \\rotatebox[origin=c]{90}{$E_{som}$}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(22.75,3.4){\\color{red}{\\line(1,0){0.3}}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(23.05,3.4){\\color{red}{\\line(0,1){10.1}}} "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(23.1,12.7){\\color{red}{\\small ",FIREemiss_som_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(23.1,12.1){\\color{red}{\\scriptsize $\\frac{",FIREemiss_som_gCm2yr[1],"}{",FIREemiss_som_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Fire emissions total"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(24.1,13.5){$E_{total}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(11.85,13.5){\\color{red}{\\vector(1,0){12.1}}} "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(24.0,12.95){\\color{red}{\\small ",fire_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(25.0,13.00){\\color{red}{\\scriptsize $\\frac{",fire_gCm2yr[1],"}{",fire_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Litter C pool"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20.25,5.2){\\framebox(2.5,2)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20.75,6.45){$C_{litter}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(20.4,5.6){\\small ",litter_gCm2[2],"}               % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(21.6,5.6){\\scriptsize $\\frac{",litter_gCm2[1],"}{",litter_gCm2[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Som C pool"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20.25,1.0){\\framebox(2.5,2.75)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(20.65,2.9){$C_{som}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(20.55,2.1){\\small ",som_gCm2[2],"}                 % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(20.75,1.45){\\scriptsize $\\frac{",som_gCm2[1],"}{",som_gCm2[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Decomposition - natural and fire"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(21.5,5.2){\\vector(0,-1){1.40}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Decomposition - natural"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(21.7,4.65){\\small ",litter_to_som_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(21.65,4.1){\\scriptsize $\\frac{",litter_to_som_gCm2yr[1],"}{",litter_to_som_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Decomposition - combusted litter to som"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(20.4,4.65){\\color{red}{\\small ",FIRElitter_litter_gCm2yr[2],"}}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(20.5,4.1){\\color{red}{\\scriptsize $\\frac{",FIRElitter_litter_gCm2yr[1],"}{",FIRElitter_litter_gCm2yr[3],"}$}} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Heterotrophic respiration of litter"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(22.75,6.0){\\line(1,0){2.25}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(23.65,6.3){\\small $Rh_{litter}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(23.65,5.45){\\small ",rhet_litter_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(25.25,5.50){\\scriptsize $\\frac{",rhet_litter_gCm2yr[1],"}{",rhet_litter_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(25.0,6.0){\\vector(1,-2){0.70}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Heterotrophic respiration of som"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(22.75,2.75){\\line(1,0){2.25}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(23.65,3.0){\\small $Rh_{som}$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(23.65,2.2){\\small ",rhet_som_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(25.25,2.25){\\scriptsize $\\frac{",rhet_som_gCm2yr[1],"}{",rhet_som_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(25.0,2.75){\\vector(1,1){0.85}}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Total heterotrophic respiration"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(26.0,4.6){$Rh$}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(25.2,4.0){\\small ",rhet_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(26.9,4.0){\\scriptsize $\\frac{",rhet_gCm2yr[1],"}{",rhet_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Net Biome Productivity "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(2.95,11.5){\\dashbox{0.2}(2.2,2.2)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(3.45,13.05){NBP}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(3.45,12.4){\\small ",nbp_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(3.45,11.8){\\scriptsize $\\frac{",nbp_gCm2yr[1],"}{",nbp_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         % Net Ecosystem Exchange "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(5.95,11.5){\\dashbox{0.2}(2.2,2.2)}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\put(6.45,13.05){NEE}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(6.45,12.4){\\small ",nee_gCm2yr[2],"}                % Median",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(paste("         \\put(6.45,11.8){\\scriptsize $\\frac{",nee_gCm2yr[1],"}{",nee_gCm2yr[3],"}$} % CI",sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(" "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("         \\end{picture}"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("        }"), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(paste("   \\caption{",figure_caption,"}", sep="")), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c(paste("   \\label{",figure_label,"}"), sep=""), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\end{figure}   "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)
write(    c("\\end{document}   "), file = outfilename, ncolumns = nos_cols, sep=col_sep, append = TRUE)


# Convert the latex document into a pdf 
system(paste("pdflatex ",outfilename,sep=""))
# Convert the pdf into a png
system(paste("pdftoppm -png ",outfilename," ",outfilename_prefix,sep=""))

# DONE!
