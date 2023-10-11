   
###
## Generic script to generate summary information on the simulated C budget
## and compare CARDAMOM estimates with its assimilated and independent observations
## NOTE: we assume this is a gridded run
###

## Sections of code you may want to change are flagged
## "PointsOfChange"

#### TO DO
# Add extractions for specific domains (recaap2)- this should be based on a loop - this should also have time sometime series information
# What happens to SS in the other pools
# What is the overlap for the Csom prior
# Does concistency between variables show an association i.e. is the consistency with LAI correlated with consistency with Cwood etc.
# Comparison between CARDAMOM and EO soil moisture

###
## Analysis specific information and generic creation
###

# Set the working directory for your CARDAMOM code base
setwd("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM")

###
## Load analysis

# PointsOfChange
#load("/exports/csce/datastore/geos/users/lsmallma/CARDAMOM_R_OUTPUT/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/reccap2_permafrost_1deg_dalec2_isimip3a_agb_lca_nbe_CsomPriorNCSDC3m/infofile.RData")
#load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/Miombo_0.5deg_allWood/infofile.RData")
#load("/exports/csce/datastore/geos/users/lsmallma/CARDAMOM_R_OUTPUT/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_2_2.5deg_C7_GCP_AGB/infofile.RData")
#load("/exports/csce/datastore/geos/users/lsmallma/CARDAMOM_R_OUTPUT/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_2_2.5deg_C7_GCP_oneAGB/infofile.RData")
load("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/CARDAMOM_OUTPUTS/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/global_1deg_dalec4_trendyv12_LCA_AGB_NBE/infofile.RData")

# Set output path for figures and tables
out_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/ESSD_update/figures/"
#out_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/RECCAP2/figures/"
#out_dir = "/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/LTSS_CARBON_INTEGRATION/InternationalScience/figures_africa/"
#out_dir = "~/WORK/GREENHOUSE/models/CARDAMOM/SECO/figures/"

#
# Load the CARDAMOM files
load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))

# Specify the position within the stored ensemble for the median estimate and the desired uncertainty bands
mid_quant = 4 ; low_quant = 2 ; high_quant = 6
wanted_quant = c(low_quant,3,mid_quant,5,high_quant)

# Extract timing information
run_years = as.numeric(PROJECT$start_year) : as.numeric(PROJECT$end_year)
nos_years = length(as.numeric(PROJECT$start_year) : as.numeric(PROJECT$end_year))
steps_per_year = length(PROJECT$model$timestep_days) / nos_years

###
## Load needed libraries and framework functions
###

# Read library
library(fields)
library(compiler)
library(terra)
library(RColorBrewer)
library(plotrix)
library(zoo)
library(ncdf4)
library(abind)
library(corrplot)

# Load any CARDAMOM functions which might be useful
source("./R_functions/generate_wgs_grid.r")
source("./R_functions/calc_pixel_area.r")
source("./R_functions/read_binary_file_format.r")
source("./R_functions/function_closest2d.r")
source("./R_functions/plotconfidence.r")
source("./R_functions/read_src_model_priors.r")
source("./R_functions/regrid_functions.r")

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
## Do parameter correlation plots here

par_tmp = array(grid_output$parameters[,,,mid_quant], dim=c(prod(dim(grid_output$parameters)[1:2]),dim(grid_output$parameters)[3]))
filter = which(is.na(par_tmp[,1]))
par_tmp = par_tmp[-filter,]
par_tmp = cor(par_tmp)

# In order of p1-33         
par_labels = c("Decomp efficiency","RaGPP","GPPfol","GPProot","Leaf lifespan","Wood turnover","Fine root turnover","Litter mineralisation"
              ,"Soil mineralisation","Rhet coefficient","Ceff","Leaf growth day","GPPlab","Leaf growth period","Leaf fall day"
              ,"Leaf fall period","LMA","Initial labile","Initial foliage","Initial fine root","Initial wood","Initial litter"
              ,"Initial soil","Initial soil water","Coarse root fraction","Root depth coefficient","Max rooting depth","Fire resilience"
              ,"Foliage CC","Wood+root CC","Soil CC"
              ,"litter CC","log-likelihood")         

colnames(par_tmp) <- par_labels
rownames(par_tmp) <- par_labels

col <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_between_pixel_parameter_correlations.png",sep=""), height = 5000, width = 5000, res = 300)
par(mfrow=c(1,1), mar=c(0,0,0,0), omi=c(0.01,0.01,0.01,0.01))
corrplot(par_tmp, method="color", col=col(200),  
         type="upper", tl.cex=1.2, number.cex=0.8, cl.cex=1.2,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
         )
dev.off()

par_tmp = array(grid_output$parameters[,,,mid_quant], dim=c(prod(dim(grid_output$parameters)[1:2]),dim(grid_output$parameters)[3]))
filter = which(is.na(par_tmp[,1]))
par_tmp = par_tmp[-filter,]
var_tmp = as.vector(grid_output$mean_nee_gCm2day[,,mid_quant])
var_tmp = cbind(var_tmp,as.vector(grid_output$mean_gpp_gCm2day[,,mid_quant]))
var_tmp = cbind(var_tmp,as.vector(grid_output$mean_rauto_gCm2day[,,mid_quant]))
var_tmp = cbind(var_tmp,as.vector(grid_output$mean_rhet_gCm2day[,,mid_quant]))
var_tmp = cbind(var_tmp,as.vector(grid_output$mean_fire_gCm2day[,,mid_quant]))
var_tmp = var_tmp[-filter,]
par_tmp = t(cor(var_tmp,par_tmp))

colnames(par_tmp) <- c("NEE","GPP","Rauto","Rhet","fire")
rownames(par_tmp) <- par_labels

png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_between_pixel_parameter_correlations_fluxes.png",sep=""), height = 5000, width = 2000, res = 300)
par(mfrow=c(1,1), mar=c(0,0,0,0), omi=c(0.01,0.01,0.01,0.01))
corrplot(par_tmp, method="color", col=col(200), 
         type="full", tl.cex=1.2, number.cex=0.8, cl.cex=1.2, cl.ratio = 0.3,
         addCoef.col = "black", # Add coefficient of correlation
         tl.col="black", tl.srt=45, #Text label color and rotation
         # hide correlation coefficient on the principal diagonal
         diag=FALSE 
         )
dev.off()

###
## Determine needed spatial information

# generate the lat / long grid again
output = generate_wgs84_grid(PROJECT$latitude,PROJECT$longitude,PROJECT$resolution)
grid_lat = array(output$lat, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
grid_long = array(output$long,dim=c(PROJECT$long_dim,PROJECT$lat_dim))
cardamom_ext = output$cardamom_ext
# then generate the area estimates for each pixel (m)
area = calc_pixel_area(grid_long,grid_lat)

###
## Allow for optional replacement or creation of alternate clustering analyses
## CARDAMOM default is affinity propogation but the most commonly used is k-means

# PointsOfChange

# Now generate cluster map using all parameters including the initial conditions
par_array_median_normalised = grid_output$parameters[,,,mid_quant]
# now normalise the parameter values
for (i in seq(1,dim(par_array_median_normalised)[3])) {
     min_par_val=min(par_array_median_normalised[,,i],na.rm=TRUE)
     max_par_val=max(par_array_median_normalised[,,i],na.rm=TRUE)
     par_array_median_normalised[,,i]=((par_array_median_normalised[,,i]-min_par_val)/(max_par_val-min_par_val))
}
# Create temporary arrays needed to allow removing of NAs and convert array into (space,par)
par_array_tmp=array(NA,dim=c(prod(dim(grid_output$parameters)[1:2]),dim(par_array_median_normalised)[3]))
par_array_tmp[1:prod(dim(par_array_median_normalised)[1:2]),1:dim(par_array_median_normalised)[3]]=par_array_median_normalised
actual_forests=which(is.na(par_array_tmp[,1]) == FALSE)
par_array_tmp=par_array_tmp[actual_forests,]
par_array_tmp=array(par_array_tmp,dim=c((length(par_array_tmp)/dim(par_array_median_normalised)[3]),dim(par_array_median_normalised)[3]))

# K-means
nos_desired_clusters = 3 # specify the number of clusters wanted
cluster = kmeans(par_array_tmp, centers = nos_desired_clusters, iter.max = 10, nstart = 50)

# Extract the cluster information which we can then use in maps / aggregation
grid_output$nos_pars_clusters_kmeans=dim(cluster$centers)[1] ; grid_output$clusters_exemplars_kmeans=cluster$centers
grid_output$pars_clusters_kmeans=array(NA,dim=c(dim(par_array_median_normalised)[1:2]))
for (i in seq(1, grid_output$nos_pars_clusters_kmeans)) {
     grid_output$clusters_kmeans[actual_forests[which(cluster$cluster == i)]] = i
}
grid_output$clusters_kmeans=array(grid_output$clusters_kmeans,dim=c(dim(par_array_median_normalised)[1:2]))

#par(mfrow=c(2,2)) ; image.plot(grid_output$pars_clusters) ; image.plot(grid_output$pars_clusters_kmeans)
image.plot(grid_output$clusters_kmeans, main=3)

### 
## Create land mask / boundary overlays needed

# Load any map information you might want to overlay, such as biome or national boundaries
# if you don't want to have any set as false

# PointsOfChange

# Load a baseline land mask
# load global shape file for land sea mask
landmask = vect("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/R_functions/global_map/national_boundaries/ne_10m_admin_0_countries.shx")
# subset by continent (could also do by country)
#landmask = subset(landmask, subset=landmask$CONTINENT == "South America") # Change continent to target area or comment out if spanning zones
#landmask = subset(landmask, subset=landmask$CONTINENT == "Africa") # Change continent to target area or comment out if spanning zones
# Clip and/or extend to the extent of the CARDAMOM analysis
landmask = crop(landmask, cardamom_ext)
landmask_area = rasterize(landmask, cardamom_ext)
landmask_area = extend(landmask_area, cardamom_ext)

# Create an updated area object for the landmask region only
# extract the lat / long information needed
#long = crds(crop(cardamom_ext,landmask),df=TRUE, na.rm=FALSE)
long = crds(landmask_area,df=TRUE, na.rm=FALSE)
tmp_lat = long$y ; tmp_long = long$x ; tmp_lat = tmp_lat[length(tmp_lat):1]
tmp = dim(landmask_area)[c(2,1)]
tmp_lat = array(tmp_lat, dim=c(tmp[1],tmp[2]))
tmp_long = array(tmp_long,dim=c(tmp[1],tmp[2]))
# then generate the area estimates for each pixel (m)
landmask_area = calc_pixel_area(tmp_long,tmp_lat)

rm(tmp)

add_biomes = " "
#add_biomes = "ssa_wwf"
#add_biomes = "wwf_ecoregions"
#add_biomes = "reccap2_permafrost"
if (add_biomes == "ssa_wwf") {
    # Read in shape file for boundaries
    biomes = shapefile("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/SECO/analysis/ssa_wwf_dissolve/ssa_wwf_dissolve.shp")
    # just to be sure enforce the projection to WGS-84
    biomes = spTransform(biomes, crs(cardamom_ext))
    # Clip to the extent of the CARDAMOM analysis
    biomes = crop(biomes, cardamom_ext)
    # Clip to the area covered bythe landmask in case of geographic restriction
    biomes = crop(biomes, landmask)
    # Overwrite the existing landmask
    landmask = biomes
    # Extract the current biome names for raster code
    biome_names = levels(factor(biomes$ECO_NAME))
    # Turn biomes object into a raster now in case we want to use it later
    biomes = rasterize(biomes, cardamom_ext, factor(biomes$ECO_NAME, fun="last"))
} else if (add_biomes == "wwf_ecoregions") {
    # Read in shape file for boundaries
    biomes = shapefile("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/LTSS_CARBON_INTEGRATION/wwf_biomes_map/Ecoregions2017.shp")
    # just to be sure enforce the projection to WGS-84
    biomes = spTransform(biomes, crs(cardamom_ext))
    # Clip to the extent of the CARDAMOM analysis
    biomes = crop(biomes, cardamom_ext)
    # Clip to the area covered by the landmask in case of geographic restriction
    biomes = crop(biomes, landmask)
    # Aggregate to the biome regions rather than the more complex ecoregions
    biomes <- aggregate(biomes, by='BIOME_NAME', dissolve = TRUE)

    # Overwrite the existing landmask
    landmask = biomes

    # Extract the current biome names for raster code
    biome_names = levels(factor(biomes$BIOME_NAME))
    # Turn the Biomes object into a raster now
    biomes = rasterize(biomes,cardamom_ext, factor(biomes$BIOME_NAME), fun="last")
} else if (add_biomes == "reccap2_permafrost") {

    # open reccap2-permafrost region mask
    input_file_1=paste("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/ESSD_update/RECCAP2_regions/RECCAP2_permafrost_regions_isimip3.nc",sep="")
    data1=nc_open(input_file_1)

    # extract location variables
    biomes_lat=ncvar_get(data1, "latitude") ; biomes_long=ncvar_get(data1, "longitude")
    # Create lat / long grid from vectors
    idim = length(biomes_long) ; jdim = length(biomes_lat)
    biomes_lat = array(biomes_lat, dim=c(jdim,idim)) ; biomes_lat = t(biomes_lat)
    biomes_long = array(biomes_long, dim=c(idim,jdim))
    # read the mask
    biomes=ncvar_get(data1, "permafrost_region_mask") 
    # For the moment convert each different sub-region into 1 and the missing number flag into NA
    biomes[biomes > 0 & biomes < 115] = 1 ; biomes[biomes > 1] = NA
    biome_names = NA
    # tidy up
    nc_close(data1)

    # Convert to a raster, assuming standad WGS84 grid
    biomes = data.frame(x = as.vector(biomes_long), y = as.vector(biomes_lat), z = as.vector(biomes))
    biomes = rast(biomes, crs = ("+init=epsg:4326"), type="xyz")
              
    # Create raster with the target crs
    target = rast(crs = ("+init=epsg:4326"), ext = ext(biomes), resolution = res(biomes))
    # Check whether the target and actual analyses have the same CRS
    if (compareCRS(biomes,target) == FALSE) {
        # Resample to correct grid
        biomes = resample(biomes, target, method="ngb") ; gc() ; removeTmpFiles()
    }
    # Trim the extent of the overall grid to the analysis domain
    biomes = crop(biomes,cardamom_ext) 
    # Match resolutions if they differ
    spatial_type = "grid"
    if (spatial_type == "grid") {
        if (res(biomes)[1] != res(cardamom_ext)[1] | res(biomes)[2] != res(cardamom_ext)[2]) {

            # Create raster with the target resolution
            target = rast(crs = crs(cardamom_ext), ext = ext(cardamom_ext), resolution = res(cardamom_ext))
            # Resample to correct grid
            biomes = resample(biomes, target, method="bilinear") ; gc() ; removeTmpFiles()
        } # Aggrgeate to resolution
    } # spatial_type == "grid"

    # Convert to polygons
    biomes = rasterToPolygons(biomes, n=16, na.rm=TRUE, digits=12, dissolve=TRUE)

    # Overwrite the existing landmask
    landmask = biomes

} else {
    # DO NOTHING
    biomes = NA ; biome_names = NA
} 

# Write biome names to a file
write.table(data.frame(BiomeCode = c(1:length(biome_names)),BiomeNames = biome_names), 
            file = paste(out_dir,"/",PROJECT$name,"_biome_names.csv",sep=""), row.names=FALSE, sep=",",append=FALSE)

# PointsOfChange

# This will be used to filter the analysis to include specific locations only
use_filter = TRUE
if (use_filter) {
    #  Design a user created / loaded filter 
    landfilter = array(NA,dim=dim(grid_output$assimilated_wood_mean_gCm2))
    landfilter[which(grid_output$assimilated_wood_mean_gCm2 > 0)] = 1
} else { 
    # Use this option if you don't want to filter
    landfilter = array(1,dim=dim(grid_output$assimilated_wood_mean_gCm2)) 
}

# Update the land filter with the mask information
landfilter = rast(vals = t(landfilter[,dim(landfilter)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
landfilter = mask(landfilter, landmask, updatevalue = NA)
# Reconstruct back into an array
landfilter = (array(as.vector(landfilter), dim=c(dim(grid_output$assimilated_wood_mean_gCm2)[1],dim(grid_output$assimilated_wood_mean_gCm2)[2])))
landfilter = landfilter[,dim(landfilter)[2]:1]
# As a final process remove anywhere which is NaN in the actual analysis
landfilter[which(is.na(grid_output$mean_gpp_gCm2day[,,mid_quant]))] = NA

# Some variables which need masking by landmask right now
if (is.na(max(biome_names))) {
    grid_output$clusters[which(is.na(landfilter))] = NA
} else {
    tmp = dim(grid_output$clusters)
    tmp = array(values(biomes), dim=tmp)[,tmp[2]:1]
#    tmp[which(is.na(grid_output$clusters))] = NA
    grid_output$clusters = tmp
    grid_output$nos_clusters = length(biome_names)  

}

###
## Create plotting colour schemes needed

# set up colour scheme
smoothScatter_colours=colorRampPalette(c("white",rep(rev(brewer.pal(11,"Spectral")),each=3)))
colour_choices_default = colorRampPalette(brewer.pal(11,"Spectral")) 
colour_choices_sign = colorRampPalette(brewer.pal(11,"PRGn"))
colour_choices_gain = colorRampPalette(brewer.pal(9,"YlGnBu"))
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
LAIobs_unc = array(NA, dim=c(dims[1],dims[2],nos_years))
LAIcount = array(NA, dim=c(dims[1],dims[2],nos_years))
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
WoodC_lowerCI = array(NA, dim=c(dims[1], dims[2],length(PROJECT$model$timestep_days)))
WoodC_upperCI = array(NA, dim=c(dims[1], dims[2],length(PROJECT$model$timestep_days)))
wood_trend_map = array(NA, dim=c(dims[1], dims[2]))
wood_trend = rep(NA, PROJECT$nosites)
mean_wood = rep(NA, PROJECT$nosites)
# Initialise variablesfor aggregate time series
cumarea = 0
lai_grid = array(NA,dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2],nos_years))
lai_m2m2 = rep(0,nos_years) ; lai_lower_m2m2 = rep(0,nos_years) ; lai_upper_m2m2 = rep(0,nos_years)
cica_ratio = rep(0, nos_years) ; cica_lower_ratio = rep(0, nos_years) ; cica_upper_ratio = rep(0, nos_years)  
SurfWater_mm = rep(0, nos_years) ; SurfWater_lower_mm = rep(0, nos_years) ; SurfWater_upper_mm = rep(0, nos_years)  
wSWP_MPa = rep(0, nos_years) ; wSWP_lower_MPa = rep(0, nos_years) ; wSWP_upper_MPa = rep(0, nos_years)  
mean_et_PgH2Oyr = rep(0,nos_years) ; et_lower_PgH2Oyr = rep(0,nos_years) ; et_upper_PgH2Oyr = rep(0,nos_years)
mean_gpp_TgCyr = rep(0,nos_years) ; gpp_lower_TgCyr = rep(0,nos_years) ; gpp_upper_TgCyr = rep(0,nos_years)
mean_npp_TgCyr = rep(0,nos_years) ; npp_lower_TgCyr = rep(0,nos_years) ; npp_upper_TgCyr = rep(0,nos_years)
mean_rauto_TgCyr = rep(0,nos_years) ; rauto_lower_TgCyr = rep(0,nos_years) ; rauto_upper_TgCyr = rep(0,nos_years)
mean_rhet_TgCyr = rep(0,nos_years) ; rhet_lower_TgCyr = rep(0,nos_years) ; rhet_upper_TgCyr = rep(0,nos_years)
mean_nee_TgCyr = rep(0,nos_years) ; nee_lower_TgCyr = rep(0,nos_years) ; nee_upper_TgCyr = rep(0,nos_years)
mean_nbe_TgCyr = rep(0,nos_years) ; nbe_lower_TgCyr = rep(0,nos_years) ; nbe_upper_TgCyr = rep(0,nos_years)
mean_fire_TgCyr = rep(0,nos_years) ; fire_lower_TgCyr = rep(0,nos_years) ; fire_upper_TgCyr = rep(0,nos_years)
mean_harvest_TgCyr = rep(0,nos_years) ; harvest_lower_TgCyr = rep(0,nos_years) ; harvest_upper_TgCyr = rep(0,nos_years)
# Pool totals
wood_TgC = rep(0,nos_years) ; wood_lower_TgC = rep(0,nos_years) ; wood_upper_TgC = rep(0,nos_years)
litter_TgC = rep(0,nos_years) ; litter_lower_TgC = rep(0,nos_years) ; litter_upper_TgC = rep(0,nos_years)
woodlitter_TgC = rep(0,nos_years) ; woodlitter_lower_TgC = rep(0,nos_years) ; woodlitter_upper_TgC = rep(0,nos_years)
soil_TgC = rep(0,nos_years) ; soil_lower_TgC = rep(0,nos_years) ; soil_upper_TgC = rep(0,nos_years)
# Flux trends
gpp_trend = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
rauto_trend = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
rhet_trend = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
lai_trend = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
et_trend = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
# Flux normalised by first year trends
#gpp_trend_normalised = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
#rauto_trend_normalised = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
#rhet_trend_normalised = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
#lai_trend_normalised = array(NA, dim=c(dim(grid_output$mean_nee_gCm2day)[1],dim(grid_output$mean_nee_gCm2day)[2]))
# Timing variable needed
time_vector = seq(0,nos_years, length.out = dim(grid_output$nee_gCm2day)[3])

# Define counting function
counting<-function(var) {return(length(which(is.na(var) == FALSE)))}

# Loop through all sites
nos_sites_inc = 0
# Loop through every site to generate TIME VARYING estimates
for (n in seq(1, PROJECT$nosites)) {

     # Extract each sites location within the grid
     i_loc = grid_output$i_location[n] ; j_loc = grid_output$j_location[n]
     
     # Check that location has run
     if (is.na(i_loc) == FALSE & is.na(j_loc) == FALSE & is.na(landfilter[i_loc,j_loc]) == FALSE) {

         ###
         ## Aggregating data
         
         nos_sites_inc = nos_sites_inc + 1
         # Estimate pixel level trends
         tmp1 = rollapply(grid_output$gpp_gCm2day[n,mid_quant,], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)*365.25
         tmp2 = rollapply(grid_output$rauto_gCm2day[n,mid_quant,], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)*365.25
         tmp3 = rollapply(grid_output$rhet_gCm2day[n,mid_quant,], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)*365.25
         tmp4 = rollapply(grid_output$lai_m2m2[n,mid_quant,], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)
         tmp5 = rollapply(grid_output$ET_kgH2Om2day[n,mid_quant,], by = steps_per_year, width = steps_per_year, FUN = mean, na.rm=TRUE)*365.25
         #gpp_trend[i_loc,j_loc]   = coef(lm(grid_output$gpp_gCm2day[n,mid_quant,] ~ time_vector))[2]   # median selected
         #rauto_trend[i_loc,j_loc] = coef(lm(grid_output$rauto_gCm2day[n,mid_quant,] ~ time_vector))[2] # median selected
         #rhet_trend[i_loc,j_loc]  = coef(lm(grid_output$rhet_gCm2day[n,mid_quant,] ~ time_vector))[2]  # median selected
         #lai_trend[i_loc,j_loc]   = coef(lm(grid_output$lai_m2m2[n,mid_quant,] ~ time_vector))[2]      # median selected
         gpp_trend[i_loc,j_loc]   = coef(lm(tmp1 ~ c(1:nos_years)))[2] 
         rauto_trend[i_loc,j_loc] = coef(lm(tmp2 ~ c(1:nos_years)))[2] 
         rhet_trend[i_loc,j_loc]  = coef(lm(tmp3 ~ c(1:nos_years)))[2] 
         lai_trend[i_loc,j_loc]   = coef(lm(tmp4 ~ c(1:nos_years)))[2] 
         et_trend[i_loc,j_loc]    = coef(lm(tmp5 ~ c(1:nos_years)))[2] 
         # Estimate pixel level trends normalised by first year
         #gpp_trend_normalised[i_loc,j_loc]   = coef(lm((tmp1/tmp1[1]) ~ c(1:nos_years)))[2] 
         #rauto_trend_normalised[i_loc,j_loc] = coef(lm((tmp2/tmp2[1]) ~ c(1:nos_years)))[2] 
         #rhet_trend_normalised[i_loc,j_loc]  = coef(lm((tmp3/tmp3[1]) ~ c(1:nos_years)))[2] 
         #lai_trend_normalised[i_loc,j_loc]   = coef(lm((tmp4/tmp4[1]) ~ c(1:nos_years)))[2] 
         # Cumulate the total area actually used in the analysis
         cumarea = cumarea + area[i_loc,j_loc]
         lai_grid[i_loc,j_loc,] = rollapply(grid_output$lai_m2m2[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
         lai_m2m2               = lai_m2m2       + lai_grid[i_loc,j_loc,]
         lai_lower_m2m2         = lai_lower_m2m2 + rollapply(grid_output$lai_m2m2[n,low_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
         lai_upper_m2m2         = lai_upper_m2m2 + rollapply(grid_output$lai_m2m2[n,high_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
         # Area averaged states
         cica_ratio             = cica_ratio         + rollapply(grid_output$CiCa[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
         cica_lower_ratio       = cica_lower_ratio   + rollapply(grid_output$CiCa[n,low_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
         cica_upper_ratio       = cica_upper_ratio   + rollapply(grid_output$CiCa[n,high_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
         SurfWater_mm           = SurfWater_mm       + rollapply(grid_output$SurfWater_kgH2Om2[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
         SurfWater_lower_mm     = SurfWater_lower_mm + rollapply(grid_output$SurfWater_kgH2Om2[n,low_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
         SurfWater_upper_mm     = SurfWater_upper_mm + rollapply(grid_output$SurfWater_kgH2Om2[n,high_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
         wSWP_MPa               = wSWP_MPa           + rollapply(grid_output$wSWP_MPa[n,mid_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
         wSWP_lower_MPa         = wSWP_lower_MPa     + rollapply(grid_output$wSWP_MPa[n,low_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
         wSWP_upper_MPa         = wSWP_upper_MPa     + rollapply(grid_output$wSWP_MPa[n,high_quant,], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)         
         # Stocks
         wood_TgC               = wood_TgC          + rollapply(grid_output$wood_gCm2[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
         wood_lower_TgC         = wood_lower_TgC    + rollapply(grid_output$wood_gCm2[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
         wood_upper_TgC         = wood_upper_TgC    + rollapply(grid_output$wood_gCm2[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
         litter_TgC             = litter_TgC        + rollapply(grid_output$litter_gCm2[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
         litter_lower_TgC       = litter_lower_TgC  + rollapply(grid_output$litter_gCm2[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
         litter_upper_TgC       = litter_upper_TgC  + rollapply(grid_output$litter_gCm2[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
         if (length(which(names(grid_output) == "woodlitter_gCm2")) > 0) {
             woodlitter_TgC        = woodlitter_TgC       + rollapply(grid_output$woodlitter_gCm2[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
             woodlitter_lower_TgC  = woodlitter_lower_TgC + rollapply(grid_output$woodlitter_gCm2[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
             woodlitter_upper_TgC  = woodlitter_upper_TgC + rollapply(grid_output$woodlitter_gCm2[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
         }
         soil_TgC               = soil_TgC          + rollapply(grid_output$som_gCm2[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)
         soil_lower_TgC         = soil_lower_TgC    + rollapply(grid_output$som_gCm2[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
         soil_upper_TgC         = soil_upper_TgC    + rollapply(grid_output$som_gCm2[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean)         
         # Fluxes
         mean_et_PgH2Oyr        = mean_et_PgH2Oyr     + (rollapply(grid_output$ET_kgH2Om2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         et_lower_PgH2Oyr       = et_lower_PgH2Oyr    + (rollapply(grid_output$ET_kgH2Om2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         et_upper_PgH2Oyr       = et_upper_PgH2Oyr    + (rollapply(grid_output$ET_kgH2Om2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         mean_gpp_TgCyr         = mean_gpp_TgCyr      + (rollapply(grid_output$gpp_gCm2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         gpp_lower_TgCyr        = gpp_lower_TgCyr     + (rollapply(grid_output$gpp_gCm2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         gpp_upper_TgCyr        = gpp_upper_TgCyr     + (rollapply(grid_output$gpp_gCm2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         mean_npp_TgCyr         = mean_npp_TgCyr      + (rollapply(grid_output$npp_gCm2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         npp_lower_TgCyr        = npp_lower_TgCyr     + (rollapply(grid_output$npp_gCm2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         npp_upper_TgCyr        = npp_upper_TgCyr     + (rollapply(grid_output$npp_gCm2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         mean_rauto_TgCyr       = mean_rauto_TgCyr    + (rollapply(grid_output$rauto_gCm2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         rauto_lower_TgCyr      = rauto_lower_TgCyr   + (rollapply(grid_output$rauto_gCm2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         rauto_upper_TgCyr      = rauto_upper_TgCyr   + (rollapply(grid_output$rauto_gCm2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)         
         mean_rhet_TgCyr        = mean_rhet_TgCyr     + (rollapply(grid_output$rhet_gCm2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         rhet_lower_TgCyr       = rhet_lower_TgCyr    + (rollapply(grid_output$rhet_gCm2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         rhet_upper_TgCyr       = rhet_upper_TgCyr    + (rollapply(grid_output$rhet_gCm2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         mean_nee_TgCyr         = mean_nee_TgCyr      + (rollapply(grid_output$nee_gCm2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         nee_lower_TgCyr        = nee_lower_TgCyr     + (rollapply(grid_output$nee_gCm2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         nee_upper_TgCyr        = nee_upper_TgCyr     + (rollapply(grid_output$nee_gCm2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         mean_nbe_TgCyr         = mean_nbe_TgCyr      + (rollapply(grid_output$nbe_gCm2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         nbe_lower_TgCyr        = nbe_lower_TgCyr     + (rollapply(grid_output$nbe_gCm2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         nbe_upper_TgCyr        = nbe_upper_TgCyr     + (rollapply(grid_output$nbe_gCm2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         mean_fire_TgCyr        = mean_fire_TgCyr     + (rollapply(grid_output$fire_gCm2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         fire_lower_TgCyr       = fire_lower_TgCyr    + (rollapply(grid_output$fire_gCm2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         fire_upper_TgCyr       = fire_upper_TgCyr    + (rollapply(grid_output$fire_gCm2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         mean_harvest_TgCyr     = mean_harvest_TgCyr  + (rollapply(grid_output$harvest_gCm2day[n,mid_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         harvest_lower_TgCyr    = harvest_lower_TgCyr + (rollapply(grid_output$harvest_gCm2day[n,low_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         harvest_upper_TgCyr    = harvest_upper_TgCyr + (rollapply(grid_output$harvest_gCm2day[n,high_quant,]*grid_output$land_fraction[i_loc,j_loc]*area[i_loc,j_loc], width = steps_per_year, by = steps_per_year, mean) * 365.25)
         
         ###
         ## Determining Drivers and trends

         # Read in pixel driving data
         drivers = read_binary_file_format(paste(PROJECT$datapath,PROJECT$name,"_",PROJECT$sites[n],".bin",sep=""))
         # Determine forest harvest intensity
         HarvestFraction[i_loc,j_loc] = sum(drivers$met[,8]) / nos_years
         # Determine mean annual fire intensity and frequency
         BurnedFraction[i_loc,j_loc] = sum(drivers$met[,9]) / nos_years
         FireFreq[i_loc,j_loc] = length(which(drivers$met[,9] > 0)) / nos_years
         # Load any priors
         SoilCPrior[i_loc,j_loc] = drivers$parpriors[23] ; if (SoilCPrior[i_loc,j_loc] == -9999) {SoilCPrior[i_loc,j_loc] = NA}
         # Clear missing data from and extract observed LAI
         drivers$obs[which(drivers$obs[,3] == -9999),3] = NA ; drivers$obs[which(drivers$obs[,4] == -9999),4] = NA
         LAIcount[i_loc,j_loc,] = rollapply(drivers$obs[,3], width = steps_per_year, by = steps_per_year, counting)
         LAIobs[i_loc,j_loc,] = rollapply(drivers$obs[,3], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)
         LAIobs_unc[i_loc,j_loc,] = rollapply(drivers$obs[,4], width = steps_per_year, by = steps_per_year, mean, na.rm=TRUE)
         # If wood stock estimate available get that too
         tmp = which(drivers$obs[,13] > 0)
         if (length(tmp) > 0) {
             for (t in seq(1, length(tmp))) {
                  # Observational constraint
                  WoodCobs[i_loc,j_loc,tmp[t]] = drivers$obs[tmp[t],13]
                  WoodCobs_CI[i_loc,j_loc,tmp[t]] = drivers$obs[tmp[t],14]
                  # Corresponding model output
                  WoodC[i_loc,j_loc,tmp[t]] = grid_output$wood_gCm2[n,mid_quant,tmp[t]]
                  WoodC_lowerCI[i_loc,j_loc,tmp[t]] = grid_output$wood_gCm2[n,low_quant,tmp[t]]
                  WoodC_upperCI[i_loc,j_loc,tmp[t]] = grid_output$wood_gCm2[n,high_quant,tmp[t]]
             } # loop time steps with obs
             # Wood stock trends
             obs_period_start = tmp[1] ; obs_period_end = tmp[length(tmp)] ; obs_period_years = length(c(obs_period_start:obs_period_end))      
             WoodCobs_trend[n] = (coef(lm(WoodCobs[i_loc,j_loc,obs_period_start:obs_period_end] ~ c(1:obs_period_years)))[2] * 12) # *12 is month to yr adjustment
             WoodCobs_trend_map[i_loc,j_loc] = WoodCobs_trend[n]
             wood_trend[n] = (coef(lm(grid_output$wood_gCm2[n,mid_quant,obs_period_start:obs_period_end] ~ c(1:obs_period_years)))[2] * 12)
             wood_trend_map[i_loc,j_loc] = wood_trend[n]
             mean_obs_wood[n] = mean(WoodCobs[i_loc,j_loc,], na.rm=TRUE)
             WoodCobs_mean_CI[n] = mean(drivers$obs[tmp,14], na.rm=TRUE)
             mean_wood[n] = mean(grid_output$mean_wood_gCm2[i_loc,j_loc,mid_quant])
         } # we have more than zero obs
     } # Did this location run
} # Site loop

# LAI averaging
lai_m2m2 = lai_m2m2 / nos_sites_inc
lai_lower_m2m2 = lai_lower_m2m2 / nos_sites_inc
lai_upper_m2m2 = lai_upper_m2m2 / nos_sites_inc
# Area averaging
cica_ratio = cica_ratio / nos_sites_inc
cica_lower_ratio = cica_lower_ratio / nos_sites_inc
cica_upper_ratio = cica_upper_ratio / nos_sites_inc
SurfWater_mm = SurfWater_mm / nos_sites_inc
SurfWater_lower_mm = SurfWater_lower_mm / nos_sites_inc
SurfWater_upper_mm = SurfWater_upper_mm / nos_sites_inc
wSWP_MPa = wSWP_MPa / nos_sites_inc
wSWP_lower_MPa = wSWP_lower_MPa  / nos_sites_inc
wSWP_upper_MPa = wSWP_upper_MPa / nos_sites_inc
# Now adjust units gC/yr -> TgC/yr
mean_gpp_TgCyr      = mean_gpp_TgCyr * 1e-12
mean_npp_TgCyr      = mean_npp_TgCyr * 1e-12
mean_rauto_TgCyr    = mean_rauto_TgCyr * 1e-12
mean_rhet_TgCyr     = mean_rhet_TgCyr * 1e-12
mean_nee_TgCyr      = mean_nee_TgCyr * 1e-12
mean_nbe_TgCyr      = mean_nbe_TgCyr * 1e-12
mean_fire_TgCyr     = mean_fire_TgCyr * 1e-12
mean_harvest_TgCyr  = mean_harvest_TgCyr * 1e-12
litter_TgC          = litter_TgC * 1e-12
woodlitter_TgC      = woodlitter_TgC * 1e-12
wood_TgC            = wood_TgC * 1e-12
soil_TgC            = soil_TgC * 1e-12
# lower
gpp_lower_TgCyr      = gpp_lower_TgCyr * 1e-12
npp_lower_TgCyr      = npp_lower_TgCyr * 1e-12
rauto_lower_TgCyr    = rauto_lower_TgCyr * 1e-12
rhet_lower_TgCyr     = rhet_lower_TgCyr * 1e-12
nee_lower_TgCyr      = nee_lower_TgCyr * 1e-12
nbe_lower_TgCyr      = nbe_lower_TgCyr * 1e-12
fire_lower_TgCyr     = fire_lower_TgCyr * 1e-12
harvest_lower_TgCyr  = harvest_lower_TgCyr * 1e-12
litter_lower_TgC     = litter_lower_TgC * 1e-12
woodlitter_lower_TgC = woodlitter_lower_TgC * 1e-12
wood_lower_TgC       = wood_lower_TgC * 1e-12
soil_lower_TgC       = soil_lower_TgC * 1e-12
# upper
gpp_upper_TgCyr      = gpp_upper_TgCyr * 1e-12
npp_upper_TgCyr      = npp_upper_TgCyr * 1e-12
rauto_upper_TgCyr    = rauto_upper_TgCyr * 1e-12
rhet_upper_TgCyr     = rhet_upper_TgCyr * 1e-12
nee_upper_TgCyr      = nee_upper_TgCyr * 1e-12
nbe_upper_TgCyr      = nbe_upper_TgCyr * 1e-12
fire_upper_TgCyr     = fire_upper_TgCyr * 1e-12
harvest_upper_TgCyr  = harvest_upper_TgCyr * 1e-12
litter_upper_TgC     = litter_upper_TgC * 1e-12
woodlitter_upper_TgC = woodlitter_upper_TgC * 1e-12
wood_upper_TgC       = wood_upper_TgC * 1e-12
soil_upper_TgC       = soil_upper_TgC * 1e-12
# Now adjust units kgH2O -> PgH2O
mean_et_PgH2Oyr      = mean_et_PgH2Oyr * 1e-12
et_lower_PgH2Oyr     = et_lower_PgH2Oyr * 1e-12
et_upper_PgH2Oyr     = et_upper_PgH2Oyr * 1e-12

# Return some information to user
SignalNoise = length(which(abs(WoodCobs_trend*length(run_years)) > as.vector(WoodCobs_mean_CI) & abs(WoodCobs_trend*length(run_years)) > 1)) 
SignalNoise = SignalNoise / length(which(WoodCobs_mean_CI > 0 & abs(WoodCobs_trend*length(run_years)) > 1))
print(paste("Percentage of locations where observed change is greater than CI = ",round(SignalNoise*1e2, digits = 3)," %", sep=""))

###
## C-budgets by cluster (TgC/yr)
## Clustering carried out using both initial condition and process parameters
###

# Store original land filter 
landfilter_keep = landfilter
for (c in seq(1, grid_output$nos_clusters)) {
    # Add further filtering based on the cluster
    landfilter[which(is.na(grid_output$clusters) | grid_output$clusters != c)] = 0
    # Summary C budgets for output to table, NOTE the use of landfilter removes areas outside of the target area
    dims = dim(grid_output$mean_gpp_gCm2day)
    cluster_area = sum(area * landfilter, na.rm=TRUE) * 1e-4
    # Combine output into dataframe
    output = data.frame(Quantile = grid_output$num_quantiles, area_ha = rep(cluster_area, length(grid_output$num_quantiles)))
    # Ecosystem gross fluxes
    output$mean_et_PgH2Oyr                     = apply(grid_output$mean_ET_kgH2Om2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_gpp_TgCyr                      = apply(grid_output$mean_gpp_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_npp_TgCyr                      = apply(grid_output$mean_npp_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_rauto_TgCyr                    = apply(grid_output$mean_rauto_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_rhet_TgCyr                     = apply(grid_output$mean_rhet_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    if (exists("mean_rhet_litter_gCm2day", where = grid_output)) {output$mean_rhet_litter_TgCyr = apply(grid_output$mean_rhet_litter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)}
    if (exists("mean_rhet_woodlitter_gCm2day", where = grid_output)) {output$mean_rhet_woodlitter_TgCyr = apply(grid_output$mean_rhet_woodlitter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)}
    if (exists("mean_rhet_som_gCm2day", where = grid_output)) {output$mean_rhet_som_TgCyr = apply(grid_output$mean_rhet_som_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)}    
    output$mean_nbp_TgCyr                      = apply(grid_output$mean_nbp_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_nbe_TgCyr                      = apply(grid_output$mean_nbe_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_nee_TgCyr                      = apply(grid_output$mean_nee_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_fire_TgCyr                     = apply(grid_output$mean_fire_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_harvest_TgCyr                  = apply(grid_output$mean_harvest_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    # Ecosystem pools
    output$mean_biomass_TgC                    = apply(grid_output$mean_biomass_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
    output$mean_dom_TgC                        = apply(grid_output$mean_dom_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
    output$mean_labile_TgC                     = apply(grid_output$mean_labile_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
    output$mean_foliage_TgC                    = apply(grid_output$mean_foliage_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
    output$mean_roots_TgC                      = apply(grid_output$mean_roots_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
    output$mean_wood_TgC                       = apply(grid_output$mean_wood_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
    output$mean_litter_TgC                     = apply(grid_output$mean_litter_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
    if (exists(x = "mean_woodlitter_gCm2", where = grid_output)) {output$mean_woodlitter_TgC = apply(grid_output$mean_woodlitter_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)}
    output$mean_som_TgC                        = apply(grid_output$mean_som_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
    # Natural C internal allocation fluxes
    output$mean_alloc_labile_TgCyr             = apply(grid_output$mean_alloc_labile_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_combined_alloc_foliage_TgCyr   = apply(grid_output$mean_combined_alloc_foliage_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_alloc_foliage_TgCyr            = apply(grid_output$mean_alloc_foliage_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_alloc_roots_TgCyr              = apply(grid_output$mean_alloc_roots_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_alloc_wood_TgCyr               = apply(grid_output$mean_alloc_wood_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)   
    output$mean_foliage_to_litter_TgCyr        = apply(grid_output$mean_foliage_to_litter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_roots_to_litter_TgCyr          = apply(grid_output$mean_roots_to_litter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_wood_to_litter_TgCyr           = apply(grid_output$mean_wood_to_litter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_litter_to_som_TgCyr            = apply(grid_output$mean_litter_to_som_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    if (exists(x = "mean_woodlitter_to_som_gCm2day", where = grid_output)) {output$mean_woodlitter_to_som_gCm2day = apply(grid_output$mean_woodlitter_to_som_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)}
    # Fire fluxes
    output$mean_FIRElitter_labile_TgCyr        = apply(grid_output$mean_FIRElitter_labile_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_FIRElitter_foliage_TgCyr       = apply(grid_output$mean_FIRElitter_foliage_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_FIRElitter_roots_TgCyr         = apply(grid_output$mean_FIRElitter_roots_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_FIRElitter_wood_TgCyr          = apply(grid_output$mean_FIRElitter_wood_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_FIRElitter_litter_TgCyr        = apply(grid_output$mean_FIRElitter_litter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    if (exists(x = "mean_FIRElitter_woodlitter_gCm2day", where = grid_output)) {agg_output$mean_FIRElitter_woodlitter_TgCyr = apply(grid_output$mean_FIRElitter_woodlitter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)}
    output$mean_FIREemiss_labile_TgCyr         = apply(grid_output$mean_FIREemiss_labile_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_FIREemiss_foliage_TgCyr        = apply(grid_output$mean_FIREemiss_foliage_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_FIREemiss_roots_TgCyr          = apply(grid_output$mean_FIREemiss_roots_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_FIREemiss_wood_TgCyr           = apply(grid_output$mean_FIREemiss_wood_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_FIREemiss_litter_TgCyr         = apply(grid_output$mean_FIREemiss_litter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    if (exists(x = "mean_FIREemiss_woodlitter_gCm2day", where = grid_output)) {agg_output$mean_FIREemiss_woodlitter_TgCyr = apply(grid_output$mean_FIREemiss_woodlitter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)}
    output$mean_FIREemiss_som_TgCyr            = apply(grid_output$mean_FIREemiss_som_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    # Harvest fluxes
    output$mean_HARVESTlitter_labile_TgCyr     = apply(grid_output$mean_HARVESTlitter_labile_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_HARVESTlitter_foliage_TgCyr    = apply(grid_output$mean_HARVESTlitter_foliage_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_HARVESTlitter_roots_TgCyr      = apply(grid_output$mean_HARVESTlitter_roots_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_HARVESTlitter_wood_TgCyr       = apply(grid_output$mean_HARVESTlitter_wood_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_HARVESTextracted_labile_TgCyr  = apply(grid_output$mean_HARVESTextracted_labile_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_HARVESTextracted_foliage_TgCyr = apply(grid_output$mean_HARVESTextracted_foliage_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_HARVESTextracted_roots_TgCyr   = apply(grid_output$mean_HARVESTextracted_roots_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_HARVESTextracted_wood_TgCyr    = apply(grid_output$mean_HARVESTextracted_wood_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    output$mean_HARVESTextracted_litter_TgCyr  = apply(grid_output$mean_HARVESTextracted_litter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    if (exists(x = "mean_HARVESTextracted_woodlitter_gCm2day", where = grid_output)) {grid_output$mean_HARVESTextracted_woodlitter_TgCyr = apply(grid_output$mean_HARVESTextracted_woodlitter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)}
    output$mean_HARVESTextracted_som_TgCyr     = apply(grid_output$mean_HARVESTextracted_som_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
    # Mean transist / residence times
    output$labile_MTT_yr  = apply(array(landfilter,dim = dims)*grid_output$MTT_labile_years, 3, mean, na.rm=TRUE)
    output$foliage_MTT_yr = apply(array(landfilter,dim = dims)*grid_output$MTT_foliage_years, 3, mean, na.rm=TRUE)
    output$roots_MTT_yr   = apply(array(landfilter,dim = dims)*grid_output$MTT_roots_years, 3, mean, na.rm=TRUE)
    output$wood_MTT_yr    = apply(array(landfilter,dim = dims)*grid_output$MTT_wood_years, 3, mean, na.rm=TRUE)
    output$litter_MTT_yr  = apply(array(landfilter,dim = dims)*grid_output$MTT_litter_years, 3, mean, na.rm=TRUE)
    if (exists(x = "MTT_annual_woodlitter_years", where = grid_output)) {output$woodlitter_MTT_yr = apply(array(landfilter,dim = dims)*grid_output$MTT_woodlitter_years, 3, mean, na.rm=TRUE)}
    output$som_MTT_yr     = apply(array(landfilter,dim = dims)*grid_output$MTT_som_years,3,mean, na.rm=TRUE)
    # Mean stock changes per year
    output$dCbiomass_gCm2yr  = apply(array(landfilter,dim = dims)*(grid_output$final_dCbiomass_gCm2/nos_years),3,mean, na.rm=TRUE)
    output$dCdom_gCm2yr      = apply(array(landfilter,dim = dims)*(grid_output$final_dCdom_gCm2/nos_years),3,mean, na.rm=TRUE)
    output$dClabile_gCm2yr   = apply(array(landfilter,dim = dims)*(grid_output$final_dClabile_gCm2/nos_years),3,mean, na.rm=TRUE)
    output$dCfoliage_gCm2yr  = apply(array(landfilter,dim = dims)*(grid_output$final_dCfoliage_gCm2/nos_years),3,mean, na.rm=TRUE)
    output$dCroots_gCm2yr    = apply(array(landfilter,dim = dims)*(grid_output$final_dCroots_gCm2/nos_years),3,mean, na.rm=TRUE)
    output$dCwood_gCm2yr     = apply(array(landfilter,dim = dims)*(grid_output$final_dCwood_gCm2/nos_years),3,mean, na.rm=TRUE)
    output$dClitter_gCm2yr   = apply(array(landfilter,dim = dims)*(grid_output$final_dClitter_gCm2/nos_years),3,mean, na.rm=TRUE)
    if (exists(x = "final_dCwoodlitter_gCm2", where = grid_output)) {output$dCwoodlitter_gCm2yr = apply(array(landfilter,dim = dims)*(grid_output$final_dCwoodlitter_gCm2/nos_years),3,mean, na.rm=TRUE)}
    output$dCsom_gCm2yr = apply(array(landfilter,dim = dims)*(grid_output$final_dCsom_gCm2/nos_years),3,mean, na.rm=TRUE)

    # Fractional partitioning of tunover to different drivers - should they exist
    output$NaturalFractionOfTurnover_biomass = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_biomass,3,mean, na.rm=TRUE)
    output$FireFractionOfTurnover_biomass = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_biomass,3,mean, na.rm=TRUE)
    output$HarvestFractionOfTurnover_biomass = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_biomass,3,mean, na.rm=TRUE)
    output$NaturalFractionOfTurnover_dom = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_dom,3,mean, na.rm=TRUE)
    output$FireFractionOfTurnover_dom = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_dom,3,mean, na.rm=TRUE)
    output$HarvestFractionOfTurnover_dom = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_dom,3,mean, na.rm=TRUE)    
    output$NaturalFractionOfTurnover_labile = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_labile,3,mean, na.rm=TRUE)
    output$FireFractionOfTurnover_labile = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_labile,3,mean, na.rm=TRUE)
    output$HarvestFractionOfTurnover_labile = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_labile,3,mean, na.rm=TRUE)    
    output$NaturalFractionOfTurnover_foliage = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_foliage,3,mean, na.rm=TRUE)
    output$FireFractionOfTurnover_foliage = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_foliage,3,mean, na.rm=TRUE)
    output$HarvestFractionOfTurnover_foliage = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_foliage,3,mean, na.rm=TRUE)
    output$NaturalFractionOfTurnover_roots = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_roots,3,mean, na.rm=TRUE)
    output$FireFractionOfTurnover_roots = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_roots,3,mean, na.rm=TRUE)
    output$HarvestFractionOfTurnover_roots = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_roots,3,mean, na.rm=TRUE)
    output$NaturalFractionOfTurnover_wood = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_wood,3,mean, na.rm=TRUE)
    output$FireFractionOfTurnover_wood = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_wood,3,mean, na.rm=TRUE)
    output$HarvestFractionOfTurnover_wood = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_wood,3,mean, na.rm=TRUE)
    output$NaturalFractionOfTurnover_litter = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_litter,3,mean, na.rm=TRUE)
    output$FireFractionOfTurnover_litter = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_litter,3,mean, na.rm=TRUE)
    output$HarvestFractionOfTurnover_litter = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_litter,3,mean, na.rm=TRUE)
    if (exists(x = "NaturalFractionOfTurnover_woodlitter", where = grid_output)) {
        output$NaturalFractionOfTurnover_woodlitter = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_woodlitter,3,mean, na.rm=TRUE)
        output$FireFractionOfTurnover_woodlitter = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_woodlitter,3,mean, na.rm=TRUE)
        output$HarvestFractionOfTurnover_woodlitter = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_woodlitter,3,mean, na.rm=TRUE)
    } 
    output$NaturalFractionOfTurnover_som = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_som,3,mean, na.rm=TRUE)
    output$FireFractionOfTurnover_som = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_som,3,mean, na.rm=TRUE)
    output$HarvestFractionOfTurnover_som = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_som,3,mean, na.rm=TRUE)

    # Write out C budget
    write.table(output, file = paste(out_dir,"/",PROJECT$name,"_cluster_",c,"_C_budget.csv",sep=""), row.names=FALSE, sep=",",append=FALSE)
    
    # Reset land filter after each cluster to return back to the correct map
    landfilter = landfilter_keep

} # cluster loop

# Reset land filter after each cluster to return back to the correct map
landfilter = landfilter_keep

###
## C - Budget (TgC/yr)

# Summary C budgets for output to table, NOTE the use of landfilter removes areas outside of the target area
# These are not time varying (unlike the loop a few sections above)
dims = dim(grid_output$mean_gpp_gCm2day)
agg_output = data.frame(area_ha = rep(sum(grid_output$land_fraction*area * landfilter, na.rm=TRUE) * 1e-4,length(grid_output$num_quantiles))) # convertion m2->ha
# Ecosystem gross fluxes
agg_output$mean_et_PgH2Oyr  = apply(grid_output$mean_ET_kgH2Om2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_gpp_TgCyr  = apply(grid_output$mean_gpp_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_npp_TgCyr  = apply(grid_output$mean_npp_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_rauto_TgCyr= apply(grid_output$mean_rauto_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_rhet_TgCyr = apply(grid_output$mean_rhet_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
if (exists("mean_rhet_litter_gCm2day", where = grid_output)) {agg_output$mean_rhet_litter_TgCyr = apply(grid_output$mean_rhet_litter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)}
if (exists("mean_rhet_woodlitter_gCm2day", where = grid_output)) {agg_output$mean_rhet_woodlitter_TgCyr = apply(grid_output$mean_rhet_woodlitter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)}
if (exists("mean_rhet_som_gCm2day", where = grid_output)) {agg_output$mean_rhet_som_TgCyr = apply(grid_output$mean_rhet_som_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)}
agg_output$mean_nbp_TgCyr  = apply(grid_output$mean_nbp_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_nbe_TgCyr  = apply(grid_output$mean_nbe_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_nee_TgCyr  = apply(grid_output$mean_nee_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_fire_TgCyr = apply(grid_output$mean_fire_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_harvest_TgCyr  = apply(grid_output$mean_harvest_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
# Ecosystem pools
agg_output$mean_biomass_TgC= apply(grid_output$mean_biomass_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
agg_output$mean_dom_TgC= apply(grid_output$mean_dom_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
agg_output$mean_labile_TgC = apply(grid_output$mean_labile_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
agg_output$mean_foliage_TgC= apply(grid_output$mean_foliage_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
agg_output$mean_roots_TgC  = apply(grid_output$mean_roots_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
agg_output$mean_wood_TgC   = apply(grid_output$mean_wood_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
agg_output$mean_litter_TgC = apply(grid_output$mean_litter_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
if (exists(x = "mean_woodlitter_gCm2", where = grid_output)) {agg_output$mean_woodlitter_TgC = apply(grid_output$mean_woodlitter_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)}
agg_output$mean_som_TgC= apply(grid_output$mean_som_gCm2*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12,3,sum, na.rm=TRUE)
# Natural C internal allocation fluxes
agg_output$mean_alloc_labile_TgCyr = apply(grid_output$mean_alloc_labile_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_combined_alloc_foliage_TgCyr   = apply(grid_output$mean_combined_alloc_foliage_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_alloc_foliage_TgCyr= apply(grid_output$mean_alloc_foliage_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_alloc_roots_TgCyr  = apply(grid_output$mean_alloc_roots_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_alloc_wood_TgCyr   = apply(grid_output$mean_alloc_wood_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)   
agg_output$mean_foliage_to_litter_TgCyr= apply(grid_output$mean_foliage_to_litter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_roots_to_litter_TgCyr  = apply(grid_output$mean_roots_to_litter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_wood_to_litter_TgCyr   = apply(grid_output$mean_wood_to_litter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_litter_to_som_TgCyr= apply(grid_output$mean_litter_to_som_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
if (exists(x = "mean_woodlitter_to_som_gCm2day", where = grid_output)) {agg_output$mean_woodlitter_to_som_gCm2day = apply(grid_output$mean_woodlitter_to_som_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)}
# Fire fluxes
agg_output$mean_FIRElitter_labile_TgCyr= apply(grid_output$mean_FIRElitter_labile_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_FIRElitter_foliage_TgCyr   = apply(grid_output$mean_FIRElitter_foliage_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_FIRElitter_roots_TgCyr = apply(grid_output$mean_FIRElitter_roots_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_FIRElitter_wood_TgCyr  = apply(grid_output$mean_FIRElitter_wood_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_FIRElitter_litter_TgCyr= apply(grid_output$mean_FIRElitter_litter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
if (exists(x = "mean_FIRElitter_woodlitter_gCm2day", where = grid_output)) {agg_output$mean_FIRElitter_woodlitter_TgCyr = apply(grid_output$mean_FIRElitter_woodlitter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)}
agg_output$mean_FIREemiss_labile_TgCyr = apply(grid_output$mean_FIREemiss_labile_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_FIREemiss_foliage_TgCyr= apply(grid_output$mean_FIREemiss_foliage_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_FIREemiss_roots_TgCyr  = apply(grid_output$mean_FIREemiss_roots_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_FIREemiss_wood_TgCyr   = apply(grid_output$mean_FIREemiss_wood_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_FIREemiss_litter_TgCyr = apply(grid_output$mean_FIREemiss_litter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
if (exists(x = "mean_FIREemiss_woodlitter_gCm2day", where = grid_output)) {agg_output$mean_FIREemiss_woodlitter_TgCyr = apply(grid_output$mean_FIREemiss_woodlitter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)}
agg_output$mean_FIREemiss_som_TgCyr= apply(grid_output$mean_FIREemiss_som_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
# Harvest fluxes
agg_output$mean_HARVESTlitter_labile_TgCyr = apply(grid_output$mean_HARVESTlitter_labile_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_HARVESTlitter_foliage_TgCyr= apply(grid_output$mean_HARVESTlitter_foliage_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_HARVESTlitter_roots_TgCyr  = apply(grid_output$mean_HARVESTlitter_roots_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_HARVESTlitter_wood_TgCyr   = apply(grid_output$mean_HARVESTlitter_wood_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_HARVESTextracted_labile_TgCyr  = apply(grid_output$mean_HARVESTextracted_labile_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_HARVESTextracted_foliage_TgCyr = apply(grid_output$mean_HARVESTextracted_foliage_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_HARVESTextracted_roots_TgCyr   = apply(grid_output$mean_HARVESTextracted_roots_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_HARVESTextracted_wood_TgCyr= apply(grid_output$mean_HARVESTextracted_wood_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
agg_output$mean_HARVESTextracted_litter_TgCyr  = apply(grid_output$mean_HARVESTextracted_litter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
if (exists(x = "mean_HARVESTextracted_woodlitter_gCm2day", where = grid_output)) {agg_output$mean_HARVESTextracted_woodlitter_TgCyr = apply(grid_output$mean_HARVESTextracted_woodlitter_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)}
agg_output$mean_HARVESTextracted_som_TgCyr = apply(grid_output$mean_HARVESTextracted_som_gCm2day*array(grid_output$land_fraction*landfilter*area,dim = dims)*1e-12*365.25,3,sum, na.rm=TRUE)
# Mean transist / residence times
agg_output$labile_MTT_yr  = apply(array(landfilter,dim = dims)*grid_output$MTT_labile_years, 3, mean, na.rm=TRUE)
agg_output$foliage_MTT_yr = apply(array(landfilter,dim = dims)*grid_output$MTT_foliage_years, 3, mean, na.rm=TRUE)
agg_output$roots_MTT_yr   = apply(array(landfilter,dim = dims)*grid_output$MTT_roots_years, 3, mean, na.rm=TRUE)
agg_output$wood_MTT_yr    = apply(array(landfilter,dim = dims)*grid_output$MTT_wood_years, 3, mean, na.rm=TRUE)
agg_output$litter_MTT_yr  = apply(array(landfilter,dim = dims)*grid_output$MTT_litter_years, 3, mean, na.rm=TRUE)
if (exists(x = "MTT_annual_woodlitter_years", where = grid_output)) {agg_output$woodlitter_MTT_yr = apply(array(landfilter,dim = dims)*grid_output$MTT_woodlitter_years, 3, mean, na.rm=TRUE)}
agg_output$som_MTT_yr     = apply(array(landfilter,dim = dims)*grid_output$MTT_som_years,3,mean, na.rm=TRUE)
# Mean stock changes per year
agg_output$dCbiomass_gCm2yr  = apply(array(landfilter,dim = dims)*(grid_output$final_dCbiomass_gCm2/nos_years),3,mean, na.rm=TRUE)
agg_output$dCdom_gCm2yr  = apply(array(landfilter,dim = dims)*(grid_output$final_dCdom_gCm2/nos_years),3,mean, na.rm=TRUE)
agg_output$dClabile_gCm2yr   = apply(array(landfilter,dim = dims)*(grid_output$final_dClabile_gCm2/nos_years),3,mean, na.rm=TRUE)
agg_output$dCfoliage_gCm2yr  = apply(array(landfilter,dim = dims)*(grid_output$final_dCfoliage_gCm2/nos_years),3,mean, na.rm=TRUE)
agg_output$dCroots_gCm2yr= apply(array(landfilter,dim = dims)*(grid_output$final_dCroots_gCm2/nos_years),3,mean, na.rm=TRUE)
agg_output$dCwood_gCm2yr = apply(array(landfilter,dim = dims)*(grid_output$final_dCwood_gCm2/nos_years),3,mean, na.rm=TRUE)
agg_output$dClitter_gCm2yr   = apply(array(landfilter,dim = dims)*(grid_output$final_dClitter_gCm2/nos_years),3,mean, na.rm=TRUE)
if (exists(x = "final_dCwoodlitter_gCm2", where = grid_output)) {agg_output$dCwoodlitter_gCm2yr = apply(array(landfilter,dim = dims)*(grid_output$final_dCwoodlitter_gCm2/nos_years),3,mean, na.rm=TRUE)}
agg_output$dCsom_gCm2yr = apply(array(landfilter,dim = dims)*(grid_output$final_dCsom_gCm2/nos_years),3,mean, na.rm=TRUE)

# Fractional partitioning of tunover to different drivers - should they exist
agg_output$NaturalFractionOfTurnover_biomass = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_biomass,3,mean, na.rm=TRUE)
agg_output$FireFractionOfTurnover_biomass = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_biomass,3,mean, na.rm=TRUE)
agg_output$HarvestFractionOfTurnover_biomass = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_biomass,3,mean, na.rm=TRUE)
agg_output$NaturalFractionOfTurnover_dom = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_dom,3,mean, na.rm=TRUE)
agg_output$FireFractionOfTurnover_dom = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_dom,3,mean, na.rm=TRUE)
agg_output$HarvestFractionOfTurnover_dom = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_dom,3,mean, na.rm=TRUE)
agg_output$NaturalFractionOfTurnover_labile = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_labile,3,mean, na.rm=TRUE)
agg_output$FireFractionOfTurnover_labile = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_labile,3,mean, na.rm=TRUE)
agg_output$HarvestFractionOfTurnover_labile = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_labile,3,mean, na.rm=TRUE)
agg_output$NaturalFractionOfTurnover_foliage = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_foliage,3,mean, na.rm=TRUE)
agg_output$FireFractionOfTurnover_foliage = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_foliage,3,mean, na.rm=TRUE)
agg_output$HarvestFractionOfTurnover_foliage = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_foliage,3,mean, na.rm=TRUE)
agg_output$NaturalFractionOfTurnover_roots = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_roots,3,mean, na.rm=TRUE)
agg_output$FireFractionOfTurnover_roots = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_roots,3,mean, na.rm=TRUE)
agg_output$HarvestFractionOfTurnover_roots = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_roots,3,mean, na.rm=TRUE)
agg_output$NaturalFractionOfTurnover_wood = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_wood,3,mean, na.rm=TRUE)
agg_output$FireFractionOfTurnover_wood = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_wood,3,mean, na.rm=TRUE)
agg_output$HarvestFractionOfTurnover_wood = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_wood,3,mean, na.rm=TRUE)
agg_output$NaturalFractionOfTurnover_litter = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_litter,3,mean, na.rm=TRUE)
agg_output$FireFractionOfTurnover_litter = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_litter,3,mean, na.rm=TRUE)
agg_output$HarvestFractionOfTurnover_litter = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_litter,3,mean, na.rm=TRUE)
if (exists(x = "NaturalFractionOfTurnover_woodlitter", where = grid_output)) {
    agg_output$NaturalFractionOfTurnover_woodlitter = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_woodlitter,3,mean, na.rm=TRUE)
    agg_output$FireFractionOfTurnover_woodlitter = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_woodlitter,3,mean, na.rm=TRUE)
    agg_output$HarvestFractionOfTurnover_woodlitter = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_woodlitter,3,mean, na.rm=TRUE)
} 
agg_output$NaturalFractionOfTurnover_som = apply(array(landfilter,dim = dims)*grid_output$NaturalFractionOfTurnover_som,3,mean, na.rm=TRUE)
agg_output$FireFractionOfTurnover_som = apply(array(landfilter,dim = dims)*grid_output$FireFractionOfTurnover_som,3,mean, na.rm=TRUE)
agg_output$HarvestFractionOfTurnover_som = apply(array(landfilter,dim = dims)*grid_output$HarvestFractionOfTurnover_som,3,mean, na.rm=TRUE)
    
# Write out C budget
write.table(agg_output, file = paste(out_dir,"/",PROJECT$name,"_C_budget.csv",sep=""), row.names=FALSE, sep=",",append=FALSE)

###
## Plot overlayed PDFs for cluster specific and grid varying parameters and trait maps
## Clustering carried out using initial conitions and process parameters
###

### Maps of absolute harvest emission

# Then create new colours with high 'alpha', i.e. transparency
#c_colours = colorRampPalette(brewer.pal(9,"Set1"))
c_colours = colorRampPalette(brewer.pal(8,"Accent"))
c_colours = c_colours(grid_output$nos_clusters)
nbins = 30 # desired number of catagories, you might not get this many

skip_clusters = FALSE
if (skip_clusters == FALSE) {

png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_parameter_PDFs_by_cluster.png",sep=""), height = 2000, width = 3000, res = 300)
par(mfrow=c(6,6), mar = c(2,2,2,1))
# Loop parameters
for (p in seq(1, dim(grid_output$parameters)[3]-1)) {
     # Set to local variables
     tmp = as.vector(grid_output$parameters[,,p,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - (mean(b,e)*0.01) ; e = e + (mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          }
     } # loop clusters first time
     # Now plot each of them
     plot(cluster_var[1,]~x_axis, type="l", lwd=2, col = c_colours[1], main=paste("Parameter = ",p,sep=""), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax)) # Start with first cluster
     for (c in seq(2,grid_output$nos_clusters)) {
          lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) # Add next cluster
     } # loop clusters again
} # loop parameters
dev.off()

png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_C_budget_PDFs_by_cluster.png",sep=""), width = 3000, height = 1800, res = 300)
par(mfrow=c(3,4), mar=c(2,2,2,1), omi=c(0.1,0.1,0.14,0.1))

     ## GPP

     # Set to local variables
     tmp = as.vector(grid_output$mean_gpp_gCm2day[,,mid_quant]*365.25*1e-2)
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("GPP (MgC h",a^-1,y^-1,")",sep="")), 
                       xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## Rauto

     # Set to local variables
     tmp = as.vector(grid_output$mean_rauto_gCm2day[,,mid_quant]*365.25*1e-2)
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste(R[auto]," (MgC h",a^-1,y^-1,")",sep="")), 
                       xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## Rhet

     # Set to local variables
     tmp = as.vector(grid_output$mean_rhet_gCm2day[,,mid_quant]*365.25*1e-2)
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste(R[het]," (MgC h",a^-1,y^-1,")",sep="")), 
                       xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## Fire

     # Set to local variables
     tmp = as.vector(grid_output$mean_fire_gCm2day[,,mid_quant]*365.25*1e-2)
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("Fire (MgC h",a^-1,y^-1,")",sep="")), 
                       xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## labile C stocks

     # Set to local variables
     tmp = as.vector(grid_output$mean_labile_gCm2[,,mid_quant]*1e-2)
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("Mean labile (MgC h",a^-1,")",sep="")), 
                       xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## Foliage C stocks

     # Set to local variables
     tmp = as.vector(grid_output$mean_foliage_gCm2[,,mid_quant]*1e-2)
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("Mean foliage (MgC h",a^-1,")",sep="")), 
                       xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## Fine roots C stocks

     # Set to local variables
     tmp = as.vector(grid_output$mean_roots_gCm2[,,mid_quant]*1e-2)
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("Mean roots (MgC h",a^-1,")",sep="")), 
                       xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## Wood stock

     # Set to local variables
     tmp = as.vector(grid_output$mean_wood_gCm2[,,mid_quant]*1e-2)
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("Mean wood (MgC h",a^-1,")",sep="")), 
                       xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## Litter C stocks

     # Set to local variables
     tmp = as.vector(grid_output$mean_litter_gCm2[,,mid_quant]*1e-2)
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("Mean litter (MgC h",a^-1,")",sep="")), 
                       xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## Soil C stocks

     # Set to local variables
     tmp = as.vector(grid_output$mean_som_gCm2[,,mid_quant]*1e-2)
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("Mean soil (MgC h",a^-1,")",sep="")), 
                       xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
dev.off()

png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_trait_PDFs_by_cluster.png",sep=""), width = 3000, height = 1800, res = 300)
par(mfrow=c(3,4), mar=c(2,2,2,1), omi=c(0.1,0.1,0.14,0.1))

     ## CUE

     # Set to local variables
     tmp = 1-as.vector(grid_output$mean_rauto_gCm2day[,,mid_quant] / grid_output$mean_gpp_gCm2day[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("CUE (0-1)",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## NPP fraction to foliage

     # Set to local variables
     tmp = as.vector(grid_output$NPP_foliage_fraction[,,mid_quant])
     tmp[which(tmp > 1)] = NA # prevents against precision error in codes not picking up on very small fl allocations but turn into large fractional ones
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("NP",P[foliar]," (0-1)",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## NPP fraction to fine roots

     # Set to local variables
     tmp = as.vector(grid_output$NPP_roots_fraction[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("NP",P[root]," (0-1)",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## NPP fraction to wood

     # Set to local variables
     tmp = as.vector(grid_output$NPP_wood_fraction[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("NP",P[wood]," (0-1)",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## MRT foliage (years)

     # Set to local variables
     tmp = as.vector(grid_output$MTT_foliage_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("MR",T[foliar]," (y)",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## MRT fine root (years)

     # Set to local variables
     tmp = as.vector(grid_output$MTT_roots_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("MR",T[root]," (y)",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## MRT wood (years)

     # Set to local variables
     tmp = as.vector(grid_output$MTT_wood_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("MR",T[wood]," (y)",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## MRT litter (foliage + fine root)

     # Set to local variables
     tmp = as.vector(grid_output$MTT_litter_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("MR",T[litter]," (y)",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## MRT soil

     # Set to local variables
     tmp = as.vector(grid_output$MTT_som_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("MR",T[som]," (y)",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## Leaf Carbon per unit leaf Area (gC/m2)

     # Set to local variables
     tmp = as.vector(grid_output$parameters[,,17,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("LCA (gC",m^-2,")",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## Canopy Photosynthetic Efficiency (gC/m2/day)

     # Set to local variables
     tmp = as.vector(grid_output$parameters[,,11,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE) # must do this before plotting to get the correct ymax
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("Ceff (gC",m^-2,d^-1,")",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
dev.off()


png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_trait_PDFs_by_cluster_alt.png",sep=""), width = 3000, height = 1800, res = 300)
par(mfrow=c(3,4), mar=c(2,2,2,1), omi=c(0.1,0.1,0.14,0.1))

     ## CUE

     # Set to local variables
     tmp = 1-as.vector(grid_output$mean_rauto_gCm2day[,,mid_quant] / grid_output$mean_gpp_gCm2day[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("CUE (0-1)",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## NPP (gCm2day) to foliage

     # Set to local variables
     tmp = as.vector(grid_output$mean_combined_alloc_foliage_gCm2day[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("NP",P[foliar]," (gC",m^-2,d^-1,")",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## NPP fraction to fine roots

     # Set to local variables
     tmp = as.vector(grid_output$mean_alloc_roots_gCm2day[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("NP",P[root]," (gC",m^-2,d^-1,")",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## NPP fraction to wood

     # Set to local variables
     tmp = as.vector(grid_output$mean_alloc_wood_gCm2day[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("NP",P[wood]," (gC",m^-2,d^-1,")",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## MRT foliage (years)

     # Set to local variables
     tmp = as.vector(grid_output$MTT_foliage_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("MR",T[foliar]," (y)",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## MRT fine root (years)

     # Set to local variables
     tmp = as.vector(grid_output$MTT_roots_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("MR",T[root]," (y)",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## MRT wood (years)

     # Set to local variables
     tmp = as.vector(grid_output$MTT_wood_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("MR",T[wood]," (y)",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## MRT litter (foliage + fine roots)

     # Set to local variables
     tmp = as.vector(grid_output$MTT_litter_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("MR",T[litter]," (y)",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## MRT soil

     # Set to local variables
     tmp = as.vector(grid_output$MTT_som_years[,,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("MR",T[som]," (y)",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## Leaf Carbon per unit leaf Area (gC/m2)

     # Set to local variables
     tmp = as.vector(grid_output$parameters[,,17,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("LCA (gC",m^-2,")",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## Canopy Photosynthetic Efficiency (gC/m2/day)

     # Set to local variables
     tmp = as.vector(grid_output$parameters[,,11,mid_quant])
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE) # must do this before plotting to get the correct ymax
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("Ceff (gC",m^-2,d^-1,")",sep="")), xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
dev.off()

png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_abiotic_PDFs_by_cluster.png",sep=""), width = 3000, height = 1800, res = 300)
par(mfrow=c(2,3), mar=c(2,2,2,1), omi=c(0.1,0.1,0.14,0.1))

     ## Mean air temperature (C)

     # Set to local variables
     tmp = as.vector(grid_output$mean_temperature_C)
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("Air temperature (C)",sep="")), 
                       xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## Annual precipitation (kg/m2/yr)

     # Set to local variables
     tmp = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("Precipitation (mm ",y^-1,")",sep="")), 
                       xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## Vapour pressure deficit

     # Set to local variables
     tmp = as.vector(grid_output$mean_vpd_Pa)
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("VPD (Pa)",sep="")), 
                       xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## Annual fire Frequency

     # Set to local variables
     tmp = as.vector(FireFreq)
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("Fire Frequency (",y^-1,")",sep="")), 
                       xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
     ## Annual burnt fraction

     # Set to local variables
     tmp = as.vector(BurnedFraction)
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("Burnt Fraction (",y^-1,")",sep="")), 
                       xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop

     ## Harvest Fraction

     # Set to local variables
     tmp = as.vector(HarvestFraction)
     # Determine the x axis range and breakpoints
     b <- min(c(tmp), na.rm=TRUE) # Set the minimum for the breakpoints
     e <- max(c(tmp), na.rm=TRUE) # Set the maximum for the breakpoints
     b = b - abs(mean(b,e)*0.01) ; e = e + abs(mean(b,e)*0.01) # add a buffer
     ax <- pretty(c(b,e), n = nbins) # Make a neat vector for the breakpoints
     # Reset ymax for update across clusters
     ymax = 0
     # Create fresh cluster array
     cluster_var = array(NA, dim=c(grid_output$nos_clusters,length(ax)-1))
     # Loop through clusters
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              tmp1 = tmp[filter]
              # Plot the seperate histograms and store them in an object, do not save them yet
              tmp1 <- hist(tmp1, breaks = ax, plot = FALSE) # Save first histogram data
              cluster_var[c,] <- tmp1$counts / length(filter) ; x_axis = tmp1$mids
              # Now plot them together
              ymax = max(c(ymax,cluster_var[c,]), na.rm=TRUE)
          } # CARDAMOM analysis exists for this cluster / biome?
     } # loop clusters first time
     create_plot = TRUE
     for (c in seq(1, grid_output$nos_clusters)) {
          # Extact specific cluster
          filter = which(grid_output$clusters == c)
          if (length(filter) > 0) {
              if (create_plot) {
                  plot(cluster_var[c,]~x_axis, type="l", lwd=2, col = c_colours[c], main=expression(paste("Harvest Fraction (",y^-1,")",sep="")), 
                       xlab="", cex.main=1.3, cex.axis=1.2, ylab="", ylim=c(0,ymax))
                  create_plot = FALSE
              } else {
                  lines(cluster_var[c,]~x_axis, col = c_colours[c], lwd=2) 
              }
          } # does information for this plot exist
     } # cluster loop
     
dev.off()

png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_cluster_map.png",sep=""), height = 2000, width = 3000, res = 300)
par(mfrow=c(1,1), mar=c(0.01,1.5,0.3,7),omi=c(0.01,0.1,0.01,0.1))
var1 = grid_output$clusters 
var1 = rast(vals = t(var1[,dim(var1)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
plot(var1, main="", col=c_colours, xaxt = "n", yaxt = "n", mar=NA, bty = "n",
     cex.lab=2, cex.main=2.0, cex.axis = 2, axes = FALSE, type="continuous", pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
plot(landmask, add=TRUE, lwd=0.5)
if (is.na(biome_names[1])) {
    mtext(expression('Parameter based clustering'), side = 2, cex = 1.6, padj = -0.25, adj = 0.5)
} else {
    mtext(expression('Biome map'), side = 2, cex = 1.6, padj = -0.25, adj = 0.5)
}
dev.off()

# Combine some well known ecosystem traits with the cluster maps
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_LCA_CUE_wNPP_wMRT_clusters.png",sep=""), height = 2700, width = 4900, res = 300)
# Plot differences
par(mfrow=c(2,3), mar=c(0.6,0.4,2.9,7),omi=c(0.1,0.4,0.18,0.2))
var1 = grid_output$parameters[,,17,mid_quant]*landfilter
var2 = (1-(grid_output$mean_rauto_gCm2day[,,mid_quant] / grid_output$mean_gpp_gCm2day[,,mid_quant]))*landfilter
var3 = grid_output$clusters*landfilter
var4 = (grid_output$mean_alloc_wood_gCm2day[,,mid_quant])*landfilter*365.25*1e-2
var5 = grid_output$MTT_wood_years[,,mid_quant]*landfilter
var5[which(var5 > 60)] = 60 ; print("Not that MRT in plot *LCA_CUE_wNPP_wMRT_clusters is capped at 60 years")
var6 = (grid_output$mean_SurfWater_kgH2Om2[,,mid_quant])*landfilter
var1 = rast(vals = t(var1[,dim(var1)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = rast(vals = t(var2[,dim(var2)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t(var3[,dim(var3)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = rast(vals = t(var4[,dim(var4)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = rast(vals = t(var5[,dim(var5)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = rast(vals = t(var6[,dim(var6)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Correct spatial area and mask
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask) ; var4 = crop(var4, landmask) ; var5 = crop(var5, landmask) ; var6 = crop(var6, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask) ; var4 = mask(var4, landmask) ; var5 = mask(var5, landmask) ; var6 = mask(var6, landmask)
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# create axis
zrange1 = c(0,max(values(var1),na.rm=TRUE))
zrange2 = c(0.2,0.8)
zrange3 = c(1,max(values(var3), na.rm=TRUE))
zrange4 = c(0,max(values(var4),na.rm=TRUE))
zrange5 = c(0,max(values(var5),na.rm=TRUE))
zrange6 = range(values(var6),na.rm=TRUE)
plot(var1, main="",col = colour_choices_gain, range=zrange1, type="continuous", xaxt = "n", yaxt = "n", mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("LCA (gC ",m^-2,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, main="",col = colour_choices_default, range=zrange2, type="continuous", xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("CUE (0-1)",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, main="",col = c_colours, range=zrange3, xaxt = "n", type="continuous", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Parameter derived clusters",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var4, main="",col = colour_choices_gain, range=zrange4, xaxt = "n", type="continuous", yaxt = "n", mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Wood NPP (MgC h",a^-1,y^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var5, main="",col = colour_choices_gain, range=zrange5, xaxt = "n", type="continuous", yaxt = "n", mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Wood MRT (years)",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var6, main="",col = colour_choices_gain, range=zrange6, type="continuous", xaxt = "n", yaxt = "n", mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Mean surface water (0-30cm; kg ",m^-2,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

###
## Determine 1-posterior:prior ratio, i.e. how much have we learned?
###

# Extract parameter prior ranges from source code
prior_ranges = read_src_model_priors(PROJECT)

# Create ratio array
posterior_prior = array(NA, dim=c(dim(grid_output$parameters)[1:2],length(prior_ranges$parmin)))
for (n in seq(1, PROJECT$nosites)) {

     # Check that location has run
     if (is.na(grid_output$i_location[n]) == FALSE & is.na(grid_output$j_location[n]) == FALSE &
         is.na(landfilter[grid_output$i_location[n],grid_output$j_location[n]]) == FALSE) {
         for (p in seq(1,length(prior_ranges$parmin))) {
              tmp = grid_output$parameters[grid_output$i_location[n],grid_output$j_location[n],p,high_quant] 
              tmp = tmp - grid_output$parameters[grid_output$i_location[n],grid_output$j_location[n],p,low_quant] 
              posterior_prior[grid_output$i_location[n],grid_output$j_location[n],p] = tmp / (prior_ranges$parmax[p]-prior_ranges$parmin[p])
         } # Loop parameters
     } # Does site exist

} # Loop sites

# Generate some summary statistics
print("===All parameters===")
print(summary(apply(1-posterior_prior, 3, mean, na.rm=TRUE)))
# Print posterior parameter reductions per cluster / biome
for (c in seq(1, grid_output$nos_clusters)) {
     # Extact specific cluster
     filter = which(grid_output$clusters == c)
     if (length(filter) > 0) {
         print(paste("===All parameters, cluster = ",c,"===",sep=""))
         tmp = apply(1-posterior_prior, c(1,2), mean, na.rm=TRUE)
         print(summary(tmp[filter]))
     } # 
} # current cluster

# Write table of parameter specific values
write(c("Parameter No.","1-posterior:prior"), file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_posterior_prior_reductions.csv",sep=""),
      ncolumns = 2, append = FALSE, sep = ",")
for (p in seq(1, dim(posterior_prior)[3])) {
write(c(p,1-mean(posterior_prior[,,p],na.rm=TRUE)), file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_posterior_prior_reductions.csv",sep=""),
      ncolumns = 2, append = TRUE, sep = ",")
}
# Then write the same table but for each cluster
# Write table of parameter specific values
for (c in seq(1, grid_output$nos_clusters)) {
     # Write file header
     write(c("Parameter No.","1-posterior:prior"), file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_cluster_",c,"_posterior_prior_reductions.csv",sep=""),
           ncolumns = 2, append = FALSE, sep = ",")
     for (p in seq(1, dim(posterior_prior)[3])) {
     write(c(p,1-mean(posterior_prior[,,p][which(grid_output$clusters == c)],na.rm=TRUE)), file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_cluster_",c,"_posterior_prior_reductions.csv",sep=""),
           ncolumns = 2, append = TRUE, sep = ",")
     } # parameter loop
} # cluster loop


# Generate some plots

png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_posterior_prior_reductions.png",sep=""), height = 2000, width = 3000, res = 300)
par(mfrow=c(1,1), mar=c(0.01,1.5,0.3,7),omi=c(0.01,0.1,0.01,0.1))
tmp = area
var1 = apply(1-posterior_prior,c(1,2),mean,na.rm=TRUE)
var1 = rast(vals = t(var1[,dim(var1)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
plot(var1, main="", range=c(0,1), col=colour_choices_default, xaxt = "n", yaxt = "n", mar=NA, bty = "n",
     cex.lab=2, cex.main=2.0, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
plot(landmask, add=TRUE, lwd=0.5)
mtext(expression('Mean posterior reduction (0-1)'), side = 2, cex = 1.6, padj = -0.25, adj = 0.5)
dev.off()
} # skip_cluster

###
## Compare fraction of model ensemble members consistent with the assimilated observations

# Print summary information to the user for each dataset
print("===Assimilated LAI overlap (0-1)===")
print(summary(as.vector(landfilter*grid_output$lai_assim_data_overlap_fraction)))
print("===Assimilated GPP overlap (0-1)===")
print(summary(as.vector(landfilter*grid_output$gpp_assim_data_overlap_fraction)))
print("===Assimilated wood overlap (0-1)===")
print(summary(as.vector(landfilter*grid_output$wood_assim_data_overlap_fraction)))
print("===Assimilated NBE overlap (0-1)===")
print(summary(as.vector(landfilter*grid_output$nbe_assim_data_overlap_fraction)))
print("===Assimilated NEE overlap (0-1)===")
print(summary(as.vector(landfilter*grid_output$nee_assim_data_overlap_fraction)))
print("===Assimilated ET overlap (0-1)===")
print(summary(as.vector(landfilter*grid_output$et_assim_data_overlap_fraction)))
print("===Assimilate fire overlap (0-1)===")
print(summary(as.vector(landfilter*grid_output$fire_assim_data_overlap_fraction)))
# Are CARDAMOM models consistent with their assimilated observations
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_assimilated_observations_fraction_overlap.png",sep=""), height = 2700, width = 4900, res = 300)
# Plot differences
par(mfrow=c(3,3), mar=c(0.6,0.4,2.9,7),omi=c(0.1,0.4,0.18,0.2))
# Convert to raster
var1 = rast(vals = t(landfilter[,dim(area)[2]:1]*grid_output$lai_assim_data_overlap_fraction[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = rast(vals = t(landfilter[,dim(area)[2]:1]*grid_output$gpp_assim_data_overlap_fraction[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t(landfilter[,dim(area)[2]:1]*grid_output$wood_assim_data_overlap_fraction[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = rast(vals = t(landfilter[,dim(area)[2]:1]*grid_output$nbe_assim_data_overlap_fraction[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = rast(vals = t(landfilter[,dim(area)[2]:1]*grid_output$nee_assim_data_overlap_fraction[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = rast(vals = t(landfilter[,dim(area)[2]:1]*grid_output$et_assim_data_overlap_fraction[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = rast(vals = t(landfilter[,dim(area)[2]:1]*grid_output$fire_assim_data_overlap_fraction[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Correct spatial area and mask
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask) ; var4 = crop(var4, landmask) ; var5 = crop(var5, landmask) ; var6 = crop(var6, landmask) ; var7 = crop(var7, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask) ; var4 = mask(var4, landmask) ; var5 = mask(var5, landmask) ; var6 = mask(var6, landmask) ; var7 = mask(var7, landmask)
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# create axis
zrange1 = c(0,1)
zrange2 = c(0,1)
zrange3 = c(0,1)
zrange4 = c(0,1)
zrange5 = c(0,1)
zrange6 = c(0,1)
zrange7 = c(0,1)
# Begin plotting
plot(var1, main="",col = colour_choices_gain, range=zrange1, xaxt = "n", yaxt = "n", mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Assimilated LAI overlap (0-1)",sep="")), side = 3, cex = 1.8, padj = +0.05, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, main="",col = colour_choices_gain, range=zrange2, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Assimilated GPP overlap (0-1)",sep="")), side = 3, cex = 1.8, padj = +0.05, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, main="",col = (colour_choices_gain), range=zrange3, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Assimilated wood overlap (0-1)",sep="")), side = 3, cex = 1.8, padj = +0.05, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var4, main="",col = (colour_choices_gain), range=zrange4, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Assimilated NBE overlap (0-1)",sep="")), side = 3, cex = 1.8, padj = +0.05, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var5, main="",col = (colour_choices_gain), range=zrange5, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Assimilated NEE overlap (0-1)",sep="")), side = 3, cex = 1.8, padj = +0.05, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var6, main="",col = (colour_choices_gain), range=zrange6, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Assimilated ET overlap (0-1)",sep="")), side = 3, cex = 1.8, padj = +0.05, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var7, main="",col = (colour_choices_gain), range=zrange7, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Assimilated Fire overlap (0-1)",sep="")), side = 3, cex = 1.8, padj = +0.05, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

###
## Loading and processing of independent observations
###

# The overall scheme here is to load multiple datasets for each observation type (NBE, GPP, Fire)
# combine them together to provide a mean estimate and an estimate of uncertainty based on the maximum and minimum values for each case

###
## OCO2-MIPv10 ensemble mean and variance of atmospheric inversions

# check which file prefix we are using today
# list all available files which we will then search
avail_files = list.files("/exports/csce/datastore/geos/groups/gcel/AtmosphericInversions/OCO2v10_MIP/LNLG/",full.names=TRUE)
prefix_est = "EnsMean_gridded" # (.)* wildcard characters for unix standard c_gls*_
prefix_unc = "EnsStd_gridded" # (.)* wildcard characters for unix standard c_gls*_

# Only a single file is provided here, from this the correct time period
# must be extracted

# Exsure that both files exist
est_file = avail_files[grepl(prefix_est, avail_files)]
unc_file = avail_files[grepl(prefix_unc, avail_files)]
if (length(est_file) != 1 | length(unc_file) != 1) {
    print(paste("Incorrect number of NBE estimte and uncertainty files found ", sep=""))
    print(paste("est_file = ",est_file, sep=""))
    print(paste("unc_file = ",unc_file, sep=""))
    stop()
}

# Open both files
data_est = nc_open(est_file)
data_unc = nc_open(unc_file)
# Get model list
model_list = ncatt_get(data_est, varid=0, "group_names")
model_list = unlist(strsplit(model_list$value, ";"))
# Extract datetime information
time_in = ncvar_get(data_est, "start_date")
# Extract estimate and uncertainty information
nbe_in = ncvar_get(data_est, "land") # gC/m2/year
nbe_unc_in = ncvar_get(data_unc, "land") # gC/m2/year

# Extract latitude and longitude
lat_in = ncvar_get(data_est, "latitude")
long_in = ncvar_get(data_est, "longitude")
# Turn lat_in / long_in from vectors to arrays
lat_in = t(array(lat_in, dim=c(dim(nbe_in)[2],dim(nbe_in)[1])))
long_in = array(long_in, dim=c(dim(nbe_in)[1],dim(nbe_in)[2]))

# Close both files
nc_close(data_est) ; nc_close(data_unc)

# timing information on the number of day in a month
month_days = rep(31,length.out=12)
month_days[2] = 28 ; month_days[c(4,6,9,11)] = 30
# Determine time overlap between the model analysis and the NBE estimates
nbe_years = unique(as.numeric(time_in[1,]))

lat_done = FALSE ; missing_years = 0 ; keepers = 0 ; yrs = 1
# loop for year here
for (yr in seq(1, length(run_years))) {
     # Inform the user
     #print(paste("... ",round((yr/length(run_years))*100,0),"% completed ",Sys.time(),sep=""))

     # first check how many files we have
     if (yr == 1) {
         for (yrr in seq(1, length(run_years))) {
              # Check whether the desired year is in the file
              this_year = which(time_in[1,] == run_years[yrr] & time_in[2,] == 1)
              # There should be a single value
              if (length(this_year) > 0) {
                  keepers = keepers+1
              } else {
                  missing_years = append(missing_years,run_years[yrr])
              }
          } # loop through possible years
          rm(yrr)
     } # first year?

     # Check where the start point and end points are for the desired year
     year_start = which(as.numeric(time_in[1,]) == run_years[yr] & as.numeric(time_in[2,]) == 1)
     year_end = which(as.numeric(time_in[1,]) == run_years[yr] & as.numeric(time_in[2,]) == 12)

     # Assuming we have right year, begin running
     if (length(year_start) > 0) {

         # Ensure these are the first and final values to cover all time steps
         year_start = year_start[1] ; year_end = year_end[length(year_end)]

         # Now loop through the available time steps
         for (t in seq(year_start, year_end)) {

              # Determine the day of year variable for each time step
              month = time_in[2,t] ; doy_in = time_in[3,t]
              # January is correct already, so only adjust if month is >= February
              if (month > 1) {
                  doy_in = doy_in + sum(month_days[1:(month-1)])
              }

              # read the NBE observations
              var1 = nbe_in[,,t] # Land based net biome exchange of CO2 (gC/m2/yr)
              # check for error variable
              var2 = nbe_unc_in[,,t] # NBE error estimate(gC/m2/yr)
              # Convert units into gC/m2/day
              #var1 = var1 / 365.25 ; var2 = var2 / 365.25

              # Convert to a raster, assuming standad WGS84 grid
              var1 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var1))
              var1 = rast(var1, crs = ("+init=epsg:4326"), type="xyz")
              var2 = data.frame(x = as.vector(long_in), y = as.vector(lat_in), z = as.vector(var2))
              var2 = rast(var2, crs = ("+init=epsg:4326"), type="xyz")

              # Create raster with the target crs (technically this bit is not required)
              target = rast(crs = ("+init=epsg:4326"), ext = ext(var1), resolution = res(var1))
              # Check whether the target and actual analyses have the same CRS
              if (compareGeom(var1,target) == FALSE) {
                  # Resample to correct grid
                  var1 = resample(var1, target, method="ngb") ; gc()
                  var2 = resample(var2, target, method="ngb") ; gc()
              }
              # Extend the extent of the overall grid to the analysis domain
              var1 = extend(var1,cardamom_ext) ; var2 = extend(var2,cardamom_ext)
              # Trim the extent of the overall grid to the analysis domain
              var1 = crop(var1,cardamom_ext) ; var2 = crop(var2,cardamom_ext)
              var1[which(as.vector(var1) == -9999)] = NA ; var2[which(as.vector(var2) == -9999)] = NA
              # If this is a gridded analysis and the desired CARDAMOM resolution is coarser than the currently provided then aggregate here.
              # Despite creation of a cardamom_ext for a site run do not allow aggragation here as tis will damage the fine resolution datasets
              if (res(var1)[1] != res(cardamom_ext)[1] | res(var1)[2] != res(cardamom_ext)[2]) {

                  # Create raster with the target resolution
                  target = rast(crs = crs(cardamom_ext), ext = ext(cardamom_ext), resolution = res(cardamom_ext))
                  # Resample to correct grid
                  var1 = resample(var1, target, method="bilinear") ; gc() 
                  var2 = resample(var2, target, method="bilinear") ; gc() 

              } # Aggrgeate to resolution

              if (lat_done == FALSE) {
                  # extract dimension information for the grid, note the axis switching between raster and actual array
                  xdim = dim(var1)[2] ; ydim = dim(var1)[1]
                  # extract the lat / long information needed
                  long = crds(var1,df=TRUE, na.rm=FALSE)
                  lat  = long$y ; long = long$x
                  # restructure into correct orientation
                  long = array(long, dim=c(xdim,ydim))
                  lat = array(lat, dim=c(xdim,ydim))                  
              }
              # break out from the rasters into arrays which we can manipulate
              var1 = array(as.vector(unlist(var1)), dim=c(xdim,ydim))
              var2 = array(as.vector(unlist(var2)), dim=c(xdim,ydim))

              # remove additional spatial information
              if (lat_done == FALSE) {
                  # create holding arrays for the nbe information...
                  nbe_hold = as.vector(var1)
                  # ...and its uncertainty information...
                  nbe_unc_hold = as.vector(var2)
                  # ...and timing
                  doy_obs = doy_in
               } else {
                  # begin populating the various outputs
                  nbe_hold = append(nbe_hold,as.vector(var1))
                  nbe_unc_hold = append(nbe_unc_hold,as.vector(var2))
                  doy_obs = append(doy_obs,doy_in)
               }

               # update flag for lat / long load
               if (lat_done == FALSE) {lat_done = TRUE}

         } # loop through available time steps in the current year       

         # keep track of years actually ran
         yrs = yrs + 1
         # clean up allocated memeory
         gc()

     } # is there information for the current year?

} # year loop

# Sanity check for NBE
if (lat_done == FALSE) {stop('No NBE information could be found...')}

# remove initial value
missing_years = missing_years[-1]

# enforce minimum uncertainty value
nbe_unc_hold[abs(nbe_unc_hold) < 0.01*365.25] = 0.01*365.25 # for consistency with the daily approach

# return spatial structure to data
nbe_out = array(as.vector(nbe_hold), dim=c(xdim,ydim,length(doy_obs)))
nbe_unc_out = array(as.vector(nbe_unc_hold), dim=c(xdim,ydim,length(doy_obs)))
# Convert standard deviation into standard error
nbe_unc_out = nbe_unc_out / sqrt(length(model_list))

# Now apply annual averaging
obs_nbe_mean_gCm2yr = array(NA, dim=c(dim(nbe_out)[1:2],keepers))
obs_nbe_unc_gCm2yr = array(NA, dim=c(dim(nbe_out)[1:2],keepers))
diff_doy = diff(doy_obs) ; diff_doy = c(0,diff_doy) ; diff_doy = which(diff_doy < 0)
a = 1 
for (yr in seq(1, length(diff_doy))) {
     b = diff_doy[yr] - 1
     obs_nbe_mean_gCm2yr[,,yr] = apply(nbe_out[,,a:b], c(1,2), mean)
     obs_nbe_unc_gCm2yr[,,yr] = sqrt((apply((nbe_unc_out[,,a:b])**2, c(1,2), sum))/length(c(a:b)))
     a = diff_doy[yr]
}

# Flip latitude axis
obs_nbe_mean_gCm2yr = obs_nbe_mean_gCm2yr[,dim(obs_nbe_mean_gCm2yr)[2]:1,]
obs_nbe_unc_gCm2yr = obs_nbe_unc_gCm2yr[,dim(obs_nbe_unc_gCm2yr)[2]:1,]

# Extract the time step mean / min / max for each of these fluxes now
obs_nbe_min_gCm2yr = obs_nbe_mean_gCm2yr - obs_nbe_unc_gCm2yr
obs_nbe_max_gCm2yr = obs_nbe_mean_gCm2yr + obs_nbe_unc_gCm2yr
# Filter out the Inf values to NaN
obs_nbe_min_gCm2yr[is.infinite(obs_nbe_min_gCm2yr) == TRUE] = NA
obs_nbe_max_gCm2yr[is.infinite(obs_nbe_max_gCm2yr) == TRUE] = NA

# Ensure that the timeseries length is consistent between the observed variable and the model analysis
# This assumes that only the timesteps that overlap the model period have been read in the first place,
# so we should only be needing to add extra empty variable space.
tmp = intersect(run_years,nbe_years)
if (length(tmp) != length(run_years)) {
    # How many years before the observations need to be added?
    add_beginning = nbe_years[1]-run_years[1]
    # How many years after the observations
    add_afterward = run_years[length(run_years)] - nbe_years[length(nbe_years)]
    if (add_beginning > 0) {
        # Convert these into arrays of the correct shape but empty
        add_beginning = array(NA, dim=c(dim(obs_nbe_min_gCm2yr)[1:2],add_beginning))
        # Add the extra years 
        obs_nbe_mean_gCm2yr = abind(add_beginning,obs_nbe_mean_gCm2yr, along=3)
        obs_nbe_min_gCm2yr = abind(add_beginning,obs_nbe_min_gCm2yr, along=3)
        obs_nbe_max_gCm2yr = abind(add_beginning,obs_nbe_max_gCm2yr, along=3)
    } 
    if (add_afterward > 0) {
        # Convert these into arrays of the correct shape but empty
        add_afterward = array(NA, dim=c(dim(obs_nbe_min_gCm2yr)[1:2],add_afterward))
        # Add the extra years 
        obs_nbe_mean_gCm2yr = abind(obs_nbe_mean_gCm2yr,add_afterward, along=3)
        obs_nbe_min_gCm2yr = abind(obs_nbe_min_gCm2yr,add_afterward, along=3)
        obs_nbe_max_gCm2yr = abind(obs_nbe_max_gCm2yr,add_afterward, along=3)
    }
} # extra years needed

# Generate aggregate values at the domain level
obs_nbe_mean_domain_TgCyr = apply(obs_nbe_mean_gCm2yr*array(grid_output$land_fraction*area, dim=c(dim(area)[1:2],dim(obs_nbe_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
obs_nbe_min_domain_TgCyr = apply(obs_nbe_min_gCm2yr*array(grid_output$land_fraction*area, dim=c(dim(area)[1:2],dim(obs_nbe_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
obs_nbe_max_domain_TgCyr = apply(obs_nbe_max_gCm2yr*array(grid_output$land_fraction*area, dim=c(dim(area)[1:2],dim(obs_nbe_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)

###
## Extract GPP estimates from Copernicus, FLUXCOM, FLUXSATv2 & MODIS

## Read from already prepared combined maps

# Read first file to get additional information
gpp_years = c(2001:2017)
gpp_years = intersect(gpp_years,run_years)

for (t in seq(1, length(gpp_years))) {
     input = nc_open(paste("/exports/csce/datastore/geos/groups/gcel/GPP_ESTIMATES/combined_gpp/global_0.5deg_monthly/Combined_GPP_OBS_",gpp_years[t],".nc",sep=""))
     #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
     input_data = ncvar_get(input, "GPP_annual")
     # Must go in as a 3D array, so check that is the case
     if (length(dim(input_data)) == 2) {input_data = array(input_data, dim=c(dim(input_data),1))}
     # Extract data source information
     input_source = ncvar_get(input,"DataSource")
     # Extract latitude / longitude information
     input_lat = ncvar_get(input,"lat_axis") ; input_long = ncvar_get(input, "long_axis")
     # Turn lat_in / long_in from vectors to arrays
     input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
     input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
     # Check for lat / long in -90 / 90, -180 / 180 repectively
     check_long = which(input_long > 180) 
     if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
     # Begin regridding
#     input_data = regrid_gdal_func(out_dir, input_data,input_lat,input_long,cardamom_ext,landmask)
     input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
     # If this is the first year, define the output object
     if (t == 1) {
         obs_gpp_ensemble_gCm2yr = array(NA, dim=c(dim(input_data$var)[1:2],length(gpp_years),length(input_source)))
         obs_gpp_mean_gCm2yr = array(NA, dim=c(dim(input_data$var)[1:2],length(gpp_years)))
         obs_gpp_min_gCm2yr = array(NA, dim=c(dim(input_data$var)[1:2],length(gpp_years)))
         obs_gpp_max_gCm2yr = array(NA, dim=c(dim(input_data$var)[1:2],length(gpp_years)))
     }
     # Assign to output variable
     # Unit convertion (gC/m2/d -> gC/m2/yr)
     obs_gpp_mean_gCm2yr[,,t] = input_data$var * 365.25

#     # GPP min
#     input_data = ncvar_get(input, "GPP_annual_min")    
#     # Must go in as a 3D array, so check that is the case
#     if (length(dim(input_data)) == 2) {input_data = array(input_data, dim=c(dim(input_data),1))}      
#     # Begin regridding
#     input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext,landmask)
#     # Assign to output variable
#     # Unit convertion (gC/m2/d -> gC/m2/yr)
#     obs_gpp_min_gCm2yr[,,t] = input_data$var * 365.25

#     # GPP max
#     input_data = ncvar_get(input, "GPP_annual_max")
#     # Must go in as a 3D array, so check that is the case
#     if (length(dim(input_data)) == 2) {input_data = array(input_data, dim=c(dim(input_data),1))}     
#     # Begin regridding
#     input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext,landmask)
#     # Assign to output variable
#     # Unit convertion (gC/m2/d -> gC/m2/yr)
#     obs_gpp_max_gCm2yr[,,t] = input_data$var * 365.25

     # GPP ensemble
     input_data = ncvar_get(input, "GPP_annual_ensemble")
     # Must go in as a 3D array, so check that is the case
     if (length(dim(input_data)) == 2) {input_data = array(input_data, dim=c(dim(input_data),1))}     
     # Begin regridding
#     input_data = regrid_gdal_func(out_dir,input_data,input_lat,input_long,cardamom_ext,landmask)
     input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
     # Assign to output variable
     # Unit convertion (gC/m2/d -> gC/m2/yr)
     obs_gpp_ensemble_gCm2yr[,,t,] = input_data$var * 365.25
     
     # Tidy
     nc_close(input) ; rm(input_data) ; gc()
}

# Ensure the spatial orientation of the processed variable matches that of CARDAMOM
obs_gpp_ensemble_gCm2yr = obs_gpp_ensemble_gCm2yr[,dim(obs_gpp_ensemble_gCm2yr)[2]:1,,]*array(landfilter, dim=c(dim(landmask_area)[1:2],dim(obs_gpp_ensemble_gCm2yr)[3],dim(obs_gpp_ensemble_gCm2yr)[4]))
obs_gpp_mean_gCm2yr = obs_gpp_mean_gCm2yr[,dim(obs_gpp_mean_gCm2yr)[2]:1,]*array(landfilter, dim=c(dim(landmask_area)[1:2],dim(obs_gpp_mean_gCm2yr)[3]))
obs_gpp_min_gCm2yr = apply(obs_gpp_ensemble_gCm2yr,c(1,2,3),min, na.rm=TRUE)*array(landfilter, dim=c(dim(landmask_area)[1:2],dim(obs_gpp_ensemble_gCm2yr)[3]))
obs_gpp_max_gCm2yr = apply(obs_gpp_ensemble_gCm2yr,c(1,2,3),max, na.rm=TRUE)*array(landfilter, dim=c(dim(landmask_area)[1:2],dim(obs_gpp_ensemble_gCm2yr)[3]))

# Ensure that the timeseries length is consistent between the observed variable and the model analysis
# This assumes that only the timesteps that overlap the model period have been read in the first place,
# so we should only be needing to add extra empty variable space.
tmp = intersect(run_years,gpp_years)
if (length(tmp) != length(run_years)) {
    # How many years before the observations need to be added?
    nos_add_beginning = gpp_years[1]-run_years[1]
    # How many years after the observations
    nos_add_afterward = run_years[length(run_years)] - gpp_years[length(gpp_years)]
    if (nos_add_beginning > 0) {
        # Convert these into arrays of the correct shape but empty
        add_beginning = array(NA, dim=c(dim(obs_gpp_min_gCm2yr)[1:2],nos_add_beginning))
        # Add the extra years 
        obs_gpp_mean_gCm2yr = abind(add_beginning,obs_gpp_mean_gCm2yr, along=3)
        obs_gpp_min_gCm2yr = abind(add_beginning,obs_gpp_min_gCm2yr, along=3)
        obs_gpp_max_gCm2yr = abind(add_beginning,obs_gpp_max_gCm2yr, along=3)
        # Convert these into arrays of the correct shape but empty
        add_beginning = array(NA, dim=c(dim(obs_gpp_ensemble_gCm2yr)[1:2],nos_add_beginning,dim(obs_gpp_ensemble_gCm2yr)[4]))
        # Add the extra years 
        obs_gpp_ensemble_gCm2yr = abind(add_beginning,obs_gpp_ensemble_gCm2yr, along=3)
    } 
    if (nos_add_afterward > 0) {
        # Convert these into arrays of the correct shape but empty
        add_afterward = array(NA, dim=c(dim(obs_gpp_min_gCm2yr)[1:2],nos_add_afterward))
        # Add the extra years 
        obs_gpp_mean_gCm2yr = abind(obs_gpp_mean_gCm2yr,add_afterward, along=3)
        obs_gpp_min_gCm2yr = abind(obs_gpp_min_gCm2yr,add_afterward, along=3)
        obs_gpp_max_gCm2yr = abind(obs_gpp_max_gCm2yr,add_afterward, along=3)
        # Convert these into arrays of the correct shape but empty
        add_afterward = array(NA, dim=c(dim(obs_gpp_ensemble_gCm2yr)[1:2],nos_add_afterward,dim(obs_gpp_ensemble_gCm2yr)[4]))
        # Add the extra years 
        obs_gpp_ensemble_gCm2yr = abind(obs_gpp_ensemble_gCm2yr,add_afterward, along=3)
    }
} # extra years needed

# Create domain averaged values for each year and data source, note that aggregation MUST happen within product type before across products
tmp = apply(obs_gpp_ensemble_gCm2yr*array(landmask_area*grid_output$land_fraction*landfilter, dim=c(dim(landmask_area)[1:2],dim(obs_gpp_ensemble_gCm2yr)[3],dim(obs_gpp_ensemble_gCm2yr)[4]))*1e-12,c(3,4),sum, na.rm=TRUE)
# where the whole grid is zero can lead to zero being introduced - remove these
tmp[which(tmp == 0)] = NA 
# Generate aggregate values at the domain level - these must come from the raw product specific variables
obs_gpp_mean_domain_TgCyr = apply(tmp,1,mean, na.rm=TRUE)
obs_gpp_min_domain_TgCyr = apply(tmp,1,min, na.rm=TRUE)
obs_gpp_max_domain_TgCyr = apply(tmp,1,max, na.rm=TRUE)
# Check for introduced Inf values
obs_gpp_min_domain_TgCyr[which(is.infinite(obs_gpp_min_domain_TgCyr))] = NA
obs_gpp_max_domain_TgCyr[which(is.infinite(obs_gpp_max_domain_TgCyr))] = NA

###
## Independent fire emissions estimate
## Two estimates available, gfed and gfas
###

## Read from already prepared combined maps

# Read first file to get additional information
fire_years = c(2004:2016)
fire_years = intersect(fire_years,run_years)

for (t in seq(1, length(fire_years))) {
     input = nc_open(paste("/exports/csce/datastore/geos/groups/gcel/FIRE_ESTIMATES/combined_fire/global_1deg_monthly/Combined_FIRE_OBS_",fire_years[t],".nc",sep=""))
     #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
     input_data = ncvar_get(input, "Fire_annual")
     # Unit convertion (gC/m2/d -> gC/m2/yr)
     input_data = input_data * 365.25
     # Must go in as a 3D array, so check that is the case
     if (length(dim(input_data)) == 2) {input_data = array(input_data, dim=c(dim(input_data),1))}
     # Extract latitude / longitude information   
     input_lat = ncvar_get(input,"lat_axis") ; input_long = ncvar_get(input, "long_axis")
     # Turn lat_in / long_in from vectors to arrays
     input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
     input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
     # Check for lat / long in -90 / 90, -180 / 180 repectively
     check_long = which(input_long > 180) 
     if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
     # Begin regridding
     input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
     # If this is the first year, define the output object
     if (t == 1) {
         obs_fire_mean_gCm2yr = array(NA, dim=c(dim(input_data$var)[1:2],length(fire_years)))
         obs_fire_min_gCm2yr = array(NA, dim=c(dim(input_data$var)[1:2],length(fire_years)))
         obs_fire_max_gCm2yr = array(NA, dim=c(dim(input_data$var)[1:2],length(fire_years)))                
     }
     # Assign to output variable
     obs_fire_mean_gCm2yr[,,t] = input_data$var

     # Fire min
     input_data = ncvar_get(input, "Fire_annual_min")    
     # Unit convertion (gC/m2/d -> gC/m2/yr)
     input_data = input_data * 365.25
     # Must go in as a 3D array, so check that is the case
     if (length(dim(input_data)) == 2) {input_data = array(input_data, dim=c(dim(input_data),1))}      
     # Begin regridding
     input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
     # Assign to output variable
     obs_fire_min_gCm2yr[,,t] = input_data$var

     # Fire max
     input_data = ncvar_get(input, "Fire_annual_max")
     # Unit convertion (gC/m2/d -> gC/m2/yr)
     input_data = input_data * 365.25
     # Must go in as a 3D array, so check that is the case
     if (length(dim(input_data)) == 2) {input_data = array(input_data, dim=c(dim(input_data),1))}     
     # Begin regridding
     input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
     # Assign to output variable
     obs_fire_max_gCm2yr[,,t] = input_data$var

     # Tidy
     nc_close(input) ; rm(input_data) ; gc()
}

# Ensure the spatial orientation of the processed variable matches that of CARDAMOM
obs_fire_mean_gCm2yr = obs_fire_mean_gCm2yr[,dim(obs_fire_mean_gCm2yr)[2]:1,]*array(landfilter, dim=c(dim(landmask_area)[1:2],dim(obs_fire_mean_gCm2yr)[3]))
obs_fire_min_gCm2yr = obs_fire_min_gCm2yr[,dim(obs_fire_min_gCm2yr)[2]:1,]*array(landfilter, dim=c(dim(landmask_area)[1:2],dim(obs_fire_mean_gCm2yr)[3]))
obs_fire_max_gCm2yr = obs_fire_max_gCm2yr[,dim(obs_fire_max_gCm2yr)[2]:1,]*array(landfilter, dim=c(dim(landmask_area)[1:2],dim(obs_fire_mean_gCm2yr)[3]))

# Ensure that the timeseries length is consistent between the observed variable and the model analysis
# This assumes that only the timesteps that overlap the model period have been read in the first place,
# so we should only be needing to add extra empty variable space.
tmp = intersect(run_years,fire_years)
if (length(tmp) != length(run_years)) {
    # How many years before the observations need to be added?
    add_beginning = fire_years[1]-run_years[1]
    # How many years after the observations
    add_afterward = run_years[length(run_years)] - fire_years[length(fire_years)]
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
obs_fire_mean_domain_TgCyr = apply(obs_fire_mean_gCm2yr*array(grid_output$land_fraction*landmask_area, dim=c(dim(landmask_area)[1:2],dim(obs_fire_mean_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
obs_fire_min_domain_TgCyr = apply(obs_fire_min_gCm2yr*array(grid_output$land_fraction*landmask_area, dim=c(dim(landmask_area)[1:2],dim(obs_fire_min_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
obs_fire_max_domain_TgCyr = apply(obs_fire_max_gCm2yr*array(grid_output$land_fraction*landmask_area, dim=c(dim(landmask_area)[1:2],dim(obs_fire_max_gCm2yr)[3]))*1e-12,c(3),sum, na.rm=TRUE)
# where the whole grid is zero can lead to zero being introduced - remove these
obs_fire_mean_domain_TgCyr[which(obs_fire_mean_domain_TgCyr == 0)] = NA 
obs_fire_min_domain_TgCyr[which(obs_fire_min_domain_TgCyr == 0)] = NA
obs_fire_max_domain_TgCyr[which(obs_fire_max_domain_TgCyr == 0)] = NA

###
## Extract ET estimates from FLUXCOM, GLEAMv3.7b & MODIS

## Read from already prepared combined maps

# Read first file to get additional information
et_years = c(2003:2016)
et_years = intersect(et_years,run_years)

for (t in seq(1, length(et_years))) {
     input = nc_open(paste("/exports/csce/datastore/geos/groups/gcel/ET_ESTIMATES/combined_et/global_0.5deg_monthly/Combined_ET_OBS_",et_years[t],".nc",sep=""))
     #input_time = ncvar_get(input, "time") # days since 1700-12-31, 365 day years
     input_data = ncvar_get(input, "ET_annual")
     # Must go in as a 3D array, so check that is the case
     if (length(dim(input_data)) == 2) {input_data = array(input_data, dim=c(dim(input_data),1))}
     # Extract data source information
     input_source = ncvar_get(input,"DataSource")
     # Extract latitude / longitude information
     input_lat = ncvar_get(input,"lat_axis") ; input_long = ncvar_get(input, "long_axis")
     # Turn lat_in / long_in from vectors to arrays
     input_lat = t(array(input_lat, dim=c(dim(input_data)[2],dim(input_data)[1])))
     input_long = array(input_long, dim=c(dim(input_data)[1],dim(input_data)[2]))
     # Check for lat / long in -90 / 90, -180 / 180 repectively
     check_long = which(input_long > 180) 
     if (length(check_long) > 0) {input_long[check_long] = input_long[check_long] - 360}
     # Begin regridding
#     input_data = regrid_gdal_func(out_dir, input_data,input_lat,input_long,cardamom_ext,landmask)
     input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
     # If this is the first year, define the output object
     if (t == 1) {
         obs_et_ensemble_kgH2Om2yr = array(NA, dim=c(dim(input_data$var)[1:2],length(et_years),length(input_source)))
         obs_et_mean_kgH2Om2yr = array(NA, dim=c(dim(input_data$var)[1:2],length(et_years)))
         obs_et_min_kgH2Om2yr = array(NA, dim=c(dim(input_data$var)[1:2],length(et_years)))
         obs_et_max_kgH2Om2yr = array(NA, dim=c(dim(input_data$var)[1:2],length(et_years)))
     }
     # Assign to output variable
     # Unit convertion (gC/m2/d -> gC/m2/yr)
     obs_et_mean_kgH2Om2yr[,,t] = input_data$var * 365.25

#     # ET min
#     input_data = ncvar_get(input, "ET_annual_min")    
#     # Must go in as a 3D array, so check that is the case
#     if (length(dim(input_data)) == 2) {input_data = array(input_data, dim=c(dim(input_data),1))}      
#     # Begin regridding
#     input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext,landmask)
#     # Assign to output variable
#     # Unit convertion (gC/m2/d -> gC/m2/yr)
#     obs_et_min_kgH2Om2yr[,,t] = input_data$var * 365.25

#     # ET max
#     input_data = ncvar_get(input, "ET_annual_max")
#     # Must go in as a 3D array, so check that is the case
#     if (length(dim(input_data)) == 2) {input_data = array(input_data, dim=c(dim(input_data),1))}     
#     # Begin regridding
#     input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext,landmask)
#     # Assign to output variable
#     # Unit convertion (gC/m2/d -> gC/m2/yr)
#     obs_et_max_kgH2Om2yr[,,t] = input_data$var * 365.25

     # ET ensemble
     input_data = ncvar_get(input, "ET_annual_ensemble")
     # Must go in as a 3D array, so check that is the case
     if (length(dim(input_data)) == 2) {input_data = array(input_data, dim=c(dim(input_data),1))}     
     # Begin regridding
#     input_data = regrid_gdal_func(out_dir,input_data,input_lat,input_long,cardamom_ext,landmask)
     input_data = regrid_func(input_data,input_lat,input_long,cardamom_ext)
     # Assign to output variable
     # Unit convertion (gC/m2/d -> gC/m2/yr)
     obs_et_ensemble_kgH2Om2yr[,,t,] = input_data$var * 365.25
     
     # Tidy
     nc_close(input) ; rm(input_data) ; gc()
}

# Ensure the spatial orientation of the processed variable matches that of CARDAMOM
obs_et_ensemble_kgH2Om2yr = obs_et_ensemble_kgH2Om2yr[,dim(obs_et_ensemble_kgH2Om2yr)[2]:1,,]*array(landfilter, dim=c(dim(landmask_area)[1:2],dim(obs_et_ensemble_kgH2Om2yr)[3],dim(obs_et_ensemble_kgH2Om2yr)[4]))
obs_et_mean_kgH2Om2yr = obs_et_mean_kgH2Om2yr[,dim(obs_et_mean_kgH2Om2yr)[2]:1,]*array(landfilter, dim=c(dim(landmask_area)[1:2],dim(obs_et_mean_kgH2Om2yr)[3]))
obs_et_min_kgH2Om2yr = apply(obs_et_ensemble_kgH2Om2yr,c(1,2,3),min, na.rm=TRUE)*array(landfilter, dim=c(dim(landmask_area)[1:2],dim(obs_et_ensemble_kgH2Om2yr)[3]))
obs_et_max_kgH2Om2yr = apply(obs_et_ensemble_kgH2Om2yr,c(1,2,3),max, na.rm=TRUE)*array(landfilter, dim=c(dim(landmask_area)[1:2],dim(obs_et_ensemble_kgH2Om2yr)[3]))

# Ensure that the timeseries length is consistent between the observed variable and the model analysis
# This assumes that only the timesteps that overlap the model period have been read in the first place,
# so we should only be needing to add extra empty variable space.
tmp = intersect(run_years,et_years)
if (length(tmp) != length(run_years)) {
    # How many years before the observations need to be added?
    nos_add_beginning = et_years[1]-run_years[1]
    # How many years after the observations
    nos_add_afterward = run_years[length(run_years)] - et_years[length(et_years)]
    if (nos_add_beginning > 0) {
        # Convert these into arrays of the correct shape but empty
        add_beginning = array(NA, dim=c(dim(obs_et_min_kgH2Om2yr)[1:2],nos_add_beginning))
        # Add the extra years 
        obs_et_mean_kgH2Om2yr = abind(add_beginning,obs_et_mean_kgH2Om2yr, along=3)
        obs_et_min_kgH2Om2yr = abind(add_beginning,obs_et_min_kgH2Om2yr, along=3)
        obs_et_max_kgH2Om2yr = abind(add_beginning,obs_et_max_kgH2Om2yr, along=3)
        # Convert these into arrays of the correct shape but empty
        add_beginning = array(NA, dim=c(dim(obs_et_ensemble_kgH2Om2yr)[1:2],nos_add_beginning,dim(obs_et_ensemble_kgH2Om2yr)[4]))
        # Add the extra years 
        obs_et_ensemble_kgH2Om2yr = abind(add_beginning,obs_et_ensemble_kgH2Om2yr, along=3)
    } 
    if (nos_add_afterward > 0) {
        # Convert these into arrays of the correct shape but empty
        add_afterward = array(NA, dim=c(dim(obs_et_min_kgH2Om2yr)[1:2],nos_add_afterward))
        # Add the extra years 
        obs_et_mean_kgH2Om2yr = abind(obs_et_mean_kgH2Om2yr,add_afterward, along=3)
        obs_et_min_kgH2Om2yr = abind(obs_et_min_kgH2Om2yr,add_afterward, along=3)
        obs_et_max_kgH2Om2yr = abind(obs_et_max_kgH2Om2yr,add_afterward, along=3)
        # Convert these into arrays of the correct shape but empty
        add_afterward = array(NA, dim=c(dim(obs_et_ensemble_kgH2Om2yr)[1:2],nos_add_afterward,dim(obs_et_ensemble_kgH2Om2yr)[4]))
        # Add the extra years 
        obs_et_ensemble_kgH2Om2yr = abind(obs_et_ensemble_kgH2Om2yr,add_afterward, along=3)
    }
} # extra years needed

# Create domain averaged values for each year and data source, note that aggregation MUST happen within product type before across products
tmp = apply(obs_et_ensemble_kgH2Om2yr*array(landmask_area*grid_output$land_fraction*landfilter, dim=c(dim(landmask_area)[1:2],dim(obs_et_ensemble_kgH2Om2yr)[3],dim(obs_et_ensemble_kgH2Om2yr)[4]))*1e-12,c(3,4),sum, na.rm=TRUE)
# where the whole grid is zero can lead to zero being introduced - remove these
tmp[which(tmp == 0)] = NA 
# Generate aggregate values at the domain level - these must come from the raw product specific variables
obs_et_mean_domain_PgH2Oyr = apply(tmp,1,mean, na.rm=TRUE)
obs_et_min_domain_PgH2Oyr = apply(tmp,1,min, na.rm=TRUE)
obs_et_max_domain_PgH2Oyr = apply(tmp,1,max, na.rm=TRUE)
# Check for introduced Inf values
obs_et_min_domain_PgH2Oyr[which(is.infinite(obs_et_min_domain_PgH2Oyr))] = NA
obs_et_max_domain_PgH2Oyr[which(is.infinite(obs_et_max_domain_PgH2Oyr))] = NA

###
## Estimate spatial consistency of DALEC analyses and independent estimates

# Determine locations where we are consistent with independent evaluation
# DALEC consistency with range described by independent estimates
nbe_sig_latitude = 0   ; nbe_sig_longitude = 0
gpp_sig_latitude = 0   ; gpp_sig_longitude = 0
fire_sig_latitude = 0  ; fire_sig_longitude = 0
et_sig_latitude = 0    ; et_sig_longitude = 0
# Create objects for mean CARDAMOM observational overlap 
grid_output$gpp_obs_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
grid_output$nbe_obs_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
grid_output$fire_obs_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
grid_output$et_obs_overlap_fraction = array(NA, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
# Specify the overlap threshold required for consistency
sig_threshold = 0.50
nos_site_inc = 0
for (i in seq(1, PROJECT$long_dim)) {
     for (j in seq(1,PROJECT$lat_dim)) {
          if (is.na(grid_output$mean_nbe_gCm2day[i,j,high_quant]) == FALSE & is.na(landfilter[i,j]) == FALSE) {
              nos_site_inc = nos_site_inc + 1
              # Determine correct site location for extracting time period specific information
              n = which(grid_output$i_location == i & grid_output$j_location == j)
              ## Does DALEC ensemble and independent estimates overlap?

              ## Where are we consistent with ET observations
              tmp = array(NA, dim=c(length(wanted_quant),nos_years))
              # Determine mean annual flux per quantile
              for (q in seq(1, length(wanted_quant))) {
                   tmp[q,] = rollapply(grid_output$ET_kgH2Om2day[n,wanted_quant[q],], by = steps_per_year, width = steps_per_year, FUN = mean)
              }
              # scale to annual value
              tmp = tmp * 365.25 ; nobs = 0 ; npdf = 0 ; grid_output$et_obs_overlap_fraction[i,j] = 0
              # Loop through time to assess model overlap with observations
              nobs = 0 ; grid_output$et_obs_overlap_fraction[i,j] = 0
              for (t in seq(1, nos_years)) {
                   if (is.na(obs_et_mean_kgH2Om2yr[i,j,t]) == FALSE) {
                       if ((obs_et_min_kgH2Om2yr[i,j,t] - obs_et_max_kgH2Om2yr[i,j,t]) != 0 ) {
                           # Create list object containing each observations distributions
                           hist_list = list(o = c(obs_et_min_kgH2Om2yr[i,j,t],obs_et_max_kgH2Om2yr[i,j,t]),
                                            m = tmp[,t])
                           # Estimate average model ensemble within observated range
                           tmp2 = (ensemble_within_range(hist_list$o,hist_list$m))
                           if (tmp2 > 0) {npdf = npdf + 1}
                           grid_output$et_obs_overlap_fraction[i,j] = grid_output$et_obs_overlap_fraction[i,j] + tmp2
                           nobs = nobs + 1
                       } 
                   }
              } # looping for cal period
              if (nobs > 0) {
                  grid_output$et_obs_overlap_fraction[i,j] = grid_output$et_obs_overlap_fraction[i,j] / nobs
                  npdf = npdf / nobs
              } else {
                  grid_output$et_obs_overlap_fraction[i,j] = 0
              }
              # Where are we consistent with independent estimates
              if (npdf > sig_threshold) {
                  et_sig_latitude = append(et_sig_latitude,j-0.5)
                  et_sig_longitude = append(et_sig_longitude,i-0.5)
              } # ET

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
              # Where are we consistent with independent estimates
              if (npdf > sig_threshold) {
                  gpp_sig_latitude = append(gpp_sig_latitude,j-0.5)
                  gpp_sig_longitude = append(gpp_sig_longitude,i-0.5)
              } # GPP

              ## Where are we consistent with NBE observations
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
              # Where are we consistent with independent estimates
              if (npdf > sig_threshold) {
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
              # Where are we consistent with independent estimates
              if (npdf > sig_threshold) {
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
et_sig_latitude     = et_sig_latitude[-1]
et_sig_longitude    = et_sig_longitude[-1]

# Print % of pixels that are consistent
print(paste("NBE pixel consistent  =",round(100*(length(nbe_sig_latitude) / nos_site_inc), digits=3)," %",sep=" "))
print(paste("GPP pixel consistent  =",round(100*(length(gpp_sig_latitude) / nos_site_inc), digits=3)," %",sep=" "))
print(paste("Fire pixel consistent =",round(100*(length(fire_sig_latitude) / nos_site_inc), digits=3)," %",sep=" "))
print(paste("ET pixel consistent   =",round(100*(length(et_sig_latitude) / nos_site_inc), digits=3)," %",sep=" "))

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
              if ((grid_output$final_dCwood_gCm2[i,j,high_quant] > 0 & grid_output$final_dCwood_gCm2[i,j,low_quant] > 0) | 
                  (grid_output$final_dCwood_gCm2[i,j,high_quant] < 0 & grid_output$final_dCwood_gCm2[i,j,low_quant] < 0)) {
                  dCwood_sig_latitude = append(dCwood_sig_latitude,j-0.5)
                  dCwood_sig_longitude = append(dCwood_sig_longitude,i-0.5)
              }
              # Is soil confidently a source or sink?
              if ((grid_output$final_dCsom_gCm2[i,j,high_quant] > 0 & grid_output$final_dCsom_gCm2[i,j,low_quant] > 0) | 
                  (grid_output$final_dCsom_gCm2[i,j,high_quant] < 0 & grid_output$final_dCsom_gCm2[i,j,low_quant] < 0)) {
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
print(paste("NBP ~ dCwood R2 = ",round(summary(lm(as.vector(-grid_output$mean_nbe_gCm2[,,mid_quant]) ~ as.vector(grid_output$final_dCwood_gCm2[,,mid_quant])))$adj.r.squared,digits=3),sep=""))
# Statisical correlation between NBE and som change
print(paste("NBP ~ dCsom R2  = ",round(summary(lm(as.vector(-grid_output$mean_nbe_gCm2[,,mid_quant]) ~ as.vector(grid_output$final_dCsom_gCm2[,,mid_quant])))$adj.r.squared,digits=3),sep=""))

###
## Plot Observations

# Climate variables
# Assign variables
var1 = grid_output$mean_temperature_C
var2 = grid_output$mean_precipitation_kgH2Om2yr
var3 = grid_output$mean_vpd_Pa*1e-3
var4 = grid_output$mean_radiation_MJm2day
# Apply filter
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var4[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Convert to raster
var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = rast(vals = t((var4)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Trim to data area
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3) ; var4 = trim(var4) 
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# Determine ranges
zrange1 = c(-1,1)*max(abs(range(values(var1),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(values(var2),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(values(var3),na.rm=TRUE)))
zrange4 = c(0,1)*max(abs(range(values(var4),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_mean_meteorology.png",sep=""), height = 2100, width = 5000*0.55, res = 300)
par(mfrow=c(2,2), mar=c(0.5,0.5,2.8,7),omi=c(0.1,0.4,0.12,0.2))
# Mean annual median estimates
plot(var1, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=(colour_choices_loss))
mtext(expression(paste('Air temperature (C)',sep="")), side=3, cex = 1.9, padj = 0.9)     
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, range=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext(expression(paste('Precipitation (mm ',y^-1,')',sep="")), side=3, cex = 1.9, padj = 0.9)     
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, range=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste('VPD (kPa)',sep="")), side=3, cex = 1.9, padj = 0.9)     
plot(landmask, add=TRUE, lwd=0.5)
plot(var4, range=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=1.9, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste('SW Radiation (MJ ',m^-2,d^-1,')',sep="")), side=3, cex = 1.9, padj = 0.9)     
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

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
var3 = apply(LAIobs*array(landfilter,dim=dim(LAIobs)),3,mean,na.rm=TRUE)
var3_unc  = apply(LAIobs_unc*array(landfilter,dim=dim(LAIobs)),3,mean,na.rm=TRUE)
var4  = lai_m2m2
zrange = range(c(var3,var4,var3+var3_unc,var3-var3_unc), na.rm=TRUE) * c(0.8,1.2)
plotCI(y = var3, x = run_years, uiw = var3_unc, main="", cex.lab=2.4, cex.main=2, cex.axis=2.4, ylim=zrange,
      col="black", pch=16, cex=0.5, lwd=4, ylab="", xlab="")
lines(var3~run_years, col="black", lwd=3, lty = 2) 
lines(var4~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var4~run_years, col=model_colours[1], pch=16)
legend("topleft", legend = c("EO","CARDAMOM"), col = c("black",model_colours[1]), lty = c(1,2), pch=c(NA,NA), horiz = FALSE, bty = "n", cex=2.1, lwd=3, ncol = 2)
mtext(expression(paste('Year',sep="")), side = 1, cex = 2.4, padj = 1.85)
mtext(expression(paste('Analysis-wide LAI (',m^2,'/',m^2,')',sep="")), side = 2, cex = 2.4, padj = -1.05)
abline(0,1, col="grey", lwd=3)
# Now plot initial soil
plot(as.vector(1e-2*grid_output$parameters[,,23,mid_quant]) ~ as.vector(1e-2*SoilCPrior), pch=1, cex = 1.6, cex.lab=2.4, cex.axis = 2.4, cex.main=2.0, ylab="", xlab="", main="", col=model_colours[1])
mtext(expression(paste('CARDAMOM',sep="")), side = 1, cex = 2.4, padj = 1.85)
mtext(expression(paste('Initial soil C (MgC/ha)',sep="")), side = 2, cex = 2.4, padj = -1.00)
abline(0,1, col="grey", lwd=3)
dev.off()

png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_observed_modelled_wood_timeseries.png",sep=""), height = 2000, width = 3500, res = 300)
par(mfrow=c(1,1), mar=c(4.2,4.7,2.8,2),omi=c(0.01,0.01,0.01,0.01))
## Plot modelled\ median and CI timeseries with corresponding observation and uncertainty, if available
var1 = NA ; var2 = NA ; var3 = NA ; var4 = NA ; var5 = NA
# Modelled wood
var1 = wood_TgC ; var2 = wood_lower_TgC ; var3 = wood_upper_TgC
# Observed wood
var4 = apply(WoodCobs*array(grid_output$land_fraction*landfilter*area,dim=dim(WoodCobs)),3,sum,na.rm=TRUE) ; var4[which(var4 == 0)] = NA
var4 = rollapply(var4, FUN = mean, by = 12, width = 12, na.rm=TRUE)*1e-12
var5 = apply(WoodCobs_CI**2*array(grid_output$land_fraction*landfilter*area,dim=dim(WoodCobs)),3,sum,na.rm=TRUE)
var5 = sqrt(rollapply(var5, FUN = mean, by = 12, width = 12, na.rm=TRUE)*1e-12)
# Begin plotting
zrange = range(c(var1,var2,var3,var4,var5), na.rm=TRUE)*c(0.9,1.1)
plot(var1~run_years, ylim=zrange, ylab="", xlab="", pch=16, cex.lab=2.3, cex.axis = 2.2, cex=2, col="white")
mtext(expression(paste('Year',sep="")), side = 1, cex = 2.4, padj = 1.85)
mtext(expression(paste('Wood stocks (TgC)',sep="")), side = 2, cex = 2.4, padj = -1.3)
plotconfidence(t(rbind(var1,var2,var3)),run_years,2,model_colours[1])
if (length(which(is.na(var4) == FALSE)) > 0) {
    plotCI(x = run_years, y = var4, uiw = var5, main="", cex.lab=2.4, cex.main=2, cex.axis=2.4, ylim=zrange,
           col="black", lwd=4, ylab="", xlab="", add=TRUE)
}
legend("topleft", legend = c("Obs","CARDAMOM"), col = c("black",model_colours[1]), lty = c(1,1), pch=c(NA,NA), 
        horiz = FALSE, bty = "n", cex=2.1, lwd=3, ncol = 2)
dev.off()

# Domain wide NBE (yaxis) model (xaxis), include independent estimates
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_CiCa_soilwater_wSWP_timeseries_comparison_plusCI.png",sep=""), height=3800, width=2500, res=300)
par(mfrow=c(3,1),mai=c(0.3,0.65,0.3,0.2),omi=c(0.2,0.2,0.3,0.005))
# CiCa ratio
var1 = cica_ratio ; var2 = cica_lower_ratio ; var3 = cica_upper_ratio
zrange = range(c(var1,var2,var3), na.rm=TRUE)
plot(var1~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, ylab="", xlab="", lty=1)
#plotconfidence(var2,run_years,2,obs_colours[1])
lines(var1~run_years, col=model_colours[1], lwd=3, lty = 1) ; points(var1~run_years, col=model_colours[1], pch=16)
lines(var2~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var2~run_years, col=model_colours[1], pch=16)
lines(var3~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var3~run_years, col=model_colours[1], pch=16)
mtext(expression(paste("CiCa (0-1)",sep="")), side=2, padj=-2.05,cex=1.5)

# Surface water content (0-30cm; KgH2O/m2)
var1 = SurfWater_mm ; var2 = SurfWater_lower_mm ; var3 = SurfWater_upper_mm
zrange = range(c(var1,var2,var3), na.rm=TRUE)
plot(var1~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, ylab="", xlab="", lty=1)
#plotconfidence(var2,run_years,2,obs_colours[1])
lines(var1~run_years, col=model_colours[1], lwd=3, lty = 1) ; points(var1~run_years, col=model_colours[1], pch=16)
lines(var2~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var2~run_years, col=model_colours[1], pch=16)
lines(var3~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var3~run_years, col=model_colours[1], pch=16)
mtext(expression(paste("Surface Water (0-30cm; mm)",sep="")), side=2, padj=-2.05, cex=1.5)

# Plant access weighted soil water potential (MPa)
var1 = wSWP_MPa ; var2 = wSWP_lower_MPa ; var3 = wSWP_upper_MPa
zrange = range(c(var1,var2,var3), na.rm=TRUE)
plot(var1~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, ylab="", xlab="", lty=1)
#plotconfidence(var2,run_years,2,obs_colours[1])
lines(var1~run_years, col=model_colours[1], lwd=3, lty = 1) ; points(var1~run_years, col=model_colours[1], pch=16)
lines(var2~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var2~run_years, col=model_colours[1], pch=16)
lines(var3~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var3~run_years, col=model_colours[1], pch=16)
mtext("Year", side=1, padj=2.0,cex=1.6)
mtext(expression(paste("wSWP (MPa)",sep="")), side=2, padj=-2.05,cex=1.5)
dev.off()

# Domain wide NBE (yaxis) model (xaxis), include independent estimates
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_wood_litter_soil_timeseries_comparison_plusCI.png",sep=""), height=3800, width=2500, res=300)
par(mfrow=c(3,1),mai=c(0.3,0.65,0.3,0.2),omi=c(0.2,0.2,0.3,0.005))
# Now plot wood stocks
var1 = wood_TgC ; var2 = wood_lower_TgC ; var3 = wood_upper_TgC
zrange = range(c(var1,var2,var3), na.rm=TRUE)
plot(var1~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, ylab="", xlab="", lty=1)
#plotconfidence(var2,run_years,2,obs_colours[1])
lines(var1~run_years, col=model_colours[1], lwd=3, lty = 1) ; points(var1~run_years, col=model_colours[1], pch=16)
lines(var2~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var2~run_years, col=model_colours[1], pch=16)
lines(var3~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var3~run_years, col=model_colours[1], pch=16)
mtext(expression(paste("Wood (TgC)",sep="")), side=2, padj=-2.05,cex=1.5)

# Now plot foliar and fine root litter stocks
var1 = litter_TgC ; var2 = litter_lower_TgC ; var3 = litter_upper_TgC
zrange = range(c(var1,var2,var3), na.rm=TRUE)
plot(var1~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, ylab="", xlab="", lty=1)
#plotconfidence(var2,run_years,2,obs_colours[1])
lines(var1~run_years, col=model_colours[1], lwd=3, lty = 1) ; points(var1~run_years, col=model_colours[1], pch=16)
lines(var2~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var2~run_years, col=model_colours[1], pch=16)
lines(var3~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var3~run_years, col=model_colours[1], pch=16)
mtext(expression(paste("Litter (TgC)",sep="")), side=2, padj=-2.05, cex=1.5)

# Now plot soil stocks
var1 = soil_TgC ; var2 = soil_lower_TgC ; var3 = soil_upper_TgC
zrange = range(c(var1,var2,var3), na.rm=TRUE)
plot(var1~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, ylab="", xlab="", lty=1)
#plotconfidence(var2,run_years,2,obs_colours[1])
lines(var1~run_years, col=model_colours[1], lwd=3, lty = 1) ; points(var1~run_years, col=model_colours[1], pch=16)
lines(var2~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var2~run_years, col=model_colours[1], pch=16)
lines(var3~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var3~run_years, col=model_colours[1], pch=16)
mtext("Year", side=1, padj=2.0,cex=1.6)
mtext(expression(paste("Soil (TgC)",sep="")), side=2, padj=-2.05,cex=1.5)
dev.off()
                      
# Determine whether we have any observed wood trend information, i.e. do we have more than 1 wood stock
if (length(which(is.na(WoodCobs_trend) == FALSE)) > 0) {

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
    var2 = rast(vals = t(landfilter[,dim(area)[2]:1]*var2[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    var3 = abs((WoodCobs_trend_map*length(run_years)) / apply(WoodCobs_CI,c(1,2),mean,na.rm=TRUE)) # signal:uncertainty
    filter = quantile(var3, prob=c(0.025, 0.975), na.rm=TRUE) 
    var3[var3 < filter[1]] = filter[1] ; var3[var3 > filter[2]] = filter[2]
    var4 = var3 ; var4[var4 > 1] = 1 ; var4[var4 < 1] = 0
    var3 = rast(vals = t(landfilter[,dim(area)[2]:1]*var3[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    var4 = rast(vals = t(landfilter[,dim(area)[2]:1]*var4[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
    # Correct spatial area and mask
    var2 = crop(var2, landmask) ; var3 = crop(var3, landmask) ; var4 = crop(var4, landmask)
    var2 = mask(var2, landmask) ; var3 = mask(var3, landmask) ; var4 = mask(var4, landmask)
    # legend position
    ee = ext(var2) ; e = rep(NA, 4)
    e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
    e[3] = ee[3] ; e[4] = ee[4]
    # create axis
    zrange2 = c(-1,1) * max(abs(range(values(var2),na.rm=TRUE)), na.rm=TRUE)
    zrange3 = c(0,1) * max(abs(range(values(var3),na.rm=TRUE)), na.rm=TRUE)
    plot(var2, main="",col = colour_choices_sign, range=zrange2, xaxt = "n", yaxt = "n", mar=NA, bty = "n",
         cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
    mtext(expression(paste("Obs ",Delta,"wood (MgC h",a^-1,"y",r^-1,")",sep="")), side = 3, cex = 1.6, padj = +0.3, adj = 0.5)
    plot(landmask, add=TRUE, lwd=0.5)
    plot(var3, main="",col = colour_choices_CI, range=zrange3, xaxt = "n", yaxt = "n", mar=NA, bty = "n",
         cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
    mtext(expression(paste("Obs ",Delta,"wood:CI",sep="")), side = 3, cex = 1.6, padj = +0.3, adj = 0.5)
    plot(landmask, add=TRUE, lwd=0.5)
    # No ext = option in legend for catagorical label, especially when only a single value present.
    plot(var4, main="",col = colour_choices_default, range=c(0,1), xaxt = "n", yaxt = "n", mar=NA, bty = "n",
         cex.lab=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(cex=1.0))
    mtext(expression(paste("Signal:Noise Catagorical",sep="")), side = 3, cex = 1.6, padj = +0.3, adj = 0.5)
    plot(landmask, add=TRUE, lwd=0.5)
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

png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_Csom_prior.png",sep=""), height = 2000, width = 3000, res = 300)
par(mfrow=c(1,1), mar=c(0.01,1.5,0.3,7),omi=c(0.01,0.1,0.01,0.1))
tmp = area
var1 = rast(vals = t(SoilCPrior[,dim(SoilCPrior)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var1 = var1 * 1e-2 # Convert gC/m2 to MgC/ha
zrange = c(0,quantile(as.vector(var1), prob=c(0.975), na.rm=TRUE))
#zrange = c(quantile(as.vector(var1), prob=c(0,0.99), na.rm=TRUE))
var1 = clamp(var1, upper = zrange[2], values = TRUE)
print(paste("SoilCPrior in ",gsub("%","_",PROJECT$name),"_Csom_prior.png restricted to max value = ",zrange[2],sep=""))
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
plot(var1, main="", range=zrange, col=colour_choices_gain, xaxt = "n", yaxt = "n", mar=NA, bty = "n",
     cex.lab=2, cex.main=2.0, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
plot(landmask, add=TRUE, lwd=0.5)
mtext(expression('Soil C prior (MgC ha)'), side = 2, cex = 1.6, padj = -0.25, adj = 0.5)
dev.off()

png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_mean_wSWP.png",sep=""), height = 2000, width = 3000, res = 300)
par(mfrow=c(1,1), mar=c(0.01,1.5,0.3,7),omi=c(0.01,0.1,0.01,0.1))
tmp = area
var1 = rast(vals = t(grid_output$mean_wSWP_MPa[,dim(grid_output$mean_wSWP_MPa)[2]:1,mid_quant]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
zrange = c(quantile(as.vector(var1), prob=c(0.01), na.rm=TRUE),0)
var1 = clamp(var1, lower = zrange[1], values = TRUE)
print(paste("mean_wSWP_MPa in ",gsub("%","_",PROJECT$name),"_mean_wSWP.png restricted to min value = ",zrange[1],sep=""))
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
plot(var1, main="", range=zrange, col=colour_choices_gain, xaxt = "n", yaxt = "n", mar=NA, bty = "n",
     cex.lab=2, cex.main=2.0, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
plot(landmask, add=TRUE, lwd=0.5)
mtext(expression('Mean wSWP (MPa)'), side = 2, cex = 1.6, padj = -0.25, adj = 0.5)
dev.off()

###
## Independent evaluation plots

# Are CARDAMOM models consistent with the range described by CTE NBE ensemble, FC GPP ensemble and GFED / GFAS Fire products
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NBE_GPP_FIRE_no_stippling.png",sep=""), height = 1000, width = 4000, res = 300)
# Plot differences
par(mfrow=c(1,3), mar=c(0.05,0.9,1.0,7.2), omi = c(0.01,0.2,0.3,0.1))
var1 = rast(vals = t(landfilter[,dim(area)[2]:1]*365.25*1e-2*grid_output$mean_nbe_gCm2day[,dim(area)[2]:1,mid_quant]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = rast(vals = t(landfilter[,dim(area)[2]:1]*365.25*1e-2*grid_output$mean_gpp_gCm2day[,dim(area)[2]:1,mid_quant]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t(landfilter[,dim(area)[2]:1]*365.25*1e-2*grid_output$mean_fire_gCm2day[,dim(area)[2]:1,mid_quant]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Correct spatial area and mask
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask)
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3)
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# create axis
zrange1 = c(-1,1) * max(abs(range(values(var1),na.rm=TRUE)), na.rm=TRUE)
zrange2 = c(0,max(values(var2), na.rm=TRUE))
zrange3 = c(0,max(values(var3), na.rm=TRUE))
plot(var1, main="",col = rev(colour_choices_default), range=zrange1, xaxt = "n", yaxt = "n", mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("NBE (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, main="",col = colour_choices_gain, range=zrange2, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("GPP (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, main="",col = (colour_choices_loss), range=zrange3, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Fire (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Are CARDAMOM models consistent with the range described by CTE NBE ensemble, FC GPP ensemble and GFED / GFAS Fire products
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NBE_GPP_FIRE_evaluation_stippling.png",sep=""), height = 1000, width = 4000, res = 300)
# Plot differences
par(mfrow=c(1,3), mar=c(0.05,0.9,1.0,7.2), omi = c(0.01,0.2,0.3,0.1))
var1 = rast(vals = t(landfilter[,dim(area)[2]:1]*365.25*1e-2*grid_output$mean_nbe_gCm2day[,dim(area)[2]:1,mid_quant]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = rast(vals = t(landfilter[,dim(area)[2]:1]*365.25*1e-2*grid_output$mean_gpp_gCm2day[,dim(area)[2]:1,mid_quant]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t(landfilter[,dim(area)[2]:1]*365.25*1e-2*grid_output$mean_fire_gCm2day[,dim(area)[2]:1,mid_quant]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Correct spatial area and mask
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask)
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3)
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# create axis
zrange1 = c(-1,1) * max(abs(range(values(var1),na.rm=TRUE)), na.rm=TRUE)
zrange2 = c(0,max(values(var2), na.rm=TRUE))
zrange3 = c(0,max(values(var3), na.rm=TRUE))
plot(var1, main="",col = rev(colour_choices_default), range=zrange1, xaxt = "n", yaxt = "n", mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("NBE (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
points(grid_long[nbe_sig_longitude+0.5,1],grid_lat[1,nbe_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, main="",col = colour_choices_gain, range=zrange2, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("GPP (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
points(grid_long[gpp_sig_longitude+0.5,1],grid_lat[1,gpp_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, main="",col = (colour_choices_loss), range=zrange3, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Fire (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
points(grid_long[fire_sig_longitude+0.5,1],grid_lat[1,fire_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Are CARDAMOM models consistent with the range described by ET ensemble, GPP ensemble and GFED / GFAS Fire products
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_ET_GPP_FIRE_evaluation_bias_stippling.png",sep=""), height = 1000, width = 4000, res = 300)
# Plot differences
par(mfrow=c(1,3), mar=c(0.05,0.9,1,7.2), omi = c(0.01,0.2,0.3,0.1))
var1 = (365.25*grid_output$mean_ET_kgH2Om2day[,,mid_quant]) - apply(obs_et_mean_kgH2Om2yr, c(1,2), mean, na.rm=TRUE)
var1 = rast(vals = t(landfilter[,dim(area)[2]:1]*var1[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = (365.25*grid_output$mean_gpp_gCm2day[,,mid_quant]) - apply(obs_gpp_mean_gCm2yr, c(1,2), mean, na.rm=TRUE)
var2 = rast(vals = t(landfilter[,dim(area)[2]:1]*1e-2*var2[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = (365.25*grid_output$mean_fire_gCm2day[,,mid_quant]) - apply(obs_fire_mean_gCm2yr, c(1,2), mean, na.rm=TRUE)
var3 = rast(vals = t(landfilter[,dim(area)[2]:1]*1e-2*var3[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Correct spatial area and mask
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask)
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3)
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# create axis
zrange1 = c(-1,1) * max(abs(range(values(var1),na.rm=TRUE)), na.rm=TRUE)
zrange2 = c(-1,1) * max(abs(range(values(var2),na.rm=TRUE)), na.rm=TRUE)
zrange3 = c(-1,1) * max(abs(range(values(var3),na.rm=TRUE)), na.rm=TRUE)
plot(var1, main="",col = colour_choices_sign, range=zrange1, xaxt = "n", yaxt = "n", mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("ET (kgH2O ",m^-2," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
points(grid_long[et_sig_longitude+0.5,1],grid_lat[1,et_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, main="",col = colour_choices_sign, range=zrange2, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("GPP (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
points(grid_long[gpp_sig_longitude+0.5,1],grid_lat[1,gpp_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, main="",col = (colour_choices_sign), range=zrange3, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Fire (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
points(grid_long[fire_sig_longitude+0.5,1],grid_lat[1,fire_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Are CARDAMOM models consistent with the range described by CTE NBE ensemble, FC GPP ensemble and GFED / GFAS Fire products
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NBE_GPP_FIRE_evaluation_bias_stippling.png",sep=""), height = 1000, width = 4000, res = 300)
# Plot differences
par(mfrow=c(1,3), mar=c(0.05,0.9,1,7.2), omi = c(0.01,0.2,0.3,0.1))
var1 = (365.25*grid_output$mean_nbe_gCm2day[,,mid_quant]) - apply(obs_nbe_mean_gCm2yr, c(1,2), mean, na.rm=TRUE)
var1 = rast(vals = t(landfilter[,dim(area)[2]:1]*1e-2*var1[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = (365.25*grid_output$mean_gpp_gCm2day[,,mid_quant]) - apply(obs_gpp_mean_gCm2yr, c(1,2), mean, na.rm=TRUE)
var2 = rast(vals = t(landfilter[,dim(area)[2]:1]*1e-2*var2[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = (365.25*grid_output$mean_fire_gCm2day[,,mid_quant]) - apply(obs_fire_mean_gCm2yr, c(1,2), mean, na.rm=TRUE)
var3 = rast(vals = t(landfilter[,dim(area)[2]:1]*1e-2*var3[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Correct spatial area and mask
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask)
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3)
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# create axis
zrange1 = c(-1,1) * max(abs(range(values(var1),na.rm=TRUE)), na.rm=TRUE)
zrange2 = c(-1,1) * max(abs(range(values(var2),na.rm=TRUE)), na.rm=TRUE)
zrange3 = c(-1,1) * max(abs(range(values(var3),na.rm=TRUE)), na.rm=TRUE)
plot(var1, main="",col = colour_choices_sign, range=zrange1, xaxt = "n", yaxt = "n", mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("NBE (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
points(grid_long[nbe_sig_longitude+0.5,1],grid_lat[1,nbe_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, main="",col = colour_choices_sign, range=zrange2, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("GPP (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
points(grid_long[gpp_sig_longitude+0.5,1],grid_lat[1,gpp_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, main="",col = (colour_choices_sign), range=zrange3, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Fire (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
points(grid_long[fire_sig_longitude+0.5,1],grid_lat[1,fire_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Are CARDAMOM models consistent with the range described by CTE NBE ensemble, FC GPP ensemble and GFED / GFAS Fire products
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NBE_GPP_FIRE_evaluation_fraction_overlap.png",sep=""), height = 1000, width = 4000, res = 300)
# Plot differences
par(mfrow=c(1,3), mar=c(0.05,0.9,0.9,7.2), omi = c(0.01,0.2,0.3,0.1))
# C1
var1 = rast(vals = t(landfilter[,dim(area)[2]:1]*grid_output$nbe_obs_overlap_fraction[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = rast(vals = t(landfilter[,dim(area)[2]:1]*grid_output$gpp_obs_overlap_fraction[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t(landfilter[,dim(area)[2]:1]*grid_output$fire_obs_overlap_fraction[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Correct spatial area and mask
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask)
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3)
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# create axis
zrange1 = c(0,1)
zrange2 = c(0,1)
zrange3 = c(0,1)
plot(var1, main="",col = colour_choices_gain, range=zrange1, xaxt = "n", yaxt = "n", mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("NBE overlap fraction",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, main="",col = colour_choices_gain, range=zrange2, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("GPP overlap fraction",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, main="",col = (colour_choices_gain), range=zrange3, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Fire overlap fraction",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Are CARDAMOM models consistent with the range described by ET ensemble, FC GPP ensemble and GFED / GFAS Fire products
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_ET_GPP_FIRE_no_stippling.png",sep=""), height = 1000, width = 4000, res = 300)
# Plot differences
par(mfrow=c(1,3), mar=c(0.05,0.9,1.0,7.2), omi = c(0.01,0.2,0.3,0.1))
var1 = rast(vals = t(landfilter[,dim(area)[2]:1]*365.25*grid_output$mean_ET_kgH2Om2day[,dim(area)[2]:1,mid_quant]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = rast(vals = t(landfilter[,dim(area)[2]:1]*365.25*1e-2*grid_output$mean_gpp_gCm2day[,dim(area)[2]:1,mid_quant]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t(landfilter[,dim(area)[2]:1]*365.25*1e-2*grid_output$mean_fire_gCm2day[,dim(area)[2]:1,mid_quant]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Correct spatial area and mask
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask)
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3)
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# create axis
zrange1 = c(0,1) * max(abs(range(values(var1),na.rm=TRUE)), na.rm=TRUE)
zrange2 = c(0,max(values(var2), na.rm=TRUE))
zrange3 = c(0,max(values(var3), na.rm=TRUE))
plot(var1, main="",col = colour_choices_gain, range=zrange1, xaxt = "n", yaxt = "n", mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("ET (kgH2O ",m^-2," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, main="",col = colour_choices_gain, range=zrange2, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("GPP (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, main="",col = (colour_choices_loss), range=zrange3, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Fire (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Are CARDAMOM models consistent with the range described by ET ensemble, FC GPP ensemble and GFED / GFAS Fire products
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_ET_GPP_FIRE_evaluation_stippling.png",sep=""), height = 1000, width = 4000, res = 300)
# Plot differences
par(mfrow=c(1,3), mar=c(0.05,0.9,1.0,7.2), omi = c(0.01,0.2,0.3,0.1))
var1 = rast(vals = t(landfilter[,dim(area)[2]:1]*365.25*grid_output$mean_ET_kgH2Om2day[,dim(area)[2]:1,mid_quant]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = rast(vals = t(landfilter[,dim(area)[2]:1]*365.25*1e-2*grid_output$mean_gpp_gCm2day[,dim(area)[2]:1,mid_quant]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t(landfilter[,dim(area)[2]:1]*365.25*1e-2*grid_output$mean_fire_gCm2day[,dim(area)[2]:1,mid_quant]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Correct spatial area and mask
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask)
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3)
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# create axis
zrange1 = c(0,1) * max(abs(range(values(var1),na.rm=TRUE)), na.rm=TRUE)
zrange2 = c(0,max(values(var2), na.rm=TRUE))
zrange3 = c(0,max(values(var3), na.rm=TRUE))
plot(var1, main="",col = rev(colour_choices_gain), range=zrange1, xaxt = "n", yaxt = "n", mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("ET (kgH2O ",m^-2," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
points(grid_long[et_sig_longitude+0.5,1],grid_lat[1,et_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, main="",col = colour_choices_gain, range=zrange2, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("GPP (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
points(grid_long[gpp_sig_longitude+0.5,1],grid_lat[1,gpp_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, main="",col = (colour_choices_loss), range=zrange3, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Fire (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
points(grid_long[fire_sig_longitude+0.5,1],grid_lat[1,fire_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Are CARDAMOM models consistent with the range described by ET ensemble, GPP ensemble and GFED / GFAS Fire products
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_ET_GPP_FIRE_evaluation_fraction_overlap.png",sep=""), height = 1000, width = 4000, res = 300)
# Plot differences
par(mfrow=c(1,3), mar=c(0.05,0.9,0.9,7.2), omi = c(0.01,0.2,0.3,0.1))
# C1
var1 = rast(vals = t(landfilter[,dim(area)[2]:1]*grid_output$et_obs_overlap_fraction[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = rast(vals = t(landfilter[,dim(area)[2]:1]*grid_output$gpp_obs_overlap_fraction[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t(landfilter[,dim(area)[2]:1]*grid_output$fire_obs_overlap_fraction[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Correct spatial area and mask
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask)
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3)
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# create axis
zrange1 = c(0,1)
zrange2 = c(0,1)
zrange3 = c(0,1)
plot(var1, main="",col = colour_choices_gain, range=zrange1, xaxt = "n", yaxt = "n", mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("ET overlap fraction",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, main="",col = colour_choices_gain, range=zrange2, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("GPP overlap fraction",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, main="",col = (colour_choices_gain), range=zrange3, xaxt = "n", yaxt = "n",  mar=NA, bty = "n",
           cex.lab=2.6, cex.main=2.6, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
mtext(expression(paste("Fire overlap fraction",sep="")), side = 3, cex = 1.8, padj = -0.1, adj = 0.5)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Domain wide NBE (yaxis) model (xaxis), include independent estimates
model_flags=c("CARDAMOM")
obs_flags=c("OCO2-MIPv10","FC/MODIS/Copernicus/FluxSatv2","GFEDv4.1s / GFAS")
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NBE_GPP_Fire_timeseries_comparison_plusCI.png",sep=""), height=3800, width=2500, res=300)
par(mfrow=c(3,1),mai=c(0.3,0.65,0.3,0.2),omi=c(0.2,0.2,0.3,0.005))
# Now plot NBE, annual time series TgC/yr
var1  = obs_nbe_mean_domain_TgCyr
var2  = cbind(cbind(c(obs_nbe_mean_domain_TgCyr),c(obs_nbe_min_domain_TgCyr)),c(obs_nbe_max_domain_TgCyr))
var3  = mean_nbe_TgCyr ; var4  = nbe_lower_TgCyr ; var5  = nbe_upper_TgCyr
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
mtext(expression(paste("Net Biome Exchange (TgC y",r^-1,")",sep="")), side=2, padj=-1.6,cex=1.5)
#mtext("Year", side=1, padj=2.0,cex=1.6)

# Now plot GPP
var3  = cbind(cbind(c(obs_gpp_mean_domain_TgCyr),c(obs_gpp_min_domain_TgCyr)),c(obs_gpp_max_domain_TgCyr))
var4  = mean_gpp_TgCyr ; var5  = gpp_lower_TgCyr ; var6  = gpp_upper_TgCyr   
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
mtext(expression(paste("Gross Primary Productivity (TgC y",r^-1,")",sep="")), side=2, padj=-1.6, cex=1.5)

# Now plot fire
var3  = cbind(cbind(c(obs_fire_mean_domain_TgCyr),c(obs_fire_min_domain_TgCyr)),c(obs_fire_max_domain_TgCyr))
var4  = mean_fire_TgCyr  ; var5  = fire_lower_TgCyr ; var6  = fire_upper_TgCyr
zrange = range(c(var3,var4,var5,var6), na.rm=TRUE)*c(0.9,1.1)
plot(var4~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, lty=2, ylab="", xlab="")
plotconfidence(var3,run_years,2,obs_colours[3])
lines(var4~run_years, col=model_colours[1], lwd=4, lty = 1) ; points(var4~run_years, col=model_colours[1], pch=16)
lines(var5~run_years, col=model_colours[1], lwd=4, lty = 2) ; points(var5~run_years, col=model_colours[1], pch=16)
lines(var6~run_years, col=model_colours[1], lwd=4, lty = 2) ; points(var5~run_years, col=model_colours[1], pch=16)
mtext("Year", side=1, padj=2.0,cex=1.6)
mtext(expression(paste("Fire Emissions (TgC y",r^-1,")",sep="")), side=2, padj=-1.6,cex=1.5)
dev.off()

# Domain wide ET (yaxis) model (xaxis), include independent estimates
model_flags=c("CARDAMOM")
obs_flags=c("FC/GLEAMv3.7b/MODIS","FC/MODIS/Copernicus/FluxSatv2","GFEDv4.1s / GFAS")
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_ET_GPP_Fire_timeseries_comparison_plusCI.png",sep=""), height=3800, width=2500, res=300)
par(mfrow=c(3,1),mai=c(0.3,0.65,0.3,0.2),omi=c(0.2,0.2,0.3,0.005))
# Now plot ET, annual time series PgH2O/yr
var1  = obs_et_mean_domain_PgH2Oyr
var2  = cbind(cbind(c(obs_et_mean_domain_PgH2Oyr),c(obs_et_min_domain_PgH2Oyr)),c(obs_et_max_domain_PgH2Oyr))
var3  = mean_et_PgH2Oyr ; var4 = et_lower_PgH2Oyr ; var5 = et_upper_PgH2Oyr
zrange = range(c(var1,var2,var3,var4,var5), na.rm=TRUE)
zrange[2] = zrange[2] + 50
plot(var3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, ylab="", xlab="", lty=1)
plotconfidence(var2,run_years,2,obs_colours[1])
lines(var3~run_years, col=model_colours[1], lwd=3, lty = 1) ; points(var3~run_years, col=model_colours[1], pch=16)
lines(var4~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var4~run_years, col=model_colours[1], pch=16)
lines(var5~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var5~run_years, col=model_colours[1], pch=16)
abline(0,0,col="grey", lwd=2)
legend("topleft", legend = c(obs_flags,model_flags), col = c(obs_colours[1:3],model_colours), 
       lty = c(rep(1,length(obs_flags)),rep(1,length(model_flags))), pch=rep(NA,length(c(obs_flags,model_flags))), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)
mtext(expression(paste("Evapotranspiration (PgH2O y",r^-1,")",sep="")), side=2, padj=-1.6,cex=1.5)
#mtext("Year", side=1, padj=2.0,cex=1.6)

# Now plot GPP
var3  = cbind(cbind(c(obs_gpp_mean_domain_TgCyr),c(obs_gpp_min_domain_TgCyr)),c(obs_gpp_max_domain_TgCyr))
var4  = mean_gpp_TgCyr ; var5  = gpp_lower_TgCyr ; var6  = gpp_upper_TgCyr   
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
mtext(expression(paste("Gross Primary Productivity (TgC y",r^-1,")",sep="")), side=2, padj=-1.6, cex=1.5)

# Now plot fire
var3  = cbind(cbind(c(obs_fire_mean_domain_TgCyr),c(obs_fire_min_domain_TgCyr)),c(obs_fire_max_domain_TgCyr))
var4  = mean_fire_TgCyr  ; var5  = fire_lower_TgCyr ; var6  = fire_upper_TgCyr
zrange = range(c(var3,var4,var5,var6), na.rm=TRUE)*c(0.9,1.1)
plot(var4~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, lty=2, ylab="", xlab="")
plotconfidence(var3,run_years,2,obs_colours[3])
lines(var4~run_years, col=model_colours[1], lwd=4, lty = 1) ; points(var4~run_years, col=model_colours[1], pch=16)
lines(var5~run_years, col=model_colours[1], lwd=4, lty = 2) ; points(var5~run_years, col=model_colours[1], pch=16)
lines(var6~run_years, col=model_colours[1], lwd=4, lty = 2) ; points(var5~run_years, col=model_colours[1], pch=16)
mtext("Year", side=1, padj=2.0,cex=1.6)
mtext(expression(paste("Fire Emissions (TgC y",r^-1,")",sep="")), side=2, padj=-1.6,cex=1.5)
dev.off()

# Domain wide ET (yaxis) model (xaxis), include independent estimates
model_flags=c("CARDAMOM")
obs_flags=c("OCO2-MIPv10","FC/GLEAMv3.7b/MODIS","FC/MODIS/Copernicus/FluxSatv2","GFEDv4.1s / GFAS")
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NBE_ET_GPP_Fire_timeseries_comparison_plusCI.png",sep=""), height=3800, width=2500, res=300)
par(mfrow=c(3,1),mai=c(0.3,0.65,0.3,0.2),omi=c(0.2,0.2,0.3,0.005))
# Now plot NBE, annual time series TgC/yr
var1  = obs_nbe_mean_domain_TgCyr
var2  = cbind(cbind(c(obs_nbe_mean_domain_TgCyr),c(obs_nbe_min_domain_TgCyr)),c(obs_nbe_max_domain_TgCyr))
var3  = mean_nbe_TgCyr ; var4  = nbe_lower_TgCyr ; var5  = nbe_upper_TgCyr
zrange = range(c(var1,var2,var3,var4,var5), na.rm=TRUE)
zrange[2] = zrange[2] + 500
plot(var3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, ylab="", xlab="", lty=1)
plotconfidence(var2,run_years,2,obs_colours[1])
lines(var3~run_years, col=model_colours[1], lwd=3, lty = 1) ; points(var3~run_years, col=model_colours[1], pch=16)
lines(var4~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var4~run_years, col=model_colours[1], pch=16)
lines(var5~run_years, col=model_colours[1], lwd=3, lty = 2) ; points(var5~run_years, col=model_colours[1], pch=16)
abline(0,0,col="grey", lwd=2)
legend("topleft", legend = c(obs_flags,model_flags), col = c(obs_colours[1:4],model_colours), 
       lty = c(rep(1,length(obs_flags)),rep(1,length(model_flags))), pch=rep(NA,length(c(obs_flags,model_flags))), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)
mtext(expression(paste("Net Biome Exchange (TgC y",r^-1,")",sep="")), side=2, padj=-1.6,cex=1.5)
#mtext("Year", side=1, padj=2.0,cex=1.6)

# Now plot ET, annual time series PgH2O/yr
var1  = obs_et_mean_domain_PgH2Oyr
var2  = cbind(cbind(c(obs_et_mean_domain_PgH2Oyr),c(obs_et_min_domain_PgH2Oyr)),c(obs_et_max_domain_PgH2Oyr))
var3  = mean_et_PgH2Oyr ; var4 = et_lower_PgH2Oyr ; var5 = et_upper_PgH2Oyr
zrange = range(c(var1,var2,var3,var4,var5), na.rm=TRUE)
zrange[2] = zrange[2] + 50
plot(var3~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, ylab="", xlab="", lty=1)
plotconfidence(var2,run_years,2,obs_colours[1])
lines(var3~run_years, col=model_colours[2], lwd=3, lty = 1) ; points(var3~run_years, col=model_colours[2], pch=16)
lines(var4~run_years, col=model_colours[2], lwd=3, lty = 2) ; points(var4~run_years, col=model_colours[2], pch=16)
lines(var5~run_years, col=model_colours[2], lwd=3, lty = 2) ; points(var5~run_years, col=model_colours[2], pch=16)
abline(0,0,col="grey", lwd=2)
legend("topleft", legend = c(obs_flags,model_flags), col = c(obs_colours[1:3],model_colours), 
       lty = c(rep(1,length(obs_flags)),rep(1,length(model_flags))), pch=rep(NA,length(c(obs_flags,model_flags))), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)
mtext(expression(paste("Evapotranspiration (PgH2O y",r^-1,")",sep="")), side=2, padj=-1.6,cex=1.5)
#mtext("Year", side=1, padj=2.0,cex=1.6)

# Now plot GPP
var3  = cbind(cbind(c(obs_gpp_mean_domain_TgCyr),c(obs_gpp_min_domain_TgCyr)),c(obs_gpp_max_domain_TgCyr))
var4  = mean_gpp_TgCyr ; var5  = gpp_lower_TgCyr ; var6  = gpp_upper_TgCyr   
zrange = range(c(var3,var4,var5,var6), na.rm=TRUE)*c(0.9,1.0)
plot(var4~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd = 4, ylab="", xlab="", lty = 2)
plotconfidence(var3,run_years,2,obs_colours[2])
lines(var4~run_years, col=model_colours[3], lwd = 4, lty = 1) ; points(var4~run_years, col=model_colours[3], pch=16)
lines(var5~run_years, col=model_colours[3], lwd = 4, lty = 2) ; points(var5~run_years, col=model_colours[3], pch=16)
lines(var6~run_years, col=model_colours[3], lwd = 4, lty = 2) ; points(var6~run_years, col=model_colours[3], pch=16)
#legend("bottomright", legend = c(obs_flags[-5],model_flags), col = c(obs_colours[1:4],model_colours), 
#       lty = c(rep(1,length(obs_flags[-5])),rep(2,length(model_flags))), pch=rep(NA,length(c(obs_flags[-5],model_flags))), horiz = FALSE, bty = "n", cex=1.8, lwd=3, ncol = 2)
#mtext("Year", side=1, padj=2.0,cex=1.6)
mtext(expression(paste("Gross Primary Productivity (TgC y",r^-1,")",sep="")), side=2, padj=-1.6, cex=1.5)

# Now plot fire
var3  = cbind(cbind(c(obs_fire_mean_domain_TgCyr),c(obs_fire_min_domain_TgCyr)),c(obs_fire_max_domain_TgCyr))
var4  = mean_fire_TgCyr  ; var5  = fire_lower_TgCyr ; var6  = fire_upper_TgCyr
zrange = range(c(var3,var4,var5,var6), na.rm=TRUE)*c(0.9,1.1)
plot(var4~run_years, main="", cex.lab=2, cex.main=2, cex.axis=1.8, ylim=zrange,
      col=model_colours[1], type="l", lwd=4, lty=2, ylab="", xlab="")
plotconfidence(var3,run_years,2,obs_colours[3])
lines(var4~run_years, col=model_colours[4], lwd=4, lty = 1) ; points(var4~run_years, col=model_colours[4], pch=16)
lines(var5~run_years, col=model_colours[4], lwd=4, lty = 2) ; points(var5~run_years, col=model_colours[4], pch=16)
lines(var6~run_years, col=model_colours[4], lwd=4, lty = 2) ; points(var5~run_years, col=model_colours[4], pch=16)
mtext("Year", side=1, padj=2.0,cex=1.6)
mtext(expression(paste("Fire Emissions (TgC y",r^-1,")",sep="")), side=2, padj=-1.6,cex=1.5)
dev.off()

###
## Statistical significance / trend maps for C-budget terms

# Comparison of NBE, wood and soil stock change over the analysis period by model
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NBP_dCwood_dCsom.png",sep=""), height = 1000, width = 3500, res = 300)
# Plot differences
par(mfrow=c(1,3), mar=c(0.05,1,0.05,7.0), omi = c(0.01,0.4,0.3,0.05))
# Create raster
var1 = rast(vals = t(365.25*-grid_output$mean_nbe_gCm2day[,dim(area)[2]:1,mid_quant]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = rast(vals = t((1/nos_years)*grid_output$final_dCwood_gCm2[,dim(area)[2]:1,mid_quant]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((1/nos_years)*grid_output$final_dCsom_gCm2[,dim(area)[2]:1,mid_quant]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Crop to size
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask)
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3)
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
## Restrict parameter range to +/- 800 gC/m2/yr
#var1[var1 > 800] = 800 ; var1[var1 < -800] = -800
#var2[var2 > 800] = 800 ; var2[var2 < -800] = -800
#var3[var3 > 800] = 800 ; var3[var3 < -800] = -800
# Convert Units gC/m2/yr -> MgC/ha/yr
var1 = var1 * 1e-2   ; var2 = var2 * 1e-2 ; var3 = var3 * 1e-2 
tmp = c(range(values(var1), na.rm=TRUE),range(values(var2), na.rm=TRUE))
tmp1 = c(range(values(var3), na.rm=TRUE))
zrange = max(abs(range(tmp, na.rm=TRUE))) * c(-1,1)
zrange1 = max(abs(c(range(tmp1, na.rm=TRUE)))) * c(-1,1)
# C1 Mean annual NBP, dCwood, dCsom
plot(var1, ylab="", xlab="", main="",  mar=NA, bty = "n",
     xaxt = "n", yaxt = "n", range=zrange,
     col=colour_choices_sign, cex.lab=2, cex.main=2.2, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
plot(landmask, add=TRUE, lwd=0.5)
points(grid_long[nbp_sig_longitude+0.5,1],grid_lat[1,nbp_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
mtext(expression(paste("NBP (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.6, padj = +0.5, adj = 0.5)
plot(var2, ylab="", xlab="", main="",  mar=NA, bty = "n",
     xaxt = "n", yaxt = "n", range=zrange,
     col=colour_choices_sign, cex.lab=2, cex.main=2.2, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
plot(landmask, add=TRUE, lwd=0.5)
points(grid_long[dCwood_sig_longitude+0.5,1],grid_lat[1,dCwood_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
mtext(expression(paste("Wood Change (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.6, padj = +0.5, adj = 0.5)
plot(var3, ylab="", xlab="", main="",  mar=NA, bty = "n",
     xaxt = "n", yaxt = "n", range=zrange1,
     col=colour_choices_sign, cex.lab=2, cex.main=2.2, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
plot(landmask, add=TRUE, lwd=0.5)
points(grid_long[dCsom_sig_longitude+0.5,1],grid_lat[1,dCsom_sig_latitude+0.5], xlab="", ylab="", pch=16,cex=0.4, col="cyan")
mtext(expression(paste("Soil Change (MgC h",a^-1," y",r^-1,")",sep="")), side = 3, cex = 1.6, padj = +0.5, adj = 0.5)
dev.off()

# GPP, Rauto, Rhet trend
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_GPP_Rauto_Rhet_LAI_trend.png",sep=""), height = 2800, width = 3300, res = 300)
# Plot differences
par(mfrow=c(2,2), mar=c(0.1,0.1,1.0,8),omi=c(0.05,0.1,0.2,0.15))
# C1 - GPP. Rauto, Rhet and LAI trends
var1 = rast(vals = t(1e-2*gpp_trend[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var2 = rast(vals = t(1e-2*rauto_trend[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t(1e-2*rhet_trend[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var16 = rast(vals = t(lai_trend[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Crop to size
var1 = crop(var1, landmask) ; var2 = crop(var2, landmask) ; var3 = crop(var3, landmask)
var16 = crop(var16, landmask)
var1 = mask(var1, landmask) ; var2 = mask(var2, landmask) ; var3 = mask(var3, landmask)
var16 = mask(var16, landmask) 
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3) ; var16 = trim(var16)
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# Determine gradients
zrange = max(abs(quantile(c(values(var1),values(var2),values(var3)), prob = c(0.001,0.999), na.rm=TRUE))) * c(-1,1)
zrange1 = max(abs(quantile(values(var16), prob = c(0.001,0.999), na.rm=TRUE))) * c(-1,1)
# Restrict to sensible bounds
#var1[var1 < zrange[1]] = zrange[1] ; var1[var1 > zrange[2]] = zrange[2]
#var2[var2 < zrange[1]] = zrange[1] ; var2[var2 > zrange[2]] = zrange[2]
#var3[var3 < zrange[1]] = zrange[1] ; var3[var3 > zrange[2]] = zrange[2]
# C1 GPP, Rauto, Rhet trend (MgC/ha/yr2)
plot(var1, ylab="", xlab="", main="", mar=NA, bty = "n",
     xaxt = "n", yaxt = "n", range=zrange,
     col=colour_choices_sign, cex.lab=2, cex.main=2.2, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
plot(landmask, add=TRUE, lwd=0.5)
mtext(expression(paste('GPP Trend (MgC h',a^-1,'y',r^-2,')',sep="")), cex=1.8, padj = +0.25)
plot(var2, ylab="", xlab="", main="",  mar=NA, bty = "n",
     xaxt = "n", yaxt = "n", range=zrange,
     col=colour_choices_sign, cex.lab=2, cex.main=2.2, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
plot(landmask, add=TRUE, lwd=0.5)
mtext(expression(paste(R[auto],' Trend (MgC h',a^-1,'y',r^-2,')',sep="")), cex=1.8, padj = +0.25)
plot(var3, ylab="", xlab="", main="",  mar=NA, bty = "n",
           xaxt = "n", yaxt = "n", range=zrange,
           col=colour_choices_sign, cex.lab=2, cex.main=2.2, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
plot(landmask, add=TRUE, lwd=0.5)
mtext(expression(paste(R[het],' Trend (MgC h',a^-1,'y',r^-2,')',sep="")), cex=1.8, padj = +0.25)
plot(var16, ylab="", xlab="", main="",  mar=NA, bty = "n",
     xaxt = "n", yaxt = "n", range=zrange1,
     col=colour_choices_sign, cex.lab=2, cex.main=2.2, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
plot(landmask, add=TRUE, lwd=0.5)
mtext(expression(paste('LAI Trend (',m^2,'',m^-2,' y',r^-2,')',sep="")), cex=1.8, padj = +0.25)
dev.off()

###
## Relative growth rate

png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NPP_biomass_ratio.png",sep=""), height = 2000, width = 3000, res = 300)
par(mfrow=c(1,1), mar=c(0.01,1.5,0.3,7),omi=c(0.01,0.1,0.01,0.1))
var1 = (grid_output$mean_npp_gCm2day[,,mid_quant]*365.25) / grid_output$mean_biomass_gCm2[,,mid_quant]
var1 = rast(vals = t(var1[,dim(var1)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
plot(var1, main="", range = c(0,max(values(var1), na.rm=TRUE)), col=colour_choices_loss, xaxt = "n", yaxt = "n", mar=NA, bty = "n",
     cex.lab=2, cex.axis = 2, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
plot(landmask, add=TRUE, lwd=0.5)
mtext(expression('Annual NPP:Biomass'), side = 2, cex = 1.6, padj = -0.25, adj = 0.5)
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
var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = rast(vals = t((var4)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = rast(vals = t((var5)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var6 = rast(vals = t((var6)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = rast(vals = t((var7)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = rast(vals = t((var8)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Trim to data area
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3) ; var4 = trim(var4) 
var5 = trim(var5) ; var6 = trim(var6) ; var7 = trim(var7) ; var8 = trim(var8) 
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
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
plot(var1, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=rev(colour_choices_default))
mtext(expression(paste('NBE (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)     
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, range=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext(expression(paste('GPP (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, range=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste('Reco (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var4, range=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste('Fire (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
# Mean annual estimates uncertainty
plot(var5, range=zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_CI)
mtext(expression(paste('NBE CI (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var6, range=zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_CI)
mtext(expression(paste('GPP CI (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var7, range=zrange7, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_CI)
mtext(expression(paste('Reco CI (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var8, range=zrange8, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_CI)
mtext(expression(paste('Fire CI (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# C fluxes
# Assign variables
var1 = grid_output$mean_nbe_gCm2day[,,mid_quant]*1e-2*365.25
var2 = grid_output$mean_gpp_gCm2day[,,mid_quant]*1e-2*365.25
var3 = grid_output$mean_reco_gCm2day[,,mid_quant]*1e-2*365.25 
var4 = grid_output$mean_fire_gCm2day[,,mid_quant]*1e-2*365.25 
# Apply filter
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var4[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Convert to raster
var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = rast(vals = t((var4)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Trim to data area
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3) ; var4 = trim(var4) 
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# Determine ranges
zrange1 = c(-1,1)*max(abs(range(values(var1),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(c(values(var2),values(var3)),na.rm=TRUE)))
zrange3 = zrange2 #c(0,1)*max(abs(range(values(var3),na.rm=TRUE)))
zrange4 = c(0,1)*max(abs(range(values(var4),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_C_fluxes_median_v2.png",sep=""), height = 2100, width = 5000*0.55, res = 300)
par(mfrow=c(2,2), mar=c(0.5,0.5,2.8,7),omi=c(0.1,0.4,0.12,0.2))
# Mean annual median estimates
plot(var1, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.0, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=rev(colour_choices_default))
mtext(expression(paste('NBE (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, range=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.0, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext(expression(paste('GPP (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, range=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.0, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste('Reco (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var4, range=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.0, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste('Fire (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Correlation between:
# GPP~Biomass, GPP~Fire, Biomass~Fire
var1 = as.vector(grid_output$mean_gpp_gCm2day[,,mid_quant]) * 365.25 * 1e-2
var2 = as.vector(grid_output$mean_biomass_gCm2[,,mid_quant]) * 1e-2
var3 = as.vector(grid_output$mean_fire_gCm2day[,,mid_quant]) * 365.25 * 1e-2
gppbiomass = lm(var2~var1)
gppfire = lm(var3~var1)
firebiomass = lm(var2~var3)
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_GPP_Biomass_Fire_correlation.png",sep=""), height = 2700, width = 4900, res = 300)
par(mfrow=c(2,3), mar=c(4.0,5.0,2.5,2.5),omi=c(0.2,0.2,0.18,0.1))
plot(var1,var2, pch=16, xlab = expression(paste('GPP (MgC h',a^-1,' y',r^-1,')',sep="")), ylab = expression(paste('Biomass (MgC h',a^-1,')',sep="")),
     cex.axis = 1.8, cex.lab = 1.8, cex.main=1.8, main=paste('R2 = ',round(summary(gppbiomass)$adj.r.squared,digits = 2),sep=""), cex = 1.2)
abline(gppbiomass, col="red", lwd=2)
plot(var1,var3, pch=16, xlab = expression(paste('GPP (MgC h',a^-1,' y',r^-1,')',sep="")), ylab = expression(paste('Fire (MgC h',a^-1,' y',r^-1,')',sep="")),
     cex.axis = 1.8, cex.lab = 1.8, cex.main=1.8, main=paste('R2 = ',round(summary(gppfire)$adj.r.squared,digits = 2),sep=""), cex = 1.2)
abline(gppfire, col="red", lwd=2)
plot(var3,var2, pch=16, xlab = expression(paste('Fire (MgC h',a^-1,' y',r^-1,')',sep="")), ylab = expression(paste('Biomass (MgC h',a^-1,')',sep="")),
     cex.axis = 1.8, cex.lab = 1.8, cex.main=1.8, main=paste('R2 = ',round(summary(firebiomass)$adj.r.squared,digits = 2),sep=""), cex = 1.2)
abline(firebiomass, col="red", lwd=2)
hist(var1, main="", ylab="No. pixels", xlab=expression(paste('GPP (MgC h',a^-1,' y',r^-1,')',sep="")),
     cex.lab=1.8, cex.axis = 1.8)
hist(var2, main="", ylab=" ", xlab=expression(paste('Biomass (MgC h',a^-1,')',sep="")),
     cex.lab=1.8, cex.axis = 1.8)
hist(var3, main="", ylab="", xlab=expression(paste('Fire (MgC h',a^-1,' y',r^-1,')',sep="")),
     cex.lab=1.8, cex.axis = 1.8)
dev.off()

# C fluxes
# Assign variables
var1 = grid_output$mean_nbe_gCm2day[,,mid_quant]*1e-2*365.25
var2 = grid_output$mean_gpp_gCm2day[,,mid_quant]*1e-2*365.25
var3 = grid_output$mean_reco_gCm2day[,,mid_quant]*1e-2*365.25 
var4 = grid_output$mean_fire_gCm2day[,,mid_quant]*1e-2*365.25 
# Apply filter
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var4[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Convert to raster
var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = rast(vals = t((var4)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Trim to data area
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3) ; var4 = trim(var4) 
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# Determine ranges
zrange1 = c(-1,1)*max(abs(range(values(var1),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(values(var2),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(values(var3),na.rm=TRUE)))
zrange4 = c(0,1)*max(abs(range(values(var4),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_C_fluxes_median.png",sep=""), height = 2100*0.5, width = 5000, res = 300)
par(mfrow=c(1,4), mar=c(0.5,0.5,2.8,7),omi=c(0.1,0.4,0.2,0.2))
# Mean annual median estimates
plot(var1, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=rev(colour_choices_default))
mtext(expression(paste('NBE (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, range=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext(expression(paste('GPP (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, range=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste('Reco (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var4, range=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste('Fire (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
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
var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = rast(vals = t((var4)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = rast(vals = t((var5)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var6 = rast(vals = t((var6)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var7 = rast(vals = t((var7)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var8 = rast(vals = t((var8)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Trim to data area
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3) ; var4 = trim(var4) 
var5 = trim(var5) ; var6 = trim(var6) ; var7 = trim(var7) ; var8 = trim(var8) 
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
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
plot(var1, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=rev(colour_choices_default))
mtext(expression(paste('NBE (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, range=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext(expression(paste('GPP (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, range=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste('Reco (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var4, range=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste('Fire (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
# Mean annual estimates uncertainty
plot(var5, range=zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_CI)
mtext(expression(paste('NBE CI (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var6, range=zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_CI)
mtext(expression(paste('GPP CI (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var7, range=zrange7, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_CI)
mtext(expression(paste('Reco CI (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var8, range=zrange8, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_CI)
mtext(expression(paste('Fire (MgC h',a^-1,' y',r^-1,')',sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

###
## Plot the final C stocks, change and uncertainty

# Final stocks
# Assign variables
var1 = grid_output$final_Ctotal_gCm2[,,mid_quant]*1e-2
var2 = grid_output$final_biomass_gCm2[,,mid_quant]*1e-2 
var3 = grid_output$final_dom_gCm2[,,mid_quant]*1e-2 
var4 = (grid_output$final_Ctotal_gCm2[,,high_quant]-grid_output$final_Ctotal_gCm2[,,low_quant])*1e-2 
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
var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = rast(vals = t((var4)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = rast(vals = t((var5)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = rast(vals = t((var6)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Trim to data area
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3) ; var4 = trim(var4) 
var5 = trim(var5) ; var6 = trim(var6) 
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
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
plot(var1, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext(expression(paste("Total (MgC h",a^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, range=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext(expression(paste("Biomass (MgC h",a^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, range=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext(expression(paste("DOM (MgC h",a^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
# Final C stocks, confidence interval
plot(var4, range=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_CI)
mtext(expression(paste("Total CI (MgC h",a^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var5, range=zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_CI)
mtext(expression(paste("Biomass CI (MgC h",a^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var6, range=zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_CI)
mtext(expression(paste("DOM CI (MgC h",a^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Final stocks
# Assign variables
var1 = grid_output$final_Ctotal_gCm2[,,mid_quant]*1e-2
var2 = grid_output$final_biomass_gCm2[,,mid_quant]*1e-2 
var3 = grid_output$final_dom_gCm2[,,mid_quant]*1e-2 
# Apply filter
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Convert to raster
var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Trim to data area
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3) 
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# ranges
zrange1 = c(0,1)*max(abs(range(values(var1),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(values(var2),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(values(var3),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_Final_stock_median.png",sep=""), height = 2700*0.5, width = 4900, res = 300)
par(mfrow=c(1,3), mar=c(0.6,0.4,2.9,7),omi=c(0.1,0.4,0.18,0.2))
# Final C stocks, median estimate
plot(var1, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext(expression(paste("Total (MgC h",a^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, range=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext(expression(paste("Biomass (MgC h",a^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, range=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext(expression(paste("DOM (MgC h",a^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Final stocks wood & DOM, wood MRT ~ burnt fraction, DOM ~ burnt fraction
# Assign variables
var1 = grid_output$mean_wood_gCm2[,,mid_quant]*1e-2 
var2 = grid_output$mean_som_gCm2[,,mid_quant]*1e-2 
var3 = grid_output$MTT_wood_years[,,mid_quant]
var4 = grid_output$MTT_som_years[,,mid_quant]
# Apply filter
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var4[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Convert to raster
var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = rast(vals = t((var4)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Trim to data area
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3) ; var4 = trim(var4)
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# ranges
zrange1 = c(0,1)*max(abs(range(values(var1),na.rm=TRUE)))
zrange2 = c(0,1)*max(abs(range(values(var2),na.rm=TRUE)))
zrange3 = c(0,1)*max(abs(range(values(var3),na.rm=TRUE)))
zrange4 = c(0,1)*max(abs(range(values(var4),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_mean_wood_som_stock_woodMRT_somMRT_fire_correlation_median.png",sep=""), height = 1600, width = 3000, res = 300)
par(mfrow=c(2,3), mar=c(4.0,3.0,2.9,1.0),omi=c(0.01,0.10,0.10,0.35))
# Correlation between Wood MRT and fire
plot(as.vector(grid_output$MTT_wood_years[,,mid_quant]) ~ as.vector(BurnedFraction), pch=16,
     cex.axis = 1.5, cex.lab = 1.5, cex = 1.2, xlab="", ylab="")
mtext(side = 2, text = "Wood MRT (years)", cex = 1.0, padj = -2.50)
# Mean C stocks, median estimate
plot(var1, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.0, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext(expression(paste("Wood (MgC h",a^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, range=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.0, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext(expression(paste("Soil (MgC h",a^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
# Correlation between Soil MRT and fire
plot(as.vector(grid_output$MTT_som_years[,,mid_quant]) ~ as.vector(BurnedFraction), pch=16,
     cex.axis = 1.5, cex.lab = 1.5, cex = 1.2, xlab="Annual burnt Fraction", ylab="")
mtext(side = 2, text = "Soil MRT (years)", cex = 1.0, padj = -2.50)
# MRTs for wood and soil
plot(var3, range=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.0, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext(expression(paste("Wood MRT(years)",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var4, range=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.0, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext(expression(paste("Soil MRT (years)",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Change stocks
# Assign variables
var1 = grid_output$final_dCtotal_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var2 = grid_output$final_dCbiomass_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var3 = grid_output$final_dCdom_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var4 = (grid_output$final_dCtotal_gCm2[,,high_quant]-grid_output$final_dCtotal_gCm2[,,low_quant])*1e-2*(1/nos_years)
var5 = (grid_output$final_dCbiomass_gCm2[,,high_quant]-grid_output$final_dCbiomass_gCm2[,,low_quant])*1e-2*(1/nos_years)
var6 = (grid_output$final_dCdom_gCm2[,,high_quant]-grid_output$final_dCdom_gCm2[,,low_quant])*1e-2*(1/nos_years)
# Apply filter
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var4[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var5[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var6[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Convert to raster
var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = rast(vals = t((var4)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = rast(vals = t((var5)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = rast(vals = t((var6)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Trim to data area
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3) ; var4 = trim(var4) 
var5 = trim(var5) ; var6 = trim(var6) 
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
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
plot(var1, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=(colour_choices_default))
mtext(expression(paste(Delta,"Total (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, range=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=(colour_choices_default))
mtext(expression(paste(Delta,"Biomass (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, range=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=(colour_choices_default))
mtext(expression(paste(Delta,"DOM (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
# Final stock changes, confidence interval
plot(var4, range=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_CI)
mtext(expression(paste(Delta,"Total CI (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var5, range=zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_CI)
mtext(expression(paste(Delta,"Biomass CI (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var6, range=zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_CI)
mtext(expression(paste(Delta,"DOM CI (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# NPP and litter fluxes from plant tissues
# Assign variables
var1 = grid_output$mean_combined_alloc_foliage_gCm2day[,,mid_quant]*1e-2*365.25
var2 = grid_output$mean_alloc_roots_gCm2day[,,mid_quant]*1e-2*365.25
var3 = grid_output$mean_alloc_wood_gCm2day[,,mid_quant]*1e-2*365.25
var4 = grid_output$mean_foliage_gCm2[,,mid_quant]*1e-2*(1/grid_output$MTT_foliage_years[,,mid_quant])
var5 = grid_output$mean_roots_gCm2[,,mid_quant]*1e-2*(1/grid_output$MTT_roots_years[,,mid_quant])
var6 = grid_output$mean_wood_gCm2[,,mid_quant]*1e-2*(1/grid_output$MTT_wood_years[,,mid_quant])
# Apply filter
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var4[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var5[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var6[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Convert to raster
var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = rast(vals = t((var4)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = rast(vals = t((var5)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = rast(vals = t((var6)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Trim to data area
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3) ; var4 = trim(var4) 
var5 = trim(var5) ; var6 = trim(var6) 
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# ranges
zrange1 = c(0,1)*max(abs(range(c(values(var1),values(var2),values(var3),values(var4),values(var5),values(var6)),na.rm=TRUE)))
zrange2 = zrange1
zrange3 = zrange1
zrange4 = zrange1
zrange5 = zrange1
zrange6 = zrange1
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NPP_litter_fluxes.png",sep=""), height = 2700, width = 4900, res = 300)
par(mfrow=c(2,3), mar=c(0.5,0.4,2.8,7),omi=c(0.1,0.4,0.2,0.2))
# Median NPP fluxes to plant tissues
plot(var1, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=(colour_choices_gain))
mtext(expression(paste("Foliar NPP (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, range=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=(colour_choices_gain))
mtext(expression(paste("Root NPP (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, range=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=(colour_choices_gain))
mtext(expression(paste("Wood NPP (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
# Median litter fluxes from plant tissues
plot(var4, range=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste("Foliar litter (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var5, range=zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste("Root litter (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var6, range=zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste("Wood litter (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Fire emission and litter fluxes from plant tissues
# Assign variables
var1 = grid_output$mean_FIREemiss_foliage_gCm2day[,,mid_quant]*1e-2*365.25
var2 = grid_output$mean_FIREemiss_roots_gCm2day[,,mid_quant]*1e-2*365.25
var3 = grid_output$mean_FIREemiss_wood_gCm2day[,,mid_quant]*1e-2*365.25
var4 = grid_output$mean_FIRElitter_foliage_gCm2day[,,mid_quant]*1e-2*365.25
var5 = grid_output$mean_FIRElitter_roots_gCm2day[,,mid_quant]*1e-2*365.25
var6 = grid_output$mean_FIRElitter_wood_gCm2day[,,mid_quant]*1e-2*365.25
# Apply filter
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var4[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var5[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var6[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Convert to raster
var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = rast(vals = t((var4)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = rast(vals = t((var5)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = rast(vals = t((var6)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Trim to data area
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3) ; var4 = trim(var4) 
var5 = trim(var5) ; var6 = trim(var6)
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# ranges
zrange1 = c(0,1)*sum(c(max(values(var1), na.rm=TRUE),max(values(var2), na.rm=TRUE),max(values(var3),na.rm=TRUE)))
zrange2 = zrange1
zrange3 = zrange1
zrange4 = c(0,1)*max(range(c(values(var4),values(var5),values(var6)),na.rm=TRUE))
zrange5 = zrange4
zrange6 = zrange4
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_Fire_emission_litter_fluxes.png",sep=""), height = 2700, width = 4900, res = 300)
par(mfrow=c(2,3), mar=c(0.5,0.4,2.8,7),omi=c(0.1,0.4,0.2,0.2))
# Median NPP fluxes to plant tissues
plot(var1, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=(colour_choices_loss))
mtext(expression(paste("Foliar fire emission (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, range=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=(colour_choices_loss))
mtext(expression(paste("Root fire emission (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, range=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=(colour_choices_loss))
mtext(expression(paste("Wood fire emission (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
# Median litter fluxes from plant tissues
plot(var4, range=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste("Foliar fire litter (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var5, range=zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste("Root fire litter (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var6, range=zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste("Wood fire litter (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Change stocks
# Assign variables
var1 = grid_output$final_dCtotal_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var2 = grid_output$final_dCbiomass_gCm2[,,mid_quant]*1e-2*(1/nos_years)
var3 = grid_output$final_dCdom_gCm2[,,mid_quant]*1e-2*(1/nos_years)
# Apply filter
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Convert to raster
var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Trim to data area
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3) 
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# ranges
zrange1 = c(-1,1)*max(abs(range(values(var1),na.rm=TRUE)))
zrange2 = c(-1,1)*max(abs(range(values(var2),na.rm=TRUE)))
zrange3 = c(-1,1)*max(abs(range(values(var3),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_Final_stock_change_median.png",sep=""), height = 2700*0.5, width = 4900, res = 300)
par(mfrow=c(1,3), mar=c(0.5,0.4,2.8,7),omi=c(0.1,0.4,0.2,0.2))
# Final stock changes, median estimates
plot(var1, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=(colour_choices_default))
mtext(expression(paste(Delta,"Total (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, range=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=(colour_choices_default))
mtext(expression(paste(Delta,"Biomass (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, range=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=(colour_choices_default))
mtext(expression(paste(Delta,"DOM (MgC h",a^-1,"y",r^-1,")",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

###
## Plot the MRTwood and NPPwood

# Traits
# Assign variables
wood_mrt_limit = 60
print(paste("Wood MRT plotting range has been limited to ",wood_mrt_limit," years",sep=""))
var1 = pmin(wood_mrt_limit,grid_output$MTT_wood_years[,,mid_quant])
var1 = array(var1, dim = dim(grid_output$MTT_wood_years)[1:2])
var2 = grid_output$NPP_wood_fraction[,,mid_quant]
var3 = grid_output$SS_wood_gCm2[,,mid_quant]*1e-2
var4 = pmin(wood_mrt_limit*2,grid_output$MTT_wood_years[,,high_quant]) - pmin(wood_mrt_limit*2,grid_output$MTT_wood_years[,,low_quant])
var4 = array(var4, dim = dim(grid_output$MTT_wood_years)[1:2])
var5 = grid_output$NPP_wood_fraction[,,high_quant] - grid_output$NPP_wood_fraction[,,low_quant]
var6 = (grid_output$SS_wood_gCm2[,,high_quant]*1e-2) - (grid_output$SS_wood_gCm2[,,low_quant]*1e-2)
# Apply filter
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var4[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var5[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var6[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Convert to raster
var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var4 = rast(vals = t((var4)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var5 = rast(vals = t((var5)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var6 = rast(vals = t((var6)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Trim to data area
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3) ; var4 = trim(var4) 
var5 = trim(var5) ; var6 = trim(var6) 
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
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
plot(var1, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext("MRT Wood (years)", side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, range=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext("NPP wood (0-1)", side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, range=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext("SS wood (MgC/ha)", side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
# Ecosystem traits, confidence intervals
plot(var4, range=zrange4, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext("MRT Wood CI (years)", side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var5, range=zrange5, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext("NPP wood CI (0-1)", side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var6, range=zrange6, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext("SS wood CI (MgC/ha)", side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Traits
# Assign variables
wood_mrt_limit = 60
print(paste("Wood MRT plotting range has been limited to ",wood_mrt_limit," years",sep=""))
var1 = pmin(wood_mrt_limit,grid_output$MTT_wood_years[,,mid_quant])
var1 = array(var1, dim = dim(grid_output$MTT_wood_years)[1:2])
var2 = grid_output$NPP_wood_fraction[,,mid_quant]
var3 = grid_output$SS_wood_gCm2[,,mid_quant]*1e-2
var4 = pmin(wood_mrt_limit*2,grid_output$MTT_wood_years[,,high_quant]) - pmin(wood_mrt_limit*2,grid_output$MTT_wood_years[,,low_quant])
# Apply filter
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Convert to raster
var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Trim to data area
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3) 
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# ranges
zrange1 = c(0,1)*max(abs(range(values(var1),na.rm=TRUE)))
zrange2 = c(0,1)
zrange3 = c(0,1)*max(abs(range(values(var3),na.rm=TRUE)))
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NPP_MRT_SS_median.png",sep=""), height = 2700*0.5, width = 4900, res = 300)
par(mfrow=c(1,3), mar=c(0.5,0.3,2.8,8),omi=c(0.1,0.3,0.2,0.2))
# Ecosystem traits, median estimates
plot(var1, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext("MRT Wood (years)", side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, range=zrange2, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext("NPP wood (0-1)", side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, range=zrange3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_gain)
mtext("SS wood (MgC/ha)", side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

###
## Partition the importance of disturbance on residence times

# Extract the proportions         
var1 = grid_output$NaturalFractionOfTurnover_wood[,,mid_quant]
var2 = grid_output$FireFractionOfTurnover_wood[,,mid_quant]
var3 = grid_output$HarvestFractionOfTurnover_wood[,,mid_quant]
# Filter for the miombo AGB map locations
var1[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var2[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
var3[which(landfilter == 0 | is.na(landfilter) == TRUE)] = NA
# Filter for the miombo AGB map locations
var1[which(landfilter == 1 & is.na(var1) == TRUE)] = 0
var2[which(landfilter == 1 & is.na(var2) == TRUE)] = 0
var3[which(landfilter == 1 & is.na(var3) == TRUE)] = 0
# Convert to raster
var1 = rast(vals = t((var1)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = rast(vals = t((var2)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((var3)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# Trim to data area
var1 = trim(var1) ; var2 = trim(var2) ; var3 = trim(var3)
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
# specify ranges
zrange1 = c(0,1)
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_Wood_turnover_contribution.png",sep=""), height = 1300, width = 4900, res = 300)
par(mfrow=c(1,3), mar=c(0.5,0.4,3.0,7),omi=c(0.1,0.3,0.1,0.2))
# Partitioning of wood turnover, median estimate
plot(var1, range=zrange1, xaxt = "n", yaxt = "n", type="continuous", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext("Natural MRT comp (0-1)", side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, range=zrange1, xaxt = "n", yaxt = "n", type="continuous", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext("Fire MRT comp (0-1)", side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
legend_args = list(ext=e, cex=1.0) ; if (length(unlist(unique(var3))) == 1) {legend_args = list(cex=1.0)}
plot(var3, range=zrange1, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = legend_args,
     main = "", col=colour_choices_loss)
mtext("Biomass removal MRT comp (0-1)", side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Convert to raster
var1 = rast(vals = t((HarvestFraction)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext)) 
var2 = rast(vals = t((BurnedFraction)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
var3 = rast(vals = t((FireFreq)[,dim(area)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_disturbance_drivers.png",sep=""), height = 1300, width = 4900, res = 300)
par(mfrow=c(1,3), mar=c(0.5,0.4,3.0,7),omi=c(0.1,0.3,0.1,0.2))
# Partitioning of wood turnover, median estimate
plot(var1, range=c(0,1), xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste("Deforested fraction (0-1)",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var2, range=c(0,1), xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste("Burnt fraction (0-1)",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
plot(var3, xaxt = "n", yaxt = "n", cex.lab=2, cex.main=2.5, mar=NA, bty = "n",
     cex.axis = 2.5, axes = FALSE, pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0),
     main = "", col=colour_choices_loss)
mtext(expression(paste("No. fires per year",sep="")), side=3, cex = 1.8, padj = 0.9)
plot(landmask, add=TRUE, lwd=0.5)
dev.off()

# Create a mask of pixels which contain minimum disturbance
threshold = 0.99 ; loess_span = 0.25 ; n_interp = 50 ; se_ci_scalar = 1.98
natural_dom = grid_output$NaturalFractionOfTurnover_wood[,,mid_quant]
natural_dom[which(natural_dom >= threshold)] = 1 ; natural_dom[which(natural_dom < threshold)] = NA
# Plot Foliage, fine root, wood, litter(foliar+fine root+wood?), soil mean residence times against main meteorology
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_MRT_meteorology_association_no_disturbance.png",sep=""), height = 2200, width = 4500, res = 300)
par(mfrow=c(3,5), mar=c(4,2,1.4,1), omi = c(0.1,0.2,0.12,0.1))
# Temperature
vary = as.vector(natural_dom*grid_output$MTT_foliage_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
mtext(expression(paste("Foliar MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
#mtext(expression('C1'), side = 2, cex = 1.6, padj = -2.5, adj = 0.5)
vary = as.vector(natural_dom*grid_output$MTT_roots_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
mtext(expression(paste("Root MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
vary = as.vector(natural_dom*grid_output$MTT_wood_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab="Mean Temperature (C)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
mtext(expression(paste("Wood MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
vary = as.vector(natural_dom*grid_output$MTT_litter_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
mtext(expression(paste("Litter MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
vary = as.vector(natural_dom*grid_output$MTT_som_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
mtext(expression(paste("Soil MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
# Precipitation
vary = as.vector(natural_dom*grid_output$MTT_foliage_years[,,mid_quant])
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
vary = as.vector(natural_dom*grid_output$MTT_roots_years[,,mid_quant])
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
vary = as.vector(natural_dom*grid_output$MTT_wood_years[,,mid_quant])
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab="Mean precipitation (mm/yr)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
vary = as.vector(natural_dom*grid_output$MTT_litter_years[,,mid_quant])
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
vary = as.vector(natural_dom*grid_output$MTT_som_years[,,mid_quant])
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
# Vapour pressure deficit
vary = as.vector(natural_dom*grid_output$MTT_foliage_years[,,mid_quant])
varx = as.vector(grid_output$mean_vpd_Pa)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
vary = as.vector(natural_dom*grid_output$MTT_roots_years[,,mid_quant])
varx = as.vector(grid_output$mean_vpd_Pa)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
vary = as.vector(natural_dom*grid_output$MTT_wood_years[,,mid_quant])
varx = as.vector(grid_output$mean_vpd_Pa)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab="Mean VPD (Pa)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
vary = as.vector(natural_dom*grid_output$MTT_litter_years[,,mid_quant])
varx = as.vector(grid_output$mean_vpd_Pa)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
vary = as.vector(natural_dom*grid_output$MTT_som_years[,,mid_quant])
varx = as.vector(grid_output$mean_vpd_Pa)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
dev.off()

# Create a mask of pixels which contain minimum disturbance
threshold = 0 ; loess_span = 0.25 ; n_interp = 50 ; se_ci_scalar = 1.98
natural_dom = grid_output$NaturalFractionOfTurnover_wood[,,mid_quant]
natural_dom[which(natural_dom >= threshold)] = 1 ; natural_dom[which(natural_dom < threshold)] = NA
# Plot Foliage, fine root, wood, litter(foliar+fine root+wood?), soil mean residence times against main meteorology
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_MRT_meteorology_association.png",sep=""), height = 2200, width = 4500, res = 300)
par(mfrow=c(3,5), mar=c(4,2,1.4,1), omi = c(0.1,0.2,0.12,0.1))
# Temperature
vary = as.vector(natural_dom*grid_output$MTT_foliage_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
mtext(expression(paste("Foliar MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
#mtext(expression('C1'), side = 2, cex = 1.6, padj = -2.5, adj = 0.5)
vary = as.vector(natural_dom*grid_output$MTT_roots_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
mtext(expression(paste("Root MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
vary = as.vector(natural_dom*grid_output$MTT_wood_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab="Mean Temperature (C)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
mtext(expression(paste("Wood MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
vary = as.vector(natural_dom*grid_output$MTT_litter_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
mtext(expression(paste("Litter MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
vary = as.vector(natural_dom*grid_output$MTT_som_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
mtext(expression(paste("Soil MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
# Precipitation
vary = as.vector(natural_dom*grid_output$MTT_foliage_years[,,mid_quant])
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
vary = as.vector(natural_dom*grid_output$MTT_roots_years[,,mid_quant])
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
vary = as.vector(natural_dom*grid_output$MTT_wood_years[,,mid_quant])
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab="Mean precipitation (mm/yr)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
vary = as.vector(natural_dom*grid_output$MTT_litter_years[,,mid_quant])
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
vary = as.vector(natural_dom*grid_output$MTT_som_years[,,mid_quant])
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
# Vapour pressure deficit
vary = as.vector(natural_dom*grid_output$MTT_foliage_years[,,mid_quant])
varx = as.vector(grid_output$mean_vpd_Pa)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
vary = as.vector(natural_dom*grid_output$MTT_roots_years[,,mid_quant])
varx = as.vector(grid_output$mean_vpd_Pa)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
vary = as.vector(natural_dom*grid_output$MTT_wood_years[,,mid_quant])
varx = as.vector(grid_output$mean_vpd_Pa)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab="Mean VPD (Pa)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
vary = as.vector(natural_dom*grid_output$MTT_litter_years[,,mid_quant])
varx = as.vector(grid_output$mean_vpd_Pa)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
vary = as.vector(natural_dom*grid_output$MTT_som_years[,,mid_quant])
varx = as.vector(grid_output$mean_vpd_Pa)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
dev.off()

# Plot Foliage, fine root, wood, litter(foliar+fine root+wood?), soil mean residence times against main meteorology
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_wood_soil_MRT_MAT_MAP_association.png",sep=""), height = 2500, width = 3000, res = 300)
par(mfrow=c(2,2), mar=c(4,2,1.4,1), omi = c(0.1,0.4,0.12,0.1))
## MRT Wood - Temperature
# Assign whole dataset for plotting
vary = as.vector(grid_output$MTT_wood_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.0, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
# Create a mask of pixels which highly disturbed only
threshold = 0.95 ; loess_span = 0.25 ; n_interp = 50 ; se_ci_scalar = 1.98
natural_dom = grid_output$NaturalFractionOfTurnover_wood[,,mid_quant]
natural_dom[which(natural_dom >= threshold)] = NA ; natural_dom[which(natural_dom < threshold)] = 1
vary = as.vector(natural_dom*grid_output$MTT_wood_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
# Overlay loess of whole dataset including disturbance
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="red", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="red")
# Create a mask of pixels which contain minimum disturbance
threshold = 0.95 ; loess_span = 0.25 ; n_interp = 50 ; se_ci_scalar = 1.98
natural_dom = grid_output$NaturalFractionOfTurnover_wood[,,mid_quant]
natural_dom[which(natural_dom >= threshold)] = 1 ; natural_dom[which(natural_dom < threshold)] = NA
# Overlay loess of whole dataset minimum disturbance
vary = as.vector(natural_dom*grid_output$MTT_wood_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
# Overlay loess of minimum disturbance 
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
mtext(expression(paste("Wood MRT (",y^-1,")",sep="")), cex = 1.3, padj = -1.5, side = 2)

## MRT Wood - Precipitation
# Assign whole dataset for plotting
vary = as.vector(grid_output$MTT_wood_years[,,mid_quant])
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
plot(vary ~ varx, main="", ylab="", xlab=" ", 
     pch=16, cex=1.0, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
# Create a mask of pixels which contain highly disturbed only
threshold = 0.95 ; loess_span = 0.25 ; n_interp = 50 ; se_ci_scalar = 1.98
natural_dom = grid_output$NaturalFractionOfTurnover_wood[,,mid_quant]
natural_dom[which(natural_dom >= threshold)] = NA ; natural_dom[which(natural_dom < threshold)] = 1
# Assign whole dataset for plotting
vary = as.vector(natural_dom*grid_output$MTT_wood_years[,,mid_quant])
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
# Overlay loess of whole dataset including disturbance
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="red", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="red")
# Create a mask of pixels which contain minimum disturbance
threshold = 0.95 ; loess_span = 0.25 ; n_interp = 50 ; se_ci_scalar = 1.98
natural_dom = grid_output$NaturalFractionOfTurnover_wood[,,mid_quant]
natural_dom[which(natural_dom >= threshold)] = 1 ; natural_dom[which(natural_dom < threshold)] = NA
# Overlay loess of whole dataset minimum disturbance
vary = as.vector(natural_dom*grid_output$MTT_wood_years[,,mid_quant])
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
# Overlay loess of minimum disturbance 
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
#mtext(expression(paste("Wood MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)

## MRT Soil - Temperature
# Assign whole dataset for plotting
vary = as.vector(grid_output$MTT_som_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
plot(vary ~ varx, main="", ylab="", xlab=expression("Mean Temperature ("*degree*"C)"), 
     pch=16, cex=1.0, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
# Create a mask of pixels which contains highly disturbed only
threshold = 0.95 ; loess_span = 0.25 ; n_interp = 50 ; se_ci_scalar = 1.98
natural_dom = grid_output$NaturalFractionOfTurnover_wood[,,mid_quant]
natural_dom[which(natural_dom >= threshold)] = NA ; natural_dom[which(natural_dom < threshold)] = 1
# Assign whole dataset for plotting
vary = as.vector(natural_dom*grid_output$MTT_som_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
# Overlay loess of whole dataset including disturbance
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="red", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="red")
# Create a mask of pixels which contain minimum disturbance
threshold = 0.95 ; loess_span = 0.25 ; n_interp = 50 ; se_ci_scalar = 1.98
natural_dom = grid_output$NaturalFractionOfTurnover_wood[,,mid_quant]
natural_dom[which(natural_dom >= threshold)] = 1 ; natural_dom[which(natural_dom < threshold)] = NA
# Overlay loess of whole dataset minimum disturbance
vary = as.vector(natural_dom*grid_output$MTT_som_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
# Overlay loess of minimum disturbance 
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
mtext(expression(paste("Soil MRT (",y^-1,")",sep="")), cex = 1.3, padj = -1.5, side = 2)
## MRT Soil - Precipitation
# Assign whole dataset for plotting
vary = as.vector(grid_output$MTT_som_years[,,mid_quant])
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=expression(paste("Mean precipitation (kg",H[2],"O ",m^{-2},y^{-1},")",sep="")), 
     pch=16, cex=1.0, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
# Create a mask of pixels which contain highly disturbed only
threshold = 0.95 ; loess_span = 0.25 ; n_interp = 50 ; se_ci_scalar = 1.98
natural_dom = grid_output$NaturalFractionOfTurnover_wood[,,mid_quant]
natural_dom[which(natural_dom >= threshold)] = NA ; natural_dom[which(natural_dom < threshold)] = 1
# Assign whole dataset for plotting
vary = as.vector(natural_dom*grid_output$MTT_som_years[,,mid_quant])
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
# Overlay loess of whole dataset including disturbance
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="red", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="red")
# Create a mask of pixels which contain minimum disturbance
threshold = 0.95 ; loess_span = 0.25 ; n_interp = 50 ; se_ci_scalar = 1.98
natural_dom = grid_output$NaturalFractionOfTurnover_wood[,,mid_quant]
natural_dom[which(natural_dom >= threshold)] = 1 ; natural_dom[which(natural_dom < threshold)] = NA
# Overlay loess of whole dataset minimum disturbance
vary = as.vector(natural_dom*grid_output$MTT_som_years[,,mid_quant])
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
# Overlay loess of minimum disturbance 
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
#mtext(expression(paste("Soil MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
legend("topright", legend=c("Highly disturbed", "Minimal disturbed"), col=c("red","cyan"), pch=c(NA,NA), lty=c(1,1), lwd=1.4, bty = "n", cex=1.2)
dev.off()

# Plot Foliage, fine root, wood, litter(foliar+fine root+wood?), soil mean residence times against main meteorology
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_wood_soil_MRT_ratio_MAT_MAP_association.png",sep=""), height = 1400, width = 3000, res = 300)
par(mfrow=c(1,2), mar=c(4,2,1.4,1.5), omi = c(0.1,0.5,0.12,0.1))
## MRT Wood : MRT Soil - Temperature
# Assign whole dataset for plotting
vary = as.vector(grid_output$MTT_wood_years[,,mid_quant] / grid_output$MTT_som_years[,,mid_quant])
varx = as.vector(grid_output$mean_temperature_C)
plot(vary ~ varx, main="", ylab="", xlab=expression("Mean Temperature ("*degree*"C)"), 
     pch=16, cex=1.0, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
abline(1,0,col="grey", lwd=1.5)
# Create a mask of pixels which contains highly disturbed only
threshold = 0.95 ; loess_span = 0.25 ; n_interp = 50 ; se_ci_scalar = 1.98
natural_dom = grid_output$NaturalFractionOfTurnover_wood[,,mid_quant]
natural_dom[which(natural_dom >= threshold)] = NA ; natural_dom[which(natural_dom < threshold)] = 1
# Assign whole dataset for plotting
vary = as.vector(natural_dom*(grid_output$MTT_wood_years[,,mid_quant] / grid_output$MTT_som_years[,,mid_quant]))
varx = as.vector(grid_output$mean_temperature_C)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
# Overlay loess of whole dataset including disturbance
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="red", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="red")
# Create a mask of pixels which contain minimum disturbance
threshold = 0.95 ; loess_span = 0.25 ; n_interp = 50 ; se_ci_scalar = 1.98
natural_dom = grid_output$NaturalFractionOfTurnover_wood[,,mid_quant]
natural_dom[which(natural_dom >= threshold)] = 1 ; natural_dom[which(natural_dom < threshold)] = NA
# Overlay loess of whole dataset minimum disturbance
vary = as.vector(natural_dom*(grid_output$MTT_wood_years[,,mid_quant] / grid_output$MTT_som_years[,,mid_quant]))
varx = as.vector(grid_output$mean_temperature_C)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
# Overlay loess of minimum disturbance 
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
mtext(expression(paste("Wood MRT : Soil MRT",sep="")), cex = 1.8, padj = -2.3, side = 2)
## MRT Soil - Precipitation
# Assign whole dataset for plotting
vary = as.vector((grid_output$MTT_wood_years[,,mid_quant] / grid_output$MTT_som_years[,,mid_quant]))
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
plot(vary ~ varx, main="", ylab="", xlab=expression(paste("Mean precipitation (kg",H[2],"O ",m^{-2},y^{-1},")",sep="")), 
     pch=16, cex=1.0, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
abline(1,0,col="grey", lwd=1.5)
# Create a mask of pixels which contain highly disturbed only
threshold = 0.95 ; loess_span = 0.25 ; n_interp = 50 ; se_ci_scalar = 1.98
natural_dom = grid_output$NaturalFractionOfTurnover_wood[,,mid_quant]
natural_dom[which(natural_dom >= threshold)] = NA ; natural_dom[which(natural_dom < threshold)] = 1
# Assign whole dataset for plotting
vary = as.vector(natural_dom*(grid_output$MTT_wood_years[,,mid_quant] / grid_output$MTT_som_years[,,mid_quant]))
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
# Overlay loess of whole dataset including disturbance
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="red", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="red")
# Create a mask of pixels which contain minimum disturbance
threshold = 0.95 ; loess_span = 0.25 ; n_interp = 50 ; se_ci_scalar = 1.98
natural_dom = grid_output$NaturalFractionOfTurnover_wood[,,mid_quant]
natural_dom[which(natural_dom >= threshold)] = 1 ; natural_dom[which(natural_dom < threshold)] = NA
# Overlay loess of whole dataset minimum disturbance
vary = as.vector(natural_dom*(grid_output$MTT_wood_years[,,mid_quant] / grid_output$MTT_som_years[,,mid_quant]))
varx = as.vector(grid_output$mean_precipitation_kgH2Om2yr)
newvarx = seq(min(varx, na.rm=TRUE),max(varx, na.rm=TRUE), length.out=n_interp)
# Overlay loess of minimum disturbance 
tmp_lm = loess(vary ~ varx, span=loess_span)
tmp_lm = predict(tmp_lm, newdata = data.frame(varx = newvarx), se=TRUE)
plotCI(y = tmp_lm$fit, x = newvarx, uiw = tmp_lm$se.fit*se_ci_scalar, lwd=1, col="cyan", add=TRUE, pch=NA)
lines(y = tmp_lm$fit, x = newvarx, lwd=1, col="cyan")
#mtext(expression(paste("Soil MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
legend("topright", legend=c("Highly disturbed", "Minimal disturbed"), col=c("red","cyan"), pch=c(NA,NA), lty=c(1,1), lwd=1.4, bty = "n", cex=1.2)
dev.off()

png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_wood_soil_MRT_ratio_map.png",sep=""), height = 2000, width = 3500, res = 300)
par(mfrow=c(2,2), mar=c(0.01,1.5,0.3,6),omi=c(0.01,0.1,0.01,0.15))
var1 = grid_output$MTT_wood_years[,,mid_quant] / grid_output$MTT_som_years[,,mid_quant]
var1 = rast(vals = t(var1[,dim(var1)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
zrange = c(0,max(as.vector(var1), na.rm=TRUE))
#zbreaks = c(zrange[1],(1-zrange[1])*0.25,(1-zrange[1])*0.5,(1-zrange[1])*0.75,1,(zrange[2]-1)*0.25,(zrange[2]-1)*0.5,(zrange[2]-1)*0.75,zrange[2])
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
plot(var1, main="", col=colour_choices_loss,mar=NA, range = zrange, xaxt = "n", yaxt = "n", bty = "n",
     cex.lab=2, cex.main=2.0, cex.axis = 2, axes = FALSE, type="continuous", pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
plot(landmask, add=TRUE, lwd=0.5)
mtext(expression('MTTwood : MTTsom'), side = 3, cex = 1.6, padj = 1.35, adj = 0.5)
var1 = grid_output$mean_wood_gCm2[,,mid_quant] / grid_output$mean_som_gCm2[,,mid_quant]
var1 = rast(vals = t(var1[,dim(var1)[2]:1]), ext = ext(cardamom_ext), crs = crs(cardamom_ext), res=res(cardamom_ext))
zrange = c(0,max(as.vector(var1), na.rm=TRUE))
# legend position
ee = ext(var1) ; e = rep(NA, 4)
e[1] = ee[2] + (abs(diff(ee[1:2]))* 0.027) ; e[2] = e[1] + (abs(diff(ee[1:2]))* 0.027)
e[3] = ee[3] ; e[4] = ee[4]
plot(var1, main="", col=colour_choices_loss,mar=NA, range = zrange, xaxt = "n", yaxt = "n", bty = "n",
     cex.lab=2, cex.main=2.0, cex.axis = 2, axes = FALSE, type="continuous", pax=list(cex.axis=2.0,hadj=0.1), plg = list(ext=e, cex=1.0))
plot(landmask, add=TRUE, lwd=0.5)
mtext(expression('Wood : SOM'), side = 3, cex = 1.6, padj = 1.35, adj = 0.5)

par(mar=c(4,4.5,0.3,3.0))
vary = as.vector(grid_output$MTT_wood_years[,,mid_quant] / grid_output$MTT_som_years[,,mid_quant])
varx = as.vector(grid_output$mean_wood_gCm2[,,mid_quant])
plot(vary ~ varx, main="", ylab="", xlab=expression(paste("Wood (gC",m^{-2},")",sep="")), 
     pch=16, cex=1.0, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
mtext(expression('MTTwood : MTTsom'), side = 2, cex = 1.6, padj = -2.10, adj = 0.5)
abline(1,0,col="grey", lwd=1.5)
vary = as.vector(grid_output$mean_wood_gCm2[,,mid_quant] / grid_output$mean_som_gCm2[,,mid_quant])
varx = as.vector(grid_output$mean_wood_gCm2[,,mid_quant])
plot(vary ~ varx, main="", ylab="", xlab=expression(paste("Wood (gC",m^{-2},")",sep="")), 
     pch=16, cex=1.0, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
mtext(expression('Wood : SOM'), side = 2, cex = 1.6, padj = -2.10, adj = 0.5)
abline(1,0,col="grey", lwd=1.5)
dev.off()

# Temperature
summary(lm(as.vector(grid_output$MTT_foliage_years[,,mid_quant])~as.vector(grid_output$mean_temperature_C)))$adj.r.squared
summary(lm(as.vector(grid_output$MTT_roots_years[,,mid_quant])~as.vector(grid_output$mean_temperature_C)))$adj.r.squared
summary(lm(as.vector(grid_output$MTT_wood_years[,,mid_quant])~as.vector(grid_output$mean_temperature_C)))$adj.r.squared
summary(lm(as.vector(grid_output$MTT_litter_years[,,mid_quant])~as.vector(grid_output$mean_temperature_C)))$adj.r.squared
summary(lm(as.vector(grid_output$MTT_som_years[,,mid_quant])~as.vector(grid_output$mean_temperature_C)))$adj.r.squared
# Precipitation
summary(lm(as.vector(grid_output$MTT_foliage_years[,,mid_quant])~as.vector(grid_output$mean_precipitation_kgH2Om2yr)))$adj.r.squared
summary(lm(as.vector(grid_output$MTT_roots_years[,,mid_quant])~as.vector(grid_output$mean_precipitation_kgH2Om2yr)))$adj.r.squared
summary(lm(as.vector(grid_output$MTT_wood_years[,,mid_quant])~as.vector(grid_output$mean_precipitation_kgH2Om2yr)))$adj.r.squared
summary(lm(as.vector(grid_output$MTT_litter_years[,,mid_quant])~as.vector(grid_output$mean_precipitation_kgH2Om2yr)))$adj.r.squared
summary(lm(as.vector(grid_output$MTT_som_years[,,mid_quant])~as.vector(grid_output$mean_precipitation_kgH2Om2yr)))$adj.r.squared
# Vapour pressure deficit
summary(lm(as.vector(grid_output$MTT_foliage_years[,,mid_quant])~as.vector(grid_output$mean_vpd_Pa)))$adj.r.squared
summary(lm(as.vector(grid_output$MTT_roots_years[,,mid_quant])~as.vector(grid_output$mean_vpd_Pa)))$adj.r.squared
summary(lm(as.vector(grid_output$MTT_wood_years[,,mid_quant])~as.vector(grid_output$mean_vpd_Pa)))$adj.r.squared
summary(lm(as.vector(grid_output$MTT_litter_years[,,mid_quant])~as.vector(grid_output$mean_vpd_Pa)))$adj.r.squared
summary(lm(as.vector(grid_output$MTT_som_years[,,mid_quant])~as.vector(grid_output$mean_vpd_Pa)))$adj.r.squared
# Temperature, precipitation, VPD, harvest fraction, burnt fraction and annual number of fires
lm_fMTT = lm(as.vector(grid_output$MTT_foliage_years[,,mid_quant]) ~ 
           as.vector(grid_output$mean_temperature_C) + 
           as.vector(grid_output$mean_precipitation_kgH2Om2yr) + 
           as.vector(grid_output$mean_vpd_Pa) + 
           as.vector(HarvestFraction) + 
           as.vector(BurnedFraction) + 
           as.vector(FireFreq))
lm_rMTT = lm(as.vector(grid_output$MTT_roots_years[,,mid_quant]) ~ 
           as.vector(grid_output$mean_temperature_C) + 
           as.vector(grid_output$mean_precipitation_kgH2Om2yr) + 
           as.vector(grid_output$mean_vpd_Pa) + 
           as.vector(HarvestFraction) + 
           as.vector(BurnedFraction) + 
           as.vector(FireFreq))
lm_wMTT = lm(as.vector(grid_output$MTT_wood_years[,,mid_quant]) ~ 
           as.vector(grid_output$mean_temperature_C) + 
           as.vector(grid_output$mean_precipitation_kgH2Om2yr) + 
           as.vector(grid_output$mean_vpd_Pa) + 
           as.vector(HarvestFraction) + 
           as.vector(BurnedFraction) + 
           as.vector(FireFreq))
lm_lMTT = lm(as.vector(grid_output$MTT_litter_years[,,mid_quant]) ~ 
           as.vector(grid_output$mean_temperature_C) + 
           as.vector(grid_output$mean_precipitation_kgH2Om2yr) + 
           as.vector(grid_output$mean_vpd_Pa) + 
           as.vector(HarvestFraction) + 
           as.vector(BurnedFraction) + 
           as.vector(FireFreq))
lm_sMTT = lm(as.vector(grid_output$MTT_som_years[,,mid_quant]) ~ 
           as.vector(grid_output$mean_temperature_C) + 
           as.vector(grid_output$mean_precipitation_kgH2Om2yr) + 
           as.vector(grid_output$mean_vpd_Pa) + 
           as.vector(HarvestFraction) + 
           as.vector(BurnedFraction) + 
           as.vector(FireFreq))
summary(lm_fMTT) ; summary(lm_rMTT) ; summary(lm_wMTT) ; summary(lm_lMTT) ; summary(lm_sMTT)
slm_fMTT = step(lm_fMTT, direction = "both") 
slm_rMTT = step(lm_rMTT, direction = "both") 
slm_wMTT = step(lm_wMTT, direction = "both")
slm_lMTT = step(lm_lMTT, direction = "both") 
slm_sMTT = step(lm_sMTT, direction = "both")
summary(slm_fMTT) ; summary(slm_rMTT) ; summary(slm_wMTT) ; summary(slm_lMTT) ; summary(slm_sMTT)
r2_fMTT = summary(slm_fMTT)$adj.r.squared 
r2_rMTT = summary(slm_rMTT)$adj.r.squared 
r2_wMTT = summary(slm_wMTT)$adj.r.squared 
r2_lMTT = summary(slm_lMTT)$adj.r.squared
r2_sMTT = summary(slm_sMTT)$adj.r.squared

# Write out to a simple file
write.table(data.frame(r2_fMTT,r2_rMTT,r2_wMTT,r2_lMTT,r2_sMTT), 
            file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_best_lm_r2_MRTs.csv",sep=""), sep=",",
            row.names=FALSE)

# Plot Foliage, fine root, wood, litter(foliar+fine root+wood?), soil mean residence times against main disturbance
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_MRT_disturbance_association.png",sep=""), height = 2200, width = 4500, res = 300)
par(mfrow=c(3,5), mar=c(4,2,1.4,1), omi = c(0.1,0.2,0.12,0.1))
# Mean annual number of fires
plot(grid_output$MTT_foliage_years[,,mid_quant]~(FireFreq), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
mtext(expression(paste("Foliar MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
#mtext(expression('C1'), side = 2, cex = 1.6, padj = -2.5, adj = 0.5)
plot(grid_output$MTT_roots_years[,,mid_quant]~(FireFreq), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
mtext(expression(paste("Root MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
plot(grid_output$MTT_wood_years[,,mid_quant]~(FireFreq), main="", ylab="", xlab="No. annual fires", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
mtext(expression(paste("Wood MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
plot(grid_output$MTT_litter_years[,,mid_quant]~(FireFreq), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
mtext(expression(paste("Litter MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
plot(grid_output$MTT_som_years[,,mid_quant]~(FireFreq), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
mtext(expression(paste("Soil MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
# Mean annual burned fraction
plot(grid_output$MTT_foliage_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_output$MTT_roots_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_output$MTT_wood_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab="Annual burned fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
plot(grid_output$MTT_litter_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_output$MTT_som_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
# Mean annual forest harvest fraction
plot(grid_output$MTT_foliage_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_output$MTT_roots_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_output$MTT_wood_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab="Annual harvested fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
plot(grid_output$MTT_litter_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_output$MTT_som_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
dev.off()

# Plot Foliage, fine root, wood NPP allocation fractions main meteorology
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NPP_meteorology_association.png",sep=""), height = 2200, width = 2800, res = 300)
par(mfrow=c(3,3), mar=c(4,2,1.4,1), omi = c(0.1,0.2,0.12,0.1))
# Temperature
plot(grid_output$NPP_foliage_fraction[,,mid_quant]~(grid_output$mean_temperature_C), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
mtext(expression("Foliar NPP (0-1)"), cex = 1.3, padj = 0, side = 3)
#mtext(expression('C1'), side = 2, cex = 1.6, padj = -2.5, adj = 0.5)
plot(grid_output$NPP_roots_fraction[,,mid_quant]~(grid_output$mean_temperature_C), main="", ylab="", xlab="Mean Temperature (C)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
mtext(expression("Root NPP (0-1)"), cex = 1.3, padj = 0, side = 3)
plot(grid_output$NPP_wood_fraction[,,mid_quant]~(grid_output$mean_temperature_C), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
mtext(expression("Wood NPP (0-1)"), cex = 1.3, padj = 0, side = 3)
# Precipitation
plot(grid_output$NPP_foliage_fraction[,,mid_quant]~(grid_output$mean_precipitation_kgH2Om2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_output$NPP_roots_fraction[,,mid_quant]~(grid_output$mean_precipitation_kgH2Om2yr), main="", ylab="", xlab="Mean precipitation (mm/yr)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_output$NPP_wood_fraction[,,mid_quant]~(grid_output$mean_precipitation_kgH2Om2yr), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
# Vapour pressure deficit
plot(grid_output$NPP_foliage_fraction[,,mid_quant]~(grid_output$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_output$NPP_roots_fraction[,,mid_quant]~(grid_output$mean_vpd_Pa), main="", ylab="", xlab="Mean VPD (Pa)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_output$NPP_wood_fraction[,,mid_quant]~(grid_output$mean_vpd_Pa), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
dev.off()

# Plot Foliage, fine root, wood NPP allocation fractions against main disturbance
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NPP_disturbance_association.png",sep=""), height = 2200, width = 2800, res = 300)
par(mfrow=c(3,3), mar=c(4,2,1.4,1), omi = c(0.1,0.2,0.12,0.1))
# Temperature
plot(grid_output$NPP_foliage_fraction[,,mid_quant]~(FireFreq), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
mtext(expression("Foliar NPP (0-1)"), cex = 1.3, padj = 0, side = 3)
#mtext(expression('C1'), side = 2, cex = 1.6, padj = -2.5, adj = 0.5)
plot(grid_output$NPP_roots_fraction[,,mid_quant]~(FireFreq), main="", ylab="", xlab="No. annual fires", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
     mtext(expression("Root NPP (0-1)"), cex = 1.3, padj = 0, side = 3)
plot(grid_output$NPP_wood_fraction[,,mid_quant]~(FireFreq), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
     mtext(expression("Wood NPP (0-1)"), cex = 1.3, padj = 0, side = 3)
# Precipitation
plot(grid_output$NPP_foliage_fraction[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_output$NPP_roots_fraction[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab="Annual burned fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_output$NPP_wood_fraction[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
# Vapour pressure deficit
plot(grid_output$NPP_foliage_fraction[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_output$NPP_roots_fraction[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab="Annual harvested fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(grid_output$NPP_wood_fraction[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
dev.off()

# Plot Foliage, fine root, wood NPP allocation fractions main meteorology
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NPP_flux_meteorology_association.png",sep=""), height = 2200, width = 2800, res = 300)
par(mfrow=c(3,3), mar=c(4,2,1.4,1), omi = c(0.1,0.2,0.12,0.1))
# Temperature
plot(as.vector(grid_output$mean_combined_alloc_foliage_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_temperature_C), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
mtext(expression(paste("Foliar NPP (MgC h",a^-1,y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
#mtext(expression('C1'), side = 2, cex = 1.6, padj = -2.5, adj = 0.5)
plot(as.vector(grid_output$mean_alloc_roots_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_temperature_C), main="", ylab="", xlab="Mean Temperature (C)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
mtext(expression(paste("Root NPP (MgC h",a^-1,y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
plot(as.vector(grid_output$mean_alloc_wood_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_temperature_C), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
mtext(expression(paste("Wood NPP (MgC h",a^-1,y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
# Precipitation
plot(as.vector(grid_output$mean_combined_alloc_foliage_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_precipitation_kgH2Om2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(as.vector(grid_output$mean_alloc_roots_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_precipitation_kgH2Om2yr), main="", ylab="", xlab="Mean precipitation (mm/yr)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(as.vector(grid_output$mean_alloc_wood_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_precipitation_kgH2Om2yr), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
# Vapour pressure deficit
plot(as.vector(grid_output$mean_combined_alloc_foliage_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(as.vector(grid_output$mean_alloc_roots_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_vpd_Pa), main="", ylab="", xlab="Mean VPD (Pa)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(as.vector(grid_output$mean_alloc_wood_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_vpd_Pa), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
dev.off()

# Plot Foliage, fine root, wood NPP allocation fractions against main disturbance
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NPP_flux_disturbance_association.png",sep=""), height = 2200, width = 2800, res = 300)
par(mfrow=c(3,3), mar=c(4,2,1.4,1), omi = c(0.1,0.2,0.12,0.1))
# Temperature
plot(as.vector(grid_output$mean_combined_alloc_foliage_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(FireFreq), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
     mtext(expression(paste("Foliar NPP (MgC h",a^-1,y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
#mtext(expression('C1'), side = 2, cex = 1.6, padj = -2.5, adj = 0.5)
plot(as.vector(grid_output$mean_alloc_roots_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(FireFreq), main="", ylab="", xlab="No. annual fires", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
     mtext(expression(paste("Root NPP (MgC h",a^-1,y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
plot(as.vector(grid_output$mean_alloc_wood_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(FireFreq), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
     mtext(expression(paste("Wood NPP (MgC h",a^-1,y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
# Precipitation
plot(as.vector(grid_output$mean_combined_alloc_foliage_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(as.vector(grid_output$mean_alloc_roots_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(BurnedFraction), main="", ylab="", xlab="Annual burned fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(as.vector(grid_output$mean_alloc_wood_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(BurnedFraction), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
# Vapour pressure deficit
plot(as.vector(grid_output$mean_combined_alloc_foliage_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(as.vector(grid_output$mean_alloc_roots_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(HarvestFraction), main="", ylab="", xlab="Annual harvested fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8)
plot(as.vector(grid_output$mean_alloc_wood_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(HarvestFraction), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8)
dev.off()

# Plot Foliage, fine root, wood NPP allocation fractions main meteorology
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NPP_flux_meteorology_association_heatmap.png",sep=""), height = 2200, width = 2800, res = 300)
par(mfrow=c(3,3), mar=c(4,2,1.4,3.8), omi = c(0.1,0.2,0.12,0.1))
fudgeit.leg.lab=""
# Temperature
smoothScatter(as.vector(grid_output$mean_combined_alloc_foliage_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_temperature_C), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(min(grid_output$mean_temperature_C, na.rm=TRUE),max(grid_output$mean_temperature_C, na.rm=TRUE)))
mtext(expression(paste("Foliar NPP (MgC h",a^-1,y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
smoothScatter(as.vector(grid_output$mean_alloc_roots_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_temperature_C), main="", 
     ylab="", xlab="Mean Temperature (C)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(min(grid_output$mean_temperature_C, na.rm=TRUE),max(grid_output$mean_temperature_C, na.rm=TRUE)))
mtext(expression(paste("Root NPP (MgC h",a^-1,y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
fudgeit.leg.lab="Relative Density"
smoothScatter(as.vector(grid_output$mean_alloc_wood_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_temperature_C), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(min(grid_output$mean_temperature_C, na.rm=TRUE),max(grid_output$mean_temperature_C, na.rm=TRUE)))
mtext(expression(paste("Wood NPP (MgC h",a^-1,y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)     
# Precipitation
fudgeit.leg.lab=""
smoothScatter(as.vector(grid_output$mean_combined_alloc_foliage_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_precipitation_kgH2Om2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_output$mean_precipitation_kgH2Om2yr, na.rm=TRUE)))
smoothScatter(as.vector(grid_output$mean_alloc_roots_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_precipitation_kgH2Om2yr), main="", 
     ylab="", xlab="Mean precipitation (mm/yr)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_output$mean_precipitation_kgH2Om2yr, na.rm=TRUE)))
fudgeit.leg.lab="Relative Density"
smoothScatter(as.vector(grid_output$mean_alloc_wood_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_precipitation_kgH2Om2yr), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_output$mean_precipitation_kgH2Om2yr, na.rm=TRUE)))
# Vapour pressure deficit
fudgeit.leg.lab=""
smoothScatter(as.vector(grid_output$mean_combined_alloc_foliage_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_output$mean_vpd_Pa, na.rm=TRUE)))
smoothScatter(as.vector(grid_output$mean_alloc_roots_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_vpd_Pa), main="", ylab="", xlab="Mean VPD (Pa)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_output$mean_vpd_Pa, na.rm=TRUE)))
fudgeit.leg.lab="Relative Density"
smoothScatter(as.vector(grid_output$mean_alloc_wood_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(grid_output$mean_vpd_Pa), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_output$mean_vpd_Pa, na.rm=TRUE)))
dev.off()

# Plot Foliage, fine root, wood NPP allocation fractions against main disturbance
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NPP_flux_disturbance_association_heatmap.png",sep=""), height = 2200, width = 2800, res = 300)
par(mfrow=c(3,3), mar=c(4,2,1.4,3.8), omi = c(0.1,0.2,0.12,0.1))
fudgeit.leg.lab=""
# Temperature
smoothScatter(as.vector(grid_output$mean_combined_alloc_foliage_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(FireFreq), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
mtext(expression(paste("Foliar NPP (MgC h",a^-1,y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
smoothScatter(as.vector(grid_output$mean_alloc_roots_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(FireFreq), main="", ylab="", xlab="No. annual fires", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
mtext(expression(paste("Root NPP (MgC h",a^-1,y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
fudgeit.leg.lab="Relative Density"
smoothScatter(as.vector(grid_output$mean_alloc_wood_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(FireFreq), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
mtext(expression(paste("Wood NPP (MgC h",a^-1,y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
# Precipitation
fudgeit.leg.lab=""
smoothScatter(as.vector(grid_output$mean_combined_alloc_foliage_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
smoothScatter(as.vector(grid_output$mean_alloc_roots_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(BurnedFraction), main="", ylab="", xlab="Annual burned fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
fudgeit.leg.lab="Relative Density"
smoothScatter(as.vector(grid_output$mean_alloc_wood_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(BurnedFraction), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
# Vapour pressure deficit
fudgeit.leg.lab=""
smoothScatter(as.vector(grid_output$mean_combined_alloc_foliage_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
smoothScatter(as.vector(grid_output$mean_alloc_roots_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(HarvestFraction), main="", ylab="", xlab="Annual harvested fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
fudgeit.leg.lab="Relative Density"
smoothScatter(as.vector(grid_output$mean_alloc_wood_gCm2day[,,mid_quant]*1e-2*365.25)~as.vector(HarvestFraction), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
dev.off()

# Plot Foliage, fine root, wood, litter(foliar+fine root+wood?), soil mean residence times against main meteorology
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_MRT_meteorology_association_heatmap.png",sep=""), height = 2200, width = 4500, res = 300)
par(mfrow=c(3,5), mar=c(4,2,1.4,3.8), omi = c(0.1,0.2,0.12,0.1))
fudgeit.leg.lab=""
# Temperature
smoothScatter(grid_output$MTT_foliage_years[,,mid_quant]~(grid_output$mean_temperature_C), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
mtext(expression(paste("Foliar MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
#mtext(expression('C1'), side = 2, cex = 1.6, padj = -2.5, adj = 0.5)
smoothScatter(grid_output$MTT_roots_years[,,mid_quant]~(grid_output$mean_temperature_C), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
mtext(expression(paste("Root MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
smoothScatter(grid_output$MTT_wood_years[,,mid_quant]~(grid_output$mean_temperature_C), main="", ylab="", 
     xlab="Mean Temperature (C)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours,   
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
mtext(expression(paste("Wood MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
smoothScatter(grid_output$MTT_litter_years[,,mid_quant]~(grid_output$mean_temperature_C), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
mtext(expression(paste("Litter MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_output$MTT_som_years[,,mid_quant]~(grid_output$mean_temperature_C), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
mtext(expression(paste("Soil MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
# Precipitation
fudgeit.leg.lab=""
smoothScatter(grid_output$MTT_foliage_years[,,mid_quant]~(grid_output$mean_precipitation_kgH2Om2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
smoothScatter(grid_output$MTT_roots_years[,,mid_quant]~(grid_output$mean_precipitation_kgH2Om2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
smoothScatter(grid_output$MTT_wood_years[,,mid_quant]~(grid_output$mean_precipitation_kgH2Om2yr), main="", ylab="", 
     xlab="Mean precipitation (mm/yr)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
smoothScatter(grid_output$MTT_litter_years[,,mid_quant]~(grid_output$mean_precipitation_kgH2Om2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_output$MTT_som_years[,,mid_quant]~(grid_output$mean_precipitation_kgH2Om2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
# Vapour pressure deficit
fudgeit.leg.lab=""
smoothScatter(grid_output$MTT_foliage_years[,,mid_quant]~(grid_output$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
smoothScatter(grid_output$MTT_roots_years[,,mid_quant]~(grid_output$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
smoothScatter(grid_output$MTT_wood_years[,,mid_quant]~(grid_output$mean_vpd_Pa), main="", ylab="", xlab="Mean VPD (Pa)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
smoothScatter(grid_output$MTT_litter_years[,,mid_quant]~(grid_output$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_output$MTT_som_years[,,mid_quant]~(grid_output$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
dev.off()

# Plot Foliage, fine root, wood, litter(foliar+fine root+wood?), soil mean residence times against main disturbance
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_MRT_disturbance_association_heatmap.png",sep=""), height = 2200, width = 4500, res = 300)
par(mfrow=c(3,5), mar=c(4,2,1.4,3.8), omi = c(0.1,0.2,0.12,0.1))
# Mean annual number of fires
fudgeit.leg.lab=""
smoothScatter(grid_output$MTT_foliage_years[,,mid_quant]~(FireFreq), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500)
mtext(expression(paste("Foliar MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
#mtext(expression('C1'), side = 2, cex = 1.6, padj = -2.5, adj = 0.5)
smoothScatter(grid_output$MTT_roots_years[,,mid_quant]~(FireFreq), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
mtext(expression(paste("Root MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
smoothScatter(grid_output$MTT_wood_years[,,mid_quant]~(FireFreq), main="", ylab="", xlab="No. annual fires", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
mtext(expression(paste("Wood MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
smoothScatter(grid_output$MTT_litter_years[,,mid_quant]~(FireFreq), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
mtext(expression(paste("Litter MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_output$MTT_som_years[,,mid_quant]~(FireFreq), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
mtext(expression(paste("Soil MRT (",y^-1,")",sep="")), cex = 1.3, padj = 0, side = 3)
# Mean annual burned fraction
fudgeit.leg.lab=""
smoothScatter(grid_output$MTT_foliage_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
smoothScatter(grid_output$MTT_roots_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
smoothScatter(grid_output$MTT_wood_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab="Annual burned fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
smoothScatter(grid_output$MTT_litter_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_output$MTT_som_years[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
# Mean annual forest harvest fraction
fudgeit.leg.lab=""
smoothScatter(grid_output$MTT_foliage_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
smoothScatter(grid_output$MTT_roots_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
smoothScatter(grid_output$MTT_wood_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab="Annual harvested fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
smoothScatter(grid_output$MTT_litter_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_output$MTT_som_years[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
dev.off()

# Plot Foliage, fine root, wood NPP allocation fractions main meteorology
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NPP_meteorology_association_heatmap.png",sep=""), height = 2200, width = 2800, res = 300)
par(mfrow=c(3,3), mar=c(4,2,1.4,3.8), omi = c(0.1,0.2,0.12,0.1))
fudgeit.leg.lab=""
# Temperature
smoothScatter(grid_output$NPP_foliage_fraction[,,mid_quant]~(grid_output$mean_temperature_C), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_output$mean_temperature_C, na.rm=TRUE)))
mtext(expression(paste("Foliar NPP (0-1)",sep="")), cex = 1.3, padj = 0, side = 3)
smoothScatter(grid_output$NPP_roots_fraction[,,mid_quant]~(grid_output$mean_temperature_C), main="", 
     ylab="", xlab="Mean Temperature (C)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_output$mean_temperature_C, na.rm=TRUE)))
mtext(expression(paste("Root NPP (0-1)",sep="")), cex = 1.3, padj = 0, side = 3)
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_output$NPP_wood_fraction[,,mid_quant]~(grid_output$mean_temperature_C), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_output$mean_temperature_C, na.rm=TRUE)))
mtext(expression(paste("Wood NPP (0-1)",sep="")), cex = 1.3, padj = 0, side = 3)
# Precipitation
fudgeit.leg.lab=""
smoothScatter(grid_output$NPP_foliage_fraction[,,mid_quant]~(grid_output$mean_precipitation_kgH2Om2yr), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_output$mean_precipitation_kgH2Om2yr, na.rm=TRUE)))
smoothScatter(grid_output$NPP_roots_fraction[,,mid_quant]~(grid_output$mean_precipitation_kgH2Om2yr), main="", 
     ylab="", xlab="Mean precipitation (mm/yr)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_output$mean_precipitation_kgH2Om2yr, na.rm=TRUE)))
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_output$NPP_wood_fraction[,,mid_quant]~(grid_output$mean_precipitation_kgH2Om2yr), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_output$mean_precipitation_kgH2Om2yr, na.rm=TRUE)))
# Vapour pressure deficit
fudgeit.leg.lab=""
smoothScatter(grid_output$NPP_foliage_fraction[,,mid_quant]~(grid_output$mean_vpd_Pa), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_output$mean_vpd_Pa, na.rm=TRUE)))
smoothScatter(grid_output$NPP_roots_fraction[,,mid_quant]~(grid_output$mean_vpd_Pa), main="", ylab="", xlab="Mean VPD (Pa)", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_output$mean_vpd_Pa, na.rm=TRUE)))
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_output$NPP_wood_fraction[,,mid_quant]~(grid_output$mean_vpd_Pa), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(grid_output$mean_vpd_Pa, na.rm=TRUE)))
dev.off()

# Plot Foliage, fine root, wood NPP allocation fractions against main disturbance
png(file = paste(out_dir,"/",gsub("%","_",PROJECT$name),"_NPP_disturbance_association_heatmap.png",sep=""), height = 2200, width = 2800, res = 300)
par(mfrow=c(3,3), mar=c(4,2,1.4,3.8), omi = c(0.1,0.2,0.1,0.1))
fudgeit.leg.lab=""
# Temperature
smoothScatter(grid_output$NPP_foliage_fraction[,,mid_quant]~(FireFreq), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
mtext(expression(paste("Foliar NPP (0-1)",sep="")), cex = 1.3, padj = 0, side = 3)
smoothScatter(grid_output$NPP_roots_fraction[,,mid_quant]~(FireFreq), main="", ylab="", xlab="No. annual fires", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
fudgeit.leg.lab="Relative Density"
mtext(expression(paste("Root NPP (0-1)",sep="")), cex = 1.3, padj = 0, side = 3)
smoothScatter(grid_output$NPP_wood_fraction[,,mid_quant]~(FireFreq), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(FireFreq, na.rm=TRUE)*1.0))
mtext(expression(paste("Wood NPP (0-1)",sep="")), cex = 1.3, padj = 0, side = 3)
# Precipitation
fudgeit.leg.lab=""
smoothScatter(grid_output$NPP_foliage_fraction[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
smoothScatter(grid_output$NPP_roots_fraction[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab="Annual burned fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_output$NPP_wood_fraction[,,mid_quant]~(BurnedFraction), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(BurnedFraction, na.rm=TRUE)*1.0))
# Vapour pressure deficit
fudgeit.leg.lab=""
smoothScatter(grid_output$NPP_foliage_fraction[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab=" ", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
smoothScatter(grid_output$NPP_roots_fraction[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab="Annual harvested fraction", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8, cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
fudgeit.leg.lab="Relative Density"
smoothScatter(grid_output$NPP_wood_fraction[,,mid_quant]~(HarvestFraction), main="", ylab="", xlab="", 
     pch=16, cex=1.4, cex.lab=1.8, cex.axis = 1.8,cex.main=1.8, transformation = function(x) (x-min(x, na.rm=TRUE)) / diff(range(x, na.rm=TRUE)), 
     colramp=smoothScatter_colours, 
     nrpoints = 0, postPlotHook = fudgeit, nbin = 1500, xlim = c(0,max(HarvestFraction, na.rm=TRUE)*1.0))
dev.off()





#> # PointsOfChange
#> # Load the map object 
#> datain = nc_open("/home/lsmallma/WORK/GREENHOUSE/models/CARDAMOM/useful_landmasks/RECCAP_AfricaSplit_MASK11_Mask.nc")
#> class_longitude = ncvar_get(datain, "longitude")
#> class_latitude = ncvar_get(datain,"latitude")
#> class_map = ncvar_get(datain, "Region_Map")
#> class_names = ncvar_get(datain, "Region_Names")
#> class_values = unique(as.vector(class_map))
#> class_values
# [1] 12  2  4 11  5 10  9  1  7  8  3  6
#> class_names
# [1] "1:North_America"    "2:South_America"    "3:Europe"          
# [4] "4:Northern_Africa"  "5:Southern_Africa"  "6:North_Asia"      
# [7] "7:Central_Asia"     "8:East_Asia"        "9:South_Asia"      
#[10] "10:South_East_Asia" "11:Oceania"         "12:Oceans"         
#> strsplit(class_names,":")

## Write loop of what should be plotted...
