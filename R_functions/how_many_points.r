
###
## Function which determines how many grid cells
## are within the defined box
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

## available_countries, a function to provide a list of the countries which can be specified in the site_name
## to define the CARDAMOM analysis area
available_countries <-function() {

   # Load the shapefile CARDAMOM uses as default to define its land sea mask
   landmask = shapefile("./R_functions/global_map/national_boundaries/ne_10m_admin_0_countries.shx")
   # Extract the list of country names used in the mask
   country_match = factor(landmask$SOVEREIGNT) ; country_match = levels(country_match)
   # For consistency / allowability of using the country name in a file path,
   # remove the spaces
   country_match = gsub(" ","",country_match,fixed=TRUE)

   # return to the user
   return(country_match)

} # end function

## lcm2007_to_ctessel, a function which matches the dominant classifications of the lcm2007 to the appropriate C/D-TESSEL PFT

lcm2007_to_ctessel<- function(input_pft) {
    ## LCM2007 class in order 1-23 = C/D-TESSEL equivalent
    # 1) "Broadleaf forest"        = 5
    # 2) "Needleleaf forest"       = 3
    # 3) "Arable"                  = 1
    # 4) "Improved Grassland"      = 2
    # 5) "Rough Grassland"         = 2
    # 6) "Natural Graassland"      = 2
    # 7) "Calcareous Grassland"    = 2
    # 8) "Acid Grassland"          = 2
    # 9) "Fen Marsh Swamp"         = 13
    #10) "Heather"                 = 2
    #11) "Heather Grassland"       = 2
    #12) "Bog"                     = 13
    #13) "Montane"                 = 9
    #14) "Inland rock"             = 8
    #15) "Saltwater"               = 15
    #16) "Freshwater"              = 14
    #17) "Supra-littoral rock"     = 20
    #18) "Supra-littoral sediment" = 20
    #19) "Littoral rock"           = 20
    #20) "Littoral sediment"       = 20
    #21) "Saltmarsh"               = 20
    #22) "Urban"                   = 19
    #23) "Suburban"                = 19

    # vector of corresponding C/D-TESSEL PFTs in order of the LCM2007 types
    tessel_types=c(5,3,1,2,2,2,2,2,13,2,2,13,9,8,15,14,20,20,20,20,20,19,19)
    # use input LCM2007 cover type to select and return the ctessel PFT
    lcm2007_to_ctessel=tessel_types[input_pft]
    # if location does not have a pft in the lcm2007 make 0 and this will use default values from ECMWF
    if (input_pft == 0) { lcm2007_to_ctessel = 0 }
    # now return needef value
    return(lcm2007_to_ctessel)
}

## Use byte compile
lcm2007_to_ctessel<-cmpfun(lcm2007_to_ctessel)

# corine2006_to_ctessel, a function which matches the dominant classifications of the corine2006 to the appropriate C/D-TESSEL PFT
corine2006_to_ctessel<- function(input_pft) {
    ## Corine2006 class in order 1-43 = C/D-TESSEL equivalent
    # 1) "Continuous urban fabric"         = 0
    # 2) "Discontinuous urban fabric"      = 19
    # 3) "Industrial or urban"             = 0
    # 4) "Road or Rail + associated"       = 0
    # 5) "Port Areas"                      = 0
    # 6) "Airports"                        = 0
    # 7) "Mineral extraction"              = 0
    # 8) "Dump sites"                      = 0
    # 9) "Construction sites"              = 0
    #10) "Green urban areas"               = 19
    #11) "Sport / Leisure facilities"      = 0
    #12) "Non-irrigated arable"            = 1
    #13) "Irrigated arable"                = 10
    #14) "Rice fields"                     = 10
    #15) "Vineyards"                       = 17
    #16) "Fruit tree/berry plantation"     = 18
    #17) "Olive groves"                    = 18
    #18) "Pastures"                        = 2
    #19) "Annual crops + fixed associated" = 1
    #20) "Complex cultivation patterns"    = 1
    #21) "Agriculture with signif natural" = 1
    #22) "Agro-forest areas"               = 3
    #23) "Broadleaf forest"                = 5
    #24) "Coniferous forest"               = 3
    #25) "Mixed forest"                    = 18
    #26) "Natural grassland"               = 2
    #27) "Moors and heathland"             = 2
    #28) "Sclerophyllous veg"              = 17
    #29) "Transitional wood-shrub"         = 19
    #30) "Beaches, Dune, sands"            = 20
    #31) "Bare Rock"                       = 8
    #32) "Sparsely vegetated"              = 11
    #33) "Burnt areas"                     = 0
    #34) "Glaciers and snow"               = 12
    #35) "Inland marsh"                    = 13
    #36) "Peat bog"                        = 13
    #37) "Salt marshes"                    = 20
    #38) "Salines"                         = 20
    #39) "Inter-tidal flats"               = 20
    #40) "Water courses"                   = 14
    #41) "Water bodies"                    = 14
    #42) "Coastal lagoons"                 = 14
    #43) "Estaries"                        = 20

    # vector of corresponding C/D-TESSEL PFTs in order of the Corine2006 types
    tessel_types=c(0,19,0,0,0,0,0,0,0,19,0,1,10,10,17,18,18,2,1,1,1,3,5,3,18,2,2,17,19,20,8,11,0,152,13,13,20,20,20,14,14,14,20)
    # use input Corine2006 cover type to select and return the ctessel PFT
    corine2006_to_ctessel=tessel_types[input_pft]
    # if location does not have a pft in the corine make 0 and this will use default values from ECMWF
    if (input_pft==0) {corine2006_to_ctessel=0}
    # now return needef value
    return(corine2006_to_ctessel)
}
## Use byte compile
corine2006_to_ctessel<-cmpfun(corine2006_to_ctessel)

#lat = sites_cardamom_lat ; long = sites_cardamom_long ; resolution = cardamom_resolution ; grid_type = cardamom_grid_type ; sitename = sites_cardamom
how_many_points<- function (path_to_landsea,lat,long,resolution,grid_type,sitename) {

    # check input data
    if (length(which(long > 180)) > 0) {stop("Long should be -180 to +180")}

    # generate UK or WGS-84 lat long grid
    if (grid_type == "UK") {
        output = generate_uk_grid(lat,long,resolution)
        area = array(PROJECT$resolution**2, dim=c(PROJECT$long_dim,PROJECT$lat_dim))
    } else if (grid_type=="wgs84") {
        output = generate_wgs84_grid(lat,long,resolution)
        grid_lat = array(output$lat, dim=c(output$long_dim,output$lat_dim))
        grid_long = array(output$long,dim=c(output$long_dim,output$lat_dim))
        # then generate the area estimates for each pixel
        area = calc_pixel_area(grid_long,grid_lat)
    } else {
        stop('have selected invalid grid type, the valid options are "UK" and "wgs84"')
    }
    # Extract useful information (output re-used later)
    lat = output$lat ; long = output$long
    lat_dim = output$lat_dim ; long_dim = output$long_dim
    cardamom_ext = output$cardamom_ext

    # now work out how many of these are land points
    # determine whether using LCM 2007 or default ECMWF land cover map
    if (use_lcm == "LCM2007") {
#        data2=nc_open("/home/lsmallma/WORK/GREENHOUSE/LCM2007/LCM2007_with_lat_long.nc")
#        lcm=ncvar_get(data2,"LCM2007")
        lcm = raster("/home/lsmallma/WORK/GREENHOUSE/LCM2007/Download_lcm2007_143707/lcm-2007-1km_397874/dominant_target_class/LCM2007_GB_1K_Dominant_TargetClass.tif")
        # Reproject onto the WGS84 grid
        lcm = projectRaster(lcm, ext = cardamom_ext, crs = CRS("+init=epsg:4326"), method = "ngb")
        # Aggregate to approximately the right resolution
        if (grid_type == "UK") {
            target_ratio = max(0.1666667,(0.001*(resolution/111))) / res(lcm)
        } else {
            target_ratio = max(0.1666667,resolution) / res(lcm)
        }
        agg_fun = function(pixels, na.rm) {
           if ((length(which(pixels > 0))/length(pixels)) > 0.2) {
               return(modal(pixels[pixels > 0], na.rm=na.rm))
           } else {
               return(0)
           }
        }
        lcm = aggregate(lcm, fact = floor(target_ratio), fun = agg_fun)
        # Extract lat / long
        lat_lcm = coordinates(lcm)
        # Convert into arrays
        long_lcm = array(lat_lcm[,1], dim=c(dim(lcm)[2],dim(lcm)[1])) ; lat_lcm = array(lat_lcm[,2], dim=c(dim(lcm)[2],dim(lcm)[1]))
        lcm = array(lcm, dim=c(dim(lcm)[2],dim(lcm)[1]))
        # Now flip the lat dimension to get it the right way
        lat_lcm = lat_lcm[,dim(lat_lcm)[2]:1] ; long_lcm = long_lcm[,dim(long_lcm)[2]:1] ; lcm = lcm[,dim(lcm)[2]:1]
    } else if (use_lcm == "CORINE2006") {
        data2=nc_open("/home/lsmallma/WORK/GREENHOUSE/Corine_lcm/Corine2006_at250m_with_lat_long.nc")
        lcm=ncvar_get(data2,"Corine2006")
    } else if (use_lcm == "CORINE2006_1km") {
        data2=nc_open("/home/lsmallma/WORK/GREENHOUSE/Corine_lcm/Corine2006_at1km_with_lat_long.nc")
        lcm=ncvar_get(data2,"Corine2006")
    } else if (use_lcm == "ECMWF") {
        # load global surfclim file and info file for surfclim
        data2=nc_open("./R_functions/global_map/ECMWF/surfclim_all.nc")
        # extract high vegetation cover fraction
        hi_veg_frac=ncvar_get(data2, "cvh")
        # extract low vegetation cover fraction
        low_veg_frac=ncvar_get(data2, "cvl")
        # extract high vegetation type
        hi_veg_type=ncvar_get(data2, "tvh")
        # extract low vegetation type
        low_veg_type=ncvar_get(data2, "tvl")
        hi_veg_frac=as.vector(hi_veg_frac) ; low_veg_frac=as.vector(low_veg_frac)
        hi_veg_type=as.vector(hi_veg_type) ; low_veg_type=as.vector(low_veg_type)
        lcm=hi_veg_type ; lcm[which(low_veg_frac > hi_veg_frac)]=low_veg_type[which(low_veg_frac > hi_veg_frac)]
        lat_lcm=ncvar_get(data2, "latitude")
        long_lcm=ncvar_get(data2, "longitude")
        # restructure to 2-D array which matches the actual data structure...
        lat_tmp = length(lat_lcm) ; long_tmp = length(long_lcm)
        lat_lcm = array(rep(lat_lcm, each = long_tmp), dim=c(long_tmp,lat_tmp))
        long_lcm = array(long_lcm, dim=c(long_tmp,lat_tmp))
        long_lcm[which(long_lcm > 180)] = long_lcm[which(long_lcm > 180)]-360
        lcm = array(lcm, dim=c(dim(lat_lcm)[1],dim(lat_lcm)[2]))
    } else {
        stop("no land cover option found / set")
    }
    # download location data
    if (use_lcm != "ECMWF" & use_lcm != "LCM2007") {
        lat_lcm=ncvar_get(data2,"lat")
        long_lcm=ncvar_get(data2,"long")
    }
    # house keeping
    if (exists("data2")) {nc_close(data2)}

    # raw total pixels
    print(paste("Raw pixel total is ",length(lat)," next filter for water bodies"))

    # find locations
    if (use_parallel) {
        cl <- makeCluster(numWorkers, type = "PSOCK")
        # load R libraries in cluster
        clusterExport(cl,"load_r_libraries") ; clusterEvalQ(cl, load_r_libraries())
        output=parLapply(cl,1:length(lat),fun=closest2d_2,lat=lat_lcm,long=long_lcm,lat_in=lat,long_in=long)
        stopCluster(cl)
        # extract the i,j values seperately
        output_i=unlist(output,use.names=FALSE)[which((1:length(unlist(output, use.names = FALSE))*0.5) != floor(1:length(unlist(output, use.names=FALSE))*0.5))]
        output_j=unlist(output,use.names=FALSE)[which((1:length(unlist(output, use.names = FALSE))*0.5) == floor(1:length(unlist(output, use.names=FALSE))*0.5))]
     } else {
        output=lapply(1:length(lat),FUN=closest2d_2,lat=lat_lcm,long=long_lcm,lat_in=lat,long_in=long)
        # extract the i,j values seperately
        output_i=unlist(output, use.names=FALSE)[which((1:length(unlist(output, use.names=FALSE))*0.5) != floor(1:length(unlist(output, use.names=FALSE))*0.5))]
        output_j=unlist(output, use.names=FALSE)[which((1:length(unlist(output, use.names=FALSE))*0.5) == floor(1:length(unlist(output, use.names=FALSE))*0.5))]
    } # parallel or not

    # Inform the user
    print("Generating land sea mask")

    if (path_to_landsea == "default") {
        # load global shape file for land sea mask
        landmask = shapefile("./R_functions/global_map/national_boundaries/ne_10m_admin_0_countries.shx")
        # just to be sure enforce the projection to WGS-84
        landmask = spTransform(landmask,CRS("+init=epsg:4326"))
        # Clip to the extent of the CARDAMOM analysis
        landmask = crop(landmask, cardamom_ext)

        # create raster, passing the raster values corresponding to the sovereign state
        # NOTE: the actual value assigned is linked the factor levels
        landsea = rasterize(landmask,cardamom_ext,factor(landmask$SOVEREIGNT), fun = "last")
        landsea_frac = rasterize(landmask,cardamom_ext,factor(landmask$SOVEREIGNT), fun = "last", getCover=TRUE)

        # Sometimes we want to simulate a particular country, which we will check now...
        country_match = factor(landmask$SOVEREIGNT) ; country_match = levels(country_match)
        country_match = gsub(" ","",country_match,fixed=TRUE)
        #sitename =  "UnitedKingdom"
        # does our site name (as specified in the grid verison of analysis) correspond to a country name as
        # given in the land mask we are using...?
        if (length(which(grepl(sitename,country_match) == TRUE)) > 0 & select_country) {
            # if so then loop through the land areas which fall within the correct country
            country_match = which(grepl(sitename,country_match) == TRUE)
            keep = rep(0,length(landsea))
            for (i in seq(1, length(country_match))) {
                 keep[as.vector(landsea) == country_match[i]] = 1
            }
        } else {
            # otherwise just assume we are interested in all land areas...
            keep = rep(0,length(landsea))
            keep[is.na(as.vector(landsea)) == FALSE] = 1
        } # country or all land area filter?
        # Set non country areas to NA, and all other to 1
        landsea[keep == 0] = NA
        # Add a buffer based on the land sea fraction to avoid missing land area we want
        landsea_frac_buffer = boundaries(landsea, type="outer")*landsea_frac
        # Set all actual data to 1
        landsea[as.vector(landsea) > 0] = 1
        # set missing data to 0
        landsea[is.na(as.vector(landsea))] = 0
        # Now combine the maps, giving a complete landsea fractional map
        landsea = (landsea*landsea_frac) + landsea_frac_buffer
        # Reset any newly created NaN from the merge
        landsea[is.na(as.vector(landsea))] = 0

    } else {

        # Assume that we have been given a geotiff file where the presence of a value > 0  should be included in the masked area
        landsea = raster(path_to_landsea)
        # just to be sure enforce the projection to WGS-84
        target = raster(crs = crs(cardamom_ext), ext = extent(landsea), resolution = res(cardamom_ext))
        # Resample to correct grid
        landsea = resample(landsea, target, method="ngb", na.rm=TRUE)
        # Clip to the extent of the CARDAMOM analysis
        landsea = crop(landsea, cardamom_ext)

        # Assume all positive values are to be included
        landsea[as.vector(landsea) > 0] = 1
        # Assume everywhere else is not to be included
        landsea[as.vector(landsea) != 1] = 0
        landsea[is.na(as.vector(landsea))] = 0

    } # default landsea mask

    # trim to the actual data area
    #landsea = trim(landsea, padding = 3)
    # extract lat/long information for the raster version
    landsea_long = coordinates(landsea)[,1] ; landsea_lat = coordinates(landsea)[,2]
    # arrange them into the correct lat / long orientations
    landsea_dim = dim(landsea) ; landsea = as.vector(landsea)

    # find locations from the landsea mask which correspond with overall grid defined in the control file
    if (use_parallel) {
        cl <- makeCluster(numWorkers, type = "PSOCK")
        # load R libraries in cluster
        clusterExport(cl,"load_r_libraries") ; clusterEvalQ(cl, load_r_libraries())
        output = parLapply(cl,1:length(lat),fun=closest2d_1,lat=landsea_lat,long=landsea_long,lat_in=lat,long_in=long)
        stopCluster(cl)
        # extract the i,j values seperately
        output_k = unlist(output, use.names=FALSE)
     } else {
        output = lapply(1:length(lat),FUN=closest2d_1,lat=landsea_lat,long=landsea_long,lat_in=lat,long_in=long)
        # extract the i,j values seperately
        output_k = unlist(output, use.names=FALSE)
    }

    # selecting only these areas of the landsea mask
#    landsea = as.vector(landsea)[output_k]

    # Check against our plant functional type / land cover maps
    # and determine which location want to keep based on land cover etc.
    remove = 0 ; pft_keep = 0
    # now iterate through the sites
    for (pft in seq(1, length(lat))) {
         # update the user, but only sometimes
         if (pft%%2000 == 0 | pft < 500) {print(paste("Ocean filter ",round((pft/length(lat))*100,0),"% complete",sep=""))}
         # convert incoming pft to common values (in this case CTESSEL)
         if (use_lcm == "LCM2007") {
             new_pft = lcm2007_to_ctessel(lcm[output_i[pft],output_j[pft]])
         } else if (use_lcm == "CORINE2006") {
             new_pft = corine2006_to_ctessel(lcm[output_i[pft],output_j[pft]])
         } else if (use_lcm == "CORINE2006_1km") {
             new_pft = corine2006_to_ctessel(lcm[output_i[pft],output_j[pft]])
         } else if (use_lcm == "forestry_commission" | use_lcm == "forestry_commission_LCM2007" | use_lcm == "forestry_commission_public_private") {
             new_pft = lcm[output_i[pft],output_j[pft]]
             if (new_pft < 0 | length(new_pft) == 0) {new_pft = 0}
         } else if (use_lcm == "ECMWF") {
             new_pft = lcm[output_i[pft],output_j[pft]]
         }
         # now exclude if not a land site
         if (new_pft == 0 | new_pft == 14 | new_pft == 15 | landsea[output_k[pft]] < 0.5) {
             remove = append(remove,pft)
         } else {
             pft_keep = append(pft_keep,new_pft)
         }
    } # sites for loop

    # remove initial value
    pft_keep = pft_keep[-1] ; remove = remove[-1] ; rm(lat_lcm,long_lcm)

    # generate the site names prior to removing undesired locations to ensure consistent naming
    b = 1 ; sites = rep("NA",times=(length(lat)))
    for (n in seq(1, length(lat))) {
         if (n%%2000 == 0 | n < 500) {print(paste("Have generated ",round((b/(length(lat)+length(remove)))*100,0),"% of site IDs" ,sep=""))}
         # we want the numbers to match their location within the domain not of the land pixels.
         # this is needed for easy reconstruction later
         sites[b] = sprintf('%05i',n) ; b = b+1
    }
    # now sub-select for the ones we want
    if (length(remove) > 0) {lat = lat[-remove] ; long = long[-remove] ; sites = sites[-remove]}
    # Inform the user of the number of pixels
    print(paste("In total there are ",length(sites)," land pixels to run",sep=""))

    # re-arrange landsea mask so that it matches with the actual grid
    landsea = as.vector(array(landsea, dim=c(long_dim,lat_dim))[,lat_dim:1])

    # Combine outputs
    output = list(nosites=length(lat),waterpixels=remove,landsea=landsea,ctessel_pft=pft_keep,lat_dim=lat_dim,long_dim=long_dim,sites=sites)
    # If the grid area has been calculated we will keep this too
    if (exists("area")) {
        output$area_m2 = area
    }

    # clean up
    gc(reset=TRUE, verbose=FALSE)
    # Return back to user
    return(output)
} # end function how_many_points

## Use byte compile
how_many_points<-cmpfun(how_many_points)
