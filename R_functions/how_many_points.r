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
# Functions which determines how many grid cells
# are within the defined project domain. Also, some which translate 
# between various plant functional type classifications
# 
# Author: T. Luke Smallman (02/05/2024)
#
#########################################################################################

available_countries <-function(cardamom_dir) {

   ## Function description
   ## available_countries(), provides a list of the countries 
   ## which can be specified in the site_name variable.
   ## Can define a more specific CARDAMOM analysis area

   if (missing(cardamom_dir)) { 
       # Load the shapefile CARDAMOM uses as default to define its land sea mask
       landmask = vect("./R_functions/global_map/national_boundaries/ne_10m_admin_0_countries.shx")
   } else {
       # Load the shapefile CARDAMOM uses as default to define its land sea mask
       landmask = vect(paste(cardamom_dir,"/R_functions/global_map/national_boundaries/ne_10m_admin_0_countries.shx",sep=""))
   }
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
how_many_points<- function(path_to_landsea,lat,long,resolution,grid_type,sitename) {

    # Spatial grid
    output = generate_grid(cardamom_grid_type,lat,long,resolution)
    # Structure area (in m) as a grid
    area = output$area
    # x ~ y coordinate information as a grid
    lat = output$lat ; long = output$long
    # x ~ y dimensions
    lat_dim = output$lat_dim ; long_dim = output$long_dim
    # Empty raster template of the cardamom domain
    cardamom_ext = output$cardamom_ext
    # Tidy
    rm(output)

    # now work out how many of these are land points
    # Determine the land cover map used
    if (use_lcm == "ECMWF") {

        # load global surfclim file and info file for surfclim
        data2 = nc_open("./R_functions/global_map/ECMWF/surfclim_all.nc")
        # extract high vegetation cover fraction
        hi_veg_frac = ncvar_get(data2, "cvh")
        # extract low vegetation cover fraction
        low_veg_frac = ncvar_get(data2, "cvl")
        # extract high vegetation type
        hi_veg_type = ncvar_get(data2, "tvh")
        # extract low vegetation type
        low_veg_type = ncvar_get(data2, "tvl")
        hi_veg_frac = as.vector(hi_veg_frac) ; low_veg_frac = as.vector(low_veg_frac)
        hi_veg_type = as.vector(hi_veg_type) ; low_veg_type = as.vector(low_veg_type)
        # Assign dominant cover type
        lcm = hi_veg_type ; lcm[which(low_veg_frac > hi_veg_frac)] = low_veg_type[which(low_veg_frac > hi_veg_frac)]
        # Read in the location information
        lat_lcm = ncvar_get(data2, "latitude") ; long_lcm = ncvar_get(data2, "longitude")
        # Extract dimensional information
        ydim = length(lat_lcm) ; xdim = length(long_lcm)
        # Create 2-D array of the latitude and longitude 1-D
        lat_lcm = array(rep(lat_lcm, each = xdim), dim=c(xdim,ydim))
        long_lcm = array(long_lcm, dim=c(xdim,ydim))
        # Change longitude from 0-360 to -180->180
        long_lcm[which(long_lcm > 180)] = long_lcm[which(long_lcm > 180)]-360
        # Create 2-d array of the land cover map itself
        lcm = array(lcm, dim=c(dim(lat_lcm)[1],dim(lat_lcm)[2]))

        # We now must create a raster which we can adjust the projection grid 
        # to match that requested by the current analysis
        lcm = data.frame(x = as.vector(long_lcm), y = as.vector(lat_lcm), z = as.vector(lcm))
        lcm = rast(lcm, crs = ("epsg:4326"), type="xyz")

        # If we have an epsg then we want to know if it differs from the one desired by the analysis
        epsg = crs(lcm, describe = TRUE)$code
        if (epsg != gsub("epsg:","",cardamom_grid_type)) {
            # Ensure that the extent of the input object is consistent 
            # with the possible extent of the selected epsg
            lcm = crop(lcm, ext(unlist(crs(cardamom_grid_type, describe=TRUE)$extent)))        
            # If it does not match we need to reproject it
            lcm = project(lcm, cardamom_grid_type, method="near", align = FALSE) ; gc()
        }

        # extract dimension information for the grid, 
        # note the axis switching between raster and actual array
        xdim = dim(lcm)[2] ; ydim = dim(lcm)[1]
        # extract the lat / long information needed
        long_lcm = crds(lcm, df=TRUE, na.rm=FALSE)
        lat_lcm  = long_lcm$y ; long_lcm = long_lcm$x
        # restructure into correct orientation
        long_lcm = array(long_lcm, dim=c(xdim,ydim))
        lat_lcm = array(lat_lcm, dim=c(xdim,ydim))
        # Flip for human viewing in R
        long_lcm = long_lcm[,ydim:1] ; lat_lcm = lat_lcm[,ydim:1]

        # Extract the fractional cover, 
        # assumed to be unchanged for any given year.
        lcm = array(values(lcm), dim=c(xdim,ydim))
        # Flip the north-south axis to invert for human eyes
        lcm = lcm[,ydim:1]

    } else if (use_lcm != "ECMWF" & use_lcm != "") {
        
        # In which case we assume that we have been directed towards a geotif with 0 or NA indicating
        # areas not to include in the analysis and 1 to indicate areas to include
        lcm = rast(use_lcm)

        # Extract the epsg from the file
        epsg = crs(lcm, describe = TRUE)$code
        if (is.null(epsg) | epsg == "") { stop(paste("the use_lcm specification leads to a geotif which does not contain epsg information."))}
        # If we have an epsg then we want to know if it differs from the one desired by the analysis
        if (epsg != gsub("epsg:","",cardamom_grid_type)) {
            # Ensure that the extent of the input object is consistent 
            # with the possible extent of the selected epsg
            lcm = crop(lcm, ext(unlist(crs(cardamom_grid_type, describe=TRUE)$extent)))                
            # If it does not match we need to reproject it
            lcm = project(lcm, cardamom_grid_type, method="near", align = FALSE) ; gc()
        }
        # Extend the extent of the overall grid to the analysis domain
        lcm = extend(lcm,cardamom_ext)
        # Trim the extent of the overall grid to the analysis domain
        lcm = crop(lcm,cardamom_ext) 
        # Aggregate based on area average
        lcm = resample(lcm, cardamom_ext, method="average") ; gc() 

        # extract dimension information for the grid, 
        # note the axis switching between raster and actual array
        xdim = dim(lcm)[2] ; ydim = dim(lcm)[1]
        # extract the lat / long information needed
        long_lcm = crds(lcm, df=TRUE, na.rm=FALSE)
        lat_lcm  = long_lcm$y ; long_lcm = long_lcm$x
        # restructure into correct orientation
        long_lcm = array(long_lcm, dim=c(xdim,ydim))
        lat_lcm = array(lat_lcm, dim=c(xdim,ydim))
        long_lcm = long_lcm[,ydim:1] ; lat_lcm = lat_lcm[,ydim:1]

        # Extract the fractional cover, 
        # assumed to be unchanged for any given year.
        lcm = array(values(lcm), dim=c(xdim,ydim))
        # Flip the north-south axis to invert for human eyes
        lcm = lcm[,ydim:1]

    } else {
        stop("no land cover option found / set")
    }

    # Clear memory and the contents of the temporary files at this time
    gc() ; tmpFiles(current=TRUE, orphan=TRUE, old=FALSE, remove=TRUE) 

    # raw total pixels
    print(paste("Raw pixel total is ",length(lat)," next filter for water bodies"))

    # In either case we will now make use of terra library functions to determine
    # where in the land cover map (lcm) is the locations in the CARDAMOM domain pixels
    tiff_in = data.frame(x = as.vector(long_lcm), y = as.vector(lat_lcm), z = as.vector(lcm))
    tiff_in = rast(tiff_in, crs = cardamom_grid_type, type="xyz")

    # Find locations
    output = closest2d_tif(tiff_in,lat,long)
    output_i = output$i_loc ; output_j = output$j_loc
    rm(tiff_in,output)

    # Inform the user
    print("Generating land sea mask")

    if (path_to_landsea == "default") {

        # load global shape file for land sea mask
        landmask = vect("./R_functions/global_map/national_boundaries/ne_10m_admin_0_countries.shx")
        # Check the EPSG code
        epsg = crs(landmask, describe = TRUE)$code
        if (is.null(epsg) | epsg == "") { stop(paste("the default landmask leads to a geotif which does not contain epsg information."))}
        # If we have an epsg then we want to know if it differs from the one desired by the analysis
        if (epsg != gsub("epsg:","",cardamom_grid_type)) {
            # Ensure that the extent of the input object is consistent 
            # with the possible extent of the selected epsg
            landmask = crop(landmask, ext(unlist(crs(cardamom_grid_type, describe=TRUE)$extent)))                
            # If it does not match we need to reproject it
            landmask = project(landmask, cardamom_grid_type) ; gc()
        }

        # Clip to the extent of the CARDAMOM analysis
        landmask = crop(landmask, cardamom_ext)

        # create raster, passing the raster values corresponding to the sovereign state
        # NOTE: the actual value assigned is linked the factor levels
        landsea = rasterize(landmask,cardamom_ext,factor(landmask$SOVEREIGNT), fun = "max")
        landsea_frac = rasterize(landmask,cardamom_ext,factor(landmask$SOVEREIGNT), fun = "max", cover=TRUE)

        # Sometimes we want to simulate a particular country, which we will check now...
        country_match = factor(landmask$SOVEREIGNT) ; country_match = levels(country_match)
        country_match = gsub(" ","",country_match,fixed=TRUE)
        #sitename =  "UnitedKingdom"
        # does our site name (as specified in the grid verison of analysis) correspond to a country name as
        # given in the land mask we are using...?
        if (length(which(grepl(sitename,country_match) == TRUE)) > 0 & select_country) {
            # if so then loop through the land areas which fall within the correct country
            country_match_loc = which(grepl(sitename,country_match) == TRUE)
            for (i in seq(1, length(country_match_loc))) {
                 landsea[which(as.vector(landsea) == country_match_loc[i])] = -1
            }
            landsea[which(as.vector(landsea) > 0)] = 0
            landsea[which(as.vector(landsea) == -1)] = 1
        } else {
            # otherwise just assume we are interested in all land areas...
            landsea[landsea > 0] = 1
        } # country or all land area filter?
        # Add a buffer based on the land sea fraction to avoid missing land area we want
        landsea_frac_buffer = boundaries(landsea, inner=FALSE)*landsea_frac
        # Set all actual data to 1
        landsea[landsea > 0] = 1
        # set missing data to 0
        landsea[is.na(landsea)] = 0
        # Now combine the maps, giving a complete landsea fractional map
        landsea = (landsea*landsea_frac) + landsea_frac_buffer
        # Reset any newly created NaN from the merge
        landsea[is.na(landsea)] = 0

        # Set the threshold below which we assume that the pixel will be excluded
        #cover_threshold = 0.5
        cover_threshold = 0.1

    } else {

        # Assume that we have been given a geotiff file where 
        # the presence of a value > 0  should be included in the masked area
        landsea = rast(path_to_landsea)
        # Check that all values rage between zero and 1
        tmp = max(values(landsea), na.rm=TRUE)
        if (tmp > 1) {
            print("The path_to_landsea file is expected to be a fractional cover, however, the maximum value is greater than 1")
            print("The file will instead be treated as a binary presence map, i.e. > 0 is present and set to 1.")
            landsea[landsea > 0] = 1 ; rm(tmp)
        }
        # Extract the epsg from the file
        epsg = crs(landsea, describe = TRUE)$code
        if (is.null(epsg) | epsg == "") { stop(paste("the use_lcm specification leads to a geotif which does not contain epsg information."))}
        # If we have an epsg then we want to know if it differs from the one desired by the analysis
        if (epsg != gsub("epsg:","",cardamom_grid_type)) {
            # Ensure that the extent of the input object is consistent 
            # with the possible extent of the selected epsg
            landsea = crop(landsea, ext(unlist(crs(cardamom_grid_type, describe=TRUE)$extent)))              
            # If it does not match we need to reproject it
            landsea = project(landsea, cardamom_grid_type, method="near", align = FALSE) ; gc()
        }

        # Ensure the extents of the landsea mask matches the CARDAMOM analysis
        landsea = crop(landsea, cardamom_ext)
        landsea = extend(landsea, cardamom_ext)
        # Set all NA to 0, to allow for any averaging
        landsea[is.na(landsea)] = 0

        # Create raster with the target resolution
        target = rast(crs = crs(cardamom_ext), ext = ext(landsea), resolution = res(cardamom_ext))
        # Now depending on whether we are in the correct resolution
        if (res(landsea)[1] >= res(cardamom_ext)[1] & res(landsea)[2] >= res(cardamom_ext)[2]) {
            # If the resolution of the land mask is greater or equal to the CARDAMOM analysis we
            # can use nearest neighbour
            landsea = resample(landsea, target, method="near")
        } else { 
            # If the resolution is finer than the CARDAMOM anaysis we will have to use average
            landsea = resample(landsea, target, method="average") ; gc() 
        } # Aggrgeate to resolution

        # Set the threshold below which we assume that the pixel will be excluded
        #cover_threshold = 0.01 # currently, equal to 1 ha, assuming a 1 km grid
        cover_threshold = 0.04 # currently, equal to 4 ha, assuming a 1 km grid
        #cover_threshold = 0.20 # currently, equal to ~the largest third of Improved grassland areas, assuming a 1 km grid.
        #                       # biased, yes, but to compromise on the number of pixels being simulated.
        #cover_threshold = 0.30 # Compromise to get a geographical spread but fewer pixels

    } # default landsea mask

    # Find locations
    output = closest2d_tif(landsea,lat,long)
    output_x = output$i_loc ; output_y = output$j_loc
    rm(output)

    # extract dimension information for the grid, 
    # note the axis switching between raster and actual array
    xdim = dim(landsea)[2] ; ydim = dim(landsea)[1]
    # Turn into array consistent with the rest of the analysis
    landsea = array(values(landsea), dim=c(xdim,ydim))
    # Flip the north-south axis to invert for human eyes
    landsea = landsea[,ydim:1]

    # Check against our plant functional type / land cover maps
    # and determine which location want to keep based on land cover etc.
    remove = 0 ; pft_keep = 0
    # now iterate through the sites
    for (pft in seq(1, length(lat))) {
         # update the user, but only sometimes
         if (pft%%10000 == 0 | pft < 10) {print(paste("Ocean filter ",round((pft/length(lat))*100,0),"% complete",sep=""))}
         # convert incoming pft to common values (in this case CTESSEL)
         if (use_lcm == "ECMWF") {
             new_pft = lcm[output_i[pft],output_j[pft]]
             if (is.na(new_pft)) {new_pft = 0} # if NA values, set to zero
         } else {
             # All other cases assume we should have 0-1
             new_pft = lcm[output_i[pft],output_j[pft]]
             if (new_pft > 0.01) {new_pft = 1} else {new_pft = 0}
         }
         # now exclude if not a land site, i.e. water, rock, ice, urban
         if (new_pft == 0 | new_pft == 8 | new_pft == 12 | 
             new_pft == 14 | new_pft == 15 | 
             landsea[output_x[pft],output_y[pft]] < cover_threshold) {
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
         if (n%%10000 == 0 | n < 10) {print(paste("Have generated ",round((b/(length(lat)+length(remove)))*100,0),"% of site IDs" ,sep=""))}
         # we want the numbers to match their location within the domain not of the land pixels.
         # this is needed for easy reconstruction later
         sites[b] = sprintf('%05i',n) ; b = b+1
    }
    # now sub-select for the ones we want
    if (length(remove) > 0) {lat = lat[-remove] ; long = long[-remove] ; sites = sites[-remove]}
    # Inform the user of the number of pixels
    print(paste("In total there are ",length(sites)," land pixels to run",sep=""))

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
