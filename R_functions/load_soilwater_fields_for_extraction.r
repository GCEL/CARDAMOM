
###
## Function to load gridded dataset of prior estimates for soil moisture content
###

# This function is T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_soilwater_fields_for_extraction<-function(latlon_in,soilwater_source) {

  if (soilwater_source == "GLEAM") {

    # let the user know this might take some time
    print("Loading processed GLEAM soil water fraction for subsequent sub-setting ...")

    # open processed modis files
    input_file_1=paste(path_to_gleam,"/GLEAM_soil_moisture_prior_0.25.nc",sep="")
    data1=nc_open(input_file_1)

    # extract location variables
    lat=ncvar_get(data1, "lat") ; long=ncvar_get(data1, "long")
    # read the file
    soil_water=ncvar_get(data1, "GLEAM_soil_moisture_prior")
    soil_water_unc=ncvar_get(data1, "GLEAM_soil_moisture_prior_SD")

    # close files after use
    nc_close(data1)

    # output variables
    return(list(soil_water=soil_water,soil_water_unc = soil_water_unc,lat=lat,long=long))

  } else {
    # output variables
    return(list(soil_water=-9999,soil_water_unc = -9999,lat=-9999,long=-9999))
  }

} # function end
