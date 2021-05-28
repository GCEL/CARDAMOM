
###
## Function to load a gridded dataset of prior values for crop planting and harvest dates
###

# This function is by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_sacks_calendar_fields_for_extraction<-function(latlon_in,crop_management_source) {

  if (crop_management_source == "sacks_crop_calendar") {

    # let the user know this might take some time
    print("Loading processed 5 min Sacks crop calendar fields for subsequent sub-setting ...")

    # open processed modis files
    input_file_1=paste(path_to_crop_management,"/Wheat.Winter.crop.calendar.fill.nc",sep="")
    data1=nc_open(input_file_1)

    # extract location variables
    lat=ncvar_get(data1, "latitude") ; long=ncvar_get(data1, "longitude")
    # read the modis lai drivers
    plant=ncvar_get(data1, "plant")
    #	plant_range=ncvar_get(data1, "plant.range")
    harvest=ncvar_get(data1, "harvest")
    #	harvest_range=ncvar_get(data1, "harvest.range")

    if (length(dim(latlon_in)) > 1) {
      max_lat=max(latlon_in[,1])+0.5 ; max_long=max(latlon_in[,2])+0.5
      min_lat=min(latlon_in[,1])-0.5 ; min_long=min(latlon_in[,2])-0.5
    } else {
      max_lat=max(latlon_in[1])+0.5 ; max_long=max(latlon_in[2])+0.5
      min_lat=min(latlon_in[1])-0.5 ; min_long=min(latlon_in[2])-0.5
    }
    # determine which locations to keep form inputs
    keep_lat=which(lat > min_lat & lat < max_lat)
    keep_long=which(long > min_long & long < max_long)
    # filter the inputs
    plant=plant[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
    plant_range=-9999 #plant_range[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
    harvest=harvest[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
    harvest_range=-9999 #harvest_range[min(keep_long):max(keep_long),min(keep_lat):max(keep_lat)]
    # now filter the inputs spatial data
    lat=lat[keep_lat] ; long=long[keep_long]
    # close files after use
    nc_close(data1)

    # clean
    rm(keep_lat,keep_long,max_lat,max_long,min_lat,min_long) ; gc(reset=TRUE,verbose=FALSE)

    # output variables
    return(list(plant=plant,plant_range=plant_range,harvest=harvest,harvest_range=harvest_range,lat=lat,long=long))

  } else {
    # output variables
    return(list(plant=-9999,plant_range=-9999,harvest=-9999,harvest_range=-9999,lat=-9999,long=-9999))
  }

} # function end
