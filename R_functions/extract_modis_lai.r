
###
## Function extracts location specific information on leaf area index from a
## globally gridded MODIS dataset that has already been loaded.
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

# declare function
bintodec <- function(x) {
  sum(2^(which(rev(unlist(strsplit(as.character(x), "")) == 1))-1))
}

extract_modis_lai<- function(timestep_days,spatial_type,resolution,grid_type,latlon_in,lai_all,years_to_load) {

  # convert input data long to conform to what we need
  check1=which(lai_all$long > 180) ; if (length(check1) > 0) { lai_all$long[check1]=lai_all$long[check1]-360 }
  # find the nearest location
  output=closest2d(1,lai_all$lat,lai_all$long,latlon_in[1],latlon_in[2],2)
  i1=unlist(output)[1] ; j1=unlist(output)[2]
  print(paste("LAI data extracted for current location ",Sys.time(),sep=""))

  # return long to 0-360
  if (length(check1) > 0) { lai_all$long[check1]=lai_all$long[check1]+360 }


  # work out number of pixels to average over
  if (spatial_type == "grid") {
      # resolution of the product
      product_res = abs(lai_all$lat[1,2]-lai_all$lat[1,1])+abs(lai_all$long[2,1]-lai_all$long[1,1])
      product_res = product_res * 0.5 # NOTE: averaging needed for line above
      if (grid_type == "wgs84") {
          # radius is ceiling of the ratio of the product vs analysis ratio
          radius = floor(0.5*(resolution / product_res))
          max_radius = radius+4
      } else if (grid_type == "UK") {
          # Estimate radius for UK grid assuming radius is determine by the longitude size
          # 6371e3 = mean earth radius (m)
          radius = round(rad2deg(sqrt((resolution / 6371e3**2))) / product_res, digits=0)
          #radius = max(0,floor(1*resolution*1e-3*0.5))
          max_radius = radius+4
      } else {
         stop("have not specified the grid used in this analysis")
      }
  } else {
      radius = 0
      max_radius = 4
  }

  answer=NA
  while (is.na(answer) == TRUE) {
    # work out average areas
    average_i=max(1,(i1-radius)):min(dim(lai_all$lai_all)[1],(i1+radius)) ; average_j=max(1,(j1-radius)):min(dim(lai_all$lai_all)[2],(j1+radius))
    # carry out averaging
    lai=array(NA, dim=c(dim(lai_all$lai_all)[3]))
    for (n in seq(1, dim(lai_all$lai_all)[3])) {
      lai[n]=mean(lai_all$lai_all[average_i,average_j,n], na.rm=TRUE)
    }
    # are any of the my data points now filled
    answer=max(as.vector(lai),na.rm=TRUE)
    # what proportion of my data points are within a realistic range
    npoints=length(which(lai > 0 & lai < 12))/length(lai)
    # do I have at least 20 % of data points filled
    if (is.na(answer) | answer == -Inf | npoints < 0.2) {radius=radius+1 ; answer=NA}
    # restrict how far to look before giving up
    if (radius >= max_radius) {answer=1}
  }
  # warning to the used
  print(paste("NOTE LAI averaged over a pixel radius (i.e. centre + radius) of ",radius," points",sep=""))
  # convert missing data back to -9999
  lai[which(is.na(lai))]=-9999.0
  # next work out how many days we should have in the year
  doy_out=0
  for (i in seq(1, length(years_to_load))) {
       nos_days = nos_days_in_year(years_to_load[i])
#       # is current year a leap or not
#       nos_days = 365
#       mod=as.numeric(years_to_load[i])-round((as.numeric(years_to_load[i])/4))*4
#      if (mod == 0) {
#          nos_days = 366
#          mod=as.numeric(years_to_load[i])-round((as.numeric(years_to_load[i])/100))*100
#          if (mod == 0) {
#              nos_days  = 365
#              mod=as.numeric(years_to_load[i])-round((as.numeric(years_to_load[i])/400))*400
#              if (mod == 0) {
#                  nos_days  = 366
#              }
#          }
#      }
      # count up days needed
      doy_out=append(doy_out,1:nos_days)
  }
  doy_out=doy_out[-1]

  # just incase there is no missing data we best make sure there is a value which can be assessed
  if (length(lai_all$missing_years) == 0) { lai_all$missing_years=1066 }

  # declare output variable
  lai_out=array(-9999, dim=length(doy_out))
  # now line up the obs days with all days
  b=1 ; i=1 ; a=1 ; start_year=as.numeric(years_to_load[1])
  while (b <= length(lai_all$doy_obs)) {
    # if we are in a year which is missing then we do not allow consideration of DOY
    if (start_year != lai_all$missing_years[a]) {
      if (doy_out[i] == lai_all$doy_obs[b]) {
        lai_out[i]=lai[b] ; b=b+1
      } # end if doy matches
    } # end if missing year
    # but we do keep counting through the total vector length which we expect
    i=i+1
    # each time we come back to doy_out[i]==1 we need to count on the year
    if (doy_out[i] == 1) {
      # and if we have just been in a missing year we need to count on the missing years vector to
      if (start_year == lai_all$missing_years[a]) {a=min(length(lai_all$missing_years),a+1)}
      start_year=start_year+1
    } # end if doy_out[i] == 1
  } # end while condition

  if (length(timestep_days) == 1 & timestep_days[1] == 1) {
    # well actually we do nothing
  } else {
    # generally this now deals with time steps which are not daily.
    # However if not monthly special case
    if (length(timestep_days) == 1) {
      run_day_selector=seq(1,length(lai_out),timestep_days)
      timestep_days=rep(timestep_days, length.out=length(lai_out))
    }
    print("...calculating monthly averages for lai")
    # determine the actual daily positions
    run_day_selector=cumsum(timestep_days)
    # create needed variables
    lai_agg=array(NA,dim=length(run_day_selector))
    for (y in seq(1,length(run_day_selector))) {
      pick=lai_out[(run_day_selector[y]-timestep_days[y]):run_day_selector[y]]
      lai_agg[y]=mean(pick[which(pick != -9999)],na.rm=TRUE)
    }
    # convert missing values back to -9999
    lai_agg[which(is.na(lai_agg))]=-9999
    # update with new output information
    lai_out=lai_agg
    # clean up
    rm(lai_agg,y) ; gc()
  } # monthly aggregation etc

  # clean up
  rm(i1,j1,check1,lai,i,nos_days,doy_out,radius,n,a) ; gc(reset=TRUE,verbose=FALSE)

  # assume LAI uncertainty for MODIS is 0.25 m2/m2
  lai_unc_out = array(-9999, dim=length(lai_out))
  lai_unc_out[which(lai_out > 0)] = 0.25

  # pass the information back
  output = list(lai = lai_out, lai_unc = lai_unc_out)
  return(output)

} # end of function

# old function source code for extracting direct from the tiles as downloaded from NASA

## MODIS defaults
# Modis tile location, for lat / long conversion
#      h=17 ; v=03
# for correct file selection ".h**v**."
#      modis_tile=paste(".h17v03.",sep="")
# Radius around location to average over
#     radius=max(1,floor(1*cardamom_resolution*1e-3*0.5))
# To carry out linear interpolation TRUE or FALSE
#    interp_vals=FALSE
##

# This section determnines cartesian coordinates for modis grid
#      vlevels=seq(85,-85,by=-10)
#      hlevels=seq(-175,175,by=10)

#      vlevel=vlevels[v+1]
#      hlevel=hlevels[h+1]

#     hrvlevels=seq(vlevel+5-5/1200,vlevel-5+5/1200,by = -10/1200)
#     hrhlevels=seq(hlevel-5+5/1200,hlevel+5-5/1200, by = 10/1200)

# contruct x and y coordinates such that they are in common array
#    long=array(hrhlevels,dim=c(length(hrhlevels),length(hrvlevels)))
#    lat=array(hrvlevels[length(hrvlevels):1],dim=c(length(hrhlevels),length(hrvlevels)))
#    lat=t(lat)

# correct for sinusoidal
#      long=long/cos(lat*pi/180)
#      long[long>180]=long[long>180]-360
#      long[long< -180]=long[long< -180]+360

# find locations desired
#      output=closest2d(1,lat,long,latlon_in[1],latlon_in[2],2)
#      output_i=unlist(output)[1];output_j=unlist(output)[2]

# list all files present first
#      list=list.files(path=path_to_lai,full.names=T)

# find out which files are the type we want
#      is_it<-grepl(c(".h5"),list)
#      # now re-list based on this
#     list=list[is_it]
#      # check for tile we want only
#      is_it<-grepl(modis_tile,list)
#      # now re-list based on this
#      list=list[is_it]

#    doy_out=0
#     for (i in seq(1, length(years_to_do))) {
#	  # is current year a leap or not
#	  nos_days = 365
#	  mod=as.numeric(years_to_do[i])-round((as.numeric(years_to_do[i])/4))*4
#	  if (mod == 0) {
#	      nos_days = 366
#	      mod=as.numeric(years_to_do[i])-round((as.numeric(years_to_do[i])/100))*100
#	      if (mod == 0) {
#		  nos_days  = 365
#		  mod=as.numeric(years_to_do[i])-round((as.numeric(years_to_do[i])/400))*400
#		  if (mod == 0) {
#		      nos_days  = 366
#		  }
#	      }
#	  }
# count up days needed
#	  doy_out=append(doy_out,1:nos_days)
#     }
#      doy_out=doy_out[-1]
# fill in blanks
#      lai_out=array(-9999,dim=c(length(doy_out))) ; rm(nos_days)

#     for (yr in seq(1, length(years_to_do))) {
#	  # do only the year we want
#	  is_it<-grepl(paste("A",years_to_do[yr],sep=""),list)
#	  # now re-list based on this
#	  sublist=list[is_it]

#	  bits=unlist(strsplit(sublist, ".h"))
#	  bits=unlist(strsplit(bits, paste("A",years_to_do[yr],sep="")))
#	  if (yr == 1) {doy=bits[2]} else {doy=append(doy,bits[2])}
#	  for (i in seq(6,length(bits),by=4)) {
#	      doy=append(doy,bits[i])
#	  }
#	  doy=as.numeric(doy)

#tile grid
#0-35
#0-17
#more info on https://lpdaac.usgs.gov/products/modis_overview
#0-0 NW
#35-17 SE

# read in the first file now
#	  lai_1km=h5read(list[1], "MOD_Grid_MOD15A2/Data Fields/Lai_1km")
#	  lai_std_1km=h5read(list[1], "MOD_Grid_MOD15A2/Data Fields/LaiStdDev_1km")
#	  Fpar_1km=h5read(list[1], "MOD_Grid_MOD15A2/Data Fields/Fpar_1km")
#	  FparLai_QC=h5read(list[1], "MOD_Grid_MOD15A2/Data Fields/FparLai_QC")
#	  FparStdDev_1km=h5read(list[1], "MOD_Grid_MOD15A2/Data Fields/FparStdDev_1km")

# define filters
#	  filter=array(1,dim=dim(lai_1km))
#	  fparfilter=array(1,dim=dim(lai_1km))
#	  laifilter=array(1,dim=dim(lai_1km))

# Additional filter conditions
# FparExtra_QC 6 BITFIELDS IN 8 BITWORD
# LANDSEA PASS-THROUGH START 0 END 1 VALIDS 4
# LANDSEA   00 = 0 LAND       AggrQC(3,5)values{001}
# LANDSEA   01 = 1 SHORE      AggrQC(3,5)values{000,010,100}
# LANDSEA   10 = 2 FRESHWATER AggrQC(3,5)values{011,101}
# LANDSEA   11 = 3 OCEAN      AggrQC(3,5)values{110,111}
# SNOW_ICE (from Aggregate_QC bits) START 2 END 2 VALIDS 2
# SNOW_ICE  0 = No snow/ice detected
# SNOW_ICE  1 = Snow/ice were detected
# AEROSOL START 3 END 3 VALIDS 2
# AEROSOL   0 = No or low atmospheric aerosol levels detected
# AEROSOL   1 = Average or high aerosol levels detected
# CIRRUS (from Aggregate_QC bits {8,9} ) START 4 END 4 VALIDS 2
# CIRRUS    0 = No cirrus detected
# CIRRUS    1 = Cirrus was detected
# INTERNaNL_CLOUD_MASK START 5 END 5 VALIDS 2
# INTERNaNL_CLOUD_MASK 0 = No clouds
# INTERNaNL_CLOUD_MASK 1 = Clouds were detected
# CLOUD_SHADOW START 6 END 6 VALIDS 2
# CLOUD_SHADOW        0 = No cloud shadow detected
# CLOUD_SHADOW        1 = Cloud shadow detected
# SCF_BIOME_MASK START 7 END 7 VALIDS 2
# SCF_BIOME_MASK  0 = Biome outside interval <1,4>
# SCF_BIOME_MASK  1 = Biome in interval <1,4>

# FparLai_QC 5 BITFIELDS IN 8 BITWORD
# MODLAND_QC START 0 END 0 VALIDS 2
# MODLAND_QC   0 = Good Quality (main algorithm with or without saturation)
# MODLAND_QC   1 = Other Quality (back-up algorithm or fill value)
# SENSOR START 1 END 1 VALIDS 2
# SENSOR       0  = Terra
# SENSOR       1  = Aqua
# DEADDETECTOR START 2 END 2 VALIDS 2
# DEADDETECTOR 0 = Detectors apparently fine for up to 50% of channels 1,2
# DEADDETECTOR 1 = Dead detectors caused >50% adjacent detector retrieval
# CLOUDSTATE START 3 END 4 VALIDS 4 (this inherited from Aggregate_QC bits {0,1} cloud state)
# CLOUDSTATE   00 = 0 Significant clouds NOT present (clear)
# CLOUDSTATE   01 = 1 Significant clouds WERE present
# CLOUDSTATE   10 = 2 Mixed cloud present on pixel
# CLOUDSTATE   11 = 3 Cloud state not defined,assumed clear
# SCF_QC START 5 END 7 VALIDS 5
# SCF_QC       000=0 Main (RT) algorithm used, best result possible (no saturation)
# SCF_QC       001=1 Main (RT) algorithm used, saturation occured. Good, very usable.
# SCF_QC       010=2 Main algorithm failed due to bad geometry, empirical algorithm used
# SCF_QC       011=3 Main algorithm faild due to problems other than geometry, empirical algorithm used
# SCF_QC       100=4 Pixel not produced at all, value coudn't be retrieved (possible reasons: bad L1B data, unusable MODAGAGG data)

# convert from binary to real numbers
#	  opt1=bintodec('00000000');
#	  opt2=bintodec('00100000');
#	  opt3=bintodec('00000010');
#	  opt4=bintodec('00100010');

# filer based on options defined by Bloom,
#	  # awaiting documentation .....
#	  filter[FparLai_QC != opt1 & FparLai_QC != opt2 & FparLai_QC != opt3 & FparLai_QC != opt4]="NaN"
#	  filter=array(as.numeric(filter), dim=dim(filter))
#

# 255 = _Fillvalue, assigned when:
#     * the MODAGAGG suf. reflectance for channel VIS, NIR was assigned its _Fillvalue, or
#     * land cover pixel itself was assigned _Fillvalus 255 or 254.
# 254 = land cover assigned as perennial salt or inland fresh water.
# 253 = land cover assigned as barren, sparse vegetation (rock, tundra, desert.)
# 252 = land cover assigned as perennial snow, ice.
# 251 = land cover assigned as "permanent" wetlands/inundated marshlands.
# 250 = land cover assigned as urban/built-up.
# 249 = land cover assigned as "unclassified" or not able to determine.
# 248 = no standard deviation available, pixel produced using backup method.

# declare flags to miss out
#	  nonos=c(255,254,253,252,251,250,249,248)

#	  laifilter[lai_std_1km==nonos[1] | lai_std_1km==nonos[2] | lai_std_1km==nonos[3] | lai_std_1km==nonos[4] | lai_std_1km==nonos[5] | lai_std_1km==nonos[6] | lai_std_1km==nonos[7]]="NaN"
#	  laifilter=array(as.numeric(laifilter), dim=dim(laifilter))

#	  # ready for output with time series
#	  lai_1km=lai_1km[,dim(lai_1km)[2]:1]*0.1
#	  lai_std_1km=lai_std_1km[,dim(lai_1km)[2]:1]
#	  filter=filter[,dim(lai_1km)[2]:1]
#	  laifilter=laifilter[,dim(lai_1km)[2]:1]
#
#	  # filter away the 'bad data'
#	  lai_1km=lai_1km*laifilter*filter
#	  lai_std_1km=lai_std_1km*laifilter*filter

#	  # work out average area first
#	  average_i=(output_i-radius):(output_i+radius) ; average_j=(output_j-radius):(output_j+radius)
#	  # then crude temp before downloading all needed lai data
#	  average_i[which(average_i < 1)]=1 ; average_i[which(average_i > dim(lai_1km)[1])]=dim(lai_1km)[1]
#	  average_j[which(average_j < 1)]=1 ; average_j[which(average_j > dim(lai_1km)[2])]=dim(lai_1km)[2]
#	  print("WARNING LAI currently fudged to single tile - some values may be bollocks")

#	  if (yr == 1) {
#	      lai=mean(lai_1km[average_i,average_j], na.rm=TRUE)
#	      lai_std=mean(lai_std_1km[average_i,average_j], na.rm=TRUE)
#	  }else {# arrange the dataset with initial values for each site
#	      lai=append(lai,mean(lai_1km[average_i,average_j], na.rm=TRUE))
#	      lai_std=append(lai_std,mean(lai_std_1km[average_i,average_j], na.rm=TRUE))
#	  }

#	  # now loop through the remaining dataset
#	  for (i in seq(2, length(sublist),by=1)) {

#print(c("File ",i,"of ", length(sublist)))
#	      lai_1km=h5read(sublist[i], "MOD_Grid_MOD15A2/Data Fields/Lai_1km")
#	      lai_std_1km=h5read(sublist[i], "MOD_Grid_MOD15A2/Data Fields/LaiStdDev_1km")
#	      Fpar_1km=h5read(sublist[i], "MOD_Grid_MOD15A2/Data Fields/Fpar_1km")
#	      FparLai_QC=h5read(sublist[i], "MOD_Grid_MOD15A2/Data Fields/FparLai_QC")
#	      FparStdDev_1km=h5read(sublist[i], "MOD_Grid_MOD15A2/Data Fields/FparStdDev_1km")

# define filters
#	      filter=array(1,dim=dim(lai_1km))
#	      fparfilter=array(1,dim=dim(lai_1km))
#	      laifilter=array(1,dim=dim(lai_1km))

#	      # filer based on options defined by Bloom,
#	      # awaiting documentation .....
#	      filter[FparLai_QC != opt1 & FparLai_QC != opt2 & FparLai_QC != opt3 & FparLai_QC != opt4]="NaN"
#	      filter=array(as.numeric(filter), dim=dim(filter))

# declare flags to miss out
#	      nonos=c(255,254,253,252,251,250,249,248)

#	      laifilter[lai_std_1km==nonos[1] | lai_std_1km==nonos[2] | lai_std_1km==nonos[3] | lai_std_1km==nonos[4] | lai_std_1km==nonos[5] | lai_std_1km==nonos[6] | lai_std_1km==nonos[7]]="NaN"
#	      laifilter=array(as.numeric(laifilter), dim=dim(laifilter))

# ready for output with time series
#	      lai_1km=lai_1km[,dim(lai_1km)[2]:1]*0.1
#	      lai_std_1km=lai_std_1km[,dim(lai_1km)[2]:1]
#	      filter=filter[,dim(lai_1km)[2]:1]
#	      laifilter=laifilter[,dim(lai_1km)[2]:1]

#	      # filter away the 'bad data'
#	      lai_1km=lai_1km*laifilter*filter
#	      lai_std_1km=lai_std_1km*laifilter*filter

#	      # arrange the dataset with initial values for each site
#	      lai=append(lai,mean(lai_1km[average_i,average_j], na.rm=TRUE))
#	      lai_std=append(lai_std,mean(lai_std_1km[average_i,average_j], na.rm=TRUE))
#	  }

#	  # convert back to numeric the altered variables
#	  lai=as.numeric(lai)
#	  # make missing values -9999
#	  lai[which(is.na(lai))]=-9999

#    } # end of year loop
#    b=1 ; i=1
#   while (b <= length(doy)) {
#	if (doy_out[i] == doy[b]) {
#	    lai_out[i]=lai[b] ; b=b+1
#	}
#	i=i+1
#    }
#    # send back to answer...
#    return(lai_out)
