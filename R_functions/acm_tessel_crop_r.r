
###
## Function containing the DALEC_CDEA model in R form
###

# Code is based on the DALEC models
# Refs:
# Williams et al., (1997)
# Sus et al., (2010)
# Bloom & Williams (2015)
# R library(randomForets)
# All R coding by T. L. Smallman (t.l.smallman@ed.ac.uk, UoE)

# DRIVERS
# maxt= max daily temperature (oC)
# mint= min daily temperature (oC)
# radiation = sum shortwave radiation (MJ.day-1)
# lat = latidude in degrees
# Ceff = canopy efficiency parameter from DALEC
# yearday = day of year
# co2 = co2 (ppm) often held constant

#acm_tessel <- function (LAI,maxt,mint,co2,yearday,lat,radiation,constants) {
#  # This is an R version of the aggregated canopy model (ACMv1; Williams et al., 1997).
#  # ACMv1 is a fast GPP response function sensitive to day length, radiation, atmospheric CO2,
#  # temperature and a location specific parameter of canopy efficiency
#  # This version of the code was written by T. L. Smallman (t.l.smallman@ed.ac.uk, UoE)
#  gc=0;pp=0;qq=0;ci=0;e0=0;mult=0;dayl=0;cps=0;dec=0;nit=1
#
#  # default constants
#  #constants=c(0,0.0156935,4.22273,208.868,0.0453194,0.37836,7.19298,0.011136,2.1001,0.789798,-2,1)
#
#  # determine temperature range
#  trange=0.5*(maxt-mint)
#  # daily canopy conductance
#  gc=abs(constants[11])**(constants[10])/((constants[6]*constants[12]+trange))
#  # maximum rate of temperature and nitrogen (canopy efficiency) limited photosynthesis (gC.m-2.day-1)
#  pn=LAI*nit*constants[1]*exp(constants[8]*maxt)
#  # pp and qq represent limitation by diffusion and metabolites respecitively
#  pp=pn/gc ; qq=constants[3]-constants[4]
#  # calculate internal CO2 concentration (ppm)
#  ci=0.5*(co2+qq-pp+((co2+qq-pp)**2-4*(co2*qq-pp*constants[3]))**0.5)
#  # limit maximum quantium efficiency by leaf area
#  e0=constants[7]*LAI^2/(LAI^2+constants[9])
#  # calculate day length (hours)
#  dec = - asin( sin( 23.45 * pi / 180. ) * cos( 2. * pi * ( yearday + 10. ) / 365. ) )
#  sinld = sin( lat*(pi/180.) ) * sin( dec )
#  cosld = cos( lat*(pi/180.) ) * cos( dec )
#  aob = sinld / cosld
#  dayl = 12.0 * ( 1. + 2. * asin( aob ) / pi )
#  # calculate CO2 limited rate of photosynthesis
#  pd=gc*(co2-ci)
#  # calculate combined light and CO2 limited photosynthesis
#  cps=e0*radiation*pd/(e0*radiation+pd)
#  # correct for day length variation
#  model=cps*(constants[2]*dayl+constants[5])
#  return(model)
#}
#
#predictRegTree<-function(x,nsample,mdim,lDaughter,rDaughter,nodestatus,split,nodepred,splitVar,idx1) {
#
#    # R code version of functions from the R library(randomForest)
#    # R translation by T. L. Smallman (t.l.smallman@ed.ac.uk, UoE)
#
#    ypred=c(1,dim=c(nsample))
#    for (i in seq(1, nsample)) {
#         k = 1
#         #/* go down the tree */
#         while (nodestatus[k,idx1] != -1) {
#            m = splitVar[k,idx1] #- 1
#            if (x[m,i] <= split[k,idx1]) {
#                k = lDaughter[k,idx1] #- 1
#            } else {
#                k = rDaughter[k,idx1] #- 1
#            }
#         }
#         #/* terminal node: assign prediction and move on to next */
#         ypred[i] = nodepred[k,idx1]
#    }
#
#    # return output
#    return(ypred)
#
#} # end function predictRegTree
#
#regForest<-function(x,mdim,n,ntree,lDaughter,rDaughter,nodestatus,xsplit,avnodes,mbest) {
#
#    # R code version of functions from the R library(randomForest)
#    # R translation by T. L. Smallman (t.l.smallman@ed.ac.uk, UoE)

#    # define output variable
#    ypred = array(0, dim=c(n))
#
#    # initial conditions
#    idx1 = 1
#
#    # looks like we run each tree for each location first then move onto the next tree and keep adding things up
#    for (i in seq(1, ntree)) {
#         # zeros of output variable for this tree
#         # n = number of predictions required
#         ytree = array(0, dim=c(n))
#         # key output from predictRegTree appears to be ytree
#         ytree = predictRegTree(x,n,mdim,lDaughter,rDaughter,nodestatus,xsplit,avnodes,mbest,idx1)
#
#         # add this trees solution to the given requested values
#         ypred = ypred + ytree
#
#         #/* increment the offset */
#         idx1 = idx1 + 1
#    }
#
#    # return variable of interest
#    return(ypred/ntree)
#
#} # end function regForest
#
#lookup_gpp<-function (object,new_data,lat,doy) {
#
#    # This function will calculate the GPP based on randomForests generated in the randomForest library
#    # The function has been simplied from the original, this however means that we are dependent on the user
#    # correctly presenting the new data in rows (each point) and columns (meteorological drivers)
#    # in order of "lai","maxt","mint","RAD","dayl","CO2","lagged_precip"
#    # This function was created by T. L Smallman (t.l.smallman@ed.ac.uk, UoE)
#
#    # calculate the hours in day
#    dec = - asin( sin( 23.45 * pi / 180. ) * cos( 2. * pi * ( doy + 10. ) / 365. ) )
#    sinld = sin( lat*(pi/180.) ) * sin( dec )
#    cosld = cos( lat*(pi/180.) ) * cos( dec )
#    aob = sinld / cosld
#    dayl = 12.0 * ( 1. + 2. * asin( aob ) / pi )
#
#    # create newdata array now with dimensions of (iter,var)
#    x=array(NA,dim=c(length(new_data$lai),(length(new_data))))
#    # merge lai and met drivers into array
#    x[,1]=new_data$lai
#    x[,2]=rep(new_data$maxt, times=length(new_data$lai))
#    x[,3]=rep(new_data$mint, times=length(new_data$lai))
#    x[,4]=rep(new_data$radiation, times=length(new_data$lai))
#    x[,5]=rep(dayl, times=length(new_data$lai))
#    x[,6]=rep(new_data$co2, times=length(new_data$lai))
#    x[,7]=rep(new_data$lagged_precip, times=length(new_data$lai))
#
#    # extract dimensional information about the model
#    mdim <- ncol(x) # number of drivers
#    ntest <- nrow(x) # number of predictions requested
#    ntree <- object$ntree # how many trees
#    # transpose the input variables
#    x <- t(data.matrix(x))
#    ans=regForest(x,mdim,ntest,ntree,object$leftDaughter,
#		  object$rightDaughter,object$nodestatus,
#		  object$xbestsplit,object$nodepred,
#		  object$bestvar)
#
#    # return out solutions
#    return(ans)
#
#} # end function lookup_gpp
#
#read_binary_response_surface<- function(infile) {
#
#      # This function was created to read CARDAMOM-DALEC binary driver files
#      # and extract information needed to run a GPP emulator based on a random forest model
#      # This code was written by T. L. Smallman (t.l.smallman@ed.ac.uk)
#
#      print("Beginning read of GPP response surface binary input files...")
#      # Open and read the DALEC binary driver file
#      # open this chains binary file into R, instructing 'r' to read and 'b' for binary
#      bob=file(infile,'rb') ; nos_var=1e6
#      bd=readBin(bob, double(),nos_var)
#      # keep reading until we have read all that can be read
#      set1=NA
#      while (length(set1) > 0) {
#	      set1=readBin(bob, double(),nos_var)
#	      bd=append(bd,set1)
#      }
#      # now close this chain
#      close(bob)
#
#      # begin preparing to disagregate the different sections of the file
#      k=0
#      # extract static data (100 places)
#      static=bd[(k+1):(k+100)]
#      k=k+100
#
#      # PFT
#      md=list(pft=static[1])
#      # number of random trees
#      md$ntree=static[2]
#      # dimension 1 of random tree
#      md$dim_1=static[3]
#      # dimension 2 of random tree
#      md$dim_2=static[4]
#      # number of model inputs needed
#      md$nos_inputs=static[5]
#
#      # extract various components of the random tree
#      md$leftDaughter=bd[(k+1):(k+(md$dim_1*md$dim_2))] ; k=k+(md$dim_1*md$dim_2)
#      # restructure correctly
#      md$leftDaughter=array(md$leftDaughter,dim=c(md$dim_1,md$dim_2))
#
#      # extract various components of the random tree
#      md$rightDaughter=bd[(k+1):(k+(md$dim_1*md$dim_2))] ; k=k+(md$dim_1*md$dim_2)
#      # restructure correctly
#      md$rightDaughter=array(md$rightDaughter,dim=c(md$dim_1,md$dim_2))
#
#      # extract various components of the random tree
#      md$nodestatus=bd[(k+1):(k+(md$dim_1*md$dim_2))] ; k=k+(md$dim_1*md$dim_2)
#      # restructure correctly
#      md$nodestatus=array(md$nodestatus,dim=c(md$dim_1,md$dim_2))
#
#      # extract various components of the random tree
#      md$xbestsplit=bd[(k+1):(k+(md$dim_1*md$dim_2))] ; k=k+(md$dim_1*md$dim_2)
#      # restructure correctly
#      md$xbestsplit=array(md$xbestsplit,dim=c(md$dim_1,md$dim_2))
#
#      # extract various components of the random tree
#      md$nodepred=bd[(k+1):(k+(md$dim_1*md$dim_2))] ; k=k+(md$dim_1*md$dim_2)
#      # restructure correctly
#      md$nodepred=array(md$nodepred,dim=c(md$dim_1,md$dim_2))
#
#      # extract various components of the random tree
#      md$bestvar=bd[(k+1):(k+(md$dim_1*md$dim_2))] ; k=k+(md$dim_1*md$dim_2)
#      # restructure correctly
#      md$bestvar=array(md$bestvar,dim=c(md$dim_1,md$dim_2))
#
#      # pass back out information
#      return(md)
#
#} # end of function

crop_development_parameters<-function(infile) {

    #! subroutine reads in the fixed crop development files which are linked the
    # the development state of the crops. The development model varies between
    # which species. e.g. winter wheat and barley, spring wheat and barley
    # The original Fortran code was written by O. Sus and translated to R by
    # T. L. Smallman (t.l.smallman@ed.ac.uk, UoE)

    # open the crop development file and pass link
    crop_development <- file(infile, "r", blocking = FALSE)

    # read in the amount of carbon available (as labile) in each seed..
    stock_seed_labile = as.numeric(unlist(strsplit(readLines(crop_development, n=1, ok=FALSE),","))[2]) # empty

    # read in C partitioning/fraction data and corresponding developmental
    # stages (DS)
    # shoot
    dims=array(0,dim=9)
    temp=readLines(crop_development, n=1, ok=FALSE)
    temp=as.numeric(unlist(strsplit(readLines(crop_development, n=1, ok=FALSE),",")))
    rows=temp[1] ; columns=temp[2] ; dims[1:3]=rows
    DS_shoot=array(NA,dim=c(rows))
    fol_frac=array(NA,dim=c(rows))
    stem_frac=array(NA,dim=c(rows))
    for (i in seq(1, rows)) {
         temp=as.double(unlist(strsplit(readLines(crop_development, n=1, ok=FALSE),",")))
         DS_shoot[i] = temp[1] ; fol_frac[i] = temp[2] ; stem_frac[i] = temp[3]
    }

    # root
    temp=readLines(crop_development, n=1, ok=FALSE)
    temp=as.numeric(unlist(strsplit(readLines(crop_development, n=1, ok=FALSE),",")))
    rows=temp[1] ; columns=temp[2] ; dims[4:5]=rows
    DS_root=array(NA,dim=c(rows))
    root_frac=array(NA,dim=c(rows))
    for (i in seq(1, rows)) {
         temp=as.double(unlist(strsplit(readLines(crop_development, n=1, ok=FALSE),",")))
         DS_root[i] = temp[1] ; root_frac[i] = temp[2]
    }

    # loss rates of leaves and roots
    # leaves
    temp=readLines(crop_development, n=1, ok=FALSE)
    temp=as.numeric(unlist(strsplit(readLines(crop_development, n=1, ok=FALSE),",")))
    rows=temp[1] ; columns=temp[2] ; dims[6:7]=rows
    DS_LRLV=array(NA,dim=c(rows))
    LRLV=array(NA,dim=c(rows))
    for (i in seq(1, rows)) {
         temp=as.double(unlist(strsplit(readLines(crop_development, n=1, ok=FALSE),",")))
         DS_LRLV[i] = temp[1] ; LRLV[i] = temp[2]
    }

    # roots
    temp=readLines(crop_development, n=1, ok=FALSE)
    temp=as.numeric(unlist(strsplit(readLines(crop_development, n=1, ok=FALSE),",")))
    rows=temp[1] ; columns=temp[2] ; dims[8:9]=rows
    DS_LRRT=array(NA,dim=c(rows))
    LRRT=array(NA,dim=c(rows))
    for (i in seq(1, rows)) {
         temp=as.double(unlist(strsplit(readLines(crop_development, n=1, ok=FALSE),",")))
         DS_LRRT[i] = temp[1] ; LRRT[i] = temp[2]
    }

    # close connection for next use
    close(crop_development, type="r")

    # return output
    return(list(stock_seed_labile=stock_seed_labile,
                DS_shoot = DS_shoot,fol_frac = fol_frac,stem_frac = stem_frac,
                DS_root = DS_root, root_frac = root_frac,
                DS_LRLV = DS_LRLV, LRLV = LRLV,
                DS_LRRT = DS_LRRT, LRRT = LRRT,crop_dims=dims))

} # crop_development_parameters
