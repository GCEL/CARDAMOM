
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
## yearday = day of year
## co2 = co2 (ppm) often held constant
#
#acm_tessel <- function (LAI,maxt,mint,co2,yearday,lat,radiation,constants) {
#  # The aggregated canopy model, a GPP response function
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

#lookup_gpp<-function (lai,new_data,doy,lat,response_matrix,coefficients,nos_fixed,maximum,minimum) {
#
#      # calculate the hours in day
#      dec = - asin( sin( 23.45 * pi / 180. ) * cos( 2. * pi * ( doy + 10. ) / 365. ) )
#      sinld = sin( lat*(pi/180.) ) * sin( dec )
#      cosld = cos( lat*(pi/180.) ) * cos( dec )
#      aob = sinld / cosld
#      dayl = 12.0 * ( 1. + 2. * asin( aob ) / pi )
#      # load hours in day
#      # its actually three but we have not added position one until the line below
#      new_data[2]=dayl
#
#      # create newdata array now with dimensions of (iter,var)
#      newdata=array(NA,dim=c(length(lai),(length(new_data)+1)))
#      # merge lai and met drivers into array
#      newdata[,1]=lai
#      newdata[,2]=rep((0.5*(new_data[1]+new_data[3])), times=length(lai))
#      newdata[,3]=rep(new_data[2], times=length(lai))
#      newdata[,4]=rep(new_data[3], times=length(lai))
#      newdata[,5]=rep(new_data[4], times=length(lai))
#      newdata[,6]=rep(new_data[5], times=length(lai))
#      newdata[,7]=rep(new_data[6], times=length(lai))
#
#      # new_data vector containing inputs (lai, min temp, dayl, max temp, SW rad, CO2)
#
#      # function checks look up matrix of smoothing functions and coefficience to calculate the GAM estimate
#      spacing_par=dim(response_matrix)[1]   ## the number of values generated in the response surface for each input
#      x0 <- c(1,rep(0,times=nos_fixed))         ## intercept column, ensures that matrix multiplcation of intercept and any fixed factors is 1
#      x0 = array(x0, dim=c(length(lai),(nos_fixed+1)))
#      dx <- 1/spacing_par    ## covariate spacing in `newd'
#      for (j in seq(0,(dim(newdata)[2]-(nos_fixed+1)),1)) { ## loop through smooth terms#
#	cols <-(nos_fixed+1)+j*9 +1:9      ## relevant cols of Xp
#	i <- floor((newdata[,j+1]-minimum[j+1])/(maximum[j+1]-minimum[j+1])*spacing_par)  ## find relevant rows of Xp
#	# pmin allows multiple comparison against the same number.
#	i = pmin(spacing_par-1,pmax(1,i))
#	w1 <- 1-(((newdata[,j+1]-minimum[j+1])/(maximum[j+1]-minimum[j+1])*spacing_par)-i) ## interpolation weights
#       w1 = pmin(1,pmax(0,w1))
#	## find approx. predict matrix row portion, by interpolation
#	x0 <- cbind(x0,response_matrix[i+1,cols]*w1 + response_matrix[i,cols]*(1-w1))
#      }
#      fv <- x0%*%coefficients #+ newdata[6]*coefficients[2]     ## evaluate and add offset
#      # in our case the link logit was used so must return exp(fv)
#      fv=exp(fv)
#      return(as.numeric(fv))
#
#} # end function lookup_gpp
#
#predictRegTree<-function(x,nsample,mdim,lDaughter,rDaughter,nodestatus,split,nodepred,splitVar,idx1) {
#
#    ypred=c(1,dim=c(nsample))
#    for (i in seq(1, nsample)) {
#	k = 1
#	#/* go down the tree */
#	while (nodestatus[k,idx1] != -1) {
#	    m = splitVar[k,idx1] #- 1
#	    if (x[m,i] <= split[k,idx1]) {
#		k = lDaughter[k,idx1] #- 1
#	    } else {
#		k = rDaughter[k,idx1] #- 1
#	    }
#	}
#	#/* terminal node: assign prediction and move on to next */
#	ypred[i] = nodepred[k,idx1]
#    }
#
#    # return output
#    return(ypred)
#
#} # end function predictRegTree
#
#regForest<-function(x,mdim,n,ntree,lDaughter,rDaughter,nodestatus,xsplit,avnodes,mbest) {
#
#    # define output variable
#    ypred = array(0, dim=c(n))
#
#    # initial conditions
#    idx1 = 1
#
#    # looks like we run each tree for each location first then move onto the next tree and keep adding things up
#    for (i in seq(1, ntree)) {
#	# zeros of output variable for this tree
#	# n = number of predictions required
#	ytree = array(0, dim=c(n))
#	# key output from predictRegTree appears to be ytree
#	ytree = predictRegTree(x,n,mdim,lDaughter,rDaughter,nodestatus,xsplit,avnodes,mbest,idx1)
#
#	# add this trees solution to the given requested values
#	ypred = ypred + ytree
#
#	#/* increment the offset */
#	idx1 = idx1 + 1
#
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

#    # extract dimensional information about the model
#    mdim <- ncol(x) # number of drivers
#    ntest <- nrow(x) # number of predictions requested
#    ntree <- object$ntree # how many trees
#    # transpose the input variables
#    x <- t(data.matrix(x))
#
#    ans=regForest(x,mdim,ntest,ntree,object$leftDaughter,
#		  object$rightDaughter,object$nodestatus,
#		  object$xbestsplit,object$nodepred,
#		  object$bestvar)
#
#    # return out solutions
#    return(ans)
#
#} # end function lookup_gpp

read_binary_response_surface<- function(infile) {

      print("Beginning read of GPP response surface binary input files...")
      # Open and read the DALEC binary driver file
      # open this chains binary file into R, instructing 'r' to read and 'b' for binary
      bob=file(infile,'rb') ; nos_var=1e6
      bd=readBin(bob, double(),nos_var)
      # keep reading until we have read all that can be read
      set1=NA
      while (length(set1) > 0) {
	      set1=readBin(bob, double(),nos_var)
	      bd=append(bd,set1)
      }
      # now close this chain
      close(bob)

      # begin preparing to disagregate the different sections of the file
      k=0
      # extract static data (100 places)
      static=bd[(k+1):(k+100)]
      k=k+100

      # PFT
      md=list(pft=static[1])
      # number of random trees
      md$ntree=static[2]
     # dimension 1 of random tree
      md$dim_1=static[3]
      # dimension 2 of random tree
      md$dim_2=static[4]
      # number of model inputs needed
      md$nos_inputs=static[5]

      # extract various components of the random tree
      md$leftDaughter=bd[(k+1):(k+(md$dim_1*md$dim_2))] ; k=k+(md$dim_1*md$dim_2)
      # restructure correctly
      md$leftDaughter=array(md$leftDaughter,dim=c(md$dim_1,md$dim_2))

      # extract various components of the random tree
      md$rightDaughter=bd[(k+1):(k+(md$dim_1*md$dim_2))] ; k=k+(md$dim_1*md$dim_2)
      # restructure correctly
      md$rightDaughter=array(md$rightDaughter,dim=c(md$dim_1,md$dim_2))

      # extract various components of the random tree
      md$nodestatus=bd[(k+1):(k+(md$dim_1*md$dim_2))] ; k=k+(md$dim_1*md$dim_2)
      # restructure correctly
      md$nodestatus=array(md$nodestatus,dim=c(md$dim_1,md$dim_2))

      # extract various components of the random tree
      md$xbestsplit=bd[(k+1):(k+(md$dim_1*md$dim_2))] ; k=k+(md$dim_1*md$dim_2)
      # restructure correctly
      md$xbestsplit=array(md$xbestsplit,dim=c(md$dim_1,md$dim_2))

      # extract various components of the random tree
      md$nodepred=bd[(k+1):(k+(md$dim_1*md$dim_2))] ; k=k+(md$dim_1*md$dim_2)
      # restructure correctly
      md$nodepred=array(md$nodepred,dim=c(md$dim_1,md$dim_2))

      # extract various components of the random tree
      md$bestvar=bd[(k+1):(k+(md$dim_1*md$dim_2))] ; k=k+(md$dim_1*md$dim_2)
      # restructure correctly
      md$bestvar=array(md$bestvar,dim=c(md$dim_1,md$dim_2))

      # pass back out information
      return(md)

} # end of function

#read_binary_response_surface<- function(infile) {
#
#      print("...Beginning read of gpp emulator parameters...")
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
#      # number of coefficients (fixed and smoothed)
#      md$nos_coef=static[2]
#      # number of fixed coefficients
#      md$nos_fixed=static[3]
#      # dimension 1 of response surface (same as interpolation interval)
#      md$dim_1=static[4]
#      # dimension 2 of response surface (same as interpolation interval)
#      md$dim_2=static[5]
#      # number of model inputs needed
#      md$nos_inputs=static[6]
#
#      # extract maximum driver values
#      md$maximums=bd[(k+1):(k+md$nos_inputs)] ; k=k+md$nos_inputs
#      # extract minimum driver values
#      md$minimums=bd[(k+1):(k+md$nos_inputs)] ; k=k+md$nos_inputs
#      # extract model coefficients
#      md$coefficients=bd[(k+1):(k+md$nos_coef)] ; k=k+md$nos_coef
#      # extract response surface
#      md$response_surface=bd[(k+1):(k+(md$dim_1*md$dim_2))] ; k=k+(md$dim_1*md$dim_2)
#      # restructure correctly
#      md$response_surface=array(md$response_surface,dim=c(md$dim_1,md$dim_2))
#
#      # pass back out information
#      return(md)
#
#} # end of function

# PARAMETERS
# 17 values

# p(1) Litter to SOM conversion rate  - m_r
# p(2) Fraction of GPP respired - f_a
# p(3) Fraction of NPP allocated to foliage - f_f
# p(4) Fraction of NPP allocated to roots - f_r
# p(5) Leaf lifespan - L_f
# p(6) Turnover rate of wood - t_w
# p(7) Turnover rate of roots - t_r
# p(8) Litter turnover rate - t_l
# p(9) SOM turnover rate  - t_S
# p(10) Parameter in exponential term of temperature - \theta
# p(11) = date of Clab release - B_day
# p(12) = Fraction allocated to Clab - f_l
# p(13) = lab release duration period - R_l
# p(14) = date of leaf fall - F_day
# p(15) = leaf fall duration period - R_f
# p(16) = LMA

# C pools initial conditions
# p(17) labile C (gC.m-2)
# p(18) foliar C (gC.m-2)
# p(19) wood C (gC.m-2)
# p(20) root C (gC.m-2)
# p(21) litter C (gC.m-2)
# p(22) som C (gC.m-2)

#ospolynomial <- function(L,w){
#    # Function calculates the day offset for Labile release and leaf turnover functions
#    mxc=c(0.000023599784710,0.000332730053021,0.000901865258885,-0.005437736864888, -0.020836027517787, 0.126972018064287, -0.188459767342504)
#    # load log of leaf / labile turnovers
#    LLog=log(L-1)
#    ospolynomial=(mxc[1]*LLog**6.+ mxc[2]*LLog**5.+ mxc[3]*LLog**4.+ mxc[4]*LLog**3.+ mxc[5]*LLog**2.+ mxc[6]*LLog**1.+ mxc[7])*w
#}
#
#acm_tessel_r <- function(met,p,lat,pft,exepath) {
#
#    stime=proc.time()["elapsed"]
#
#    ###
#    ## load met drivers
#
#    mint=met[,2]
#    maxt=met[,3]
#    radiation=met[,4]
#    co2=met[,5]
#    doy=met[,6]
#    lagged_precip=met[,7]
#
#    ###
#    ## load initial conditions
#
#    C_labile=array(p[17,],dim=c(dim(p)[2],length(mint)+1))
#    C_foliar=array(p[18,],dim=c(dim(p)[2],length(mint)+1))
#    C_root=array(p[19,],dim=c(dim(p)[2],length(mint)+1))
#    C_wood=array(p[20,],dim=c(dim(p)[2],length(mint)+1))
#    C_litter=array(p[21,],dim=c(dim(p)[2],length(mint)+1))
#    C_som=array(p[22,],dim=c(dim(p)[2],length(mint)+1))
#    lai=C_foliar/p[16,]
#
#    ###
#    ## declare needed variables
#
#    gpp=array(-9999,dim=c(dim(p)[2],length(mint)+1))
#    nee=array(-9999,dim=c(dim(p)[2],length(mint)+1))
#    respiration_auto=array(-9999,dim=c(dim(p)[2],length(mint)+1))
#    respiration_het_litter=array(-9999,dim=c(dim(p)[2],length(mint)+1))
#    respiration_het_som=array(-9999,dim=c(dim(p)[2],length(mint)+1))
#    wood_production=array(-9999,dim=c(dim(p)[2],length(mint)+1))
#    woodlitter_production=array(-9999,dim=c(dim(p)[2],length(mint)+1))
#    litter2som=array(-9999,dim=c(dim(p)[2],length(mint)+1))
#    extracted_C=array(-9999,dim=c(dim(p)[2],length(mint)+1))
#
#    ###
#    ## Declare some constants
#
#    # Calculating constants for leaffall and labile release factors
#    # release period coefficient, based on duration of labile turnover or leaf
#    # fall durations
#    wf=(2.**0.5)*p[15,]/2.
#    wl=(2.**0.5)*p[13,]/2.
#
#    #! magnitude coefficient
#    ff=(log(p[5,])-log(p[5,]-1.))/2.
#    fl=(log(1+1e-3)-log(1e-3))/2.
#    # set minium labile life span to one year
#    ml=1.+1e-3
#    # offset for labile and leaf turnovers
#    osf=ospolynomial(p[5,],wf)
#    osl=ospolynomial(ml,wl)
#
#    # day in radians
#    sf = 365.25/pi
#
#    # time step of model in days
#    deltat=doy[2]-doy[1]
#
#    # acm constants
#    if (pft == 1) {
#	infile=paste(exepath,"/gpp_emulator_parameters_1.bin",sep="")
#	randomForest=read_binary_response_surface(infile)
#	emulator=read_binary_response_surface(infile)
##	constants=c(14.2924979679255,0.0561498827171238,99.821467204072
##		  ,594.028459378808,8.7949932157796e-05,1.99806271380768
##		  ,0.797783286683111,0.0562354937904277,1.38146269047142e-06
##		  ,-0.809188804123206,-2.04385402571047,1.08684780843903)
#    } else if (pft == 2){
#	infile=paste(exepath,"/gpp_emulator_parameters_2.bin",sep="")
#	randomForest=read_binary_response_surface(infile)
##	emulator=read_binary_response_surface(infile)
##	constants=c(15.381311044568,0.0183200388064739,99.6448728924969
##		  ,599.299262528604,0.000749961643586201,1.99566157991778
##		  ,4.04001469202779,0.0988695360809748,0.00131439603124759
##		  ,0.213656047166534,-1.94699523816983,1.09809814067031)
#    } else if (pft == 3){
#	infile=paste(exepath,"/gpp_emulator_parameters_3.bin",sep="")
#	randomForest=read_binary_response_surface(infile)
##	emulator=read_binary_response_surface(infile)
##	constants=c(10.2316570273157,0.00712986467252934,97.9317562828492
##		  ,562.420236360716,0.0161299065163589,1.97580972754145
##		  ,8.75027052875512,0.141273386628701,0.00258422947037047
##		  ,0.862003227647435,-2.04583637141348,1.09426341985612)
#    } else if (pft == 5){
#	infile=paste(exepath,"/gpp_emulator_parameters_5.bin",sep="")
#	randomForest=read_binary_response_surface(infile)
##	emulator=read_binary_response_surface(infile)
##	constants=c(14.3571506063787,0.00817383750120144,99.456661908902
##		  ,565.615444981272,0.0131673228020551,1.9846931558873
##		  ,8.60847959834574,0.126892411934286,0.00289230125304661
##		  ,1.47135978636525,-2.00184623053322,1.09368398889993)
#    } else if (pft == 11){
#	infile=paste(exepath,"/gpp_emulator_parameters_11.bin",sep="")
#	randomForest=read_binary_response_surface(infile)
##	emulator=read_binary_response_surface(infile)
##	constants=c(16.1247736227685,0.000617309178894364,89.512048631903
##		  ,282.023083528009,0.0301297353658013,1.44748316513338
##		  ,1.71596189896891,0.241803967406994,0.000306360169728313
##		  ,-0.554150638896367,-2.0788751390956,1.07337026956825)
#    } else if (pft == 13){
#	infile=paste(exepath,"/gpp_emulator_parameters_13.bin",sep="")
#	randomForest=read_binary_response_surface(infile)
##	emulator=read_binary_response_surface(infile)
##	constants=c(2.37323087141509,0.0433056473554607,96.0793202668459
##		  ,524.190312195035,1.34018492473325e-06,1.69940113655973
##		  ,1.41439798367338,0.0975902734200319,0.000954291392234617
##		  ,-1.9108709274161,-2.07205814716648,0.956166673229136)
#    } else if (pft == 18){
#	infile=paste(exepath,"/gpp_emulator_parameters_18.bin",sep="")
#	randomForest=read_binary_response_surface(infile)
##	emulator=read_binary_response_surface(infile)
##	constants=c(9.71571673798685,0.00808290472235316,99.3802919097048
##		  ,596.29303517113,0.0269131568567349,1.93737880793734
##		  ,7.23776549676329,0.136836539382728,0.00180138325980396
##		  ,1.25991762411738,-1.93236861131181,1.09477279931996)
#    } else if (pft == 19){
#	infile=paste(exepath,"/gpp_emulator_parameters_19.bin",sep="")
##	emulator=read_binary_response_surface(infile)
#	randomForest=read_binary_response_surface(infile)
##	constants=c(13.8225676305661,0.0100104941805291,98.7337385590505
##		  ,589.032594902438,0.00246376099756583,1.97405597521383
##		  ,8.02661955979173,0.102049963328083,3.97172806458295e-05
##		  ,1.07602065071723,-2.07307691061897,1.09610425681178)
#    }
#
#    # create cluster at the beginning
##    if (use_parallel) {cl <- makeCluster(numWorkers, type = "PSOCK")}
#
#    ###
#    ## MAIN DALEC
#
#    for (step in seq(1,dim(met)[1])) {
#
#	# calculate lai in step
#	lai[,step]=C_foliar[,step]/p[16,]
#
#	# Labile release and leaffall factors for current doy
#	leaffall_factor=(2./(pi**0.5))*(ff/wf)*exp(-((sin((doy[step]-p[14,]+osf)/sf)*sf/wf)**2.))
#	labrelease_factor=(2./(pi**0.5))*(fl/wl)*exp(-((sin((doy[step]-p[11,]+osl)/sf)*sf/wl)**2.))
#
#	# Temperature rate coefficient
#	temprate=exp(p[10,]*(maxt[step]+mint[step])*0.5)
#
#	# run acm
##	gpp[,step]=acm_tessel(lai[,step],maxt[step],mint[step],co2[step],doy[step],lat,radiation[step],constants)
#	# or use GPP loop up table from GAM
##	new_data=c(mint[step],0,maxt[step],radiation[step],co2[step],lagged_precip[step])
##        gpp[,step]=lookup_gpp(lai=lai[,step],new_data=new_data,doy=doy[step],lat=lat,response_matrix=emulator$response_surface,coefficients=emulator$coefficients,nos_fixed=emulator$nos_fixed,maximum=emulator$maximum,minimum=emulator$minimum)
#	# or use the randomForest regression trees as emulator
#	new_data=list(lai=lai[,step],maxt=maxt[step],mint=mint[step],radiation=radiation[step],dayl=0,co2=co2[step],lagged_precip=lagged_precip[step])
#	gpp[,step]=lookup_gpp(randomForest,new_data,lat,doy[step])
#
#	# fluxes / allocation gC.m-2
#	respiration_auto[,step] = p[2,]*gpp[,step]
#	leaf_production = (gpp[,step]-respiration_auto[,step])*p[3,]
#	labile_production = (gpp[,step]-respiration_auto[,step]-leaf_production)*p[12,]
#	root_production = (gpp[,step]-respiration_auto[,step]-leaf_production-labile_production)*p[4,]
#	wood_production = gpp[,step]-respiration_auto[,step]-leaf_production-root_production-labile_production
#
#	# time dependancies
#	labile_release=C_labile[,step]*(1-(1-labrelease_factor)**deltat)/deltat
#	leaflitter_production = C_foliar[,step]*(1-(1-leaffall_factor)**deltat)/deltat
#	woodlitter_production[,step] = C_wood[,step]*(1-(1-p[6,])**deltat)/deltat
#	rootlitter_production = C_root[,step]*(1-(1-p[7,])**deltat)/deltat
#
#	# those with temperature AND time dependancies
#	respiration_het_litter[,step] = C_litter[,step]*(1-(1-temprate*p[8,])**deltat)/deltat
#	respiration_het_som[,step] = C_som[,step]*(1-(1-temprate*p[9,])**deltat)/deltat
#	litter2som[,step] = C_litter[,step]*(1-(1-p[1,]*temprate)**deltat)/deltat
#
#	# Pools:
#	C_labile[,step+1] = C_labile[,step] + labile_production*deltat - labile_release*deltat
#	C_foliar[,step+1] =  C_foliar[,step] + leaf_production*deltat - leaflitter_production*deltat + labile_release*deltat
#	C_wood[,step+1] = C_wood[,step] +  wood_production*deltat - woodlitter_production[,step]*deltat
#	C_root[,step+1] = C_root[,step] + root_production*deltat - rootlitter_production*deltat
#	C_litter[,step+1] = C_litter[,step] + (leaflitter_production + rootlitter_production - respiration_het_litter[,step] - litter2som[,step])*deltat
#	C_som[,step+1] = C_som[,step] + (litter2som[,step] - respiration_het_som[,step]+woodlitter_production[,step])*deltat
#
#    } # end of time loop
#
#    # at end of the job close the parallel operations
# #   if (use_parallel) {stopCluster(cl)}
#
#    ###
#    ## Output some results
#
#    # calculate nee at the end
#    nee=(respiration_auto+respiration_het_litter+respiration_het_som)-gpp
#    # calculate total ecosystem carbon too
#    Biomass=C_foliar+C_wood+C_litter+C_root+C_labile+C_som
#
#    # collect results
#    results=list(nee=nee
#		,gpp=gpp
#		,respiration_auto=respiration_auto
#		,respiration_het_litter=respiration_het_litter
#		,respiration_het_som=respiration_het_som
#		,C_foliar=C_foliar
#		,C_wood=C_wood
#		,C_root=C_root
#		,C_litter=C_litter
#		,C_som=C_som
#		,C_labile=C_labile
#		,lai=lai
#		,Biomass=Biomass
#		,litter2som=litter2som
#		,woodlitter_production=woodlitter_production
#		,extracted_C=extracted_C)
#
#    # how long did this all take
#    print(paste("......DALEC simulation time: ",round(proc.time()["elapsed"]-stime, digits=1)," seconds",sep=""))
#
#    # return results
#    return(results)
#
#}
