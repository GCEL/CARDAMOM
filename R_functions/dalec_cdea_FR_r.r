
###
## Function containing the DALEC_CDEA model in R form
###

# This is an R version of the DALEC2 model with disturbance modifications to allow for the imposition of forest rotations management
# The original fortran code was written by A. A. Bloom (UoE, now at JPL) and T. L. Smallman (UoE).
# The R translation and forest rotation management modifications by T. L. Smallman (UoE)

# DRIVERS
# maxt= max daily temperature (oC)
# mint= min daily temperature (oC)
# radiation = sum shortwave radiation (MJ.day-1)
# lat = latidude in degrees
# Ceff = canopy efficiency parameter from DALEC
# yearday = day of year
# co2 = co2 (ppm) often held constant

acm_dalec_cdea_FR <- function (LAI,maxt,mint,co2,yearday,lat,radiation,constants,p11) {
  # The aggregated canopy model, a GPP response function
  gc=0;pp=0;qq=0;ci=0;e0=0;mult=0;dayl=0;cps=0;dec=0;nit=1

  # default constants
  #constants=c(0,0.0156935,4.22273,208.868,0.0453194,0.37836,7.19298,0.011136,2.1001,0.789798,-2,1)

  # determine temperature range
  trange=0.5*(maxt-mint)
  # daily canopy conductance
  gc=abs(constants[11])**(constants[10])/((constants[6]*constants[12]+trange))
  # maximum rate of temperature and nitrogen (canopy efficiency) limited photosynthesis (gC.m-2.day-1)
  #pn=LAI*nit*constants[1]*exp(constants[8]*maxt)
  pn=LAI*nit*p11*exp(constants[8]*maxt)
  # pp and qq represent limitation by diffusion and metabolites respecitively
  pp=pn/gc ; qq=constants[3]-constants[4]
  # calculate internal CO2 concentration (ppm)
  ci=0.5*(co2+qq-pp+((co2+qq-pp)**2-4*(co2*qq-pp*constants[3]))**0.5)
  # limit maximum quantium efficiency by leaf area
  e0=constants[7]*LAI^2/(LAI^2+constants[9])
  # calculate day length (hours)
#  dec = - asin( sin( 23.45 * pi / 180. ) * cos( 2. * pi * ( yearday + 10. ) / 365. ) )
#  sinld = sin( lat*(pi/180.) ) * sin( dec )
#  cosld = cos( lat*(pi/180.) ) * cos( dec )
#  aob = sinld / cosld
#  dayl = 12.0 * ( 1. + 2. * asin( aob ) / pi )
  dec=-23.4*cos((360*(yearday+10.)/365.)*pi/180.)*pi/180.
  mult=tan(lat*pi/180)*tan(dec)
  if (mult>=1) {
    dayl=24.0
  } else if (mult<=-1) {
    dayl=0.
  } else {
    dayl=24.*acos(-mult)/pi
  }
  # calculate CO2 limited rate of photosynthesis
  pd=gc*(co2-ci)
  # calculate combined light and CO2 limited photosynthesis
  cps=e0*radiation*pd/(e0*radiation+pd)
  # correct for day length variation
  model=cps*(constants[2]*dayl+constants[5])
  return(model)
}

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
# p(11) Canopy efficiency parameter - C_eff (part of ACM)
# p(12) = date of Clab release - B_day
# p(13) = Fraction allocated to Clab - f_l
# p(14) = lab release duration period - R_l
# p(15) = date of leaf fall - F_day
# p(16) = leaf fall duration period - R_f
# p(17) = LMA

# C pools initial conditions
# p(18) labile C (gC.m-2)
# p(19) foliar C (gC.m-2)
# p(20) wood C (gC.m-2)
# p(21) root C (gC.m-2)
# p(22) litter C (gC.m-2)
# p(23) som C (gC.m-2)

ospolynomial <- function(L,w){
    # Function calculates the day offset for Labile release and leaf turnover functions
    mxc=c(0.000023599784710,0.000332730053021,0.000901865258885,-0.005437736864888, -0.020836027517787, 0.126972018064287, -0.188459767342504)
    # load log of leaf / labile turnovers
    LLog=log(L-1)
    ospolynomial=(mxc[1]*LLog**6.+ mxc[2]*LLog**5.+ mxc[3]*LLog**4.+ mxc[4]*LLog**3.+ mxc[5]*LLog**2.+ mxc[6]*LLog**1.+ mxc[7])*w
}

dalec_cdea_fr_r <- function(met,p,lat) {

    ###
    ## load met drivers

    mint=met[,2]
    maxt=met[,3]
    radiation=met[,4]
    co2=met[,5]
    doy=met[,6]

    ###
    ## load deforestation parameters

    # variables related to deforestation
    # labile_loss = total loss from labile pool from deforestation
    # foliar_loss = total loss form foliar pool from deforestation
    # roots_loss = total loss from root pool from deforestation
    # wood_loss = total loss from wood pool from deforestation
    # labile_residue = harvested labile remaining in system as residue
    # foliar_residue = harested foliar remaining in system as residue
    # roots_residue = harvested roots remaining in system as residue
    # wood_residue = harvested wood remaining in system as residue
    # wood_pellets = harvested wood removed for pellets
    # coarse_root_residue = expected coarse woody root left in system as residue

    # parameters related to deforestation
    # labile_res_frac = fraction of labile harvest left as residue
    # foliage_res_frac = fraction of foliage harvest left as residue
    # roots_res_frac = fraction of roots harvest left as residue
    # wood_res_frac = fraction of wood harvest left as residue
    # pellet_frac = fraction of wood harvest take for pellets
    # below_wood_frac = fraction of wood pool expected to be coarse root
    # below_wood_frac_res = fraction of coarse root left as residue
    # soil_loss_frac = fraction determining Csom expected to be physically
    # removed along with coarse roots

    # initial values for deforestation variables
    labile_loss = 0.    ; foliar_loss = 0.
    roots_loss = 0.     ; wood_loss = 0.
    labile_residue = 0. ; foliar_residue = 0.
    roots_residue = 0.  ; wood_residue = 0.
    wood_pellets = 0.	; C_total = 0.
    reforest_day = 0.
    soil_loss_with_roots = 0.
    coarse_root_residue = 0.

    ## default Duke run
    # parameter values for deforestation variables
    # Assumes all non-wood is left behind.
    # wood residue is assuming all branches and 'twigs' are left
    # (Ares & Brauers 2005)
#    foliage_res_frac = 0
#    roots_res_frac = 0
#    below_wood_frac = 0.32 # Adegbidi et al (2005)
#    wood_res_frac = 0 # 0.36 # sd = +/- 4 % applies to stem wood
#    below_wood_frac_res = 0
#    pellet_frac = 0
    # assumes Csom loss due to phyical removal with roots Andy Moffat et al
    # (2011) Forestry Commission Research Note
#    soil_loss_frac = 0.02 # actually between 1-3 %
    ## scen 1
    # parameter values for deforestation variables
    # Assumes all non-wood is left behind.
    # wood residue is assuming all branches are left
    # Table 4 (Ares & Brauers 2005)
    # > 15 years 15 % branches ; <  15 years 25 % branches
    foliage_res_frac = 1.0
    roots_res_frac = 1.0
    below_wood_frac = 0.32 # Adegbidi et al (2005)
    wood_res_frac = 0.20 # 0.36 # sd = +/- 4 % applies to stem wood
    below_wood_frac_res = 1.
    pellet_frac = 0
    # assumes Csom loss due to phyical removal with roots
    # Morison et al (2012) Forestry Commission Research Note
    soil_loss_frac = 0.02 # actually between 1-3 %
    ## scen 2
    # parameter values for deforestation variables
    # Assumes all non-wood is left behind.
    # wood residue is assuming all branches and 'twigs' are left
    # (Ares & Brauers 2005)
    foliage_res_frac = 0.5
    roots_res_frac = 1.0
    below_wood_frac = 0.32 # Adegbidi et al (2005)
    wood_res_frac = 0.20 # sd = +/- 4 % applies to stem wood
    below_wood_frac_res = 1.
    pellet_frac = 1
    # assumes Csom loss due to phyical removal with roots Andy Moffat et al
    # (2011) Forestry Commission Research Note
    soil_loss_frac = 0.02 # actually between 1-3 %
    ## scen 3
#    # parameter values for deforestation variables
#    # Assumes all non-wood is left behind.
#    # wood residue is assuming all branches and 'twigs' are left
#    # (Ares & Brauers 2005)
    foliage_res_frac = 0.5
    roots_res_frac = 1.0
    below_wood_frac = 0.32 # Adegbidi et al (2005)
    wood_res_frac = 0.20 # sd = +/- 4 % applies to stem wood
    below_wood_frac_res = 0
    pellet_frac = 1
    # assumes Csom loss due to phyical removal with roots Andy Moffat et al
    # (2011) Forestry Commission Research Note
    soil_loss_frac = 0.02 # actually between 1-3 %

    ###
    ## load initial conditions

    C_labile=array(p[18,],dim=c(dim(p)[2],length(mint)+1))
    C_foliar=array(p[19,],dim=c(dim(p)[2],length(mint)+1))
    C_root=array(p[20,],dim=c(dim(p)[2],length(mint)+1))
    C_wood=array(p[21,],dim=c(dim(p)[2],length(mint)+1))
    C_litter=array(p[22,],dim=c(dim(p)[2],length(mint)+1))
    C_som=array(p[23,],dim=c(dim(p)[2],length(mint)+1))
    lai=C_foliar/p[17,]

    ###
    ## declare needed variables

    gpp=array(-9999,dim=c(dim(p)[2],length(mint)+1))
    nee=array(-9999,dim=c(dim(p)[2],length(mint)+1))
    respiration_auto=array(-9999,dim=c(dim(p)[2],length(mint)+1))
    respiration_het_litter=array(-9999,dim=c(dim(p)[2],length(mint)+1))
    respiration_het_som=array(-9999,dim=c(dim(p)[2],length(mint)+1))
    wood_production=array(-9999,dim=c(dim(p)[2],length(mint)+1))
    woodlitter_production=array(-9999,dim=c(dim(p)[2],length(mint)+1))
    litter2som=array(-9999,dim=c(dim(p)[2],length(mint)+1))
    extracted_C=array(0,dim=c(dim(p)[2],length(mint)+1))

    ###
    ## Declare some constants

    # Calculating constants for leaffall and labile release factors
    # release period coefficient, based on duration of labile turnover or leaf
    # fall durations
    wf=(2.**0.5)*p[16,]/2.
    wl=(2.**0.5)*p[14,]/2.

    #! magnitude coefficient
    ff=(log(p[5,])-log(p[5,]-1.))/2.
    fl=(log(1+1e-3)-log(1e-3))/2.
    # set minium labile life span to one year
    ml=1.+1e-3
    # offset for labile and leaf turnovers
    osf=ospolynomial(p[5,],wf)
    osl=ospolynomial(ml,wl)

    # day in radians
    sf = 365.25/pi

    # time step of model in days
    deltat=doy[2]-doy[1]

    ###
    ## MAIN DALEC

    for (step in seq(1,dim(met)[1])) {

	# Labile release and leaffall factors for current doy
	leaffall_factor=(2./(pi**0.5))*(ff/wf)*exp(-((sin((doy[step]-p[15,]+osf)/sf)*sf/wf)**2.))
	labrelease_factor=(2./(pi**0.5))*(fl/wl)*exp(-((sin((doy[step]-p[12,]+osl)/sf)*sf/wl)**2.))
#	leaffall_factor=(2./(pi**0.5))*(ff/wf)*exp(-((sin((met[step,1]-p[15]+osf)/sf)*sf/wf)**2.))
#	labrelease_factor=(2./(pi**0.5))*(fl/wl)*exp(-((sin((met[step,1]-p[12]+osl)/sf)*sf/wl)**2.))

	# Temperature rate coefficient
	temprate=exp(p[10,]*(maxt[step]+mint[step])*0.5)

	# acm constants
	constants=c(0,0.0156935,4.22273,208.868,0.0453194,0.37836,7.19298,0.011136,2.1001,0.789798,-2,1)
	gpp[,step]=acm_dalec_cdea_FR(lai[,step],maxt[step],mint[step],co2[step],doy[step],lat,radiation[step],constants,p[11,])
	# ACM generates a NaN value for any situation of lai = 0 therefore we need to check against this
	gpp[which(lai[,step] == 0),step] = 0.

	# fluxes / allocation gC.m-2
	respiration_auto[,step] = p[2,]*gpp[,step]
	leaf_production = (gpp[,step]-respiration_auto[,step])*p[3,]
	labile_production = (gpp[,step]-respiration_auto[,step]-leaf_production)*p[13,]
	root_production = (gpp[,step]-respiration_auto[,step]-leaf_production-labile_production)*p[4,]
	wood_production = gpp[,step]-respiration_auto[,step]-leaf_production-root_production-labile_production

	# time dependancies
	labile_release=C_labile[,step]*(1-(1-labrelease_factor)**deltat)/deltat
	leaflitter_production = C_foliar[,step]*(1-(1-leaffall_factor)**deltat)/deltat
	woodlitter_production[,step] = C_wood[,step]*(1-(1-p[6,])**deltat)/deltat
	rootlitter_production = C_root[,step]*(1-(1-p[7,])**deltat)/deltat

	# those with temperature AND time dependancies
	respiration_het_litter[,step] = C_litter[,step]*(1-(1-temprate*p[8,])**deltat)/deltat
	respiration_het_som[,step] = C_som[,step]*(1-(1-temprate*p[9,])**deltat)/deltat
	litter2som[,step] = C_litter[,step]*(1-(1-p[1,]*temprate)**deltat)/deltat

	# Pools:
	C_foliar[,step+1] =  C_foliar[,step] + leaf_production*deltat - leaflitter_production*deltat + labile_release*deltat
	C_wood[,step+1] = C_wood[,step] +  wood_production*deltat - woodlitter_production[,step]*deltat
	C_root[,step+1] = C_root[,step] + root_production*deltat - rootlitter_production*deltat
	C_litter[,step+1] = C_litter[,step] + (leaflitter_production + rootlitter_production - respiration_het_litter[,step] - litter2som[,step])*deltat
	C_som[,step+1] = C_som[,step] + (litter2som[,step] - respiration_het_som[,step]+woodlitter_production[,step])*deltat
	C_labile[,step+1] = C_labile[,step] + labile_production*deltat - labile_release*deltat
	lai[,step+1]=C_foliar[,step]/p[17,]

	#
	# deal first with deforestation
	#

	if (step == reforest_day) {
	    C_labile[,step+1] = p[18,]
	    C_foliar[,step+1] = p[19,]
	    C_root[,step+1] = p[20,]
	    C_wood[,step+1] = p[21,]
	}

	if (met[step,8] > 0.) {

	  # assume that labile C is proportionally distributed through biomass
	  labile_res_frac = C_foliar[,step+1] + C_root[,step+1] + C_wood[,step+1]
	  # need to deal with below / above ground wood issue here too
	  wood_residue = ( (C_wood[,step+1]*below_wood_frac) / labile_res_frac )
	  labile_res_frac = (( (C_foliar[,step+1] / labile_res_frac) * foliage_res_frac )
			  + ( (C_root[,step+1]   / labile_res_frac) * roots_res_frac )
			  + wood_residue + ( ((C_wood[,step+1]-(C_wood[,step+1]*below_wood_frac)) / labile_res_frac) * wood_res_frac ))

	  # loss of carbon from each pools
	  labile_loss = C_labile[,step+1]*met[step,8]
	  foliar_loss = C_foliar[,step+1]*met[step,8]
	  roots_loss  = C_root[,step+1]*met[step,8]
	  wood_loss   = C_wood[,step+1]*met[step,8]
	  # transfer fraction of harvest waste to litter or som pools
	  labile_residue = C_labile[,step+1]*met[step,8]*labile_res_frac
	  foliar_residue = C_foliar[,step+1]*met[step,8]*foliage_res_frac
	  roots_residue  = C_root[,step+1]*met[step,8]*roots_res_frac
	  coarse_root_residue = C_wood[,step+1]*met[step,8]*below_wood_frac*below_wood_frac_res
	  soil_loss_with_roots = C_wood[,step+1]*met[step,8]*below_wood_frac*(1-below_wood_frac_res)*soil_loss_frac
	  wood_pellets   = ((C_wood[,step+1]*met[step,8])-coarse_root_residue) * (wood_res_frac*pellet_frac)
	  wood_residue   = coarse_root_residue + ( ((C_wood[,step+1]*met[step,8])-coarse_root_residue) * (wood_res_frac*(1.-pellet_frac)) )
	  # update living pools directly
	  C_labile[,step+1] = C_labile[,step+1]-labile_loss
	  C_foliar[,step+1] = C_foliar[,step+1]-foliar_loss
	  C_root[,step+1] = C_root[,step+1]-roots_loss
	  C_wood[,step+1] = C_wood[,step+1]-wood_loss
	  # update residue pools
	  C_litter[,step+1] = C_litter[,step+1] + labile_residue+foliar_residue+roots_residue
	  C_som[,step+1] = C_som[,step+1] + wood_residue - soil_loss_with_roots

	  # total woody biomass extracted
#          extracted_C[,step+1] = (wood_pellets+ (wood_loss-wood_residue)
#		 			       + (labile_loss-labile_residue)
#                                              + (foliar_loss-foliar_residue)
#                                              + (roots_loss-roots_residue))

	  # if total clearance occured then we need to ensure some minimum
	  # values
	  if (met[step,8] > 0.99) {
	      reforest_day = min(step+365,dim(met)[1])
#	      C_labile[,step+1] = p[18,]
#	      C_foliar[,step+1] = p[19,]
#	      C_root[,step+1] = p[20,]
#	      C_wood[,step+1] = p[21,]
	  } # if total clearance

	  lai[,step+1]=C_foliar[,step+1]/p[17,]

	} # end deforestation info

	#
	# then deal with fire
	#

#        if (met[step,9] > 0.) {

#        } # end burnst area issues

    } # end of time loop

    ###
    ## Output some results

    # calculate nee at the end
    nee=(respiration_auto+respiration_het_litter+respiration_het_som)-gpp
    # calculate total ecosystem carbon too
    Biomass=C_foliar+C_wood+C_litter+C_root+C_labile+C_som

    # collect results
    results=list(nee=nee
		,gpp=gpp
		,respiration_auto=respiration_auto
		,respiration_het_litter=respiration_het_litter
		,respiration_het_som=respiration_het_som
		,C_foliar=C_foliar
		,C_wood=C_wood
		,C_root=C_root
		,C_litter=C_litter
		,C_som=C_som
		,C_labile=C_labile
		,lai=lai
		,Biomass=Biomass
		,litter2som=litter2som
		,woodlitter_production=woodlitter_production
		,extracted_C=extracted_C)
    # return results
    return(results)

} # dalec_cdea_fr_r end
