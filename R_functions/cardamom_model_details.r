
###
## Function which contains basic information about the models
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

cardamom_model_details <-function(modelname,specific_pft,ctessel_pft) {

  if (modelname == "ACM") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(2,dim=c(length(ctessel_pft)))
    nopars=array(20,dim=c(length(ctessel_pft)))
    nofluxes=array(4,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="ACM",nopools=nopools,nofluxes=nofluxes,nomet=19+4,nopars=nopars)
  } else if (modelname == "DALEC.C3.M1.#" | modelname == "DALEC.14.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(35,dim=c(length(ctessel_pft)))
    nofluxes=array(21,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.C3.M1.#",shortname="DALEC.14.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A3.C3.H2.M1.#" | modelname == "DALEC.15.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(10,dim=c(length(ctessel_pft)))
    nopars=array(38,dim=c(length(ctessel_pft)))
    nofluxes=array(41,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A3.C3.H2.M1.#",shortname="DALEC.15.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.C1.D1.F2.P1.#" | modelname == "DALEC.2.") {
    # Information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(6,dim=c(length(ctessel_pft)))
    nopars=array(28,dim=c(length(ctessel_pft)))
    nofluxes=array(39,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.C1.D1.F2.P1.#",shortname="DALEC.2.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H1.P1.#" | modelname == "DALEC.3.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(6,dim=c(length(ctessel_pft)))
    nopars=array(28,dim=c(length(ctessel_pft)))
    nofluxes=array(39,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H1.P1.#", shortname="DALEC.3.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P1.#" | modelname == "DALEC.4.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(32,dim=c(length(ctessel_pft)))
    nofluxes=array(45,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H2.P1.#",shortname="DALEC.4.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A2.C1.D2.F2.H2.P1.#" | modelname == "DALEC.20.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(32,dim=c(length(ctessel_pft)))
    nofluxes=array(45,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A2.C1.D2.F2.H2.P1.#",shortname="DALEC.20.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P2.#" | modelname == "DALEC.18.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(33,dim=c(length(ctessel_pft)))
    nofluxes=array(40,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H2.P2.#",shortname="DALEC.18.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P5.#" | modelname == "DALEC.21.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(33,dim=c(length(ctessel_pft)))
    nofluxes=array(40,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H2.P5.#",shortname="DALEC.21.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P6.#" | modelname == "DALEC.22.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(34,dim=c(length(ctessel_pft)))
    nofluxes=array(40,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H2.P6.#",shortname="DALEC.22.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P1.R1.#" | modelname == "DALEC.5.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(32,dim=c(length(ctessel_pft)))
    nofluxes=array(40,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H2.P1.R1.#",shortname="DALEC.5.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P1.R1.#" | modelname == "DALEC.6.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(35,dim=c(length(ctessel_pft)))
    nofluxes=array(45,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H2.P1.R1.#",shortname="DALEC.6.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P2.R1.#" | modelname == "DALEC.7.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(36,dim=c(length(ctessel_pft)))
    nofluxes=array(45,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H2.P2.R1.#",shortname="DALEC.7.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P2.R3.#" | modelname == "DALEC.19.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(38,dim=c(length(ctessel_pft)))
    nofluxes=array(45,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H2.P2.R3.#", shortname="DALEC.19.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.C5.D1.F2.P1.#" | modelname == "DALEC.13.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(4,dim=c(length(ctessel_pft)))
    nopars=array(21,dim=c(length(ctessel_pft)))
    nofluxes=array(32,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.C5.D1.F2.P1.#",shortname="DALEC.13.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.D1.F2.#" | modelname == "DALEC.1.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(5,dim=c(length(ctessel_pft)))
    nopars=array(22,dim=c(length(ctessel_pft)))
    nofluxes=array(35,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.D1.F2.#",shortname = "DALEC.1.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.C4.D1.F2.#" | modelname == "DALEC.12.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(3,dim=c(length(ctessel_pft)))
    nopars=array(15,dim=c(length(ctessel_pft)))
    nofluxes=array(28,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.C4.D1.F2.#",shortname = "DALEC.12.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H1.P4.R2.#" | modelname == "DALEC.10.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(48,dim=c(length(ctessel_pft)))
    nofluxes=array(25,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H1.P4.R2.#",shortname = "DALEC.10.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P4.R2.#" | modelname == "DALEC.11.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(49,dim=c(length(ctessel_pft)))
    nofluxes=array(25,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H2.P4.R2.#", shortname="DALEC.11.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P7.R2.#" | modelname == "DALEC.23.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(48,dim=c(length(ctessel_pft)))
    nofluxes=array(45,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H2.P7.R2.#",shortname="DALEC.23.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P8.R2.#" | modelname == "DALEC.24.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(51,dim=c(length(ctessel_pft)))
    nofluxes=array(45,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H2.P8.R2.#",shortname="DALEC.24.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P9.R2.#" | modelname == "DALEC.25.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(49,dim=c(length(ctessel_pft)))
    nofluxes=array(45,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H2.P9.R2.#",shortname="DALEC.25.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P10.R2.#" | modelname == "DALEC.26.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(48,dim=c(length(ctessel_pft)))
    nofluxes=array(25,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H2.P10.R2.#",shortname="DALEC.26.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P3.R1.#" | modelname == "DALEC.9.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(46,dim=c(length(ctessel_pft)))
    nofluxes=array(25,dim=c(length(ctessel_pft)))
    if (specific_pft == "pft_specific") {nopars[which(ctessel_pft == 1)]=38 ; nofluxes[which(ctessel_pft == 1)]=21 ; nopools[which(ctessel_pft == 1)]=9}
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H2.P3.R1.#",shortname = "DALEC9.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H1.P3.R1.#" | modelname == "DALEC.8.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(43,dim=c(length(ctessel_pft)))
    nofluxes=array(25,dim=c(length(ctessel_pft)))
    if (specific_pft == "pft_specific") {nopars[which(ctessel_pft == 1)]=35 ; nofluxes[which(ctessel_pft == 1)]=21 ; nopools[which(ctessel_pft == 1)]=8}
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H1.P3.R1.#",shortname = "DALEC.8.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_1005" | modelname == "DALEC.27.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(38,dim=c(length(ctessel_pft)))
    nofluxes=array(43,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_1005",shortname="DALEC.27.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC_1005a" | modelname == "DALEC.28.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(38,dim=c(length(ctessel_pft)))
    nofluxes=array(43,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_1005a",shortname="DALEC.28.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H3.P1.#" | modelname == "DALEC.29.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(33,dim=c(length(ctessel_pft)))
    nofluxes=array(45,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H3.P1.#",shortname="DALEC.29.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.A3.C1.D2.F2.H2.P1.#" | modelname == "DALEC.30.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(32+6,dim=c(length(ctessel_pft)))
    nofluxes=array(45,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A3.C1.D2.F2.H2.P1.#",shortname="DALEC.30.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars)
  } else if (modelname == "DALEC.M2.#" | modelname == "DALEC_GRASS") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(6,dim=c(length(ctessel_pft)))
    nopars=array(34,dim=c(length(ctessel_pft)))
    nofluxes=array(23,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.M2.#",shortname="DALEC_GRASS",nopools=nopools,nofluxes=nofluxes,nomet=14,nopars=nopars)
  } else {
    # No model name was matched
    stop(paste("the inputed model name ('",modelname,"') could not be matched in the available library"))
  } # If modelname == "..."


} # end function cardamom_model_details

## Use byte compile
cardamom_model_details<-cmpfun(cardamom_model_details)
