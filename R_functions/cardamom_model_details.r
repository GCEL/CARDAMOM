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
# This function contains basic information about the DALEC models
# This function is based on an original Matlab function development by A. A. Bloom 
# (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).
#
#########################################################################################

cardamom_model_details <-function(modelname,specific_pft,ctessel_pft) {

  if (modelname == "ACM") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(2,dim=c(length(ctessel_pft)))
    nopars   = array(20,dim=c(length(ctessel_pft)))
    nofluxes = array(4,dim=c(length(ctessel_pft)))
    nodiags  = array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="ACM",nopools=nopools,nofluxes=nofluxes,nomet=19+4,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.C3.M1.014" | modelname == "DALEC.14.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(9,dim=c(length(ctessel_pft)))
    nopars   = array(37,dim=c(length(ctessel_pft)))
    nofluxes = array(42,dim=c(length(ctessel_pft)))
    nodiags  = array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.C3.M1.014",shortname="DALEC.14.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A3.C3.H2.M1.015" | modelname == "DALEC.15.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(10,dim=c(length(ctessel_pft)))
    nopars   = array(38,dim=c(length(ctessel_pft)))
    nofluxes = array(48,dim=c(length(ctessel_pft)))
    nodiags  = array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A3.C3.H2.M1.015",shortname="DALEC.15.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A3.H1.M2.016" | modelname == "DALEC.16.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(5,dim=c(length(ctessel_pft)))
    nopars   = array(34,dim=c(length(ctessel_pft)))
    nofluxes = array(57,dim=c(length(ctessel_pft)))
    nodiags  = array(23,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A3.H1.M2.016",shortname="DALEC.16.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A3.H2.M2.017" | modelname == "DALEC.17.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(6,dim=c(length(ctessel_pft)))
    nopars   = array(35,dim=c(length(ctessel_pft)))
    nofluxes = array(57,dim=c(length(ctessel_pft)))
    nodiags  = array(23,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A3.H2.M2.017",shortname="DALEC.17.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.C1.D1.F2.P1.002" | modelname == "DALEC.2.") {
    # Information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(6,dim=c(length(ctessel_pft)))
    nopars   = array(28,dim=c(length(ctessel_pft)))
    nofluxes = array(39,dim=c(length(ctessel_pft)))
    nodiags  = array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.C1.D1.F2.P1.002",shortname="DALEC.2.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H1.P1.003" | modelname == "DALEC.3.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(6,dim=c(length(ctessel_pft)))
    nopars   = array(31,dim=c(length(ctessel_pft)))
    nofluxes = array(51,dim=c(length(ctessel_pft)))
    nodiags  = array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H1.P1.003", shortname="DALEC.3.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P1.004" | modelname == "DALEC.4.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(7,dim=c(length(ctessel_pft)))
    nopars   = array(32,dim=c(length(ctessel_pft)))
    nofluxes = array(51,dim=c(length(ctessel_pft)))
    nodiags  = array(24,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H2.P1.004",shortname="DALEC.4.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A2.C1.D2.F2.H2.P1.020" | modelname == "DALEC.20.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(7,dim=c(length(ctessel_pft)))
    nopars   = array(32,dim=c(length(ctessel_pft)))
    nofluxes = array(51,dim=c(length(ctessel_pft)))
    nodiags  = array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A2.C1.D2.F2.H2.P1.020",shortname="DALEC.20.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A4.C6.D2.F2.H2.P11.031" | modelname == "DALEC.31.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(7,dim=c(length(ctessel_pft)))
    nopars   = array(43,dim=c(length(ctessel_pft)))
    nofluxes = array(51,dim=c(length(ctessel_pft)))
    nodiags  = array(24,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A4.C6.D2.F2.H2.P11.031",shortname="DALEC.31.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H6.P1.R5.032" | modelname == "DALEC.32.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(7,dim=c(length(ctessel_pft)))
    nopars   = array(37,dim=c(length(ctessel_pft)))
    nofluxes = array(51,dim=c(length(ctessel_pft)))
    nodiags  = array(24,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H6.P1.R5.032",shortname="DALEC.32.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A4.C6.D2.F2.H3.P12.033" | modelname == "DALEC.33.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(7,dim=c(length(ctessel_pft)))
    nopars   = array(46,dim=c(length(ctessel_pft)))
    nofluxes = array(51,dim=c(length(ctessel_pft)))
    nodiags  = array(30,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A4.C6.D2.F2.H3.P12.033",shortname="DALEC.33.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P2.018" | modelname == "DALEC.18.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(7,dim=c(length(ctessel_pft)))
    nopars   = array(33,dim=c(length(ctessel_pft)))
    nofluxes = array(51,dim=c(length(ctessel_pft)))
    nodiags  = array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H2.P2.018",shortname="DALEC.18.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P5.021" | modelname == "DALEC.21.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(7,dim=c(length(ctessel_pft)))
    nopars   = array(33,dim=c(length(ctessel_pft)))
    nofluxes = array(51,dim=c(length(ctessel_pft)))
    nodiags  = array(24,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H2.P5.021",shortname="DALEC.21.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P6.022" | modelname == "DALEC.22.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(7,dim=c(length(ctessel_pft)))
    nopars   = array(34,dim=c(length(ctessel_pft)))
    nofluxes = array(51,dim=c(length(ctessel_pft)))
    nodiags  = array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H2.P6.022",shortname="DALEC.22.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H2.P1.R1.005" | modelname == "DALEC.5.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(7,dim=c(length(ctessel_pft)))
    nopars   = array(32,dim=c(length(ctessel_pft)))
    nofluxes = array(52,dim=c(length(ctessel_pft)))
    nodiags  = array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H2.P1.R1.005",shortname="DALEC.5.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P1.R1.006" | modelname == "DALEC.6.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(8,dim=c(length(ctessel_pft)))
    nopars   = array(35,dim=c(length(ctessel_pft)))
    nofluxes = array(57,dim=c(length(ctessel_pft)))
    nodiags  = array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H2.P1.R1.006",shortname="DALEC.6.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P2.R1.007" | modelname == "DALEC.7.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(8,dim=c(length(ctessel_pft)))
    nopars   = array(36,dim=c(length(ctessel_pft)))
    nofluxes = array(57,dim=c(length(ctessel_pft)))
    nodiags  = array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H2.P2.R1.007",shortname="DALEC.7.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P2.R3.019" | modelname == "DALEC.19.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(8,dim=c(length(ctessel_pft)))
    nopars   = array(38,dim=c(length(ctessel_pft)))
    nofluxes = array(57,dim=c(length(ctessel_pft)))
    nodiags  = array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H2.P2.R3.019", shortname="DALEC.19.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.C5.D1.F2.P1.013" | modelname == "DALEC.13.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(4,dim=c(length(ctessel_pft)))
    nopars   = array(21,dim=c(length(ctessel_pft)))
    nofluxes = array(32,dim=c(length(ctessel_pft)))
    nodiags  = array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.C5.D1.F2.P1.013",shortname="DALEC.13.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.D1.F2.001" | modelname == "DALEC.1.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(5,dim=c(length(ctessel_pft)))
    nopars=array(22,dim=c(length(ctessel_pft)))
    nofluxes=array(35,dim=c(length(ctessel_pft)))
    nodiags=array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.D1.F2.001",shortname = "DALEC.1.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.C4.D1.F2.012" | modelname == "DALEC.12.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(3,dim=c(length(ctessel_pft)))
    nopars=array(15,dim=c(length(ctessel_pft)))
    nofluxes=array(28,dim=c(length(ctessel_pft)))
    nodiags=array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.C4.D1.F2.012",shortname = "DALEC.12.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H1.P4.R2.010" | modelname == "DALEC.10.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(43,dim=c(length(ctessel_pft)))
    nofluxes=array(57,dim=c(length(ctessel_pft)))
    nodiags=array(21,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H1.P4.R2.010",shortname = "DALEC.10.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P4.R2.011" | modelname == "DALEC.11.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(43,dim=c(length(ctessel_pft)))
    nofluxes=array(57,dim=c(length(ctessel_pft)))
    nodiags=array(21,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H2.P4.R2.011", shortname="DALEC.11.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P7.R2.023" | modelname == "DALEC.23.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(48,dim=c(length(ctessel_pft)))
    nofluxes=array(57,dim=c(length(ctessel_pft)))
    nodiags=array(23,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H2.P7.R2.023",shortname="DALEC.23.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H4.P1.024" | modelname == "DALEC.24.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(7,dim=c(length(ctessel_pft)))
    nopars   = array(37,dim=c(length(ctessel_pft)))
    nofluxes = array(51,dim=c(length(ctessel_pft)))
    nodiags  = array(24,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H4.P1.024",shortname="DALEC.24.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC...025" | modelname == "DALEC.25.") {
    stop("DALEC.25. has not been assigned")
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(48,dim=c(length(ctessel_pft)))
    nofluxes=array(57,dim=c(length(ctessel_pft)))
    nodiags=array(23,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC...025",shortname="DALEC.25.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A4.C6.D2.F2.H3.P10.026" | modelname == "DALEC.26.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(49,dim=c(length(ctessel_pft)))
    nofluxes=array(53,dim=c(length(ctessel_pft)))
    nodiags=array(35,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A4.C6.D2.F2.H3.P10.026",shortname="DALEC.26.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H2.P3.R1.009" | modelname == "DALEC.9.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(40,dim=c(length(ctessel_pft)))
    nofluxes=array(57,dim=c(length(ctessel_pft)))
    nodiags=array(21,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H2.P3.R1.009",shortname = "DALEC.9.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C2.D2.F2.H1.P3.R1.008" | modelname == "DALEC.8.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(39,dim=c(length(ctessel_pft)))
    nofluxes=array(57,dim=c(length(ctessel_pft)))
    nodiags=array(21,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C2.D2.F2.H1.P3.R1.008",shortname = "DALEC.8.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC_1005" | modelname == "DALEC.27.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(38,dim=c(length(ctessel_pft)))
    nofluxes=array(43,dim=c(length(ctessel_pft)))
    nodiags=array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_1005",shortname="DALEC.27.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC_1005a" | modelname == "DALEC.28.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(8,dim=c(length(ctessel_pft)))
    nopars=array(38,dim=c(length(ctessel_pft)))
    nofluxes=array(43,dim=c(length(ctessel_pft)))
    nodiags=array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC_1005a",shortname="DALEC.28.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H3.P1.029" | modelname == "DALEC.29.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(7,dim=c(length(ctessel_pft)))
    nopars=array(33,dim=c(length(ctessel_pft)))
    nofluxes=array(49,dim=c(length(ctessel_pft)))
    nodiags=array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H3.P1.029",shortname="DALEC.29.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A3.C1.D2.F2.H2.P1.030" | modelname == "DALEC.30.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(7,dim=c(length(ctessel_pft)))
    nopars   = array(38,dim=c(length(ctessel_pft)))
    nofluxes = array(51,dim=c(length(ctessel_pft)))
    nodiags  = array(24,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A3.C1.D2.F2.H2.P1.030",shortname="DALEC.30.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C8.D2.F2.H2.P1.R4.036" | modelname == "DALEC.36.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(11,dim=c(length(ctessel_pft)))
    nopars=array(50,dim=c(length(ctessel_pft)))
    nofluxes=array(67,dim=c(length(ctessel_pft)))
    nodiags=array(20,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C8.D2.F2.H2.P1.R4.036",shortname="DALEC.36.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else if (modelname == "DALEC.A1.C1.D2.F2.H5.P1.037" | modelname == "DALEC.37.") {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools  = array(7,dim=c(length(ctessel_pft)))
    nopars   = array(34,dim=c(length(ctessel_pft)))
    nofluxes = array(51,dim=c(length(ctessel_pft)))
    nodiags  = array(24,dim=c(length(ctessel_pft)))
    cardamom_model_details=list(name="DALEC.A1.C1.D2.F2.H5.P1.037",shortname="DALEC.4.",nopools=nopools,nofluxes=nofluxes,nomet=16,nopars=nopars,nodiags=nodiags)
  } else {
    # No model name was matched
    stop(paste("the inputed model name ('",modelname,"') could not be matched in the available library"))
  } # If modelname == "..."


} # end function cardamom_model_details

## Use byte compile
cardamom_model_details<-cmpfun(cardamom_model_details)
