!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! CARbon DAta MOdel fraMework (CARDAMOM) and DALEC terrestrial ecosystem model suite
! CARDAMOM is a Bayesian model-data fusion software framework. CARDAMOM is used to 
! assimilate observations and ecological theory to retrieve parameters for the 
! DALEC suite of intermediate complexity terrestrial ecosystem models. DALEC can be
! used as a fully integrated component of CARDAMOM or independently. 
! Copyright (C) 2024  University of Edinburgh,
!                     Mathew Williams (mat.williams@ed.ac.uk), 
!                     T. Luke Smallman (t.l.smallman@ed.ac.uk)
! UoE = University of Edinburgh

! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.

! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <https://www.gnu.org/licenses/>.

!!!!!!!!!!!! File specific description !!!!!!!!!!
! Module contains uniform prior parameter information for the DALEC.A1.C2.D2.F2.H2.P7.R2 model.
!
! This code is based on the original C verion of the University of Edinburgh
! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
! All code translation into Fortran, integration into the University of
! Edinburgh CARDAMOM code and subsequent modifications by:
! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
! See function / subroutine specific comments for exceptions and contributors
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module MODEL_PARAMETERS

  implicit none

  ! make all private
  private

  ! specify explicitly the public
  public :: pars_info

  contains

  !
  !------------------------------------------------------------------
  !
  subroutine pars_info
    use MCMCOPT, only: PI
    use cardamom_structures, only: DATAin

    ! Subroutine contains a list of parameter ranges for the model.
    ! These could or possibly should go into an alternate file which can be read in.
    ! This may improve the usability when it comes to reading these information
    ! in for different PFTs

    implicit none

    !
    ! declare parameters
    !

    ! Decomposition efficiency of litter/CWD to som (fraction)
    PI%parmin(1) = 0.25d0
    PI%parmax(1) = 0.75d0

    ! Fraction of GPP respired as Rm(fol,root,wood)
    PI%parmin(2) = 0.1d0
    PI%parmax(2) = 0.7d0

    ! Initial canopy growing index value 
    PI%parmin(3) = 0d0
    PI%parmax(3) = 1d0

    ! Fraction of (1-fgpp) to roots*/
    PI%parmin(4) = 0.1d0
    PI%parmax(4) = 0.80d0

    ! Potential leaf turnover rate
    PI%parmin(5) = 0.002737851d0 ! 1 year
    PI%parmax(5) = 0.016666667d0 ! 60 days

    ! Turnover fraction of wood
    PI%parmin(6) = 0.000009d0 ! 304  years
    PI%parmax(6) = 0.001d0    ! 2.74 years

    ! Turnover fraction of roots
    ! Gill and Jackson (2000), New Phytol., 147, 13–31
    ! Fig. 6 turnover by diameter class
    PI%parmin(7) = 0.001368925d0 ! 2    years 
    PI%parmax(7) = 0.01d0        ! 0.27 years

    ! Turnover of litter (fraction; temperature adjusted)
    PI%parmin(8) = 0.0001141d0 ! 24   years at 0oC
    PI%parmax(8) = 0.02d0      ! 0.13 years at 0oC

    ! Turnover of som to Rhet (fraction; temperature adjusted)
    PI%parmin(9) = 1.368925d-06   ! 2000 years at 0oC
    PI%parmax(9) = 9.126169d-05   !   30 years at 0oC !0.0001368926d0 !   20 years at 0oC

    ! Exponential coefficient for Rhet temperature response
    ! Temp factor* = Q10 = 1.2-2.2
    PI%parmin(10) = 0.019d0
    PI%parmax(10) = 0.08d0

    ! log10 avg foliar N (gN.m-2)
    ! Kattge et al., (2011) (Quantiles 25% / 75%)
    ! and Thomas et al., (2019) (Aconite canopy paper)
    PI%parmin(11) = 0.07918125d0 ! restricted to 1.2 gN/m2leaf
    PI%parmax(11) = 0.4771213d0  ! restricted to 3 gC/m2leaf

    ! Max labile turnover fraction to foliage
    PI%parmin(12) = 0.002737851d0 !  1 years
    PI%parmax(12) = 0.025d0       ! 40 days

    ! Fraction of GPP to Clab*/
    PI%parmin(13) = 0.05d0
    PI%parmax(13) = 0.5d0

    ! Canopy CGI phenology gradient threshold
    ! THIS COULD POSSIBLY BE HARDCODED TO ZERO
    PI%parmin(14) = -0.01d0
    PI%parmax(14) =  0.01d0
    ! NCCE return on new Cfol investment (gCperGPP per gCnewfol)
    PI%parmin(15) = 0.005d0
    PI%parmax(15) = 0.05d0

    ! Turnover rate for CWD
    PI%parmin(16) = 1.368925d-05 ! 200.00 years at 0oC
    PI%parmax(16) = 0.001d0      !   2.74 years at 0oC

    ! LMA
    ! Kattge et al. 2011,
    PI%parmin(17) = 20d0
    PI%parmax(17) = 180d0

    ! fraction of Cwood which is coarse root
    PI%parmin(25) = 0.15d0
    PI%parmax(25) = 0.50d0 ! increased based on evidence of savannah system 50 % below !0.30d0

    ! BUCKET - coarse root biomass (i.e. gbio/m2 not gC/m2) needed to reach 50 %
    ! of max depth
    PI%parmin(26) = 100d0
    PI%parmax(26) = 2500d0 !500d0
    ! BUCKET - maximum rooting depth
    PI%parmin(27) = 0.35d0
    PI%parmax(27) = 20d0

    ! Resilience factor for burned but not combusted C stocks
    PI%parmin(28) = 0.01d0
    PI%parmax(28) = 0.99d0
    ! Combustion completeness factor for foliage
    PI%parmin(29) = 0.01d0
    PI%parmax(29) = 0.99d0
    ! Combustion completeness factor for fine root and wood
    PI%parmin(30) = 0.01d0
    PI%parmax(30) = 0.99d0
    ! Combustion completeness factor for soil
    PI%parmin(31) = 0.01d0
    PI%parmax(31) = 0.1d0
    ! Combustion completeness factor for foliage + fine root litter
    PI%parmin(32) = 0.01d0
    PI%parmax(32) = 0.99d0
    ! Combustion completeness factor for wood litter
    PI%parmin(33) = 0.01d0
    PI%parmax(33) = 0.99d0

    ! CGI minimum temperature (oC)
    PI%parmin(34) = -10d0
    PI%parmax(34) =  5d0
    ! CGI Optimum temperature (oC)
    PI%parmin(35) = 10d0
    PI%parmax(35) = 40d0
    ! CGI temperature maximum (oC)
    PI%parmin(36) = 20d0
    PI%parmax(36) = 60d0
    ! CGI temperature kurtosis
    ! Larger number means more peaky
    PI%parmin(37) = 0.01d0
    PI%parmax(37) = 0.3d0
 
    ! wSWP water potential (MPa) at which foliar growth is fully suppressed
    PI%parmin(38) = -5d0
    PI%parmax(38) = -2d0 ! This is approximately the currently hardcoded minLWP
    ! wSWP water potential (MPa) at which foliar growth suppression begins
    PI%parmin(39) = -5d0
    PI%parmax(39) =  0d0

    ! Nitrogen use efficiency (gC/gN per m2 per day at optimum temperature)
    ! Derived from Vcmax reported in Wullschleger (1993), Journal of
    ! Experimental Botany, Vol 44, No. 262, pp. 907-920.
    ! ~40 gC/gN/day
    ! TRY database equivalent 2.5 % = 1.648512; 97.5 % = 19.906560 gC/gN/m2/day
    ! Xu et al., (2017):
    ! Variations of leaf longevity in tropical moist forests predicted by a
    ! trait-driven carbon optimality model,
    ! Ecology Letters, doi: 10.1111/ele.12804, upper value of 82 gC/gN/day
    PI%parmin(41) =  1.6d0
    PI%parmax(41) = 80.0d0

    ! Reich et al., 2008 - Leaf N linked respiration exponential coefficient
    PI%parmin(42) = 0.935d0 ! 1.639-0.01
    PI%parmax(42) = 1.774d0 ! 1.639+0.01
    ! Reich et al., 2008 - Leaf N linked respiration intercept
    ! max/min values based on observed ranges from Reich et al (2008)
    ! Figure 1
    PI%parmin(43) = 0.10d0 !0.645
    PI%parmax(43) = 1.10d0 !0.911

    ! Parameters linking the NCCE to the CMI
    ! via a logistic function. This is the NCCE
    ! at which the CMI is suppressed by 0.5
    PI%parmin(44) = -0.5d0!0.0001d0 
    PI%parmax(44) =  0.5d0!0.5d0
    ! Parameters linking the NCCE to the CMI
    ! via a logistic function. This is the gradient term
    ! of the logistic function    
    PI%parmin(45) = 10d0
    PI%parmax(45) = 1000d0

    ! Absolute minimum temperature at which no leaf growth occurs (oC)
    PI%parmin(46) = -5d0
    PI%parmax(46) =  5d0
    ! Degree C above p46 at which this limitation is at 50%
    PI%parmin(47) = 0.01d0
    PI%parmax(47) = 2d0

    ! Parameters linking the NCCE gradient to the CMI
    ! via a logistic function. This is the gradient term
    ! of the logistic function
    PI%parmin(48) = 10d0
    PI%parmax(48) = 5000d0

    ! Initial NCCE (gC/gCleaf/day) reference value for gradient calculations
    PI%parmin(49) = -0.10d0
    PI%parmax(49) =  0.10d0

    ! Parameters linking the NCCE gradient to the CMI
    ! via a logistic function. This is the NCCE gradient 
    ! 50 % value of the logistic function
    PI%parmin(50) = -0.00001d0
    PI%parmax(50) = -0.005d0

    !
    ! INITIAL VALUES DECLARED HERE
    !

    ! C labile
    PI%parmin(18) = 1d0
    PI%parmax(18) = 2000d0

    ! C foliar
    PI%parmin(19) = 1d0
    PI%parmax(19) = 2000d0

    ! C roots
    PI%parmin(20) = 1d0
    PI%parmax(20) = 2000d0

    ! C_wood
    PI%parmin(21) = 1d0
    PI%parmax(21) = 30000d0

    ! C litter
    PI%parmin(22) = 1d0
    PI%parmax(22) = 2000d0

    ! C_som
    PI%parmin(23) = 200d0
    PI%parmax(23) = 250000d0 !90000d0

    ! C CWD
    PI%parmin(24) = 1d0
    PI%parmax(24) = 10000d0

    ! Initial soil water fraction
    PI%parmin(40) = 0.05d0
    PI%parmax(40) = 1.00d0

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
