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
! Module contains uniform prior parameter information for the DALEC.A1.C8.D2.F2.H2.P1.R4.036 model.
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
  use samplers_shared, only: PARINFO

  implicit none

  ! make all private
  private

  ! specify explicitly the public
  public :: pars_info

  contains
  !
  !------------------------------------------------------------------
  !
  subroutine pars_info(PI)

    ! Subroutine contains a list of parameter ranges for the model.
    ! These could or possibly should go into an alternate file which can be read in.
    ! This may improve the usability when it comes to reading these information
    ! in for different PFTs

    implicit none

    ! Arguments
    type(PARINFO), intent(inout):: PI
    
    PI%npars = 50
    if (.not. allocated(PI%parmin)) allocate(PI%parmin(PI%npars))
    if (.not. allocated(PI%parmax)) allocate(PI%parmax(PI%npars))

    !
    ! declare parameters
    !

    ! Decomposition efficiency of litter/CWD to som (fraction)
    PI%parmin(1) = 0.25d0
    PI%parmax(1) = 0.75d0

    ! Fraction of GPP respired as Rm(fol,root,wood)
    PI%parmin(2) = 0.1d0
    PI%parmax(2) = 0.7d0

    ! Fraction of (1-fgpp) to foliage
    PI%parmin(3) = 0.1d0
    PI%parmax(3) = 0.5d0

    ! Fraction of (1-fgpp) to roots*/
    PI%parmin(4) = 0.1d0
    PI%parmax(4) = 0.80d0

    ! Leaf Lifespan (yr)
    ! Wright et al. 2004
    PI%parmin(5) = 1.001d0
    PI%parmax(5) = 6d0 !8d0

    ! TOR wood* - 1% loss per year value
    PI%parmin(6) = 0.000009d0 ! 304  years
    PI%parmax(6) = 0.001d0    ! 2.74 years

    ! Turnover fraction of roots
    ! Gill and Jackson (2000), New Phytol., 147, 13–31
    ! Fig. 6 turnover by diameter class
    PI%parmin(7) = 0.001368925d0 ! 2    years !0.0006844627d0 ! 4 years
    PI%parmax(7) = 0.01d0        ! 0.27 years

    ! Turnover of foliar litter (fraction; temperature adjusted)
    PI%parmin(8) = 0.0001141d0 ! 24   years at 0oC
    PI%parmax(8) = 0.02d0      ! 0.13 years at 0oC

    ! Turnover of fine root litter (fraction; temperature adjusted)
    PI%parmin(9) = 0.0001141d0 ! 24   years at 0oC
    PI%parmax(9) = 0.02d0      ! 0.13 years at 0oC

    ! Temp factor* = Q10 = 1.2-2.2
    PI%parmin(10) = 0.019d0
    PI%parmax(10) = 0.08d0

    ! Canopy Efficiency
    ! NUE and avN combination give a Vcmax equivalent, the canopy efficiency.
    ! Kattge et al (2011) offers a potential prior range of 3.4 - 30.7 gC/m2leaf/day.
    ! Here, to be cautious we will expand accepted range
    ! Thus CUE = NUE * avN -> 1.64 / 42.0
    ! TLS: 27/10/2021 restricted again based now on 95 %CI (12.61 / 29.68) from TRY
    PI%parmin(11) = 10d0 !5d0
    PI%parmax(11) = 100d0 !42d0 !50d0

    ! max bud burst day
    PI%parmin(12) = 365.25d0
    PI%parmax(12) = 365.25d0*4d0

    ! Fraction to Clab*/
    PI%parmin(13) = 0.01d0
    PI%parmax(13) = 0.5d0

    ! Clab Release period
    PI%parmin(14) = 10d0
    PI%parmax(14) = 100d0

    ! max leaf fall day
    PI%parmin(15) = 365.25d0
    PI%parmax(15) = 365.25d0*4d0

    ! Leaf fall period
    PI%parmin(16) = 20d0
    PI%parmax(16) = 150d0

    ! LMA (gC.m-2)
    ! Kattge et al. 2011
    PI%parmin(17) = 20d0
    PI%parmax(17) = 180d0

    ! fraction of Cwood which is coarse root
    PI%parmin(25) = 0.15d0
    PI%parmax(25) = 0.50d0

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

    ! Turnover rate for wood litter
    PI%parmin(35) = 1.368925d-05 ! 200.00 years at 0oC
    PI%parmax(35) = 0.001d0      !   2.74 years at 0oC

    ! Initial microbial activity 
    PI%parmin(39) = 0.01d0
    PI%parmax(39) = 0.5d0
    ! Inhibition constant for C dependant microbial activity 
    PI%parmin(40) = 75d0
    PI%parmax(40) = 250d0

    ! Foliar lignin fraction 
    PI%parmin(41) = 0.01d0
    PI%parmax(41) = 0.40d0
    ! Fine root lignin fraction 
    PI%parmin(42) = 0.01d0
    PI%parmax(42) = 0.40d0
    ! Wood lignin fraction 
    ! lignin fractions based on Cornwell et al., (2009), Global Change Biology, 15: 2431-2449. https://doi.org/10.1111/j.1365-2486.2009.01916.x
    PI%parmin(43) = 0.15d0
    PI%parmax(43) = 0.40d0

    ! Efficiency of substrate uptake by microbes 
    ! original value from Xenakis & Williams (2014)
    PI%parmin(44) = 0.62d0-(0.62d0*0.75d0)
    PI%parmax(44) = 0.62d0+(0.62d0*0.75d0)

    ! Maximum microbial death rate  (fraction/day)
    ! original value from Xenakis & Williams (2014)
    PI%parmin(45) = 0.24d0-(0.24d0*0.75d0)
    PI%parmax(45) = 0.24d0+(0.24d0*0.75d0)
    ! Inhibition constant for microbial death 
    ! original value from Xenakis & Williams (2014)
    PI%parmin(46) = 0.213d0-(0.213d0*0.75d0)
    PI%parmax(46) = 0.213d0+(0.213d0*0.75d0)

    ! Microbial maintenance respiration coefficient 
    PI%parmin(47) = 0.45d0
    PI%parmax(47) = 0.55d0

    ! Microbial decomposition efficiency
    PI%parmin(48) = 0.01d0
    PI%parmax(48) = 0.04d0

    ! Turnover constant for slow som (fraction/day)
    ! Slow som turnover is assumed to be a based on this fraction
    ! applied to the microbial pool C and scaled by its activity.
    PI%parmin(49) = 0.3688761d0-(0.3688761d0*0.75d0)
    PI%parmax(49) = 0.3688761d0+(0.3688761d0*0.75d0)

    ! 2nd order rate constant for microbial uptake from fast som pool (day-1)
    PI%parmin(50) = 0.1129056d0-(0.1129056d0*0.75d0)
    PI%parmax(50) = 0.1129056d0+(0.1129056d0*0.75d0)

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
    PI%parmin(20) = 1.0d0
    PI%parmax(20) = 2000d0
    ! C_wood
    PI%parmin(21) = 1d0
    PI%parmax(21) = 30000d0

    ! C foliar litter
    PI%parmin(22) = 1d0
    PI%parmax(22) = 2000d0
    ! C som (slow)
    PI%parmin(23) = 200d0
    PI%parmax(23) = 250000d0 

    ! Initial soil water fraction
    PI%parmin(24) = 0.05d0
    PI%parmax(24) = 1.00d0

    ! C wood litter
    PI%parmin(34) = 1d0
    PI%parmax(34) = 10000d0
    ! C fine root litter
    PI%parmin(36) = 1d0
    PI%parmax(36) = 2000d0
    ! C som (fast)
    PI%parmin(37) = 1d0
    PI%parmax(37) = 100d0
    ! C microbial
    PI%parmin(38) = 1d0
    PI%parmax(38) = 100d0

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
