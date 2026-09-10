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
! Module contains uniform prior parameter information for the DALEC.A3.H2.M2.017 model.
!
! This code is based on the original C verion of the University of Edinburgh
! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
! All code translation into Fortran, integration into the University of
! Edinburgh CARDAMOM code and subsequent modifications by:
! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
! S. Zhu (University of Edinburgh, now University of Southampton)
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

    !
    ! declare parameters
    !

    type(PARINFO), intent(inout):: PI

    PI%npars = 35
    if (.not. allocated(PI%parmin)) allocate(PI%parmin(PI%npars))
    if (.not. allocated(PI%parmax)) allocate(PI%parmax(PI%npars))

    ! Decomposition rate [1e-5, 0.01]
    PI%parmin(1) = 0.001d0 
    PI%parmax(1) = 0.1d0 

    ! GPP to resp fraction [~0.54]
    PI%parmin(2) = 0.20d0 
    PI%parmax(2) = 0.70d0 

    ! Canopy GSI phenology gradient threshold
    PI%parmin(3) = -1d-3
    PI%parmax(3) =  1d-2

    ! NPP belowground allocation exponential parameter [0.01, 1.00]
    ! i.e. fine roots
    PI%parmin(4) = 0.01d0 
    PI%parmax(4) = 1.0d0 

    ! Potential leaf turnover rate
    PI%parmin(5) = 0.002737851d0 ! 1 year
    PI%parmax(5) = 0.016666667d0 ! 60 days

    ! Turnover fraction of roots
    ! Gill and Jackson (2000), New Phytol., 147, 13–31
    ! Fig. 6 turnover by diameter class
    PI%parmin(6) = 0.001368925d0 ! 2    years !0.0006844627d0 ! 4 years
    PI%parmax(6) = 0.01d0        ! 0.27 years

    ! TOR litter [0.0001, 0.01]
    PI%parmin(7) = 0.001d0 
    PI%parmax(7) = 0.1d0 

    ! Turnover of som to Rhet (fraction; temperature adjusted)
    PI%parmin(8) = 1.368925d-06   ! 2000 years at 0oC
    PI%parmax(8) = 9.126169d-05   !   30 years at 0oC !0.0001368926d0 !   20 years at 0oC

    ! Exponential coefficient for Rhet temperature response
    ! Temp factor* = Q10 = 1.2-2.2
    PI%parmin(9) = 0.019d0 
    PI%parmax(9) = 0.08d0 

    ! Potential labile turnover fraction to foliage
    PI%parmin(10) = 0.002737851d0*2d0 !  6 months
    PI%parmax(10) = 0.025d0           ! 40 days

    ! Canopy Efficiency
    ! NUE and avN combination give a Vcmax equivalent, the canopy efficiency.
    ! Kattge et al (2011) offers a potential prior range of 3.4 - 30.7 gC/m2leaf/day.
    ! Here, to be cautious we will expand accepted range
    ! Thus CUE = NUE * avN -> 1.64 / 42.0
    ! TLS: 27/10/2021 restricted again based now on 95 %CI (12.61 / 29.68) from TRY
    PI%parmin(11) = 10d0 !5d0
    PI%parmax(11) = 100d0 !42d0 !50d0

    ! GSI min T (K) [225, 330] 
    PI%parmin(12) = 235d0 
    PI%parmax(12) = 330d0!290d0 

    ! GSI max T (K) [225, 330] 
    PI%parmin(13) = 273.15d0 
    PI%parmax(13) = 330d0!300d0 

    ! GSI min photoperiod (sec) [3600, 36000]
    PI%parmin(14) = 3600d0*3d0  !  3 hours
    PI%parmax(14) = 3600d0*21d0 ! 21 hours

    ! Leaf carbon per Area [20, 60]
    PI%parmin(15) = 20d0 
    PI%parmax(15) = 60d0 

    ! GSI max photoperiod (sec) [3600, 64800]
    PI%parmin(20) = 3600d0*3d0  !  3 hours
    PI%parmax(20) = 3600d0*21d0 ! 21 hours

    ! GSI min VPD (Pa) [1, 5500] 
    PI%parmin(21) = 10d0 
    PI%parmax(21) = 5500d0 

    ! GSI max VPD (Pa) [1, 5500]
    PI%parmin(22) = 10d0 
    PI%parmax(22) = 5500d0 

    ! BUCKET - coarse root biomass (i.e. gbio/m2 not gC/m2) needed to reach 50 %
    ! of max depth
    PI%parmin(24) = 50d0
    PI%parmax(24) = 250d0 !500d0
    ! BUCKET - maximum rooting depth
    PI%parmin(25) = 0.35d0
    PI%parmax(25) = 1d0

    ! Initial canopy GSI value 
    PI%parmin(26) = 0d0
    PI%parmax(26) = 1d0 

    ! Minimum amount of DM in (above ground) labile and foliage for grazing to occur. 
    ! This value is also the minimum amount of (above ground) labile and foliage 
    ! which must remain after grazing.
    !! Note kg.DM.ha-1 converted to gC/m2 equivalent assuming 47.5 % C content
    PI%parmin(27) = 500d0*0.0475d0
    PI%parmax(27) = 1500d0*0.0475d0 

    ! Minimum amount of DM in (above ground) labile and foliage for cutting to occur. 
    !! Note kg.DM.ha-1 converted to gC/m2 equivalent assuming 47.5 % C content
    PI%parmin(28) = 1500d0*0.0475d0 
    PI%parmax(28) = 3000d0*0.0475d0 

    ! leaf:stem allocation [0.05, 0.75]
    ! NOT CURRENTLY IN USE...
    PI%parmin(29) = 0.25d0 
    PI%parmax(29) = 0.75d0 

    ! GPP return on new Cfol investment (gCperGPP per gCnewfol)
    PI%parmin(30) = 0.001d0
    PI%parmax(30) = 0.05d0

    ! livestock demand in DM (1-3% of animal weight) 
    ! NOT CURRENTLY IN USE...
    PI%parmin(31) = 0.015d0 
    PI%parmax(31) = 0.035d0 

    ! Post-grazing labile loss (fraction)
    PI%parmin(32) = 0.01d0 
    PI%parmax(32) = 0.75d0 

    ! Post-cutting labile loss (fraction)
    PI%parmin(33) = 0.01d0 
    PI%parmax(33) = 0.75d0 

    ! Minimum amount of DM which must be removed for grazing instance to occur (g.C.m-2.d-1)
    PI%parmin(34) = 0.1d0!/7d0 
    PI%parmax(34) = 1.0d0!/7d0

    !
    ! INITIAL VALUES DECLARED HERE
    !

    ! initial labile pool size, gC/m2 [1, 500]
    PI%parmin(16) = 1d0 
    PI%parmax(16) = 500d0 

    ! initial foliar pool size, gC/m2 [1, 500]
    PI%parmin(17) = 1d0 
    PI%parmax(17) = 500d0 

    ! initial root pool size, gC/m2 [1, 500]
    PI%parmin(18) = 1d0 
    PI%parmax(18) = 500d0 

    ! initial litter pool size, gC/m2 [1, 500]
    PI%parmin(19) = 1d0 
    PI%parmax(19) = 500d0 

    ! initial SOM pool size, gC/m2 [5000, 10000] (UK) 19000, 21000
    PI%parmin(23) = 200d0
    PI%parmax(23) = 250000d0 !90000d0

    ! Initial soil water fraction
    PI%parmin(35) = 0.05d0
    PI%parmax(35) = 1.00d0

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
