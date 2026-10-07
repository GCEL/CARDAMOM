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
! Module contains uniform prior parameter information for the DALEC.A4.C6.D2.F2.H3.P12.033 model.
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

    !
    ! declare parameters
    !

    type(PARINFO), intent(inout):: PI

    PI%npars = 48
    if (.not. allocated(PI%parmin)) allocate(PI%parmin(PI%npars))
    if (.not. allocated(PI%parmax)) allocate(PI%parmax(PI%npars))

    ! Decomposition efficiency of litter/CWD to som (fraction)
    PI%parmin(1) = 0.25d0
    PI%parmax(1) = 0.75d0

    ! Fraction of GPP respired as autotrophic respiration
    ! for maintenance of wood and fine roots
    PI%parmin(2) = 0.1d0
    PI%parmax(2) = 0.5d0

    ! Potential rate of labile to foliage (gC/m2/day)
    PI%parmin(3) = 1d0
    PI%parmax(3) = 20d0

    ! Potential rate of labile to fine root (gC/m2/day)
    PI%parmin(4) = 1d0
    PI%parmax(4) = 20d0

    ! Initial NCCE (gC/gCleaf/day) reference value for gradient calculations
    PI%parmin(5) = -0.10d0
    PI%parmax(5) =  0.10d0

    ! Turnover of wood (fraction / day)
    PI%parmin(6) = 0.000009d0 ! 304  years
    PI%parmax(6) = 0.001d0    ! 2.74 years

    ! Turnover of fine roots (fraction / day)
    PI%parmin(7) = 0.001368925d0 !  2 years !0.0006844627d0 ! 4 years
    PI%parmax(7) = 0.017d0       ! 60 days

    ! Turnover of litter (fraction; temperature adjusted)
    PI%parmin(8) = 0.0001141d0 ! 24   years at 0oC
    PI%parmax(8) = 0.02d0      ! 0.13 years at 0oC    

    ! Turnover of som to Rhet (fraction; temperature adjusted)
    PI%parmin(9) = 1.368925d-06   ! 2000 years at 0oC
    PI%parmax(9) = 9.126169d-05   !   30 years at 0oC !0.0001368926d0 !   20 years at 0oC
!    PI%parmin(9) = 0.0000001d0 ! 27378.0 years at 0oC
!    PI%parmax(9) = 0.001d0     !     2.7 years at 0oC

    ! Temp factor* = Q10 = 1.2-2.2
    PI%parmin(10) = 0.019d0
    PI%parmax(10) = 0.08d0

    ! Vcmax, the maximum rate of carboxylation at the canopy top
    ! umolC/m2/s
    PI%parmin(11) = 10d0
    PI%parmax(11) = 100d0

    ! Minimum leaf water potential (MPa), at which photosynthesis is suppressed
    PI%parmin(12) = -8d0
    PI%parmax(12) = -0.5d0

    ! Parameters linking the NCCE to the CMI
    ! via a Michaelis-Menten function. 
    ! This is the NCCE at which the CMI is suppressed by 0.5
    PI%parmin(13) = -0.5d0
    PI%parmax(13) = -0.0005d0
    ! Parameters linking the NCCE gradient to the CMI
    ! via a Michaelis-Menten function. This is the NCCE gradient 
    ! 50 % value
    PI%parmin(14) = -0.1d0
    PI%parmax(14) = -0.00005d0
       
    ! Foliar NCCE return for a foliar loss to progress (gC/gC/m2/d)
    ! POSSIBLY THIS PARAMETER IS NOT NEED AND CODE JUST BE CODED AS A VERY SMALL POSITIVE NUMBER?
    PI%parmin(15) = 0.005d0
    PI%parmax(15) = 0.2d0
    ! Foliar NCCE return for a foliar growth to progress (gC/gC/m2/d)   
    PI%parmin(16) = 0.005d0
    PI%parmax(16) = 0.2d0

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
    PI%parmax(26) = 1000d0 !500d0

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

    ! labile:biomass at which growth is limited by 50 %
    PI%parmin(33) = 0.0001d0 ! 0.01 %
    PI%parmax(33) = 0.01d0   ! 1 %

    ! Temperature (oC) above p36 at which foliage and fine root growth is limited by 50 %
    PI%parmin(34) = 0.1d0
    PI%parmax(34) = 5d0
    ! Temperature (oC) above p37 at which wood growth is limited by 50 %
    PI%parmin(35) = 0.1d0
    PI%parmax(35) = 5d0
    ! Temperature (oC) at which foliage and fine root growth is prevented
    PI%parmin(36) =-8d0 
    PI%parmax(36) = 8d0
    ! Temperature (oC) at which wood growth is prevented
    PI%parmin(37) = 0d0
    PI%parmax(37) = 8d0

    ! Potential growth rate of wood (gC/m2/day)
    PI%parmin(38) = 0.1d0
    PI%parmax(38) = 20d0

    ! wSWP water potential (MPa) at which wood growth is fully suppressed
    PI%parmin(39) = -5d0
    PI%parmax(39) = -0.001d0
    ! wSWP water potential (MPa) at which wood growth suppression begins
    PI%parmin(40) = -5d0
    PI%parmax(40) = -0.001d0

    ! wSWP water potential (MPa) at which leaf growth is fully suppressed
    PI%parmin(41) = -5d0
    PI%parmax(41) = -0.001d0
    ! wSWP water potential (MPa) at which leaf growth suppression begins
    PI%parmin(42) = -5d0
    PI%parmax(42) = -0.001d0

    ! Baseline leaf maintenance respiration.
    ! For details see Table S3, Heskel et al., (2016), doi: http://www.pnas.org/cgi/doi/10.1073/pnas.1520282113
    PI%parmin(43) = -4.4d0
    PI%parmax(43) = -0.6d0

    ! Potential loss rate for foliage to litter (fraction/day)
    ! based on historical NCCE (gCgCday) and current step NCCE (gCgCday)
    PI%parmin(44) = 0.0003424658 ! 8 years
    PI%parmax(44) = 0.0333333333 ! 30 days

    ! Minimum foliar loss proposal (gC/m2/day)
    ! for dNCCE calculation
    PI%parmin(45) = 0.001d0   
    PI%parmax(45) = 1d0    

    ! Intrinsic canopy water use efficiency for stomatal regulation (gC/mmolH2O-1/m2leaf/s-1)
    ! A credible iWUE range spans atleast 0.00001 -> 0.01
    PI%parmin(46) = 1d-6
    PI%parmax(46) = 1d-1

    ! Initial target allocation for autotrophic carbon
    PI%parmin(48) = 0.01d0
    PI%parmax(48) = 10d0

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

    ! C litter
    PI%parmin(22) = 1d0
    PI%parmax(22) = 2000d0

    ! C_som
    PI%parmin(23) = 200d0
    PI%parmax(23) = 250000d0 !90000d0

    ! Initial soil water fraction
    PI%parmin(24) = 0.05d0
    PI%parmax(24) = 1.00d0

    ! C autotrophic
    PI%parmin(47) = 1d0
    PI%parmax(47) = 2000d0

  end subroutine pars_info

  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
