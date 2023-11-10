module MODEL_PARAMETERS

  implicit none

  !!!!!!!!!!!
  ! Authorship contributions
  !
  ! This code is based on the original C verion of the University of Edinburgh
  ! CARDAMOM framework created by A. A. Bloom (now at the Jet Propulsion Laboratory).
  ! All code translation into Fortran, integration into the University of
  ! Edinburgh CARDAMOM code and subsequent modifications by:
  ! T. L. Smallman (t.l.smallman@ed.ac.uk, University of Edinburgh)
  ! S. Zhu (University of Edinburgh)
  ! See function / subroutine specific comments for exceptions and contributors
  !!!!!!!!!!!

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

    ! Subroutine contains a list of parameter ranges for the model.
    ! These could or possibly should go into an alternate file which can be read in.
    ! This may improve the usability when it comes to reading these information
    ! in for different PFTs

    implicit none

    !
    ! declare parameters
    !

    ! Decomposition rate [1e-5, 0.01]
    PI%parmin(1) = 0.001d0 
    PI%parmax(1) = 0.1d0 

    ! GPP to resp fraction [~0.46]
    PI%parmin(2) = 0.43d0 
    PI%parmax(2) = 0.48d0 

    ! GSI sens leaf growth [1.0, 1.025]
    PI%parmin(3) = 0.75d0 
    PI%parmax(3) = 1.5d0 

    ! NPP belowground allocation exponential parameter [0.01, 1.00]
    PI%parmin(4) = 0.1d0 
    PI%parmax(4) = 1.0d0 

    ! GSI max leaf turnover [1e-5, 0.2]
    PI%parmin(5) = 0.001d0 
    PI%parmax(5) = 2.0d0 

    ! TOR roots [0.0001, 0.01]
    PI%parmin(6) = 0.001d0 
    PI%parmax(6) = 0.1d0 

    ! TOR litter [0.0001, 0.01]
    PI%parmin(7) = 0.001d0 
    PI%parmax(7) = 0.1d0 

    ! TOR SOM [1e-7, 0.001]
    PI%parmin(8) = 0.0000001d0 
    PI%parmax(8) = 0.0001d0 

    ! T factor (Q10) [0.018,  0.08]
    PI%parmin(9) = 0.01d0 
    PI%parmax(9) = 0.2d0 

    ! PNUE [7, 20]
    PI%parmin(10) = 7d0 
    PI%parmax(10) = 25d0 

    ! GSI max labile turnover [1e-6, 0.2]
    PI%parmin(11) = 0.001d0 
    PI%parmax(11) = 1.0d0 

    ! GSI min T (K) [225, 330] 
    PI%parmin(12) = 230d0 
    PI%parmax(12) = 290d0 

    ! GSI max T (K) [225, 330] 
    PI%parmin(13) = 250d0 
    PI%parmax(13) = 300d0 

    ! GSI min photoperiod (sec) [3600, 36000]
    PI%parmin(14) = 3600d0 
    PI%parmax(14) = 20000d0 

    ! Leaf Mass per Area [20, 60]
    PI%parmin(15) = 35d0 
    PI%parmax(15) = 55d0 

    ! initial labile pool size [1, 1000]
    PI%parmin(16) = 20d0 
    PI%parmax(16) = 100d0 

    ! initial foliar pool size [1, 1000]
    PI%parmin(17) = 20d0 
    PI%parmax(17) = 100d0 

    ! initial root pool size [1, 1000]
    PI%parmin(18) = 40d0 
    PI%parmax(18) = 2000d0 

    ! initial litter pool size [1, 10000]
    PI%parmin(19) = 40d0 
    PI%parmax(19) = 2000d0 

    ! GSI max photoperiod (sec) [3600, 64800]
    PI%parmin(20) = 10000d0 
    PI%parmax(20) = 40000d0 

    ! GSI min VPD (Pa) [1, 5500] 
    PI%parmin(21) = 100d0 
    PI%parmax(21) = 3000d0 

    ! GSI max VPD (Pa) [1, 5500]
    PI%parmin(22) = 1000d0 
    PI%parmax(22) = 5000d0 

    ! critical GPP for LAI growth [1e-10, 0.30]
    PI%parmin(23) = 0.001d0 
    PI%parmax(23) = 0.5d0 

    ! GSI sens for leaf senescenece [0.96, 1.00]
    PI%parmin(24) = 0.96d0 
    PI%parmax(24) = 1.0d0 

    ! GSI growing stage/step [0.50, 1.5]
    PI%parmin(25) = 0.5d0 
    PI%parmax(25) = 3.0d0 

    ! Initial GSI [1.0, 2.0]
    PI%parmin(26) = 1.0d0 
    PI%parmax(26) = 2.0d0 

    ! DM min lim for grazing (kg.DM.ha-1)
    PI%parmin(27) = 500d0 
    PI%parmax(27) = 1500d0 

    ! DM min lim for cutting (kg.DM.ha-1)
    PI%parmin(28) = 1500d0 
    PI%parmax(28) = 3000d0 

    ! leaf:stem allocation [0.05, 0.75]
    PI%parmin(29) = 0.25d0 
    PI%parmax(29) = 0.75d0 

    ! initial SOM pool size [5000, 10000] (UK) 19000, 21000
    PI%parmin(30) = 19000d0 
    PI%parmax(30) = 21000d0 

    ! livestock demand in DM (1-3% of animal weight) 
    PI%parmin(31) = 0.015d0 
    PI%parmax(31) = 0.035d0 

    ! Post-grazing labile loss (fraction)
    PI%parmin(32) = 0.01d0 
    PI%parmax(32) = 0.1d0 

    ! Post-cutting labile loss (fraction)
    PI%parmin(33) = 0.5d0 
    PI%parmax(33) = 0.9d0 

    ! min DM removal for grazing instance to occur (g.C.m-2.w-1)
    PI%parmin(34) = 0.1d0 
    PI%parmax(34) = 1.0d0 

  end subroutine pars_info

  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
