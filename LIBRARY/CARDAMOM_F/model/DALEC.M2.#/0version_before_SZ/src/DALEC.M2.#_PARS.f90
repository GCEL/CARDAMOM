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
    PI%parmax(2) = 0.8d0

    ! Background leaf turnover rate
    ! NOT IN USE
    PI%parmin(3) = 0.0002737851d0 ! 10 years
    PI%parmax(3) = 0.0009126169d0 !  3 year

    ! Fraction of (1-fgpp) to roots*/
    PI%parmin(4) = 0.1d0
    PI%parmax(4) = 0.80d0

    ! GSI max leaf turnover
    PI%parmin(5) = 0.002737851d0 ! 1 year
    PI%parmax(5) = 0.016666667d0 ! 60 days

    ! Turnover fraction of wood
    PI%parmin(6) = 0.000009d0 ! 304  years
    PI%parmax(6) = 0.001d0    ! 2.74 years

    ! Turnover fraction of roots
    ! Gill and Jackson (2000), New Phytol., 147, 13–31
    ! Fig. 6 turnover by diameter class
    PI%parmin(7) = 0.001368925d0 ! 2    years !0.0006844627d0 ! 4 years
    PI%parmax(7) = 0.02d0        ! 0.13 years

    ! Turnover of litter (fraction; temperature adjusted)
    PI%parmin(8) = 0.0001141d0 ! 24   years at 0oC
    PI%parmax(8) = 0.02d0      ! 0.13 years at 0oC

    ! Turnover of som to Rhet (fraction; temperature adjusted)
    PI%parmin(9) = 1.368925d-06   ! 2000 years at 0oC
    PI%parmax(9) = 9.126169d-05   !   30 years at 0oC !0.0001368926d0 !   20 years at 0oC
!    PI%parmin(9) = 0.0000001d0 ! 27378.0 years at 0oC
!    PI%parmax(9) = 0.001d0     !     2.7 years at 0oC

    ! Exponential coefficient for Rhet temperature response
    PI%parmin(10) = 0.019d0
    PI%parmax(10) = 0.08d0

    ! log10 avg foliar N (gN.m-2)
    ! Kattge et al., (2011) (Quantiles 25% / 75%)
    ! and Thomas et al., (2019) (Aconite canopy paper)
    PI%parmin(11) = 0.07918125d0!0d0 !-0.2218487d0 !TLS: restricted to 1.2 gN/m2leaf
    PI%parmax(11) = 0.4771213d0 ! 0.5563025d0 ! TLS: restricted to 3 gC/m2leaf

    ! Max labile turnover fraction to foliage
    PI%parmin(12) = 0.002737851d0 !  1 years
    PI%parmax(12) = 0.025d0       ! 40 days

    ! Fraction of GPP to Clab*/
    PI%parmin(13) = 0.05d0
    PI%parmax(13) = 0.35d0

    ! GSI min temperature threshold (oC)
    PI%parmin(14) = 235d0
    PI%parmax(14) = 330d0

    ! GSI max temperature threshold (oC)
    PI%parmin(15) = 273.15d0 !243d0 !235d0
    PI%parmax(15) = 330d0

    ! GSI min photoperiod threshold (sec)
    PI%parmin(16) = 3600d0*3d0  !  3 hours
    PI%parmax(16) = 3600d0*21d0 ! 21 hours

    ! LMA
    ! Kattge et al. 2011,
    PI%parmin(17) = 20d0
    PI%parmax(17) = 180d0

    ! GSI max photoperiod threshold (sec)
    PI%parmin(24) = 3600d0*3d0   !  3 hours
    PI%parmax(24) = 3600d0*21d0  ! 21 hours

    ! GSI min VPD threshold (Pa)
    PI%parmin(25) = 10d0 !100d0
    PI%parmax(25) = 5500d0

    ! GSI max VPD threshold (Pa)
    PI%parmin(26) = 10d0 !1000d0
    PI%parmax(26) = 5500d0

    ! GPP return on new Cfol investment (gCperGPP per gCnewfol)
    PI%parmin(27) = 0.001d0
    PI%parmax(27) = 0.1d0

    ! Initial GSI value
    PI%parmin(28) = 0d0
    PI%parmax(28) = 1d0

    ! fraction of Cwood which is coarse root
    PI%parmin(29) = 0.15d0
    PI%parmax(29) = 0.50d0 ! increased based on evidence of savannah system 50 % below !0.30d0

    ! GSI senstivity for leaf senescence
    PI%parmin(34) = -1d-3
    PI%parmax(34) = -1d-4

    ! Turnover rate for CWD
    PI%parmin(35) = 1.368925d-05 ! 200.00 years at 0oC
    PI%parmax(35) = 0.001d0      !   2.74 years at 0oC

    ! Optimum nitrogen use efficiency (gC/gN per m2 at optimum temperature)
    ! Derived from Vcmax reported in Wullschleger (1993), Journal of
    ! Experimental Botany, Vol 44, No. 262, pp. 907-920.
    ! ~40 gC/gN/day
    ! TRY database equivalent 2.5 % = 1.648512; 97.5 % = 19.906560
    ! Xu et al., (2017):
    ! Variations of leaf longevity in tropical moist forests predicted by a
    ! trait-driven carbon optimality model,
    ! Ecology Letters, doi: 10.1111/ele.12804, upper value of 82 gC/gN/day
    ! Thus we will compromise on the value between these but closer to the
    ! newer estimate (i.e. 30 gC/gN/day)
    PI%parmin(36) =  1.6d0
    PI%parmax(36) = 40.0d0

    ! Resilience factor for burned but not combusted C stocks
    PI%parmin(38) = 0.1d0
    PI%parmax(38) = 0.9d0
    ! Combustion completeness factor for foliage
    PI%parmin(39) = 0.01d0
    PI%parmax(39) = 0.99d0
    ! Combustion completeness factor for fine root and wood
    PI%parmin(40) = 0.01d0
    PI%parmax(40) = 0.99d0
    ! Combustion completeness factor for soil
    PI%parmin(41) = 0.001d0
    PI%parmax(41) = 0.1d0
    ! Combustion completeness factor for foliage + fine root litter
    PI%parmin(42) = 0.01d0
    PI%parmax(42) = 0.99d0
    ! Combustion completeness factor for wood litter
    PI%parmin(43) = 0.01d0
    PI%parmax(43) = 0.99d0

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
    PI%parmin(37) = 1d0
    PI%parmax(37) = 10000d0

    !
    ! Replanting pools values
    !

    ! C labile
    PI%parmin(30) = 1.0d0
    PI%parmax(30) = 100.0d0

    ! C foliar
    PI%parmin(31) = 1.0d0
    PI%parmax(31) = 100.0d0

    ! C roots
    PI%parmin(32) = 1.0d0
    PI%parmax(32) = 100.0d0

    ! C_wood derived from forestry yield curves age = 1
    PI%parmin(33) = 1.0d0
    PI%parmax(33) = 1000.0d0

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
