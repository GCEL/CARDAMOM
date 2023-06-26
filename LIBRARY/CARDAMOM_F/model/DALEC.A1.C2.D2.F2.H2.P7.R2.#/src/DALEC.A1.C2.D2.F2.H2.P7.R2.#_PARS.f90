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

    ! generic model

    !
    ! Declare parameters
    !

    ! Decomposition efficiency of litter/CWD to som (fraction)
    PI%parmin(1) = 0.25d0
    PI%parmax(1) = 0.75d0

    ! Fraction of GPP respired as Rm(wood,root)
    PI%parmin(2) = 0.05d0
    PI%parmax(2) = 0.40d0

    ! Baseline canopy turnover
    PI%parmin(3) = 0.0002737851 ! 10 years
    PI%parmax(3) = 0.002737851  ! 1 year

    ! Fraction of (1-fgpp) to roots*/
    PI%parmin(4) = 0.10d0
    PI%parmax(4) = 0.80d0

    ! Max leaf turnover
    PI%parmin(5) = 0.0006844627d0 !  4 years
    PI%parmax(5) = 0.07142857d0   ! 14 days

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
    PI%parmin(12) = 0.0006844627d0 !  4 years
    PI%parmax(12) = 0.07142857d0   ! 14 days

    ! Fraction of GPP to Clab*/
    PI%parmin(13) = 0.15d0 ! 0.05d0
    PI%parmax(13) = 0.55d0 ! 0.35d0

    ! CGI minimum temperature (oC)
    PI%parmin(14) = -10d0
    PI%parmax(14) =  10d0

    ! CGI Optimum temperatue (oC)
    PI%parmin(15) = 10d0
    PI%parmax(15) = 40d0

    ! CGI temperature maximum (oC)
    PI%parmin(16) = 20d0
    PI%parmax(16) = 60d0

    ! LMA
    ! Kattge et al. 2011,
    PI%parmin(17) = 20d0
    PI%parmax(17) = 180d0

    ! CGI temperature kurtosis
    ! Larger number means more peaky
    PI%parmin(24) = 0.01d0
    PI%parmax(24) = 0.3d0

    ! NCCE per gCleaf at which
    ! CMI is suppressed by 0.5
    PI%parmin(25) = 0.01d0 !-0.1d0
    PI%parmax(25) = 0.3d0

    ! Number of days overwhich cmi_ncce is sensitive
    ! NOT IN USE
    PI%parmin(26) = 7d0
    PI%parmax(26) = 183d0

    ! Net Canopy Carbon Export due to new leaf (gC/m2/day)
    PI%parmin(27) = 0.0001d0 !0.1d0
    PI%parmax(27) = 0.5d0 !0.025d0

    ! Initial NCCE reference value for gradient calculations
    PI%parmin(28) = 0d0
    PI%parmax(28) = 0.30d0

    ! fraction of Cwood which is coarse root
    PI%parmin(29) = 0.15d0
    PI%parmax(29) = 0.50d0 ! increased based on evidence of savannah system 50 % below !0.30d0

    ! CGI wSWP at which stress total (MPa)
    PI%parmin(30) = -6d0
    PI%parmax(30) = -0.001d0

    ! CGI wSWP at which stress begins (MPa)
    PI%parmin(31) = -6d0
    PI%parmax(31) =  0.001d0

    ! Min leaf water potential (MPa)
    PI%parmin(32) = -5d0
    PI%parmax(32) = -1d0

    ! CGI min photoperiod threshold (sec)
    PI%parmin(33) = 3600d0*3d0  !  3 hours
    PI%parmax(33) = 3600d0*21d0 ! 21 hours

    ! CGI max photoperiod threshold (sec)
    PI%parmin(34) = 3600d0*3d0   !  3 hours
    PI%parmax(34) = 3600d0*21d0  ! 21 hours

    ! Heskel et al., (2016) intercept for leaf maintenance respiration temperature response
    ! i.e. maintenance respiration at 0oC
    PI%parmin(35) = -4.40d0 ! 1.639-0.01
    PI%parmax(35) = -0.85d0 ! 1.639+0.01

    ! Initial leaf lifespan (days)
    PI%parmin(36) =  60.d0
    PI%parmax(36) = 365.25d0 * 4d0

    ! Turnover rate for CWD
    PI%parmin(38) = 1.368925d-05 ! 200.00 years at 0oC
    PI%parmax(38) = 0.001d0       !   2.74 years at 0oC

    ! BUCKET - coarse root biomass (i.e. gbio/m2 not gC/m2) needed to reach 50 % of max depth
    PI%parmin(39) = 100d0
    PI%parmax(39) = 2000d0 !500d0

    ! BUCKET - maximum rooting depth
    PI%parmin(40) = 0.35d0
    PI%parmax(40) = 10d0

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
    PI%parmin(42) =  1.6d0
    PI%parmax(42) = 40.0d0

    ! Resilience factor for burned but not combusted C stocks
    PI%parmin(43) = 0.1d0
    PI%parmax(43) = 0.9d0
    ! Combustion completeness factor for foliage
    PI%parmin(44) = 0.01d0
    PI%parmax(44) = 0.99d0
    ! Combustion completeness factor for fine root and wood
    PI%parmin(45) = 0.01d0
    PI%parmax(45) = 0.99d0
    ! Combustion completeness factor for soil
    PI%parmin(46) = 0.001d0
    PI%parmax(46) = 0.1d0
    ! Combustion completeness factor for foliage + fine root litter
    PI%parmin(47) = 0.01d0
    PI%parmax(47) = 0.99d0
    ! Combustion completeness factor for wood litter
    PI%parmin(48) = 0.01d0
    PI%parmax(48) = 0.99d0

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

    ! C wood litter
    PI%parmin(37) = 1d0
    PI%parmax(37) = 10000d0

    ! Initial soil water - fraction of field capacity
    PI%parmin(41) = 0.50d0
    PI%parmax(41) = 1.00d0

  end subroutine pars_info
  !
  !------------------------------------------------------------------
  !
end module MODEL_PARAMETERS
