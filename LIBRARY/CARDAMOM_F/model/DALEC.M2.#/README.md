# Description of the DALEC.M2.#.f90

Using DALEC_GRASS V1

DALEC_GRASS =>
short code: “DALEC.16.”
long name:  “DALEC.M2.#”

DALEC_GRASS_BUCKET  => (ACM_2 not integrated yet)
short code: “DALEC.17.”
long name: “DALEC.A1.H2.M2.#"

Section 1: DALEC_GRASS naming system in CARDAMOM
+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CARDAMOM/LIBRARY/CARDAMOM_F/general/cardamom_io.f90
@ line 146 - 153: 
> npools = 6 (pool 4 is a placeholder, ph)
> nopars = 34
> nfluxes = 21 (11 & 17 ph)
@ lines 150 - 153 (no idea yet)

+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CARDAMOM/R_functions/binary_data.r
@lines 57 - 60: Changes no need
@line 127: Not needed for DALEC_GRASS
@line 128: Replace with 'cutting/grazing'
@line 132: 21-day avg VPD (Pa) 
@lines 134 - 149: Not needed for DALEC_GRASS
@155-203: Currently DALEC_GRASS is calibrated by LAI without unc, validated against GPP, the inverted model unc is calculated using stats.bayes_mvs
@766 add a ‘else if’ section of pars priors for DALEC_GRASS @153 - 175:
‘’’
    ! load some values
    gpppars(4) = 2.0  ! g N leaf_m-2
    gpppars(7) = lat
    gpppars(9) = -2.0 ! leafWP-soilWP
    gpppars(10) = 1.0 ! totaly hydraulic resistance
    gpppars(11) = pi

    ! assign acm parameters
    constants(1)=pars(10) 
    constants(2)=0.0156935
    constants(3)=4.22273
    constants(4)=208.868
    constants(5)=0.0453194
    constants(6)=0.37836
    constants(7)=7.19298
    constants(8)=0.011136
    constants(9)=2.1001
    constants(10)=0.789798

    ! post-removal residues and root death | 0:none 1:all
    foliage_frac_res  = 0.05  ! fraction of removed foliage that goes to litter
    labile_frac_res   = 0.05  ! fraction of removed labile that goes to litter
    roots_frac_death  = 0.01  ! fraction of roots that dies and goes to litter


‘’’
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CARDAMOM/R_functions/cardamom_model_details.r
@256 add a ‘else if’ section for DALEC_GRASS:
else if (modelname == “DALEC.M2.#” | modelname ==  “DALEC.16.”) {
    # information contains is
    # The model name
    # Number of met parameters
    # Number of model parameters to be optimised
    nopools=array(6,dim=c(length(ctessel_pft)))
    nopars=array(34,dim=c(length(ctessel_pft)))
    nofluxes=array(21,dim=c(length(ctessel_pft)))
    if (specific_pft == "pft_specific") {nopars[which(ctessel_pft == 1)]=35 ; nofluxes[which(ctessel_pft == 1)]=22 ; nopools[which(ctessel_pft == 1)]=7}
    cardamom_model_details=list(name=“DALEC.M2.#”,shortname = "DALEC.16.",nopools=nopools,nofluxes=nofluxes,nomet=14,nopars=nopars)
  } 

Section 2: Simulate and post-process R code

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CARDAMOM/R_functions/simulate_all.r
DALEC_Grass outputs: lai, gpp, nee, pools, fluxes, rem
@line 28: output_dim might not be 17
@52: output_dim might not be 11, acm() is integrated in DALEC_Grass with 21 outputs, i.e., fluxes @251 in DALEC_GRASS.f90
@2831: add a ‘else if’ section for DALEC_GRASS: 
 output_dim = TBC; MTT_dim = TBC ; SS_dim = TBC
reassign DALEC_Grass outputs to ‘states_all’ with different key names to the standalone DALEC_Grass
What is ‘aNPP_dim’, ‘MTT_dim’, and ‘SS_dim’.

++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CARDAMOM/R_functions/uncertainty_figures.r
@line 27: might need daily for DALEC_Grass as if (PROJECT$model$timestep == "daily") {timestep = mean(PROJECT$model$timestep_days)}
@107: add a ‘else if’ section for DALEC_GRASS:
else if (PROJECT$model$name == "DALEC.M2.#"){
       # These models assume rooting depth is controlled by coarse root, which is a fraction of the woody pool!
       tmp = t(states_all$wood_gCm2)*as.vector(parameters[25,,])
       var = t(states_all$roots_gCm2) + tmp
       # Now estimate the rooting depth based on the equation imbedded in DALEC.A1.C2.D2.F2.H2.P3.R1.
       var = as.vector(parameters[27,,]) * (var*2) / (as.vector(parameters[26,,]) + (var*2))
       plot_root_depth = TRUE
   } 
@122 & 148 no snow in DALEC_Grass
@740 no Wood stocks in DALEC_Grass
@830 & 854 seems no Canopy growth inde in DALEC_Grass
@878 seem no Canopy mortality index in DALEC_Grass
@1023 & 1047 replace harvest with cut/grazing
@1072 DALEC_Grass has no fire pool
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CARDAMOM/R_functions/run_mcmc_results.r
@ lines 85 - 86: DALEC_Grass might not need wood pool 
@ lines 95 - 100: might need to change
@ 113: Need to consider if to remove Wood pool for DALEC_Grass
@ 126, 143, 158: Replace harvest with cut & grazing
@ 127, 144, 159: DALEC_Grass has no fire
@ 188: remove fire
@ 189: harvest -> cut/grazing
@ 190 - 209 -> remove fire-related variables
@ 236: remove fire
@ 237: harvest -> cut/grazing
@ 238 - 256-> remove fire-related variables
@ 284: remove fire
@ 285: harvest -> cut/grazing
@ 286-309 -> remove fire-related variables
@ 338: remove fire
@ 339: harvest -> cut/grazing
@ 341 – 366-> remove fire-related variables
@ 368 DALEC_Grass has no wood
@ 442: remove fire
@ 443: harvest -> cut/grazing
@ 444 - 459-> remove fire-related variables
@ 462 DALEC_Grass has no wood
@ 526: remove fire
@ 527: harvest -> cut/grazing
@ 529 - 538-> remove fire-related variables
@ 560: remove fire
@ 561: harvest -> cut/grazing
@ 563 - 572-> remove fire-related variables
@ 576 DALEC_Grass doesn’t have ET_kgH2Om2day yet, need ACM two
@ 641 no snow component
@ 660 - 677 DALEC_Grass’s ACM might have different stomatal conductance expression - seems no boundary layer gs 
@ 723 might not need MTT_wood_years_to_NPP_wood_gCm2day_correlation for DALEC_Grass
@ 917 fire_gCm2day in if condition, can be not changed
@ 924 harvest_gCm2day note here but no need to act
@ 947 how to avoid error flags when no fire?
@ 954 - 973 -> no wood
@ 1086 - 1113 Wood pool
@ 1114 - 1169 ET
@ 1198 - 1223 fire
@ 1231 wood, but might not need to change anything
@ 1235 wood, but might not need to change anything
@ 1237 wood, but might not need to change anything
@ 1242 wood, but might not need to change anything
@ 1244 wood, but might not need to change anything
@ 1307 fire, but might not need to change anything
@ 1310 fire, but might not need to change anything
@ 1313 harvest, but might not need to change anythin
@ 1315 add cut/grazing
@ 1320 - 1329 fire and harvest, but might not need to change anything
@ 1335 - 1343 fire and harvest, but might not need to change anything
@ 1346 - 1362 wood, fire, and harvest, but might not need to change anything
@ 1373 - 1389 fire and harvest, but might not need to change anything
@ 1393 - 1398 word, fire, and harvest, but might not need to change anything
@ 1442 - 1466 assimilated_wood_mean_gCm2 and assimilated_wood_mean_unc_gCm2 
@ 1506 NPP_wood_fraction, no action is needed
@ 1511 MTT_wood_years, no action is needed
@ 1513 MTT_woodlitter_years, no action is needed
@ 1519 SS_wood_gCm2, no action is needed
@ 1521 SS_woodlitter_gCm2, no action is needed
@ 1556 harvest_gCm2day -> cut/grazing
@ 1557 mean_harvest_gCm2day -> cut/grazing
@ 1558 mean_annual_harvest_gCm2day -> cut/grazing
@ 1559 - 1561 fire, no action is needed, leave it as nan
@ 1605 - 1636, leave as nan and add for cut/grazing
@ 1609 -  1718 fire and harvest, leave as nan and add for cut/grazing
@ 1784 - 1813 fire and harvest, leave as nan and add for cut/grazing
@ 1823 fire, no action is needed, leave it as nan
@ 1825 harvest, no action is needed, leave it as nan, but add something for cut/grazing
@ 1839 - 1922, wood, no action is needed, leave it as nan
@ 1947 - 1969 fire and harvest, leave as nan and add for cut/grazing
@ 1979 fire, no action is needed, leave it as nan
@ 1981  harvest, no action is needed, leave it as nan, but add something for cut/grazing
@ 1994 - 2052 wood litter, no action is needed, leave it as nan
@ 2081 - 2096 fire and harvest, leave as nan and add for manure into soil pool
@ 2105 - 2107 fire and harvest, no action is needed, leave it as nan, but add something for cut/grazing
Add harvest for the following lines in run_mcmc_results.r
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CARDAMOM/R_functions/post_process_into_grid.r (main_dev_branch)
@ line 50 and 5: may not need ‘assimilated_wood_mean_gCm2’ and ‘assimilated_wood_mean_unc_gCm2‘’, but can leave them as NaN
@ 59: ‘NPP_wood_fraction’ is 0
@ 64: ‘MTT_wood_years’ is none for DALEC_Grass
@ 66: ‘MTT_woodlitter_years’ is none
@ 72: ‘SS_wood_gCm2’ is NaN
@ 74: ‘SS_woodlitter_gCm2’ is NaN
@ 91: ‘mean_harvest_gCm2day’ is none for DALEC_Grass but add for cut and grazing
@ 92: ‘mean_fire_gCm2day’ is none
@ 108: ‘harvest_gCm2day’ is none, add for cut and grazing
@ 109: ‘fire_gCm2day’ is none
@ 122: ‘mean_annual_harvest_gCm2day’ is none, and add for cut and grazing
@ 123: ‘mean_annual_fire_gCm2day’ is none
@ 153: ‘FireFractionOfTurnover_biomass’ is none
@ 154: ‘HarvestFractionOfTurnover_biomass’ us none, but DALEC_Grass needs cut & grazing turnover fraction
@ 155 - 172: not for DALEC_Grass but don’t need to change anything (considering manure). And repeat this process/check for @ 176 - 677
@ 684: ‘fire_parameter_correlation’ is none but add for cut & grazing 
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CARDAMOM/R_functions/post_process_dalec.r (main_dev_branch)
@ line 33: ‘rhet_woodlitter_gCm2day’, DALEC_Grass doesn’t have wood litter pool, but no action is required.
@ 48: no ‘fire_gCm2day’, but no action is required
@ 55: no ‘harvest_gCm2day’, but add for cut & grazing
@ 78 - 104: alert, no fire and wood for DALEC_Grass, but considering to add correlation for cut and grazing separately with all other parameters.
@ 232 - 259: no wood but no need to change anything
@ 344 - 369: no fire but no need to change anything
Cannot add for grazing and cut (or manure) because no observation information exists for them.

Section 3: Configure the models
Use DALEC.A1.C1.D2.F2.H1.P1.# as an example below.
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CARDAMOM/LIBRARY/CARDAMOM_F/model/DALEC.M2.#/src/DALEC.M2.#.f90
Nothing need to add for ‘DALEC_GRASS.f90’ except the file name at this stage
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CARDAMOM/LIBRARY/CARDAMOM_F/model/DALEC.M2.#/src/DALEC.M2.#_PARS.f90
DALEC_Grass doesn’t provide the prior/initial values for parameters but for their range in ‘MDF.py’ @ lines 45 - 78:
pars_lims = {
    0:  [1e-3,  0.1],     # Decomp rate [1e-5, 0.01]
    1:  [0.43,   0.48],    # GPP to resp fraction [~0.46]
    2:  [0.75,   1.5],     # GSI sens leaf growth [1.0, 1.025]
    3:  [0.10,   1.0],     # NPP belowground allocation exponential parameter [0.01, 1.00]
    4:  [1e-3,   2.0],     # GSI max leaf turnover [1e-5, 0.2]
    5:  [1e-3,   1e-1],    # TOR roots [0.0001, 0.01]
    6:  [1e-3,   1e-1],    # TOR litter [0.0001, 0.01]
    7:  [1e-7,   1e-4],    # TOR SOM [1e-7, 0.001]
    8:  [0.01,   0.20],    # T factor (Q10) [0.018,  0.08]
    9:  [7,      25],      # PNUE [7, 20]
    10: [1e-3,   1.0],     # GSI max labile turnover [1e-6, 0.2]
    11: [230,    290],     # GSI min T (K) [225, 330] 
    12: [250,    300],     # GSI max T (K) [225, 330] 
    13: [3600,   20000],   # GSI min photoperiod (sec) [3600, 36000]
    14: [35,     55],      # Leaf Mass per Area [20, 60]
    15: [20,     100],     # initial labile pool size [1, 1000]
    16: [20,     100],     # initial foliar pool size [1, 1000]
    17: [40,     2000],    # initial root pool size [1, 1000]
    18: [40,     2000],    # initial litter pool size [1, 10000]
    19: [10000,  40000],   # GSI max photoperiod (sec) [3600, 64800]
    20: [100,    3000],    # GSI min VPD (Pa) [1, 5500] 
    21: [1000,   5000],    # GSI max VPD (Pa) [1, 5500]
    22: [1e-3,   0.5],     # critical GPP for LAI growth [1e-10, 0.30]
    23: [0.96,   1.00],    # GSI sens for leaf senescenece [0.96, 1.00]
    24: [0.5,    3.0],     # GSI growing stage/step [0.50, 1.5]
    25: [1.0,    2.0],     # Initial GSI [1.0, 2.0]
    26: [500,    1500],    # DM min lim for grazing (kg.DM.ha-1)
    27: [1500,   3000],    # DM min lim for cutting (kg.DM.ha-1)
    28: [0.25,   0.75],    # leaf:stem allocation [0.05, 0.75]
    29: [19000,  21000],   # initial SOM pool size [5000, 10000] (UK) 19000, 21000
    30: [0.015,  0.035],   # livestock demand in DM (1-3% of animal weight) 
    31: [0.01,   0.10],    # Post-grazing labile loss (fraction)
    32: [0.50,   0.90],    # Post-cutting labile loss (fraction)
    33: [0.1,    1.0]     # min DM removal for grazing instance to occur (g.C.m-2.w-1)
}
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CARDAMOM/LIBRARY/CARDAMOM_F/model/DALEC.M2.#/src/DALEC.M2.#_R_interface.f90
@ lines 41 and 42, 73 and 74: ‘soil_frac_clay_in(nos_soil_layers) & ! clay in soil (%)’ and ‘soil_frac_sand_in(nos_soil_layers) & ! sand in soil (%)’ might not suitable for DALEC_Grass
@ 104 - 158: change the flux indexes of ‘FLUXES’ to be compatible with DALEC_Grass, leave fire and harvest as NaNs and add for cut & grazing. This action requires change of ‘output_dim’. 
@ 179 - 182: make wood as NaN
@ 195 - 224: check FLUXES and POOLS index oders to fit DALEC_Grass, note placeholders (assigned as NaN) in DALEC_Grass for Wood.
@ 229 - 242:  check FLUXES index oders
++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
CARDAMOM/LIBRARY/CARDAMOM_F/model/DALEC.M2.#/likelihood/MODEL_LIKELIHOOD.f90
Configure ‘carbon_model’ calling parameters for DALEC_Grass
@ lines between 485 and 830: Add EDC checks in MDF.py @ 129 - 146: 

if  (  (np.isnan(pools).any()) or (np.isnan(lai).any())
or (np.isnan(gpp).any()) or (np.isnan(fluxes).any())                          
or (np.any(pools < 0)) or (np.any(fluxes < 0))
or (np.any(lai < 0))
or (np.all(fluxes[:,[0,1,2,3,4,5,6,7,8,9,11,12,13,14,15,17,18,19,20]]==0,axis=0)).any()
### Fluxes                                          
or (np.any(gpp > 25))
or ((gpp*7).sum() < 500*self.noyears )
or ((gpp*7).sum() > 2800*self.noyears)
or (np.any((fluxes[:,12]+fluxes[:,13]+fluxes[:,2]) > 20) )
or (((fluxes[:,12]+fluxes[:,13]+fluxes[:,2])*7).sum() < 500*self.noyears)
or (((fluxes[:,12]+fluxes[:,13]+fluxes[:,2])*7).sum() > 2600*self.noyears)
### Soil C  
or (abs(pars[29] - pools[-1,5]) > pars[29]*0.05) ## soil C stable      
### Management
or ( (rem[0,:]*21/float(650*0.035) > 70).any() ) # max total LSU_ha_week
or ( int(abs(REMDF.cutsno.sum())) != int(len(REMDF[REMDF.Csim>0])) ) # all cuts in inputs are simulated
) : return [-np.inf]
