# FI-Hyy - 4 years
load("~/RDM/CARD_output/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/ch1_final/FI_Hyy_16_19_tempv99_satLAIallUNC02_MRTfolNA_NEE_p12_p15_CP/RESULTS_PROCESSED/FI_Hyy_2016_2019.RData")
FI_Hyy_LAIall_states_all <- states_all
FI_Hyy_LAIall_pars <- parameters
FI_Hyy_LAIall_drivers <- drivers
rm(states_all,parameters,drivers,site_ctessel_pft)



load("~/RDM/CARD_output/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/ch1_final/FI_Hyy_16_19_tempv99_satLAIsummerUNC02_MRTfol15unc_NEE_p12_p15_CP/RESULTS_PROCESSED/FI_Hyy_2016_2019.RData")
FI_Hyy_LAIsum_states_all <- states_all
FI_Hyy_LAIsum_pars <- parameters
FI_Hyy_LAIsum_drivers <- drivers
rm(states_all,parameters,drivers,site_ctessel_pft)


#SE-Htm 4 years
load("~/RDM/CARD_output/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/ch1_final/SE_Htm_16_19_tempvDEF_satLAIallUNC02_MRTfolNA_NEE_p12_p15_CP/RESULTS_PROCESSED/SE_Htm_2016_2019.RData")
SE_Htm_LAIall_states_all <- states_all
SE_Htm_LAIall_pars <- parameters
SE_Htm_LAIall_drivers <- drivers
rm(states_all,parameters,drivers,site_ctessel_pft)



load("~/RDM/CARD_output/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/ch1_final/SE_Htm_16_19_tempvDEF_satLAIsummerUNC02_MRTfol15unc_NEE_p12_p15_CP/RESULTS_PROCESSED/SE_Htm_2016_2019.RData")
SE_Htm_LAIsum_states_all <- states_all
SE_Htm_LAIsum_pars <- parameters
SE_Htm_LAIsum_drivers <- drivers
rm(states_all,parameters,drivers,site_ctessel_pft)


# FI-Var - 4 years
load("~/RDM/CARD_output/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/ch1_final/FI_Var_16_19_tempv99_satLAIallUNC02_MRTfolna_NEE_p12_p15_CP/RESULTS_PROCESSED/FI_Var_2016_2019.RData")
FI_Var_LAIall_states_all <- states_all
FI_Var_LAIall_pars <- parameters
FI_Var_LAIall_drivers <- drivers
rm(states_all,parameters,drivers,site_ctessel_pft)



load("~/RDM/CARD_output/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/ch1_final/FI_Var_16_19_tempv99_satLAIsummerUNC02_MRTfol15unc_NEE_p12_p15_CP/RESULTS_PROCESSED/FI_Var_2016_2019.Rdata")
FI_Var_LAIsum_states_all <- states_all
FI_Var_LAIsum_pars <- parameters
FI_Var_LAIsum_drivers <- drivers
rm(states_all,parameters,drivers,site_ctessel_pft)