#include <math.h>
#include "EDCs.c"
#include "../src/DALEC_LU_FIRES_GSI_newalloc.c" // model
#include "../src/DALEC_LU_FIRES_GSI_newalloc_PARS.c" // parameter info
#include "../../../general/cardamom_structures.c"
#include "../../../method/MHMCMC/MHMCMC.c"


/*********************************************************
 |    Calculate the likelihood according to observations |
 *********************************************************/
double likelihood(DATA D)
{
    int n,dn;
    double P=0;
    double tot_exp;
    //double CPG,PP;


    /*GPP likelihood*/
    tot_exp=0;
    if (D.ngpp>0){for (n=0;n<D.ngpp;n++){dn=D.gpppts[n];tot_exp+=pow((D.M_FLUXES[dn*D.nofluxes]-D.GPP[dn])/2,2);}
    P=P-0.5*tot_exp;}
  
    
    /*LAI log likelihood*/
    tot_exp=0;
    if (D.nlai>0){for (n=0;n<D.nlai;n++){dn=D.laipts[n];tot_exp+=pow(log(D.M_LAI[dn]/D.LAI[dn])/log(2),2);}
    P=P-0.5*tot_exp;}
  
    
    /*NEE likelihood*/
    tot_exp=0;
    if (D.nnee>0){for (n=0;n<D.nnee;n++){dn=D.neepts[n];tot_exp+=pow((D.M_NEE[dn]-D.NEE[dn])/2,2);}
    P=P-0.5*tot_exp;}

    // Reco likelihood - Reco is the sum of of Ra and Rh from both litter and soil, fluxes ID 2, 12 and 13
    tot_exp=0;
    if (D.nreco>0){for (n=0;n<D.nreco;n++){dn=D.recopts[n];tot_exp+=pow(((D.M_FLUXES[(dn*D.nofluxes)+2]+D.M_FLUXES[(dn*D.nofluxes)+12]+D.M_FLUXES[(dn*D.nofluxes)+13])-D.Reco[dn])/2,2);}
    P=P-0.5*tot_exp;}

    // WOO likelihood
    tot_exp=0;
    if (D.nwoo>0){for (n=0;n<D.nwoo;n++){dn=D.woopts[n];tot_exp+=pow((D.M_POOLS[(dn*D.nopools)+3]-D.WOO[dn])/2,2);}    
    P=P-0.5*tot_exp;}

    //printf("P: %f \n", P);
    /* log-likelihood*/
    if (isnan(P)){P=log(0);}
    
    return(P);
}


/**********************************************
    Likelihood according to prior knowledge
***********************************************/
double likelihood_p(DATA D, double PARS[])
{
//*remember - LOG likelihood
double p=0;
int n;

// looping through all priors for P
for (n=0;n<D.nopars;n++)
{   
    if (D.parpriors[n]>-9998)
    {
        p=p-0.5*pow(log(PARS[n]/D.parpriors[n])/log(D.parpriorunc[n]),2);
    }
}

// total biomass (pools 1:4) defined here - using second space of parpriors for total biomass
if (D.otherpriors[1]>-9998)
{
    p=p-0.5*pow(log((PARS[17]+PARS[18]+PARS[19]+PARS[20])/D.otherpriors[1])/log(D.otherpriorunc[1]),2);
}

// for any other priors, explicitly define functions based on values in DATA.otherpriors


return p;}


/***************************************************************
    Calculate the EDC model likelihood to find starting point
****************************************************************/
double model_likelihood(DATA D, PARAMETER_INFO PI, double *PARS)
{
    struct EDCDIAGNOSTIC EDCD;
    EDCD.DIAG=0;
    int EDC,n;
    double P=0,P_p;

//    for (n=0; n<PI.npars; n++){
  //  printf("param value passed EDCs0  %f %f\n", PI.parini[0], PARS[0]);
    // check first set of EDCs: parameter values
    EDC=pow(EDC1_LU_FIRES_GSI(PARS, &EDCD, D.meantemp, D.meanrad), D.EDC);
//    for (n=0; n<PI.npars; n++){
    //printf("param value passed EDCs1  %f %f\n", PI.parini[0], PARS[0]);
    P=P+log((double)EDC);
  //  printf("hi from model_likelihood1: %d\n", EDC);
    //if parameter values are OK 
    if (EDC==1)
    {
        // compute parameter of prior values
        P=P+likelihood_p(D,PARS);
       // printf("like_p from mod_likelihood %f\n", likelihood_p(D,PARS));
        P_p=P;
//        for (n=0; n<PI.npars; n++){
    //    printf("param value passed before EDCs2a %f %f\n", PI.parini[0], PARS[0]);
        // run model
        CARBON_MODEL(D.MET, PARS, D.deltat, D.nodays, D.LAT, D.M_LAI, D.M_NEE, D.M_FLUXES, D.M_POOLS);
       // printf("LAI: %f \n", D.M_LAI[0]);
        // storing GPP
//        for (n=0; n<PI.npars; n++){
    //    printf("param value passed before EDCs2b %f %f\n", PI.parini[0], PARS[0]);
        for (n=0;n<D.nodays;n++){D.M_GPP[n]=D.M_FLUXES[n*D.nofluxes];}
//        for (n=0; n<PI.npars; n++){
       // printf("param value passed before EDCs2c %f %f\n", PI.parini[0], PARS[0]);

        // EDC2 check
        EDC=EDC2_LU_FIRES_GSI(PARS, D.MET, D.M_LAI, D.M_NEE, D.M_POOLS, D.M_FLUXES, PI, D, &EDCD, D.meantemp);
        EDC=pow(EDC,D.EDC);
       // printf("hi from model_likelihood2: %d\n", EDC);
        // printf("EDC2: %i; ",EDC);

        // LIKELIHOOD: -inf if EDC=0
        P=P+log((double)EDC);
        if (EDC==1){P=P+likelihood(D);}
    }

    // Returning the log likelihood P
    //printf("hi from model_likelihood: %f\n", P);
  //  if (EDC == 1) {printf("param value passed after EDCs2  %f %f\n", PI.parini[0], PARS[0]);}
    return P;
}

/***************************************************************
    Calculate the EDC model likelihood to find starting point
****************************************************************/
double edc_model_likelihood(DATA D, PARAMETER_INFO PI, double *PARS)
{
    struct EDCDIAGNOSTIC EDCD;
    EDCD.DIAG=1;
    int EDC, n;
    double P=0;

    //check first set of EDCs that will return detailed EDCD.PASSFAIL
    EDC=EDC1_LU_FIRES_GSI(PARS, &EDCD, D.meantemp, D.meanrad);

    // running model
    CARBON_MODEL(D.MET, PARS, D.deltat, D.nodays, D.LAT, D.M_LAI, D.M_NEE, D.M_FLUXES, D.M_POOLS);

    // EDC2 check and return detailed EDCD.PASSFAIL
    EDC=EDC*EDC2_LU_FIRES_GSI(PARS, D.MET, D.M_LAI, D.M_NEE, D.M_POOLS, D.M_FLUXES, PI, D, &EDCD, D.meantemp);

    // LIKELIHOOD
    int tot_exp=0;
    
    for (n=0;n<EDCD.nedc;n++)
    {
        tot_exp+=1-EDCD.PASSFAIL[n];    
       // if (EDCD.PASSFAIL[n] == 0) {printf("%d ", n);}
    }
   // printf("\n");
//    for (n=0; n<16;n++) {printf("%d-%d ",n,EDCD.PASSFAIL[n]);}

    P=-0.5*((double)tot_exp*10)*D.EDC;

    // overriding if model likelihood is -inf
    double ML=model_likelihood(D,PI,PARS);
    if (D.EDC==0 && (isinf(ML)==-1 || isnan(ML))){P=P-0.5*10;}
    

    // adding an exponential decay related term to find starting point quicker!
    double PC=0,C;
    double co=-log(2)/(D.nodays*D.deltat);
    for (n=0;n<D.nopools;n++)
    {
      //  printf("%d %f\n",n, D.M_POOLS[n]);
        
        C=expdecay2(D.M_POOLS,n,D.nodays+1,D.deltat,D.nopools);//printf("%f %f ", D.deltat, co);
        if (C<co && finite(C)){PC=PC-0.5*pow((C-co)/(co*10),2);}
    }

    P=P+PC;

    //printf("LAI: %f\n",D.M_LAI[0]);
  //  printf("hi from edc_model_likelihood: %f %f\n", P, ML);
    // Log likelihood
    return P;

}

void find_edc_initial_values(PARAMETER_INFO *PI, DATA D)
{
    MCMC_OPTIONS MCOPT_init;
    MCMC_OUTPUT MCOUT_init;

    // set values in MCOUT
    MCOUT_init.best_pars = calloc(PI->npars,sizeof(double));
    
    MCOUT_init.complete = 0;

    MCOPT_init.APPEND=0;
    MCOPT_init.nADAPT=20;
    MCOPT_init.fADAPT=0.5;
    MCOPT_init.nOUT=2000;
    MCOPT_init.nPRINT=0;
    MCOPT_init.nWRITE=0;

    // force MCMC to return the best parameter values as starting point
    MCOPT_init.randparini=1;
    MCOPT_init.returnpars=1;
    MCOPT_init.fixedpars=0;

    // increase step size and set prior value
    int parid = 0;
    for (parid = 0; parid < PI->npars; parid++) 
    {
        //PI->stepsize[parid] = 0.02;
        PI->parini[parid] = D.parpriors[parid];
        PI->parfix[parid] = 0; 
        // if the prior is not missing and we have not told the edc to be random keep the value
        if (PI->parini[parid] != -9999 && D.edc_random_search <1) {PI->parfix[parid]=1;}
    }

    double PEDC=-1; // the edc_likelihood
    int counter_local=0;

    while (PEDC < 0)
    {
        printf("Beginning search of starting point\n");
        for (parid = 0; parid < PI->npars; parid++) {PI->stepsize[parid]=0.0005;}
        // call the MHMCMC directing to the appropriate likelihood
        MHMCMC(edc_model_likelihood, D, *PI, MCOPT_init, &MCOUT_init);

        // store the best parameters from that loop
        for (parid = 0; parid < PI->npars; parid++) {PI->parini[parid] = MCOUT_init.best_pars[parid];}

        // call edc likelihood function to get final edc probability
        printf("Calculate PEDC of best parameters: ");
        PEDC = edc_model_likelihood(D, *PI, MCOUT_init.best_pars);
        double proba = model_likelihood(D, *PI, MCOUT_init.best_pars);
        printf("PEDC: %f; Proba: %f\n\n", PEDC,proba);

        // keep track of attempts
        counter_local=counter_local+1;
        // periodically reset the initial conditions
        if (PEDC < 0 && counter_local%3 == 0) 
        {          
            for (parid = 0; parid < PI->npars; parid++) {PI->parini[parid] = D.parpriors[parid];}
        }
    }
    
    for (parid = 0; parid < PI->npars; parid++) {PI->parini[parid]=MCOUT_init.best_pars[parid]; PI->parfix[parid]=0;}
}


























