#pragma once
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include "../../misc/math_functions.c"
#include "../../misc/oksofar.c"

#include "MCMC_FUNCTIONS.c" // file contains functions such as "step", "adapt_step_size", etc...

#include "../../general/cardamom_structures.c"

void MHMCMC(
double (MODEL_LIKELIHOOD)(DATA, PARAMETER_INFO, double *),
DATA DATA, PARAMETER_INFO PI, MCMC_OPTIONS MCO, MCMC_OUTPUT *MCOUT)
{

    /* ***********INPUTS************
     *
     * MODEL_LIKELIHOOD: A function wholly responsible for 
     * (a) running the model given the DATA and parameters, 
     * (b) comparing it to observations,and 
     * (c) returning  the (log) LIKELIHOOD.
     * The function will be run as MODEL_LIKELIHOOD(DATA,PI,PARS);
     * To facilitate this, ALL data can be 
     * passed to the MHMCMC function as a structure (in order to avoid 
     * repeated read/write computational time). 
     *
     * DATA: All data needed for the MODEL_LIKELIHOOD. It can include
     * drivers, observations, etc.
     *
     * PARINFO: This structure contains information on 
     * (a) pmin, pmax:	parameter ranges (compulsory)
     * (b) initpars:	parameter starting values (optional/recommended).
     * (c) npars:		number of pars
     *
     * MCO: This structure contains option values for the MCMC run.
     * These will be set to default values if empty. Options include:
     * (a) number of runs
     * (b) filename for writing file with results
     * (c) step adaptation frequency
     * (d) initial step size
     * */

    /* **************OUTPUTS*************
     * 
     * RESULTS FILE: File includes (a) results (b) likelihood and (c) final step size
     *
     * */

    /*ensuring rapid seed change for random numbers (for parallel tasks, etc.)*/
    struct timeval time;
    gettimeofday(&time, NULL);
    srand((time.tv_sec * 1e6) + (time.tv_usec));

    /*srand(time(0));*/

    /* replace previous file if MCO.APPEND == 0
    if (MCO.nWRITE>0)
    {
        if (MCO.APPEND == 0) {FILE *fileout = fopen(MCO.parfile,"wb");}
        if (MCO.APPEND == 1) {FILE *fileout = fopen(MCO.parfile,"ab");}

    }
*/
    /*DECLARING*/
    double P0;
    /*initialising P as -inf */
    double P=log(0);
    int n;

    if (MCO.nWRITE > 0) 
    {
        printf("Writing MCMC output in \"%s\" and \"%s\"\n", MCO.parfile, MCO.stepfile);
    }    
    else
    {
    
        printf("MCMC will not save output into files.\n");
    }

    COUNTERS N;
    N.ACC=0;
    N.ITER=0;
    N.ACCLOC=0;
    N.ACCRATE=0;
    /*New and default parameter vectors*/
    double *PARS, *PARS0, *PARSALL, *BESTPARS;
    PARS0=calloc(PI.npars,sizeof(double));
    PARS=calloc(PI.npars,sizeof(double));
    BESTPARS=calloc(PI.npars,sizeof(double));
    /*All accepted parameters*/
    /*This is now the last N parameter vectors
     * where N is the adaptation frequency*/
    /*PARSALL is only used for adaptation*/
    PARSALL=calloc(MCO.nADAPT*PI.npars,sizeof(double));

    /*Random starting parameters if MCO.randparini*/
    printf("nopars from mhmcmc: %d\n", DATA.nopars); 
    for (n=0;n<PI.npars;n++)
    {
        // if MCO.fixedpars == 0, parameter n has not a fixed value
        if (MCO.fixedpars!=1){PI.parfix[n]=0;}

        /*ONLY assigning randompars if (a) randparini==1 or (b) PI.parini[n]=-9999*/
        /*BUG IS HERE!!!4-9-2013*/
        if (MCO.randparini == 1 && PI.parfix[n]!=1)
        {
            /*random parameter if PI.parini = -9999*/
            PI.parini[n]=nor2par((double)random()/RAND_MAX,PI.parmin[n],PI.parmax[n]);
        }
        // printing parameter values
        printf("p%d=%2.2f ",n,log10(PI.parini[n]));
        if ( (n+1) % 3==0 || n == PI.npars-1) {printf("\n");}
    }



    oksofar("Established PI.parini - begining MHMCMC now");
    
    for (n=0;n<PI.npars;n++){PARS0[n]=PI.parini[n];}//printf("%d %f\n",n,PARS0[n]);}
    memcpy(BESTPARS,PARS0,PI.npars*sizeof(double));

    
    /*STEP 1 - RUN MODEL WITH INITIAL PARAMETERS*/
    P0=MODEL_LIKELIHOOD(DATA,PI,PI.parini);
    printf("Starting likelihood = %e\n",P0);
    
    if (isinf(P0)==-1){printf("WARNING! P0=-inf - MHMCMC may get stuck - if so, please check your initial conditions\n");return;}
    


    /*STEP 2 - BEGIN MCMC*/
    printf("MCO.nOUT: %d\n", MCO.nOUT);
    while (N.ACC < MCO.nOUT && (P < 0 || MCO.nWRITE > 0 || N.ACC == 0))
    {
        // sample a parameter set = "take a step"
	    step(PARS0,PARS,PI);

	    // Calculate new LIKELIHOOD
	    P=MODEL_LIKELIHOOD(DATA,PI,PARS);
	
        // test whether the parameter set is accepted
	    if (P-P0>log((double)random()/RAND_MAX))
        {
		    // store accepted solution and check whether it's the best one so far
    	    for (n=0;n<PI.npars;n++)
            {
		        PARSALL[N.ACCLOC*PI.npars+n]=PARS[n];
		        PARS0[n]=PARS[n];
		        if (P>P0) {BESTPARS[n]=PARS[n];}
            }
		
		    N.ACC=N.ACC+1;N.ACCLOC=N.ACCLOC+1;
		    P0=P;
		    /*writing results: parameters and log LIKELIHOOD*/
                    /*writing ALL results! Changed on Sat 30 Nov 2013*/
		    if (MCO.nWRITE>0 && (N.ACC % MCO.nWRITE)==0) {write_results(PARS,P,PI,MCO);}
        }
	
	    /*Continuing in any case*/
	    N.ITER=N.ITER+1;


	    /*Adapting Step Size*/
	    if (N.ITER % MCO.nADAPT==0)
        {
		    N.ACCRATE=(double)N.ACCLOC/(double)MCO.nADAPT;
		
		    /*Adapting step size*/
		    if (MCO.fADAPT*(double)MCO.nOUT>(double)N.ACC)
            {
    		    adapt_step_size(PARSALL,PI,N,MCO);
            }
		
		    N.ACCLOC=0;
        }

	    /*Printing Info to Screen*/
	    if (MCO.nPRINT>0 && N.ITER % MCO.nPRINT==0)
        {
		    printf("Total Accepted = %d out of %d\n",N.ACC,MCO.nOUT);
		    printf("Local Acceptance rate %5.1f\n %%",N.ACCRATE*100);
		    printf("Log LIKELIHOOD %e\n",P0);
        }

	    /*END OF WHILE LOOP*/
    }
	
   

    /*filling in MCOUT details*/
    /*best parameter combination*/

    for (n=0;n<PI.npars;n++){MCOUT->best_pars[n]=BESTPARS[n];}
//    memcpy(MCOUT->best_pars,BESTPARS,PI.npars*sizeof(double));
   
    /*MCMC completed*/
    MCOUT->complete=1;
    /*done with MCMC completion*/
    for (n=0;n<PI.npars;n++){MCOUT->best_pars[n]=BESTPARS[n];}

    free(BESTPARS);
    free(PARS);
    free(PARS0);
    free(PARSALL);

    //if (MCO.nWRITE > 0) {fclose(fileout);}
    /*END OF MHMCMC*/
}


int MEMORY_CLEANUP(DATA DATA, PARAMETER_INFO PI, MCMC_OPTIONS MCOPT, MCMC_OUTPUT MCOUT)
{

    free(PI.parmin);
    free(PI.parmax);
    free(PI.parini);
    free(PI.parfix);
    free(PI.stepsize);
    free(DATA.MET);
    free(DATA.LAI);
    free(DATA.NEE);
//    free(DATA.WOO);
    free(DATA.GPP);


    free(DATA.M_FLUXES);
    free(DATA.M_LAI);
    free(DATA.M_NEE);
    free(DATA.M_POOLS);
    free(DATA.M_GPP);
    free(DATA.M_ABGB_ANN);

    if (DATA.ngpp>0){free(DATA.gpppts);}
    if (DATA.nlai>0){free(DATA.laipts);}
    if (DATA.nnee>0){free(DATA.neepts);}
    if (DATA.nCabgb_stock>0){free(DATA.Cabgb_stockpts);}
   // if (DATA.nwoo>0){free(DATA.woopts);}


    free(MCOUT.best_pars);


    return 0;

}
