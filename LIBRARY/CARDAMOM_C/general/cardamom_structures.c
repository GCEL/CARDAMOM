#pragma once

// define the structure that will contain the data in the input file
typedef struct{

// drivers
double *MET; // a vector of meteorological drivers
double meantemp, meanrad; // mean conditions used in some EDCs

// observations
double *GPP;                  // GPP (gC.m-2.d-1)
double *NEE;                  // NEE (gC.m-2.d-1)
double *LAI;                  // LAI (m2.m-2)
//double *WOO;                  // Wood increment observations (gC.m-2.y-1)
double *Reco;                 // Ecosystem respiration (gC.m-2.y-1)
double *Cfol_stock;           // Time specific estimate of foliage carbon (gC.m-2)
double *Cwood_stock;          // Time specific estimate of wood carbon (gC.m-2)
double *Croots_stock;         // Time specific estimate of roots carbon (gC.m-2)
double *Csom_stock;           // Time specific estimate of SOM carbon (gC.m-2)
double *Cagb_stock;           // Time specific estimate of AGB carbon (gC.m-2)
double *Cabgb_stock;           // Time specific estimate of ABGB carbon (gC.m-2)
double *Clit_stock;           // Time specific estimate of litter carbon (gC.m-2)
double *Cstem_stock;          // Time specific estimate of stem carbon (gC.m-2)
double *Cbranch_stock;        // Time specific estimate of branch carbon (gC.m-2)
double *Ccoarseroot_stock;    // Time specific estimate of coarse root carbon (gC.m-2)
double *Cfolmax_stock;        // maximum annual foliar stock  (gC.m-2)

double *Evap;                 // evaporation (mm d-1)

// observation uncertainties
double *GPP_unc;                  // (gC.m-2.d-1)
double *NEE_unc;                  // (gC.m-2.d-1)
double *LAI_unc;                  // log(m2.m-2)
//double *WOO_unc;                  // log(gC.m-2.y-1)
double *Reco_unc;                 // (gC.m-2.y-1)
double *Cfol_stock_unc;           // (%)
double *Cwood_stock_unc;          // (%)
double *Croots_stock_unc;         // (%)
double *Csom_stock_unc;           // (%)
double *Cagb_stock_unc;           // (%)
double *Cabgb_stock_unc;           // (%)
double *Clit_stock_unc;           // (%)
double *Cstem_stock_unc;          // (%)
double *Cbranch_stock_unc;        // (%)
double *Ccoarseroot_stock_unc;    // (%)
double *Cfolmax_stock_unc;        // (%)

double *Evap_unc;                 // mm d-1

// counters for the number of observations per data stream
int ngpp, nnee, nlai, nwoo, nreco; // C fluxes
int nCfol_stock, nCwood_stock, nCroots_stock, nCsom_stock,  nClit_stock; //stocks #1
int nCagb_stock, nCstem_stock, nCbranch_stock, nCcoarseroot_stock, nCfolmax_stock; // stocks #2
int nCabgb_stock;
int nevap;

// location of observations in the data stream
int *gpppts;                // vector of length ngpp that contains indices of valid observation days
int *neepts;                // same for nee
int *laipts;                // same for lai
//int *woopts;                // same for woody increment
int *recopts;               // same for ecosystem respiration
int *Cfol_stockpts;         // same for Cfoliage
int* Cwood_stockpts;        // same for Cwood
int *Croots_stockpts;       // same for Croots  
int *Csom_stockpts;         // same for Csom
int *Cagb_stockpts;         // same for agb C
int *Cabgb_stockpts;         // same for abgb C
int *Clit_stockpts;         // same for Clitter
int *Cstem_stockpts;        // same for Cstem
int *Cbranch_stockpts;      // same for Cbranch
int *Ccoarseroot_stockpts;  // same for Ccoarseroot
int *Cfolmax_stockpts;      // same for Cfolmax
int *Evap_pts;               // same for evaporation


/*saving computational speed by allocating memory to model output*/
double *M_GPP;
double *M_NEE;
double *M_LAI;
double *M_FLUXES;
double *M_POOLS;
double *C_POOLS;
double *M_ABGB_ANN; //annual pools

//static data
int nodays;
double deltat;
double LAT;
int ID, noobs, nomet, EDC;
int PFT, yield, age, nopars_dummy;

//binary file mcmc options (need to add all options HERE except inout files
int edc_random_search;

/*priors*/
double parpriors[50];
double parpriorunc[50];
double otherpriors[50];
double otherpriorunc[50];

// number of fluxes, pools, parameters
int nofluxes, nopools, nopars;
}DATA;


//define the structure that contains the MCMC options
typedef struct{
//Adapt step size every N iterations
int nADAPT;
//Adapt step size for XX fraction of full run
double fADAPT;
//Number of requested parameter vectors
int nOUT;
//Print info every N solutions (set to 0 for silent)
int nPRINT;
//Write every N solutions
int nWRITE;
//Append previous file with same options (1-0)
int APPEND;
//output file name
char parfile[200];
//saved step filename
char stepfile[200];
//return best fit
int returnpars;
//random ini pars
int randparini;
//int fixedpars (1-0)= fix iSnitial values that are not equal to -9999*/
int fixedpars;
}MCMC_OPTIONS;


// define the structure that contains the MCMC output
typedef struct{
/*1 or 0*/
int complete;
/*best parameters here*/
double *best_pars;
}MCMC_OUTPUT;



//define the structure that contains parameters info
typedef struct{
//number of parameters
int npars;
//maximum parameter values
double *parmax;
//minimum parameter values
double *parmin;
//initial parameter values
double *parini;
//parfix
double *parfix;
//initial step size
double *stepsize;
}PARAMETER_INFO;

typedef struct{
    int ACC;
    int ITER;
    int ACCLOC;
    double ACCRATE;
}COUNTERS;
