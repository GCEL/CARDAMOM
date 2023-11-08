#pragma once
#include <string.h>
#include <stdio.h>
#include <unistd.h>
#include <stdbool.h>
#include <stdlib.h>
#include <stddef.h>
#include <sys/stat.h>

void cardamom_model_library(DATA *DATA)
{
    if (DATA->ID == 1 || DATA->ID == 3 || DATA->ID == 5 || DATA->ID == 8)
    { // DALEC - GSI or DALEC - GSI with new allocation patterns JEFF version, or using GLEAM evap stress
        DATA->nopools = 6;
        DATA->nopars = 26;
        DATA->nofluxes = 29;
    }  
    if (DATA->ID == 2 || DATA->ID == 4 || DATA->ID == 9)
    { // DALEC CDEA LU FIRE (2) and special version to use with Yi's data (4) and using GLEAM to limit GPP (9)
        DATA->nopools = 6;
        DATA->nopars = 23;
        DATA->nofluxes = 28;
    }  
    if (DATA->ID == 6)
    { // DALEC with lumped water storage
        DATA->nopools = 7;
        DATA->nopars = 26;
        DATA->nofluxes = 30;
    }
    if (DATA->ID == 7)
    { // DALEC CDEA with CWD
        DATA->nopools = 7;
        DATA->nopars = 25;
        DATA->nofluxes = 32;
    }
    if (DATA->ID == 10)
    { // DALEC CDEA with HBV soil moisture
        DATA->nopools = 7;
        DATA->nopars = 27;
        DATA->nofluxes = 30;
    }

}


/************************************************************
    This function reads the input data and parameter info
*************************************************************/



/***********************************
    Read data in the binary file
************************************/
void read_binary_data(DATA *DATA, char *infile)
{
/*
    This opens and reads the binary data files that contains
    met drivers and parameter information. Data are loaded in the
    DATA structure

    Template for DALEC MCMC files is
    - 1-100: static elements (model ID, latitude, longitude, number of data streams
    - 101-150: parameter priors
    - 151-200: uncertainty of parameter priors
    - 201-300: other priors and uncertainties
    - 300-end: temporal drivers and data
*/
    printf("Input file binary = \"%s\"\n", infile);

    double statdat[100];
//    statdat = calloc(100, sizeof(double));
    //printf("sizeof double: %d %d\n", sizeof(double),sizeof(statdat));
    FILE *fid=fopen(infile,"rb");


    fread(statdat, sizeof(double), 100, fid);

    //printf("%f\n", statdat[0]);    

    DATA->ID = (int) statdat[0];printf("Model id = %d\n", DATA->ID);
    DATA->LAT = statdat[1];
    DATA->nodays= (int) statdat[2];
    DATA->nomet= (int) statdat[3];

    DATA->noobs= (int) statdat[4];
    DATA->EDC= (int) statdat[5];
    DATA->PFT= (int) statdat[6];
    DATA->yield = (int) statdat[7];
    DATA->age = (int)statdat[8];
    //int nopars_dummy = (int)statdat[9]; // needed for next dev stage
    // allocate case specific information
    DATA->edc_random_search = (int) statdat[10];

    // read in parameter information
    fread(DATA->parpriors, sizeof(double), 50, fid);
    fread(DATA->parpriorunc, sizeof(double), 50, fid);
    fread(DATA->otherpriors, sizeof(double), 50, fid);
    fread(DATA->otherpriorunc, sizeof(double), 50, fid);

    // allocate arrays for drivers and observations
    DATA->MET=calloc(DATA->nomet*DATA->nodays,sizeof(double));

    DATA->GPP=calloc(DATA->nodays,sizeof(double));
    DATA->GPP_unc=calloc(DATA->nodays,sizeof(double));

    DATA->NEE=calloc(DATA->nodays,sizeof(double));
    DATA->NEE_unc=calloc(DATA->nodays,sizeof(double));

    DATA->LAI=calloc(DATA->nodays,sizeof(double));
    DATA->LAI_unc=calloc(DATA->nodays,sizeof(double));

    DATA->Reco=calloc(DATA->nodays,sizeof(double));
    DATA->Reco_unc=calloc(DATA->nodays,sizeof(double));

  //  DATA->WOO=calloc(DATA->nodays,sizeof(double));
  //  DATA->WOO_unc=calloc(DATA->nodays,sizeof(double));

    DATA->Cfol_stock=calloc(DATA->nodays,sizeof(double));
    DATA->Cfol_stock_unc=calloc(DATA->nodays,sizeof(double));

    DATA->Cwood_stock=calloc(DATA->nodays,sizeof(double));
    DATA->Cwood_stock_unc=calloc(DATA->nodays,sizeof(double));

    DATA->Croots_stock=calloc(DATA->nodays,sizeof(double));
    DATA->Croots_stock_unc=calloc(DATA->nodays,sizeof(double));

    DATA->Csom_stock=calloc(DATA->nodays,sizeof(double));
    DATA->Csom_stock_unc=calloc(DATA->nodays,sizeof(double));

    DATA->Cagb_stock=calloc(DATA->nodays,sizeof(double));
    DATA->Cagb_stock_unc=calloc(DATA->nodays,sizeof(double));

    DATA->Cabgb_stock=calloc(DATA->nodays,sizeof(double));
    DATA->Cabgb_stock_unc=calloc(DATA->nodays,sizeof(double));

    DATA->Clit_stock=calloc(DATA->nodays,sizeof(double));
    DATA->Clit_stock_unc=calloc(DATA->nodays,sizeof(double));

    DATA->Cstem_stock=calloc(DATA->nodays,sizeof(double));
    DATA->Cstem_stock_unc=calloc(DATA->nodays,sizeof(double));

    DATA->Cbranch_stock=calloc(DATA->nodays,sizeof(double));
    DATA->Cbranch_stock_unc=calloc(DATA->nodays,sizeof(double));

    DATA->Ccoarseroot_stock=calloc(DATA->nodays,sizeof(double));
    DATA->Ccoarseroot_stock_unc=calloc(DATA->nodays,sizeof(double));

    DATA->Cfolmax_stock=calloc(DATA->nodays,sizeof(double));
    DATA->Cfolmax_stock_unc=calloc(DATA->nodays,sizeof(double));

    DATA->Evap=calloc(DATA->nodays,sizeof(double));

    // initialise the obs counters to 0
    DATA->ngpp=0;
    DATA->nlai=0;
    DATA->nnee=0;
  //  DATA->nwoo=0;
    DATA->nreco=0;
    DATA->nCfol_stock=0;
    DATA->nCwood_stock=0;
    DATA->nCroots_stock=0;
    DATA->nCsom_stock=0;
    DATA->nClit_stock=0;
    DATA->nCagb_stock=0;
    DATA->nCabgb_stock=0;
    DATA->nCstem_stock=0;
    DATA->nCbranch_stock=0;
    DATA->nCcoarseroot_stock=0;
    DATA->nCfolmax_stock=0;

    DATA->nevap=0;

    // loop over data to extract values
    double mettemp[DATA->nomet], obstemp[DATA->noobs];
    int tstep=0, mm;

    for (tstep = 0; tstep < DATA->nodays; tstep++)
    {
        fread(mettemp,sizeof(double),DATA->nomet,fid);
        for (mm = 0; mm < DATA->nomet; mm++)
        {
            DATA->MET[mm+tstep*DATA->nomet] = mettemp[mm];
        }

        fread(obstemp,sizeof(double),DATA->noobs,fid);

        // store values of observation streams and count the number
        // of time steps with non-missing value
        DATA->GPP[tstep] = obstemp[0]; if(obstemp[0] > -9998){DATA->ngpp = DATA->ngpp+1;}
        DATA->LAI[tstep] = obstemp[1]; if(obstemp[1] > -9998){DATA->nlai = DATA->nlai+1;}
        DATA->NEE[tstep] = obstemp[2]; if(obstemp[2] > -9998){DATA->nnee = DATA->nnee+1;}

        // create DATA->Reco vector if needed
        if (DATA->noobs > 3) 
        {
            DATA->Reco[tstep] = obstemp[3]; if(obstemp[3] > -9998){DATA->nreco = DATA->nreco+1;}
        }

        //printf("DATA->noobs %d\n", DATA->noobs);
        // create DATA->WOO if wood stocks are given
        if (DATA->noobs > 4) 
        {
            DATA->Cabgb_stock[tstep] = obstemp[4]; if(obstemp[4] > -9998){DATA->nCabgb_stock = DATA->nCabgb_stock+1;}
        }

        if (DATA->noobs > 5) 
        {
            DATA->Cwood_stock[tstep] = obstemp[5]; if(obstemp[5] > -9998){DATA->nCwood_stock = DATA->nCwood_stock+1;}
        }

        if (DATA->noobs > 6) 
        {
            DATA->Evap[tstep] = obstemp[6]; if(obstemp[6] > -9998){DATA->nevap = DATA->nevap+1;}
        }

        if (DATA->noobs > 7) 
        {
            DATA->Csom_stock[tstep] = obstemp[7]; if(obstemp[7] > -9998){DATA->nCsom_stock = DATA->nCsom_stock+1;}
        }

// JFE TO BE COMPLETED WITH ADDED OBSERVATIONS AND UNCERTAINTIES
    }

    // set constant deltat
    DATA->deltat = DATA->MET[8];
    //printf("%f", DATA->deltat);

    // allocate array to save the index of days with observations
    if (DATA->ngpp > 0)  {DATA->gpppts = calloc(DATA->ngpp,sizeof(int));}
    if (DATA->nlai > 0)  {DATA->laipts = calloc(DATA->nlai,sizeof(int));}
    if (DATA->nnee > 0)  {DATA->neepts = calloc(DATA->nnee,sizeof(int));}
    if (DATA->nreco > 0) {DATA->recopts = calloc(DATA->nreco,sizeof(int));}
    if (DATA->nCabgb_stock > 0)  {DATA->Cabgb_stockpts = calloc(DATA->nCabgb_stock,sizeof(int));}
    if (DATA->nCwood_stock > 0)  {DATA->Cwood_stockpts = calloc(DATA->nCwood_stock,sizeof(int));}  
    if (DATA->nevap > 0)  {DATA->Evap_pts = calloc(DATA->nevap,sizeof(int));}  
    if (DATA->nCsom_stock > 0)  {DATA->Csom_stockpts = calloc(DATA->nCsom_stock,sizeof(int));} 
    // counters
    int idgpp = 0, idlai = 0, idnee = 0, idreco = 0, idwoo = 0, idCabgb = 0, idCwood = 0, idEvap = 0;
    int idCsom = 0;

    // save the index of time step with observed values
    for (tstep = 0; tstep < DATA->nodays; tstep++)
    {
        if (DATA->GPP[tstep] > -9998) {DATA->gpppts[idgpp] = tstep; idgpp = idgpp+1;}
        if (DATA->LAI[tstep] > -9998) {DATA->laipts[idlai] = tstep; idlai = idlai+1;}
        if (DATA->NEE[tstep] > -9998) {DATA->neepts[idnee] = tstep; idnee = idnee+1;}

        if (DATA->nreco > 0)
        {
            if (DATA->Reco[tstep] > -9998) {DATA->recopts[idreco] = tstep; idreco = idreco+1;}
        }

        if (DATA->nCabgb_stock > 0)
        {
            if (DATA->Cabgb_stock[tstep] > -9998) {DATA->Cabgb_stockpts[idCabgb] = tstep; idCabgb = idCabgb+1;}
        }

        if (DATA->nCwood_stock > 0)
        {
            if (DATA->Cwood_stock[tstep] > -9998) {DATA->Cwood_stockpts[idCwood] = tstep; idCwood = idCwood+1;}
        }

        if (DATA->nevap > 0)
        {
            if (DATA->Evap[tstep] > -9998) {DATA->Evap_pts[idEvap] = tstep; idEvap = idEvap+1;}
        }

        if (DATA->nCsom_stock > 0)
        {
            if (DATA->Csom_stock[tstep] > -9998) {DATA->Csom_stockpts[idCsom] = tstep; idCsom = idCsom+1;}
        }

    }

    // now calculate the average temperature and radiation
    DATA->meantemp = 0; DATA->meanrad = 0;
    // first sum data
    for (tstep = 0; tstep < DATA->nodays; tstep++)
    {
        DATA->meantemp = DATA->meantemp+0.5*(DATA->MET[tstep*DATA->nomet+1]+DATA->MET[tstep*DATA->nomet+2]);
        DATA->meanrad = DATA->meanrad+DATA->MET[tstep*DATA->nomet+3];
    }
    // then average
    DATA->meantemp = DATA->meantemp/DATA->nodays;
    DATA->meanrad = DATA->meanrad/DATA->nodays;

    //printf("%d", DATA->nCabgb_stock);

    fclose(fid); // close input file
            
}

void read_pari_data(PARAMETER_INFO *PI, DATA *DATA, MCMC_OUTPUT *MCOUT, char *infile)
{
    printf("Input file = \"%s\"\n", infile);

    read_binary_data(DATA, infile);
     // choose model type to assign number of pools, fluxes and parameters
    cardamom_model_library(DATA);   

    // allocate vectors to store model output
    DATA->M_LAI = calloc(DATA->nodays,sizeof(double));
    DATA->M_GPP = calloc(DATA->nodays,sizeof(double));
    DATA->M_NEE = calloc(DATA->nodays,sizeof(double));
   // DATA->M_Reco = calloc(DATA->nodays,sizeof(double));

    DATA->M_FLUXES = calloc(DATA->nodays*DATA->nofluxes,sizeof(double));
    DATA->M_POOLS = calloc((DATA->nodays+1)*DATA->nopools,sizeof(double));

    DATA->M_ABGB_ANN = calloc(DATA->nodays/12,sizeof(double));
    

    // allocate vectors for parameters
    PI->npars = DATA->nopars;
    printf("No. pars io: %d\n", DATA->nopars);
    printf("No. pars io: %d\n", PI->npars);
    PI->parmin = calloc(PI->npars,sizeof(double));
    PI->parmax = calloc(PI->npars,sizeof(double));
    PI->parini = calloc(PI->npars,sizeof(double));
    PI->parfix = calloc(PI->npars,sizeof(double));
    PI->stepsize = calloc(PI->npars,sizeof(double));

    // pars_info is a model specific function define in models/<model_name>/src/<model_name>_PARS.c
    pars_info(PI);
    
    int parid = 0;
    for (parid = 0; parid < PI->npars; parid++) {PI->stepsize[parid] = 0.01;}

    // initialise the structure that stores the best parameter sets
    MCOUT->best_pars = calloc(PI->npars,sizeof(double));   
    MCOUT->complete = 0;
}

/************************
    read MCMC options
*************************/
void read_options(MCMC_OPTIONS *MCOPT, int solutions_wanted, int freq_print, int freq_write, char *outfile)
{
    // defining MCMC_OPTIONS structure
    MCOPT->APPEND=0;  
    MCOPT->nADAPT=100;
    MCOPT->fADAPT=0.5;

    MCOPT->randparini=0;
    MCOPT->returnpars=0;
    MCOPT->fixedpars=0;

    // set number of MCMC successful parameter sets
    if (solutions_wanted > 0)
    {
        MCOPT->nOUT = solutions_wanted;
    }
    else
    {
        MCOPT->nOUT = 1000;
        printf("Using default MCOPT->nOUT = 1000\n");
    }
    
    // set frequency of console print 
    if (freq_print >= 0)
    {
        MCOPT->nPRINT = freq_print;
    }
    else
    {
        MCOPT->nPRINT = 1000;
        printf("Using default MCOPT->nPRINT = 1000\n");
    }

    // set frequency of output
    if (freq_write >= 0)
    {
        MCOPT->nWRITE = freq_write;
    }
    else
    {
        MCOPT->nWRITE = 1000;
        printf("Using default MCOPT->nWRITE = 1000\n");
    }

    // add suffix for output file type
    char stepfile[strlen(outfile)+5], parfile[strlen(outfile)+5];

    // create file names
    strcat(strcpy(parfile,outfile),"PARS");
    strcat(strcpy(stepfile,outfile),"STEP");

    // store names in structure
    strcpy(MCOPT->parfile,parfile);
    strcpy(MCOPT->stepfile,stepfile);

    return;
}



/******************************************************************
    check for existing output files to possibly restart the run
*******************************************************************/
int check_for_existing_output_files(char *parfile, char *stepfile, PARAMETER_INFO *PI)
{
    int par_exists, step_exists;
    int restart_flag;


    /* returns restart_flag that can take the following values:
    - 0: no output files were detected, results will be written in new file
    - -1: error can be only one output files has been detected, or files are corrupted, or
    do not correspond to the current model version (different no. pars)
    if both output files have been detected and contain output consistent with the model version
    restart_flag is an integer that contains the number of parameter sets in the parfile

    We assume that restart is only used to finish an aborted run with the same number of 
    solutions and output frequency...
    */

    // check whether the output and step files exist
    par_exists = access(parfile, F_OK);
    step_exists = access(stepfile, F_OK);

    if (par_exists != -1 && step_exists !=-1)
    {
        printf("Found both output files, now checking their content\n");

        // get step sizes
        struct stat st;
        stat(stepfile, &st);
        if (st.st_size == sizeof(double)*PI->npars)
        {
            printf("Stepfile \"%s\" contains the right number of parameters.\n", stepfile);
          
            FILE *filestp = fopen(stepfile,"rb");
            fread(PI->stepsize, sizeof(double), PI->npars, filestp);
            fclose(filestp);
            restart_flag = 1 ;
        
            // now checking the parameter file    
            stat(parfile,&st);
            int noruns;
            noruns = st.st_size % (sizeof(double)*(PI->npars+1));
            if (noruns == 0)
            {
                noruns = st.st_size / (sizeof(double)*(PI->npars+1));

                printf("Parfile \"%s\" contains a total of %d runs\n", parfile, (int)noruns);
                FILE *filepar = fopen(parfile,"rb");
                fseek(filepar, sizeof(double)*(PI->npars+1)*(noruns-1), SEEK_CUR);
                fread(PI->parini, sizeof(double), PI->npars, filepar);
                fclose(filepar);   
                restart_flag = noruns;      
            }
            else
            {
                printf("Parfile does not contain a finite number of parameter sets ... it was created with another model version or crashed during writing");
                restart_flag = -1;
            }
        }
        else
        {
            printf("Stepfile does not have the correct number of parameters for this model version ... \n");
            restart_flag =-1 ;
        }
        // check number of runs

    }
    else if (par_exists == -1 && step_exists == -1)
    {
        printf("None of the output files has been found ... new output files will be created.\n");
        restart_flag = 0;
    }
    else
    {
        printf("Only one output file has been found ... please check file content.\n");
        restart_flag = -1;
    }
        
    return restart_flag;
}
