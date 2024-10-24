#pragma once
#include <math.h>
#include <stdio.h>

/*

16 Oct 2016
This version changes the way allocation of NPP is calculated. 

16 Feb 2015 - JFE:
This version adds two main components:
- updated fire module 
- economy of creating new leaves

01 Oct 2014 - JFE
GSI acts as a switch : 
- > 0.5 : labile release in foliar pool
- < 0.5 : leaf senescence

25 Sep 2014 - JFE
This version implements a Growing Season Index to scale 
the turnover rates of foliar and labile pools as a function
of minimum temperature, vapour pressure deficit and photoperiod

02 Jun 2014 - JFE
This version merges the land-use module with the fire module.

26 May 2014 - JFE
Removal is back to be done at end of time step to ease the transfer of 
its effect through time

21 May 2014 - JFE
Added a varying deltat to deal with leap years and 
MODIS at end of year 
The removal is now done at the beginning of the time-step

16 May 2014 - JFE
Added the possibility to remove some AGB 

DALEC - All Biomes: Crops, Deciduous, Evergreen,Forests
This routine is C only 


*/

double ACM(double const pars[], double const consts[]);

/* CDEA specific function
double offset(double const L, double const w) 
{
    // function is CDEA specific: 

    double mxc[7]={0.000023599784710, 0.000332730053021,    0.000901865258885,  -0.005437736864888,  -0.020836027517787,   0.126972018064287,   -0.188459767342504};

    double lf=log(L-1);
    double os=mxc[0]*pow(lf,6) + mxc[1]*pow(lf,5) + mxc[2]*pow(lf,4) + mxc[3]*pow(lf,3) + mxc[4]*pow(lf,2) + mxc[5]*lf +mxc[6];

    os=os*w;

    return os;
}
*/

void CARBON_MODEL(double const MET[], double pars[26], double const deltat, int const nr, double const lat,
double *LAI, double *NEE, double *FLUXES, double *POOLS) {

    double gpppars[11],pi;
    /*C-pools, fluxes, meteorology indices*/
    int p,f,m,nxp;
    int n=0;
    pi=3.1415927;
    double timestep; // JFE - to store timestep length
    int nn;
    /*constant gpppars terms*/
    gpppars[3]=1;
    gpppars[6]=lat;
    gpppars[8]=-2.0;
    gpppars[9]=1.0;
    gpppars[10]=pi;

    double constants[10]={pars[10],0.0156935,4.22273,208.868,0.0453194,0.37836,7.19298, 0.011136,2.1001,0.789798};

    // number of DALEC pools
    int nopools=6;
    // assigning values to pools L,F,R,W,Lit,SOM
    POOLS[0]=pars[17];
    POOLS[1]=pars[18];
    POOLS[2]=pars[19];
    POOLS[3]=pars[20];
    POOLS[4]=pars[21];
    POOLS[5]=pars[22];

    /* NOTES FOR FLUXES
    MET[:,0]: projday
    MET[:,1]: mintemp
    MET[:,2]: maxtemp
    MET[:,3]: rad
    MET[:,4]: co2
    MET[:,5]: DoY

    JFE added Jun 2014
    MET[:,6]: removal
    MET[:,7]: removal by fire (i.e. leading to emissions)
    MET[:,8]: time-step length (in days)

    JFE added Sep 2014 for GSI
    MET[:,9]:  21 day average Tmn
    MET[:,10]: 21 day average photoperiod
    MET[:,11]: 21 day average VPD
    
    number of MET drivers*/
    int nomet=12;
     
    /*fluxes - other*********
    0.GPP
    1.temprate
    2.respiration_auto
    3.leaf_production
    4.labile_production
    5.root_production
    6.wood_production
    7.labile_release
    8.leaffall_factor
    9.leaflitter_production
    10.woodlitter_production  
    11.rootlitter_production         
 	12.respiration_het_litter
  	13.respiration_het_som
  	14.litter2som
  	15.labrelease_factor

    16.release by fire - Jun 2014 / added by JFE 

    17-22: release by fire for each pool
    23-27: fluxes between pools

    28.GSI value - Sep 2014 / Added by JFE

    number of DALEC fluxes to store*/
    int nofluxes=29;

    // Fire specific variables - JFE
    double CF[6] = {0.1,0.9,0.1,0.1,0.5,0.01}; // combustion efficiencies: fraction of pool emitted when burnt
    
    // resilience factor
    double rfac=0.5; // fraction of the none-burnt pool that survives 

    // repeating loop for each timestep
    for ( n=0; n < nr; n++) {

        // JFE - replaced constant deltat by dynamic timestep

        // pool index
        p=nopools*n;
        // next pool index
        nxp=nopools*(n+1);
        // met index
        m=nomet*n;
        // flux array index
        f=nofluxes*n;
       
        timestep=MET[m+8];

        // LAI: used to calculate GPP
        LAI[n]=POOLS[p+1]/pars[16]; 
  
        // GPP
        gpppars[0]=LAI[n];
        gpppars[1]=MET[m+2];
        gpppars[2]=MET[m+1];
        gpppars[4]=MET[m+4];
        gpppars[5]=MET[m+5];
        gpppars[7]=MET[m+3];

        FLUXES[f+0]=ACM(gpppars,constants);

        // temprate - now comparable to Q10 - factor at 0C is 1
        FLUXES[f+1]=exp(pars[9]*0.5*(MET[m+2]+MET[m+1]));
        // respiration auto
        FLUXES[f+2]=pars[1]*FLUXES[f+0];
        // JFE changed allocation to leaves 16/10/15
        // FLUXES[f+3]=(FLUXES[f+0]-FLUXES[f+2])*pars[2];
        FLUXES[f+3]=((FLUXES[f+0]-FLUXES[f+2])*pars[2])*(1-pars[12]);
        // JFE allocation to labile pool 16/10/15
        //FLUXES[f+4] = (FLUXES[f+0]-FLUXES[f+2]-FLUXES[f+3])*pars[12];
        FLUXES[f+4] = ((FLUXES[f+0]-FLUXES[f+2])*pars[2])*pars[12];                            
        // allocation to roots
        FLUXES[f+5] = (FLUXES[f+0]-FLUXES[f+2]-FLUXES[f+3]-FLUXES[f+4])*pars[3];            
        // allocation to wood 
        FLUXES[f+6] = FLUXES[f+0]-FLUXES[f+2]-FLUXES[f+3]-FLUXES[f+5]-FLUXES[f+4]; 

        /* 25/09/14 - JFE
        Here we calculate the Growing Season Index based on 
        - Jolly et al. A generalized, bioclimatic index to predict foliar phenology in response to climate Global Change Biology, Volume 11, page 619-632 - 2005 (doi: 10.1111/j.1365-2486.2005.00930.x) 
        - Stoeckli, R., T. Rutishauser, I. Baker, M. A. Liniger, and A. S. Denning (2011), A global reanalysis of vegetation phenology, J. Geophys. Res., 116, G03020, doi:10.1029/2010JG001545.
        
        It is the product of 3 limiting factors for temperature, photoperiod and vapour pressure deficit 
        that grow linearly from 0 to 1 between a calibrated min and max value.
        Temperature, photoperiod and VPD are direct input
        */    

        // calculate the temperature based limiting factor
        double Tfac;
        if (MET[m+9]>pars[14]) // Tmin larger than parameter Tmmax
        {
            Tfac = 1.;
        }
        else if (MET[m+9]<(pars[13])) // Tmin lower than parameter Tmmin
        {
            Tfac = 0.;         
        }
        else
        {
            Tfac = (MET[m+9]-(pars[13]))/(pars[14]-(pars[13]));
        }

// calculate the photoperiod limiting factor
        double PHOTOfac;
        if (MET[m+10]>pars[23]) // photoperiod larger than parameter PHOTOmax
        {
            PHOTOfac = 1.;
        }
        else if (MET[m+10]<pars[15]) // photoperiod lower than parameter PHOTOmin
        {
            PHOTOfac = 0.;         
        }
        else
        {
            PHOTOfac = (MET[m+10]-pars[15])/(pars[23]-pars[15]);
        }

// calculate the VPD based limiting factor - high VPD is limiting
        double VPDfac;
        if (MET[m+11]>pars[25]) // VPD larger than parameter VPDmax
        {
            VPDfac = 0.;
        }
        else if (MET[m+11]<pars[24]) // VPD lower than parameter VPDmin
        {
            VPDfac = 1.;         
        }
        else
        {
            VPDfac = 1-(MET[m+11]-pars[24])/(pars[25]-pars[24]);
        }

// calculate the GSI and store it in the FLUXES array with ID 17
        FLUXES[f+28] = Tfac*PHOTOfac*VPDfac;
        
        FLUXES[f+8]=(1-FLUXES[f+28])*pars[4]; // leaf fall
        FLUXES[f+15]=FLUXES[f+28]*pars[11]; // lab release
       
        // Average biogenic fluxes over the time step

        // total labile release into foliar
        FLUXES[f+7] = POOLS[p+0]*(1-pow(1-FLUXES[f+15],timestep))/timestep;             
      
        // leaf to litter
        FLUXES[f+9] = POOLS[p+1]*(1-pow(1-FLUXES[f+8],timestep))/timestep;
 
        // wood to SOM
        FLUXES[f+10] = POOLS[p+3]*(1-pow(1-pars[5],timestep))/timestep;    

        // root to litter
        FLUXES[f+11] = POOLS[p+2]*(1-pow(1-pars[6],timestep))/timestep;                                    

        // respiration heterotrophic from litter
        FLUXES[f+12] = POOLS[p+4]*(1-pow(1-FLUXES[f+1]*pars[7],timestep))/timestep;
                
        // respiration heterotrophic from SOM
        FLUXES[f+13] = POOLS[p+5]*(1-pow(1-FLUXES[f+1]*pars[8],timestep))/timestep;                  

        // decomposition of litter into SOM
        FLUXES[f+14] = POOLS[p+4]*(1-pow(1-pars[0]*FLUXES[f+1],timestep))/timestep; 
       

        // Update pools

        // labile pool: +allocation of NPP -leaf production
        POOLS[nxp+0] = POOLS[p+0] + (FLUXES[f+4]-FLUXES[f+7])*timestep;

        // foliar pool: +allocation of NPP -leaf fall + production from labile 
        POOLS[nxp+1] =  POOLS[p+1] + (FLUXES[f+3] - FLUXES[f+9] + FLUXES[f+7])*timestep;

        // wood pool: +allocation of NPP -wood transfer to SOM
        POOLS[nxp+3] = POOLS[p+3] +  (FLUXES[f+6] - FLUXES[f+10])*timestep;

        // root pool: +allocation of NPP -transfer to litter
        POOLS[nxp+2] = POOLS[p+2] + (FLUXES[f+5] - FLUXES[f+11])*timestep;

        // litter pool: +leaf fall +root to litter -rh_litter -decomposition
        POOLS[nxp+4] = POOLS[p+4] + (FLUXES[f+9] + FLUXES[f+11] - FLUXES[f+12] - FLUXES[f+14])*timestep;                

        // SOM pool: +decomposition of litter -rh_som +wood_transfer
        POOLS[nxp+5]= POOLS[p+5]+ (FLUXES[f+14] - FLUXES[f+13] + FLUXES[f+10])*timestep;

        NEE[n]=-FLUXES[f+0]+FLUXES[f+2]+FLUXES[f+12]+FLUXES[f+13]; // NEE only contains biogenic fluxes

       // Reco[n] = FLUXES[f+2]+FLUXES[f+12]+FLUXES[f+13];

        // perform the deforestation removal of labile, leaf and wood if require
        if (MET[m+6] > 0.) 
        {
            POOLS[nxp+0] = POOLS[nxp+0]*(1-MET[m+6]);
            POOLS[nxp+1] = POOLS[nxp+1]*(1-MET[m+6]);
            POOLS[nxp+3] = POOLS[nxp+3]*(1-MET[m+6]);    
        } // end removal

        FLUXES[f+16]=0.; // fire emissions
        // perform the fire part and save all fluxes
        if (MET[m+7] > 0.) 
        {
        	// Calculate all fire transfers (1. combustion, and 2. litter transfer) note: all fluxes are in gC m-2 day-1
	        for (nn=0;nn<6;nn++){FLUXES[f+17+nn] = POOLS[nxp+nn]*MET[m+7]*CF[nn]/timestep;} // combustion
	        for (nn=0;nn<5;nn++){FLUXES[f+23+nn] = POOLS[nxp+nn]*MET[m+7]*(1-CF[nn])*(1-rfac)/timestep;} // litter transfer

	        // Update pools
	        // Live C pools
	        for (nn=0;nn<4;nn++) {POOLS[nxp+nn]=POOLS[nxp+nn]-(FLUXES[f+17+nn]+FLUXES[f+23+nn])*timestep;}
	        // Dead C pools Litter and SOM
	        POOLS[nxp+4]=POOLS[nxp+4]+(FLUXES[f+23]+FLUXES[f+23+1]+FLUXES[f+23+2]-FLUXES[f+17+4]-FLUXES[f+23+4])*timestep;
	        POOLS[nxp+5]=POOLS[nxp+5]+(FLUXES[f+23+3]+FLUXES[f+23+4]-FLUXES[f+17+5])*timestep;

	        // save the sum of all fluxes
            for (nn=0;nn<6;nn++){FLUXES[f+16]+=FLUXES[f+17+nn];}

        } 
        else // be sure fluxes are 0
        {
            for (nn=0;nn<6;nn++) {FLUXES[f+17+nn] = 0;}
            for (nn=0;nn<5;nn++) {FLUXES[f+23+nn] = 0;}
        }// end fire  

    } // end time steploop
} // end CARBON_MDEOl



double ACM(double const pars[], double const consts[])
{
    // function calculates GPP
    double gc,pp,qq,ci,e0,mult,dayl,cps,dec,GPP;

    gc=(double)pow(fabs(pars[8]),consts[9])/(consts[5] * pars[9] + 0.5 * ( pars[1]- pars[2]));
    pp=(double)pars[0]*pars[3]/gc*consts[0]*exp(consts[7]*pars[1]);
    qq=(double)consts[2]-consts[3];
    ci=(double)0.5*(pars[4]+qq-pp+pow(pow(pars[4]+qq-pp,2)-4*(pars[4]*qq-pp*consts[2]),0.5));
    e0=(double)consts[6]*pow(pars[0],2)/(pow(pars[0],2)+consts[8]);
    dec=(double)-23.4*cos((360.*(pars[5]+10.)/365.)*pars[10]/180.)*pars[10]/180.;
    mult=(double)tan(pars[6]*pars[10]/180)*tan(dec);


    if (mult>=1) {dayl=24.;} else if(mult<=-1) {dayl=0.;} else {dayl=(double)24.*acos(-mult) / pars[10];}
    cps=(double)e0*pars[7]*gc*(pars[4]-ci)/(e0*pars[7]+gc*(pars[4]-ci));
    GPP=cps*(consts[1]*dayl+consts[4]);
    return GPP;
}





