#pragma once
#include <math.h>
#include <stdio.h>

double ACM(double const pars[], double const consts[])
{
  /*double gc=0,pp=0,qq=0,ci=0,e0=0,mult=0,dayl=0,cps=0,dec=0;*/
  double gc,pp,qq,ci,e0,mult,dayl,cps,dec,GPP;
  /*pars= &pars;  
 *   consts= &consts;
 *   */
  gc=(double)pow(fabs(pars[8]),consts[9])/(consts[5] * pars[9] + 0.5 * ( pars[1]- pars[2]));
  pp=(double)pars[0]*pars[3]/gc*consts[0]*exp(consts[7]*pars[1]);
  qq=(double)consts[2]-consts[3];
  ci=(double)0.5*(pars[4]+qq-pp+pow(pow(pars[4]+qq-pp,2)-4*(pars[4]*qq-pp*consts[2]),0.5));
  e0=(double)consts[6]*pow(pars[0],2)/(pow(pars[0],2)+consts[8]);
  dec=(double)-23.4*cos((360.*(pars[5]+10.)/365.)*pars[10]/180.)*pars[10]/180.;
  mult=(double)tan(pars[6]*pars[10]/180)*tan(dec);
  if (mult>=1){ 
   dayl=24.;}  
  else if(mult<=-1)
  dayl=0.;
  else{
  dayl=(double)24.*acos(-mult) / pars[10];}


  cps=(double)e0*pars[7]*gc*(pars[4]-ci)/(e0*pars[7]+gc*(pars[4]-ci));
  GPP=cps*(consts[1]*dayl+consts[4]);
  return GPP;
}


double offset(double const L, double const w) {

/*solution to t = f*exp(-t^2)*/
/*see dalecstepfunction2*/

double mxc[7]={0.000023599784710, 0.000332730053021,    0.000901865258885,  -0.005437736864888,  -0.020836027517787,   0.126972018064287,   -0.188459767342504};

double lf=log(L-1);
double os=mxc[0]*pow(lf,6) + mxc[1]*pow(lf,5) + mxc[2]*pow(lf,4) + mxc[3]*pow(lf,3) + mxc[4]*pow(lf,2) + mxc[5]*lf +mxc[6];

os=os*w;

return os;
}

void CARBON_MODEL(double const *MET,double const *pars, double const deltat,int const nr, double const lat,
double *LAI, double *NEE, double *FLUXES, double *POOLS) {


    //printf("hello from DALEC_CDEA\n");
    double gpppars[11],pi;
    /*C-pools, fluxes, meteorology indices*/
    int p,f,m,nxp, i;
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

    /*assigning values to pools*/
    /*L,F,R,W,Lit,SOM*/
    POOLS[0]=pars[17];
    POOLS[1]=pars[18];
    POOLS[2]=pars[19];
    POOLS[3]=pars[20];
    POOLS[4]=pars[21];
    POOLS[5]=pars[22];

    // water storage is a fraction of wp to fc
    POOLS[6]=pars[25]*(pars[24]-pars[23])+pars[23];

    // number of DALEC pools: 6 C pools + 1 water pool
    int nopools=7;

    /* NOTES FOR POOLS AND FLUXES
    MET[:,0]: projday
    MET[:,1]: mintemp
    MET[:,2]: maxtemp
    MET[:,3]: rad
    MET[:,4]: co2
    MET[:,5]: yearday

    JFE added Jun 2014
    MET[:,6]: removal
    MET[:,7]: removal by fire (i.e. leading to emissions)
    MET[:,8]: time-step length (in days)
    
    JFE added Oct 2016 for water
    MET[:,9]: precipitation mm/d
    MET[:,10]: ET reference
    /*number of MET drivers*/
    int nomet=11;
     
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

    28: ETact
    29: perc
    

    number of DALEC fluxes to store*/
    int nofluxes=30;


    /*constants for exponents of leaffall and labrelease factors*/
    /*width*/
    double wf=pars[15]*sqrt(2)/2;
    double wl=pars[13]*sqrt(2)/2;

    /*factor*/
    double ff=(log(pars[4])-log(pars[4]-1))/2;
    double fl=(log(1.001)-log(0.001))/2;

    /*additional offset*/
    double osf=offset(pars[4],wf);
    double osl=offset(1.001,wl);

    /*scaling to biyearly sine curve*/
    double sf=365.25/pi;

    // Fire specific variables - JFE
    /*fully combusted fire flux (all pools except SOM)*/
    double CFF[5];

    /*fully combusted fire flux (all pools except SOM)*/
    double NCFF[5];

    /*combustion efficiencies*/
    double CF[6] = {0.1,0.9,0.1,0.1,0.5,0.01};

    /*resilience factor*/
    double rfac=0.5;

    // water storage parameters
    double fc = pars[24];
    double wp = pars[23];
    double opt = 0.5*(wp+fc);

    // dummy terms to store daily GPP, ET, WS, and monthly Vs and ET
    double dummyGPP,dummyET,dummyWS, dummyVs, dummyVs_sum, dummyET_sum;
    // integer to loop over number of days per time step
    int dd=0;  

    //for (n =0; n<26;n++) {printf("%i: %f\n", n, pars[n]);}


    /*repeating loop for each timestep*/
    for ( n=0; n < nr; n++) {

        /*ppol index*/
        p=nopools*n;
        /*next pool index*/
        nxp=nopools*(n+1);
        /*met index*/
        m=nomet*n;
        /*flux array index*/
        f=nofluxes*n;
        //timestep length
        timestep = MET[m+8];
        
        /*LAI*/
        LAI[n]=POOLS[p+1]/pars[16]; 
  
        /*GPP*/
        gpppars[0]=LAI[n];
        gpppars[1]=MET[m+2];
        gpppars[2]=MET[m+1];
        gpppars[4]=MET[m+4];
        gpppars[5]=MET[m+5];
        gpppars[7]=MET[m+3];

        /*GPP - daily timestep*/
        FLUXES[f+0]=ACM(gpppars,constants);

        /* added Oct 2016: calculation of water balance to scale GPP
        wp: pars[24]
        fc: pars[25]
        */

      
        // all fluxes are calculated in mm/day
        double ETpot = MET[m+10]*LAI[n]/2.88; //adjust ET from ref crop to current vegetation

        dummyWS     = POOLS[p+6];
        dummyET     = 0;
        dummyGPP    = 0;
        dummyVs     = 0;
        dummyET_sum = 0;
        dummyVs_sum = 0;
        

        for (dd=0;dd<timestep;dd++) {
            // increment a daily time step for water and GPP
            // start with ET
            if (dummyWS <= wp)
            {
                dummyET = 0.;
            }
            else if (dummyWS > wp && dummyWS <= opt)
            {
                dummyET = ETpot*(POOLS[p+6]-wp)/(opt-wp);  // et_pot[tt] * (Vt-wp)/(opt-wp)
            }
            else if (dummyWS >= opt)
            {
                dummyET = ETpot;
            }
            
            
            //save daily ET to get monthly sum
            dummyET_sum = dummyET_sum + dummyET;
            // enough water for outflow?
            double Vs = dummyWS + MET[m+9] - fc - dummyET;
            if (Vs > 0)
            {
                dummyVs = Vs*exp(-1);
            }
            else
            {   
                dummyVs = 0;
            }
            dummyVs_sum = dummyVs_sum + dummyVs;
            // update water pool: +P - ETact - Perc
            dummyWS = dummyWS + (MET[m+9] - dummyET - dummyVs);

            // correct GPP as a function of ETact/ETpot
            dummyGPP = dummyGPP + FLUXES[f+0] * ( dummyET / ETpot );

        }

        // save states and fluxes
        POOLS[nxp+6] = dummyWS;
        FLUXES[f+0] = dummyGPP/timestep;
     //   printf("GPP %f\n",dd,FLUXES[f+0]);
        FLUXES[f+28] = dummyET_sum/timestep;
        FLUXES[f+29] = dummyVs_sum/timestep;

        /*temprate - now comparable to Q10 - factor at 0C is 1*/
        FLUXES[f+1]=exp(pars[9]*0.5*(MET[m+2]+MET[m+1]));
        /*respiration auto*/
        FLUXES[f+2]=pars[1]*FLUXES[f+0];
        /*leaf production rate*/
        FLUXES[f+3]=(FLUXES[f+0]-FLUXES[f+2])*pars[2];
        /*labile production*/
        FLUXES[f+4] = (FLUXES[f+0]-FLUXES[f+2]-FLUXES[f+3])*pars[13-1];              

        /*root production*/        
        FLUXES[f+5] = (FLUXES[f+0]-FLUXES[f+2]-FLUXES[f+3]-FLUXES[f+4])*pars[4-1];            

        /*wood production rate*/       
        FLUXES[f+6] = FLUXES[f+0]-FLUXES[f+2]-FLUXES[f+3]-FLUXES[f+5]-FLUXES[f+4]; 


        /*leaf fall factor*/
        /*FLUXES[f+8]=pow((0.75+cos((MET[m+5]-pars[15-1])*(2*pi/365.25))/4),10)*pars[5-1];*/
        /*FLUXES[f+8]=1- pow((1 +0.00001 - 1/pars[4]),nf*exp(nf*sin((MET[m+0]-pars[14])/sf)*sf)/pow(1+exp(nf*sin((MET[m+0]-pars[14])/sf)*sf),2));*/
        FLUXES[f+8] = (2/sqrt(pi))*(ff/wf)*exp(-pow(sin((MET[m+0]-pars[14]+osf)/sf)*sf/wf,2));


        /*Labrelease factor*/
        FLUXES[f+15]=(2/sqrt(pi))*(fl/wl)*exp(-pow(sin((MET[m+0]-pars[11]+osl)/sf)*sf/wl,2));
         
        // JFE - replaced constant deltat by dynamic timestep / now at beginning of loop
        //timestep=MET[m+8];
       // printf("Adapting timestep %f", timestep);

        /*total labile release - mean at daily timestep */
        //FLUXES[f+7]=POOLS[p+0]*(1-pow(1-FLUXES[f+15],deltat))/deltat;
        FLUXES[f+7]=POOLS[p+0]*(1-pow(1-FLUXES[f+15],timestep))/timestep;             
      
        /*total leaf litter production*/       
        //FLUXES[f+9] = POOLS[p+1]*(1-pow(1-FLUXES[f+8],deltat))/deltat;                     
        FLUXES[f+9] = POOLS[p+1]*(1-pow(1-FLUXES[f+8],timestep))/timestep;
 
        /*total wood litter production*/       
        //FLUXES[f+10] = POOLS[p+3]*(1-pow(1-pars[6-1],deltat))/deltat;                    
        FLUXES[f+10] = POOLS[p+3]*(1-pow(1-pars[6-1],timestep))/timestep;    

        /*root litter production*/
        //FLUXES[f+11] = POOLS[p+2]*(1-pow(1-pars[7-1],deltat))/deltat;                    
        FLUXES[f+11] = POOLS[p+2]*(1-pow(1-pars[7-1],timestep))/timestep;                                    

        /*respiration heterotrophic litter*/
        //FLUXES[f+12] = POOLS[p+4]*(1-pow(1-FLUXES[f+1]*pars[8-1],deltat))/deltat;
        FLUXES[f+12] = POOLS[p+4]*(1-pow(1-FLUXES[f+1]*pars[8-1],timestep))/timestep;
                
        /*respiration heterotrophic SOM*/
        //FLUXES[f+13] = POOLS[p+5]*(1-pow(1-FLUXES[f+1]*pars[9-1],deltat))/deltat;
        FLUXES[f+13] = POOLS[p+5]*(1-pow(1-FLUXES[f+1]*pars[9-1],timestep))/timestep;                  

        /*litter to SOM*/
        //FLUXES[f+14] = POOLS[p+4]*(1-pow(1-pars[1-1]*FLUXES[f+1],deltat))/deltat; 
        FLUXES[f+14] = POOLS[p+4]*(1-pow(1-pars[1-1]*FLUXES[f+1],timestep))/timestep; 
       

        /* update pools */

        // labile pool
        //POOLS[nxp+0] = POOLS[p+0] + (FLUXES[f+4]-FLUXES[f+7])*deltat;
        POOLS[nxp+0] = POOLS[p+0] + (FLUXES[f+4]-FLUXES[f+7])*timestep;

        // foliar pool
        //POOLS[nxp+1] =  POOLS[p+1] + (FLUXES[f+3] - FLUXES[f+9] + FLUXES[f+7])*deltat;
        POOLS[nxp+1] =  POOLS[p+1] + (FLUXES[f+3] - FLUXES[f+9] + FLUXES[f+7])*timestep;

        // wood pool
        //POOLS[nxp+3] = POOLS[p+3] +  (FLUXES[f+6] - FLUXES[f+10])*deltat;
        POOLS[nxp+3] = POOLS[p+3] +  (FLUXES[f+6] - FLUXES[f+10])*timestep;

        // root pool
        //POOLS[nxp+2] = POOLS[p+2] + (FLUXES[f+5] - FLUXES[f+11])*deltat;
        POOLS[nxp+2] = POOLS[p+2] + (FLUXES[f+5] - FLUXES[f+11])*timestep;

        // litter pool
        //POOLS[nxp+4] = POOLS[p+4] + (FLUXES[f+9] + FLUXES[f+11] - FLUXES[f+12] - FLUXES[f+14])*deltat;
        POOLS[nxp+4] = POOLS[p+4] + (FLUXES[f+9] + FLUXES[f+11] - FLUXES[f+12] - FLUXES[f+14])*timestep;                

        // SOM pool
        //POOLS[nxp+5]= POOLS[p+5]+ (FLUXES[f+14] - FLUXES[f+13]+FLUXES[f+10])*deltat;
        POOLS[nxp+5]= POOLS[p+5]+ (FLUXES[f+14] - FLUXES[f+13]+FLUXES[f+10])*timestep;

        NEE[n]=-FLUXES[f+0]+FLUXES[f+2]+FLUXES[f+12]+FLUXES[f+13];

        /* perform the deforestation removal bit if required */
        if (MET[m+6] > 0.) 
        {
         //   
            POOLS[nxp+0] = POOLS[nxp+0]*(1-MET[m+6]);
            POOLS[nxp+1] = POOLS[nxp+1]*(1-MET[m+6]);
            POOLS[nxp+3] = POOLS[nxp+3]*(1-MET[m+6]);
      //      printf("End deforestation %03i\n %f\n", n, POOLS[nxp+0]);
        } // end removal

        // perform the fire part if required
        FLUXES[f+16]=0.;
        // perform the fire part if required
        if (MET[m+7] > 0.) 
        {
            /*first fluxes*/
         //   printf("Fire detected\n %f\n", POOLS[nxp+0]);
	        /*LABILE*/
        		/*Calculating all fire transfers (1. combustion, and 2. litter transfer)*/
	            /*note: all fluxes are in gC m-2 day-1*/
	        for (nn=0;nn<6;nn++){FLUXES[f+17+nn] = POOLS[nxp+nn]*MET[m+7]*CF[nn]/timestep;}
	        for (nn=0;nn<5;nn++){FLUXES[f+23+nn] = POOLS[nxp+nn]*MET[m+7]*(1-CF[nn])*(1-rfac)/timestep;}

	        /*Adding all fire pool transfers here*/
	        /*live C pools*/	
	        for (nn=0;nn<4;nn++){POOLS[nxp+nn]=POOLS[nxp+nn]-(FLUXES[f+17+nn]+FLUXES[f+23+nn])*timestep;}
	        /*dead C pools*/
	        /*litter*/
	        POOLS[nxp+4]=POOLS[nxp+4]+(FLUXES[f+23]+FLUXES[f+23+1]+FLUXES[f+23+2]-FLUXES[f+17+4]-FLUXES[f+23+4])*timestep;
	        /*som*/
	        POOLS[nxp+5]=POOLS[nxp+5]+(FLUXES[f+23+3]+FLUXES[f+23+4]-FLUXES[f+17+5])*timestep;

	        /*fires - total flux in gC m-2 day-1*/
	        /*this term is now (essentially) obsolete*/
	        /*replace in next version of DALEC_FIRES*/
	        for (nn=0;nn<6;nn++){FLUXES[f+16]+=FLUXES[f+17+nn];}

	        /*Net ecosystem exchange*/
	        /*WARNING WARNING WARNING*/
	        /*NEE should not contain fires (or should it?)*/
	        /*Always check!!!*/
        } // end fire   // end fire  
 
    } // end daily loop

} // end DALEC_CDEA_LU




