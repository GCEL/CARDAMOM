#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../../../general/cardamom_structures.c"

struct EDCDIAGNOSTIC
{
    int EDC;
    int DIAG;
    /*allocating space for 100 checks (in case future additions are made, etc.).*/
    int PASSFAIL[100];
    int nedc;
};

int EDC1_CDEA_GLEAM(double const pars[], struct EDCDIAGNOSTIC *EDCD, double meantemp, double meanrad)
{
    // ALL EDCs set as 1 initially
    EDCD->nedc=100;
    int n; for (n=0;n<EDCD->nedc;n++){EDCD->PASSFAIL[n]=1;}

    // declaring variables and constants for subsequent EDCs
    int EDC=1;
    int DIAG=EDCD->DIAG;
    double const fauto=pars[1];
    double const ffol=(1-fauto)*pars[2];
    double const flab=(1-fauto-ffol)*pars[12];
    double const froot=(1-fauto-ffol-flab)*pars[3];
    double const fwood=1-fauto-ffol-flab-froot;

    // added constant values for allocations wrt NPP
    double const ffol_npp = ffol/(1-fauto);
    double const flab_npp = flab/(1-fauto);
    double const froot_npp = froot/(1-fauto);
    double const fwood_npp = fwood/(1-fauto);
    /*fraction of GPP som under equilibrium conditions*/
    //double const fsom=fwood+(froot+flab+ffol)*pars[0]/(pars[0]+pars[7]);
 
    /*  EDC CHECK NO 1
        TOR_LIT faster than TOR_SOM */
    if ((EDC==1 || DIAG==1) && (pars[8]>pars[7])){EDC=0;EDCD->PASSFAIL[0]=0;}

    /*  EDC CHECK NO 2
        Litter2SOM greater than SOM to Atm. rate */
    if ((EDC==1 || DIAG==1) && (pars[0]<pars[8])){EDC=0;EDCD->PASSFAIL[1]=0;}

    /*  EDC CHECK NO 3
        TOR_FOL faster than TOR_WOOD */
    if ((EDC==1 || DIAG==1) && (pars[5]>pars[4])){EDC=0;EDCD->PASSFAIL[2]=0;}

    /*  EDC CHECK NO 4
        Root turnover greater than SOM turnover at meantemp */
    if ((EDC==1 || DIAG==1) && (pars[6]<pars[8]*exp(pars[9]*meantemp))){EDC=0;EDCD->PASSFAIL[3]=0;}

    /*  EDC CHECK NO 5 */
    if ((EDC==1 || DIAG==1) && (((ffol+flab)>5*froot) || ((ffol+flab)*5<froot))){EDC=0;EDCD->PASSFAIL[4]=0;}

    // check actual allocation fraction -flab,ffol,froot,fwood in that order: from EDC ID 40 onwards
    //amended to use fnpp wrt discussion w Luke and Malhi et al 10.1098/rstb.2011.0062 10.1111/nph.14189
//    if ((EDC==1 || DIAG==1) && ((flab<0.01) || (flab > 0.9))){EDC=0;EDCD->PASSFAIL[40]=0;}
//    if ((EDC==1 || DIAG==1) && ((ffol<0.01) || (ffol > 0.9))){EDC=0;EDCD->PASSFAIL[41]=0;}
//    if ((EDC==1 || DIAG==1) && ((froot<0.01) || (froot > 0.9))){EDC=0;EDCD->PASSFAIL[42]=0;}
//    if ((EDC==1 || DIAG==1) && ((fwood<0.01) || (fwood > 0.9))){EDC=0;EDCD->PASSFAIL[43]=0;}

//    if ((EDC==1 || DIAG==1) && ((flab<0.01) || (flab > 0.9))){EDC=0;EDCD->PASSFAIL[40]=0;}
    if ((EDC==1 || DIAG==1) && (ffol_npp<0.10)){EDC=0;EDCD->PASSFAIL[41]=0;}
    if ((EDC==1 || DIAG==1) && (froot_npp<0.05)){EDC=0;EDCD->PASSFAIL[42]=0;}
    if ((EDC==1 || DIAG==1) && (fwood_npp>0.6)){EDC=0;EDCD->PASSFAIL[43]=0;}
//    if ((EDC==1 || DIAG==1) && ((fwood<0.01) || (fwood > 0.6))){EDC=0;EDCD->PASSFAIL[43]=0;}
    

    /*  GSI specific checks
        Tmn,min < Tmn,max */
    //if ((EDC==1 || DIAG==1) && (pars[13]>pars[14])){EDC=0;EDCD->PASSFAIL[5]=0;}

    //  Photoperiod,min < Photoperiod,max
    //if ((EDC==1 || DIAG==1) && (pars[15]>pars[23])){EDC=0;EDCD->PASSFAIL[6]=0;}

    //  VPD,min < VPD,max
    // if ((EDC==1 || DIAG==1) && (pars[24]>pars[25])){EDC=0;EDCD->PASSFAIL[7]=0;}

 //   for (n=0;n<5;n++) {if (EDCD->PASSFAIL[n]==0) {printf("fail %d\n", n);}}


    return EDC;
}


double mean_pool(double *PA, int p, int nc, int nopools)
{
    // variable declarations
    int c;
    double meanpool=0;
    // deriving mean of pool p
    for (c=0;c<nc;c++) {meanpool=meanpool+PA[c*nopools+p]/(double)nc;}
    
    return meanpool;
}

double mean_annual_pool(double *POOLS, int year, int pool, int nopools,int deltat){
/*inputs
 * POOLS: Pools double pointer, as output from DALEC
 * year: year for which to average (first year = 0)
 * pool: the specific pool 
 * nc
declarations*/
int c;
double meanpool=0;
// deriving mean of pool p
int stday=floor(365.25*year/deltat);
int enday=floor(365.25*(year+1)/deltat);
for (c=stday;c<enday;c++){
meanpool=meanpool+POOLS[c*nopools+pool]/(enday-stday);}
// returing meanpool value*/
return meanpool;}


double expdecay2(double const *POOLS, int pool, int nodays,int deltat,int nopools)
/*only accepting pool and number of days, the rest done here*/
{
/*using 365 day averaging window for each year!*/

/*explicitly calculating number of years*/
//int noyears=floor(nodays*deltat/365);
int n;

/*yearly means and yearly means with offset Ds = 100 days*/
/*Deriving exponential decay coefficients a,b and c in equation
 * Cexp = a + b*exp(c*t)*/
/*four mean pool values to be derived*/
/* MP0 = mean pool (year 1 to year end-2)
 * MP1 = mean pool (year 2 to year end-1)
 * MP0os = mean pool (year 1+os to year end-2+os)
 * MP1os = mean pool (year 1+os to year end-2+os)*/
double MP0=0,MP1=0,MP0os=0,MP1os=0;
/*averaging window*/
int aw=floor(365./(double)deltat);
/*offset in days*/
/*OFFSET = 1 is ideal to capture rapid exponential decay without compromising quality of fit*/
int os=1;
/*deriving coefficients to *
 * in Cexp = A + B exp(C*t);*/
//double A,B,b,C;

/*mean pool within each averaging window*/
for (n=0;n<aw;n++){MP0=MP0+POOLS[n*nopools+pool];};MP0=MP0/(double)aw;
for (n=aw;n<aw*2;n++){MP1=MP1+POOLS[n*nopools+pool];};MP1=MP1/(double)aw;
for (n=os;n<aw+os;n++){MP0os=MP0os+POOLS[n*nopools+pool];};MP0os=MP0os/(double)aw;
for (n=aw+os;n<aw*2+os;n++){MP1os=MP1os+POOLS[n*nopools+pool];};MP1os=MP1os/(double)aw;

/*deriving mean gradient ratio dcdt1/dcdt0*/
/*dcdt1 is the numeric gradient between n+1 and n+365+1*/
/*dcdt0 is the numeric gradient between n and n+365*/
//double dcdtr=0,dcdt1,dcdt0;

double dcdty1=(MP1os-MP0os);
double dcdty0=(MP1-MP0);

/*denominators*/

/*using multiyear mean to determine c*/
double C=log(dcdty1/dcdty0)/((double)os*deltat);
/*deriving final exp. decay fit with startpoint = startpoint of pool*/

/*compared and validated against dalec_expdecay3.m*/
/*
printf("noyears = %d\n",noyears);
printf("MP0 = %f\n",MP0);
printf("MP1 = %f\n",MP1);
printf("MP0os = %f\n",MP0os);
printf("MP1os = %f\n",MP1os);
printf("A = %e\n",A);
printf("B = %e\n",B);
printf("C = %e\n",C);
printf("aw = %d\n",aw);
printf("dcdty1 = %f\n",dcdt1);
printf("dcdty0 = %f\n",dcdt0);
*/

/*half life must be more than noyears*/

/*if (fabs(-log(2)/C)<noyears*365.25 && C<0 && finite(C) && abs(B)>abs(A)*0.01){EDC=0;}*/


return C;


}


int EDC2_CDEA_GLEAM(double const pars[], double *MET, double *LAI, double *NEE, double *POOLS,double *FLUXES, PARAMETER_INFO PI, DATA D, struct EDCDIAGNOSTIC *EDCD, double meantemp)
{

// THESE EDC2 checks are for DALEC_CDEA
int EDC=1,n=0;
int DIAG=EDCD->DIAG; // 1

int nomet=D.nomet, nopools=D.nopools; 
int nofluxes=D.nofluxes, nodays = D.nodays; 

int deltat=(int)(MET[nomet]-MET[0]);

// deriving mean pools here
double* MPOOLS;
MPOOLS=malloc(nopools*sizeof(double));
for (n=0;n<nopools;n++){MPOOLS[n]=mean_pool(POOLS,n,D.nodays+1,nopools);};


/***********************EDCs start here****************************/



/*EDC LAI related:*/
/*now obsolete*/
/*MAX(Cfol)/clma <LAImax*/
/*double LAImax=10;
if (EDC==1 || DIAG==1){
for (n=0;n<nodays;n++){
if (LAI[n]>LAImax){EDC=0;EDCD->PASSFAIL[5-1]=0;}}}*/
/*MLAI=mean_pointer_vector(LAI,nodays);*/



/*EDC no 9*/
/*0.2*Cf < Cr < 5*Cf*/
/*Cfoliar : Croot = 5:1 or 1:5*/
if ((EDC==1 || DIAG==1) && (MPOOLS[1]>MPOOLS[2]*5 || MPOOLS[1]*5<MPOOLS[2])){EDC=0;EDCD->PASSFAIL[8]=0;}// printf("failed no. 8");}


/*EDC 10
 *G10 = 2 - growth factor in 10 years
 *for any further constraint please incorporate in likelyhood*/
int noyears=floor(nodays*deltat/365);
double G=0.1;
int y;
int poolid[2] = {3,5};
int pp;

//printf("No days %i\n", nodays);
//printf("Deltat %i\n", deltat);
//printf("No years %i\n", noyears);
/*mean annual pool array*/
double *MAPOOL=malloc(noyears*sizeof(double));

//only for wood and som
for (n=0;n<2;n++){
/*Rapid POOL Growth: cannot increase by more than factor of Gn over N yrs*/
/*note - this is GROWTH ONLY!!!*/

pp = poolid[n];

// step 1 - deriving mean annual pools
for (y=0;y<noyears;y++){MAPOOL[y]=mean_annual_pool(POOLS,y,pp,nopools,deltat);} //replaced n by pp

// step 2 - checking pool growth

//if ((EDC==1 || DIAG==1) && (MAPOOL[noyears-1]/MAPOOL[0]>(1+G*(double)noyears))){EDC=0;EDCD->PASSFAIL[9]=0;}

//JFE this pools growth factor G corrected to allow 10% change every 10 years - only for wood and som
if ((EDC==1 || DIAG==1) && (fabs(log(MAPOOL[noyears-1]/MAPOOL[0]))>log(pow(1+G,(double)noyears/10.)))){EDC=0;EDCD->PASSFAIL[9]=0;}

}




/*EDC no 11*/
/*POOL EXPONENTIAL DECAY*/
/*Performed on all seven pools*/

// exporting the decay coefficient here
/*double C;
for (n=0;n<nopools;n++)
{
    if (EDC==1 || DIAG==1) 
    {
        C=expdecay2(POOLS,n,nodays+1,deltat,nopools);
        if (fabs(-log(2)/C)<365.25*noyears && C<0 && finite(C)){EDC=0;}
        if (EDC==0){EDCD->PASSFAIL[11-1]=0;}// printf("failed no. 10");}
    }
}
*/





//double const fauto=pars[1];
//double const ffol=(1-fauto)*pars[2];
//double const flab=(1-fauto-ffol)*pars[12];
//double const froot=(1-fauto-ffol-flab)*pars[3];
//double const fwood=1-fauto-ffol-flab-froot;

/*fraction of GPP to som and lit under equilibrium conditions*/
//double const fsom=fwood+(froot+flab+ffol)*pars[0]/(pars[0]+pars[7]);
//double const flit=(froot+flab+ffol);



/*SOM attractor - must be within a factor of 2 from Csom0*/
/*half the TR is balanced by x2 Csom*/

/*equilibrium factor (in comparison to C_initial)*/
double EQF=2.0;// JFE - put 2. following chat with AAB
double meangpp=0;
for (n=0;n<nodays;n++){meangpp+=FLUXES[n*nofluxes]/(double)nodays;}




/*EDC 12 - SOM steady-state attractor proximity attractor proximity
if ((EDC==1 || DIAG==1) && meangpp*fsom/(pars[8]*exp(pars[9]*meantemp))>pars[22]*EQF){EDC=0;EDCD->PASSFAIL[12-1]=0;}
if ((EDC==1 || DIAG==1) && meangpp*fsom/(pars[8]*exp(pars[9]*meantemp))<pars[22]/EQF){EDC=0;EDCD->PASSFAIL[12-1]=0;}

EDC 13 - LIT  -steady-state attractor proximity
if ((EDC==1 || DIAG==1) && meangpp*flit/(pars[7]*exp(pars[9]*meantemp))>pars[21]*EQF){EDC=0;EDCD->PASSFAIL[13-1]=0;}
if ((EDC==1 || DIAG==1) && meangpp*flit/(pars[7]*exp(pars[9]*meantemp))<pars[21]/EQF){EDC=0;EDCD->PASSFAIL[13-1]=0;}

EDC 14 - WOO  -steady-state attractor proximity
if ((EDC==1 || DIAG==1) && meangpp*fwood/(pars[5])>pars[20]*EQF){EDC=0;EDCD->PASSFAIL[14-1]=0;}
if ((EDC==1 || DIAG==1) && meangpp*fwood/(pars[5])<pars[20]/EQF){EDC=0;EDCD->PASSFAIL[14-1]=0;}

EDC 15 - ROO  -steady-state attractor proximity
if ((EDC==1 || DIAG==1) && meangpp*froot/(pars[6])>pars[19]*EQF){EDC=0;EDCD->PASSFAIL[15-1]=0;}
if ((EDC==1 || DIAG==1) && meangpp*froot/(pars[6])<pars[19]/EQF){EDC=0;EDCD->PASSFAIL[15-1]=0;}
*/

/*Total fluxes*/
double FT[nofluxes];
int f=0;
for (f=0;f<nofluxes;f++){FT[f]=0;for (n=0;n<nodays;n++){FT[f]+=FLUXES[n*nofluxes+f];}}

// replaced expdecay steady-state proximity at t = 0, following PNAS paper
double Sproxjanmean, Sproxjan;// = malloc(1+nodays/12);
int yearid; // iterator
//printf("%d ",noyears);

double Fin[nopools];
double Fout[nopools];

/*labile*/
Fin[0]=FT[4];
Fout[0]=FT[7]+FT[17]+FT[23];
/*foliar*/
Fin[1]=FT[3]+FT[7];
Fout[1]=FT[9]+FT[18]+FT[24];
/*root*/
Fin[2]=FT[5];
Fout[2]=FT[11]+FT[19]+FT[25];
/*wood*/
Fin[3]=FT[6];
Fout[3]=FT[10]+FT[20]+FT[26];
/*litter*/
Fin[4]=FT[9]+FT[11]+FT[23]+FT[24]+FT[25];
Fout[4]=FT[12]+FT[14]+FT[21]+FT[27];
/*som*/
Fin[5]=FT[10]+FT[14]+FT[26]+FT[27];
Fout[5]=FT[13]+FT[22];


/*temporary print switch*/
int psw=0;
int nostep_peryear; // get number of time steps per year
double fin_fout_lim; // get the limit of Fin/Fout ratio // JFE added 16/11/16
if (deltat == 1) {nostep_peryear = 365;fin_fout_lim = 0.10;} else {nostep_peryear = 12; fin_fout_lim = 0.05;}
//printf("nostep_peryear %i ", nostep_peryear);
//printf("%d",nodays);

int poolid2[2] = {3,5}; // define pools to check
int pid=0;

for (pid=0;pid<6;pid++)
{if ((EDC==1 || DIAG==1) && (fabs(log(Fin[pp]/Fout[pp]))>log(EQF))) {EDC=0;EDCD->PASSFAIL[12-1]=0;}}

// in GSI variable, cfol and clab... sprox not applicable

for (pid=0;pid>2;pid++)
{
    pp = poolid2[pid];
    Sproxjanmean=0;
    for (yearid=0;yearid<noyears;yearid++)
    {
        Sproxjanmean = Sproxjanmean+POOLS[(nostep_peryear*yearid*nopools)+pp];
    //printf("step: %i ",nostep_peryear*yearid*nopools);
    }
    Sproxjanmean = Sproxjanmean/noyears;
    Sproxjan = (Fin[pp]/Fout[pp])*(Sproxjanmean/POOLS[pp]); //eq. S5
if ((EDC==1 || DIAG ==1) && (fabs((Fin[pp]/Fout[pp])-Sproxjan) > fin_fout_lim)) {EDC=0;EDCD->PASSFAIL[10]=0;EDCD->PASSFAIL[16+pp]=0;}
}
//EDC 16 - average torfol > torwood -  JFE added 26/2/15
double torfol;
int nn; 

/* commented out because sum of torfol is already stored in FT
for (nn=0;n<nodays;nn++){torfol=torfol+FLUXES[nn*nofluxes+28]);}
torfol=torfol/nodays;
*/
torfol = FT[8]/nodays;

if ((EDC==1 || DIAG ==1) && (pars[5]>torfol))
{EDC=0;EDCD->PASSFAIL[15]=0;}// printf("failed no. 15");}





/***********************EDCs done here****************************/



/*Additional faults can be stored in positions 35-40*/

/*PRIOR RANGES - ALL POOLS MUST CONFORM*/

for (n=0;n<6;n++){if ((EDC==1 || DIAG==1) && ((MPOOLS[n])>PI.parmax[n+17])){EDC=0;EDCD->PASSFAIL[35-1]=0;}}

int PEDC;
/*ensuring minimum of each pool is zero && finite*/
if (EDC==1 || DIAG==1)
{
    n=0;
    while (n<nopools && (EDC==1 || DIAG==1))
    {
        nn=0;PEDC=1;
        while (nn<nodays+1 && (PEDC==1))
        {
            if (POOLS[n+nn*nopools]<0 || isnan(POOLS[36-1])==1) {EDC=0;PEDC=0;EDCD->PASSFAIL[35+n]=0;}
            nn=nn+1;
        }
        n=n+1;
    }
}









/*FREE MEMORY*/
free(MPOOLS);
free(MAPOOL);





/*final check confirming EDC = 1 or 0*/
int Num_EDC=100;


if (DIAG==1){for (n=0;n<Num_EDC;n++){if (EDCD->PASSFAIL[n]==0){EDC=0;}}}
//printf("EDCs: %i",EDC);
/*Returning EDC */
return EDC;


}

