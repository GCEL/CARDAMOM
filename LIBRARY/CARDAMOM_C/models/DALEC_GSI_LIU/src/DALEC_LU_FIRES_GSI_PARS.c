#pragma once
#include <stdio.h>
#include "../../../general/cardamom_structures.c"


void pars_info(PARAMETER_INFO *PI)
{


PI->npars=26;

// Decomposition rate
PI->parmin[0]=0.00001;
PI->parmax[0]=0.01;

// Fraction of GPP respired
PI->parmin[1]=0.2;
PI->parmax[1]=0.8;

// Fraction of (1-fgpp) to foliage
PI->parmin[2]=0.01;
PI->parmax[2]=1.;

// Fraction of (1-fgpp) to roots
PI->parmin[3]=0.01;
PI->parmax[3]=1.;

// TOR foliar
PI->parmin[4]=0.000001;
PI->parmax[4]=0.1; // changed 0.1 to 0.99 on 03/11/14 

// TOR wood
PI->parmin[5]=0.000025;
PI->parmax[5]=0.001;

// TOR roots

PI->parmin[6]=0.0001;
PI->parmax[6]=0.01;

// TOR litter
PI->parmin[7]=0.0001;
PI->parmax[7]=0.01;

// TOR SOM
PI->parmin[8]=0.0000001;
PI->parmax[8]=0.001;

// Temp factor* = Q10 = 1.2-1.6
PI->parmin[9]=0.018;
PI->parmax[9]=0.08;

// Canopy Efficiency
PI->parmin[10]=5;
PI->parmax[10]=50;

// TORlabile
PI->parmin[11]=0.000001;
PI->parmax[11]=0.1; // changed 0.1 to 0.99 on 03/11/14 

// Fraction to Clab
PI->parmin[12]=0.01;
PI->parmax[12]=1.;

// Tmn,min for GSI
PI->parmin[13]=225.;
PI->parmax[13]=330.;

// Tmn,max for GSI
PI->parmin[14]=225.; 
PI->parmax[14]=330.; 

// Changed to min photoperiod, in sec
PI->parmin[15]=3600.; // 1h
PI->parmax[15]=82800.;// 23h

// LMA - Kattge et al. 2011
PI->parmin[16]=5;
PI->parmax[16]=200;

// INITIAL VALUES DECLARED HERE

// Clab
PI->parmin[17]=20.0;
PI->parmax[17]=2000.0;

// Cfol
PI->parmin[18]=20.0;
PI->parmax[18]=2000.0;

// Croots
PI->parmin[19]=20.0;
PI->parmax[19]=2000.0;

// Cwood
PI->parmin[20]=100.0;
PI->parmax[20]=100000.0;

// Clit
PI->parmin[21]=20.0;
PI->parmax[21]=2000.0;

// Csom
PI->parmin[22]=100.0;
PI->parmax[22]=200000.0;

// Max photoperiod for GSI
PI->parmin[23]=3600.; //6h
PI->parmax[23]=82800.; //18h

// Min VPD for GSI - in Pa
PI->parmin[24]=1.; 
PI->parmax[24]=5500.;

// Max VPD for GSI - in Pa
PI->parmin[25]=1.;
PI->parmax[25]=5500.;
}


