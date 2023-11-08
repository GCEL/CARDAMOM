#pragma once
#include <stdio.h>
#include "../../../general/cardamom_structures.c"

void pars_info(PARAMETER_INFO *PI)
{

//PI->npars=23;

/*contains 6 fields with min max log for par and par*/


/*Decomposition rate*/
PI->parmin[0]=0.00001;
PI->parmax[0]=0.01;



/*Fraction of GPP respired*/
PI->parmin[1]=0.3;
PI->parmax[1]=0.7;


/*Fraction of (1-fgpp) to foliage*/
PI->parmin[2]=0.01;
PI->parmax[2]=1.;

/*Fraction of (1-fgpp) to roots*/
PI->parmin[3]=0.01;
PI->parmax[3]=1.;


/*Leaf Lifespan*/
/*Wright et al. 2004*/
PI->parmin[4]=1.001;
PI->parmax[4]=8;


/*TOR wood* - 1% loss per year value*/
PI->parmin[5]=0.000025;
PI->parmax[5]=0.001;


/*TOR roots*/
PI->parmin[6]=0.0001;
PI->parmax[6]=0.01;

/*TOR litter*/
PI->parmin[7]=0.0001;
PI->parmax[7]=0.01;

/*TOR SOM*/
PI->parmin[8]=0.0000001;
PI->parmax[8]=0.001;


/*Temp factor* = Q10 = 1.2-1.6*/
PI->parmin[9]=0.018;
PI->parmax[9]=0.08;


/*Canopy Efficiency*/
/*set to parmin=1 for FLUXCOM only*/
/*e.g. for wetlands etc.*/
PI->parmin[10]=5;
PI->parmax[10]=50;


/*Bday*/
PI->parmin[11]=365.25;
PI->parmax[11]=365.25*4;


/*Fraction to Clab*/
PI->parmin[12]=0.01;
PI->parmax[12]=1.;


/*Clab Release period*/
PI->parmin[13]=10;
PI->parmax[13]=100;


/*Fday*/
PI->parmin[14]=365.25;
PI->parmax[14]=365.25*4;


/*Leaf fall period*/
PI->parmin[15]=20;
PI->parmax[15]=150;

/*LMA*/
/*Kattge et al. 2011*/
PI->parmin[16]=5;
PI->parmax[16]=200;




/*INITIAL VALUES DECLARED HERE*/



/*C labile*/
PI->parmin[17]=20.0;
PI->parmax[17]=2000.0;

/*C foliar*/
PI->parmin[18]=20.0;
PI->parmax[18]=2000.0;

/*C roots*/
PI->parmin[19]=20.0;
PI->parmax[19]=2000.0;

/*C_wood*/
PI->parmin[20]=100.0;
PI->parmax[20]=100000.0;

/*C litter*/
PI->parmin[21]=20.0;
PI->parmax[21]=2000.0;

/*C_som*/
PI->parmin[22]=100.0;
PI->parmax[22]=200000.0;

// added 8 May 2017 - CWD relevant parameters
//C_cwd
PI->parmin[23]=1.0;
PI->parmax[23]=10000.0;

//TORcwd
PI->parmin[24]=0.0001;
PI->parmax[24]=0.01;

}


