/* 
12 Aug 2015 - JFE
This file contains the functions specific to the MHMCMC:
    - step: perform a sampling step
    - adapt_step_size: adapt the size of a sampling step
    - write_results: append parameter sets and saves step size
*/



void step(double *pars0, double *pars, PARAMETER_INFO PI)
{
/*FIXEDPARS*/
/*ones and zeros depending on whether parameters are kept fixed*/
    int n,fp;
    double npar=0,rn=0;

    // sampling parameters
    for (n=0;n<PI.npars;n++)
    {      
	    fp=0;       
	    while (fp==0)
        {
	    // normalising parameter
		    rn=randn();
		    npar=par2nor(pars0[n],PI.parmin[n],PI.parmax[n])+rn*PI.stepsize[n]*(1-PI.parfix[n]); // this way, if parfix, the value is not changed
		    // repeat until parameter value is within 0-1 limits
		    if (npar>0 && npar<1)
            {
                fp=1;pars[n]=nor2par(npar,PI.parmin[n],PI.parmax[n]);
            }
        }
    }
}





void adapt_step_size(double *PARSALL, PARAMETER_INFO PI, COUNTERS N, MCMC_OPTIONS MCO)
{

    // first adjusting all step sizes
    int n;
    double fac=2;
    double adapfac=1.5;
    double minstepsize=10000/N.ITER;if (minstepsize>0.01){minstepsize=0.01;}

    // calculate the local acceptance rate: fraction of accepted runs since last change in step size
    N.ACCRATE=(double)N.ACCLOC/(double)MCO.nADAPT;



    // if the local acceptance rate is too low, narrow the distribution
    if (N.ACCLOC>0 && N.ACCRATE<0.23)
    {
        for (n=0;n<PI.npars;n++){PI.stepsize[n]=PI.stepsize[n]/adapfac;}
    }
    // if the acceptance rate is too high, widen the distribution
    else if(N.ACCRATE>0.44)
    {
        for (n=0;n<PI.npars;n++){PI.stepsize[n]=PI.stepsize[n]*adapfac;}
    }

    // dimension specific adjustments
    if (N.ACCLOC>3)
    {

	    /*need to remove this!!*/
	    double norparstd;
	    double *norparvec=calloc(N.ACCLOC,sizeof(double));
	    int p;
	    /*std of each parameter (for all recently N.ACCLOCepted parameters)*/
	    for (p=0;p<PI.npars;p++)
        {
		    for (n=0;n<N.ACCLOC;n++){norparvec[n]=par2nor(PARSALL[PI.npars*n+p],PI.parmin[p],PI.parmax[p]);}
        
		    /*transforming parameters and storing in pointer*/			
		    /*calculating variability*/
		    norparstd=std(norparvec,N.ACCLOC);
		    /*reducing if the stepsize is too small*/
     		if (PI.stepsize[p]<norparstd/fac && N.ACCRATE<0.23)
		    /*disabled step-dim-adapt*/
		    {PI.stepsize[p]=PI.stepsize[p]*sqrt(adapfac);}
        }

	    free(norparvec);
    }


    // overall check
    for (n=0;n<PI.npars;n++)
    {
        if (PI.stepsize[n]>1){PI.stepsize[n]=PI.stepsize[n]/adapfac;}
        if (PI.stepsize[n]<minstepsize){PI.stepsize[n]=PI.stepsize[n]*adapfac;}
        if (PI.stepsize[n]<minstepsize){PI.stepsize[n]=minstepsize;}
    }
}


void write_results(double *PARS, double PROB, PARAMETER_INFO PI, MCMC_OPTIONS MCO)
{

    int n;



    FILE *fileout=fopen(MCO.parfile,"ab");
    FILE *filestep=fopen(MCO.stepfile,"wb");
    for (n=0;n<PI.npars;n++){
            fwrite(&PARS[n],1,sizeof(double),fileout);
           	fwrite(&PI.stepsize[n],1,sizeof(double),filestep);}
        
    /*writing likelyhood*/
    fwrite(&PROB,1,sizeof(double),fileout);
    fclose(fileout);
    fclose(filestep);


    //return 0;


}


