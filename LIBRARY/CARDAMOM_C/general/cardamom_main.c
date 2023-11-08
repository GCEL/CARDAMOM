/*
7 Aug 2015 - JFE
All routines work, need to wrok on output file generation and restart...

JFE - coding started 3 Aug 2015
This is the where the main CARDAMOM program is called. The code has been
cleaned up to match with the organisation of the Fortran version and allows
easier changes.
*/

#include <stdio.h>

#include "cardamom_structures.c"
#include "cardamom_io.c"
//#include "../models/DALEC_GSI/likelihood/MODEL_LIKELIHOOD.c"
//#include "../method/MHMCMC/MHMCMC.c"

int main(int argc, char *argv[]){

/* input arguments
 * 1. file in
 * 2. file out
 * 3. number of solutions requested
 * 4. print-to-screen frequency
 * 5. write-to-file frequency */

//check that the number of argument is correct
if (argc != 6){
    printf("Error: %d argument(s) input instead of 5 required\n",argc-1);
    return 1;
}
else{
    printf("Correct number of arguments found\n");
}


// read arguments
char* infile = argv[1];
char* outfile = argv[2];
// store as integers
int solutions_wanted = atoi(argv[3]);
int freq_print = atoi(argv[4]);
int freq_write = atoi(argv[5]);

// print summary of options
printf("%d solutions wanted\n", solutions_wanted);
printf("Print status every %d accepted solutions\n", freq_print);
printf("Write every %d accepted solutions into file \"%s\"\n", freq_write, outfile);

// define data structure
DATA DATA;
// define parameter structure
PARAMETER_INFO PI;
// define MCMC_OPTIONS structure
MCMC_OPTIONS MCOPT;
// define MCMC output options
MCMC_OUTPUT MCOUT;

// read parameters and data - in cardamom_io.c 
read_pari_data(&PI, &DATA, &MCOUT, infile);

printf("No. pars: %d\n", DATA.nopars);
printf("No. pars: %d\n", PI.npars);

/* Various prints to check it works
printf("Model ID: %d\n", DATA.ID);
printf("Latitude: %f\n", DATA.LAT);
printf("No. days: %d\n", DATA.nodays);
printf("No. met: %d\n", DATA.nomet);
printf("No. obs: %d\n", DATA.noobs);
printf("EDC: %d\n", DATA.EDC);


printf("METDATA: %f\n",DATA.MET[DATA.nomet*(DATA.nodays-1)]);
printf("LAI: %f\n",DATA.LAI[DATA.nodays-1]);
printf("nLAI: %d\n", DATA.nlai);
printf("Meantemp: %f\n",DATA.meantemp);
printf("Meanrad: %f\n",DATA.meanrad);
*/

// store MCMC options in the designated structure
read_options(&MCOPT, solutions_wanted, freq_print, freq_write, outfile);

printf("MCO.parfile main %s\n", MCOPT.parfile);



// check for existing output files to allow restarting an overtime run
printf("Now checking whether this is a restart run ...\n");
int restart_flag = check_for_existing_output_files(MCOPT.parfile, MCOPT.stepfile, &PI);

if (restart_flag < 0) {return 2;}

int parid = 0;
if (restart_flag == 0)
{
    // search for a starting point
    printf("Beginning search for initial parameter values ... \n");
    find_edc_initial_values(&PI, DATA);

    // perform the main MCMC
    printf("Established starting point ... starting the MHMCMC\n");
    for (parid = 0; parid < PI.npars; parid++) {PI.stepsize[parid] = 0.001;}
}
else if (restart_flag > 0)
{
    for (parid = 0; parid < PI.npars; parid++) {printf("%d %d %f %f\n",restart_flag, parid, PI.stepsize[parid], PI.parini[parid]);}
    MCOPT.nOUT = MCOPT.nOUT - restart_flag*MCOPT.nWRITE;
    printf("Restarting the MCMC for an additional %d accepted solutions.\n", MCOPT.nOUT);
}

//open_output_files(MCOPT.parfile, MCOPT.stepfile);
if (MCOPT.nOUT > 0)
{
    MHMCMC(model_likelihood, DATA, PI, MCOPT, &MCOUT);
}
else
{
    printf("Output file already contains at least the number of requested parameter sets\n");
}

return 0;

}
