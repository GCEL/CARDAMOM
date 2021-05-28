
###
## Function to compile the model C code
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

dalec_cdea_mcmc_compile<-function (paths,PROJECT) {

    print('Compiling latest DALEC setup')
    #Global DALEC fast code runs
    system(paste("gcc ",paths$cardamom_library,"CARDAMOM_C/projects/DALEC_GLOBAL/DALEC_CDEA_ECDF_OUTPUT.c  -o ",paths$cardamom_library,"CARDAMOM_C/projects/DALEC_GLOBAL/DALEC_CDEA_ECDF_OUTPUT.exe -lm",sep=""))
    #DALEC_CDEA_TEMPLATE
    system(paste("gcc ",paths$cardamom_library,"CARDAMOM_C/projects/DALEC_CDEA_TEMPLATE/DALEC_CDEA_TEMPLATE.c -o ",PROJECT$exepath,PROJECT$name,'.exe -lm',sep=""))
    #DALEC - diagnostic
    #mex([PATHS.CARDAMOM_LIBRARY,'CARDAMOM_C/projects/DALEC_CODE/DALEC_CDEA/DALEC_CDEA_LOOP_EDC.c'],'-outdir','CARDAMOM_LOCAL/LOCALLY_COMPILED/')
    #mex([PATHS.CARDAMOM_LIBRARY,'CARDAMOM_C/projects/DALEC_CODE/DALEC_CDEA/DALEC_CDEA_PARAMETER_RANGES.c'],'-outdir','CARDAMOM_LOCAL/LOCALLY_COMPILED/')
    print('Done compiling local DALEC C code!!')

}
