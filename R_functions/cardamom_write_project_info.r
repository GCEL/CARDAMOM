
###
## Function creates a text file with all project details
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

cardamom_write_project_info <- function (PROJECT) {

    # create file name and output basic information
    file=paste(PROJECT$localpath,"/",PROJECT$name,"_run_DETAILS.txt", sep="")
    write(paste("Name: ",PROJECT$name,sep=""),sep=",",ncolumns = 1,file=file,append="F")
    write(paste("Number of sites: ",PROJECT$nosites,sep=""),sep=",",ncolumns = 1,file=file,append="T")
    write(paste("Number of MCMC chains: ",PROJECT$nochains,sep=""),sep=",",ncolumns = 1,file=file,append="T")
    write(paste("Number of accepted samples: ",PROJECT$nsamples,sep=""),sep=",",ncolumns = 1,file=file,append="T")
    write(paste("Number of accepted subsamples: ",PROJECT$nsubsamples,sep=""),sep=",",ncolumns = 1,file=file,append="T")
    write(paste("Description: ",PROJECT$description,sep=""),sep=",",ncolumns = 1,file=file,append="T")
    write(paste("Date: ",PROJECT$date,sep=""),sep=",",ncolumns = 1,file=file,append="T")

    # load and output parameter information
    pdetails = parameter_details(PROJECT$model$name,PROJECT$parameter_type,PROJECT$ctessel_pft)
    ranges = parameter_ranges(PROJECT$model$name,PROJECT$parameter_type,PROJECT$ctessel_pft)
    Pmin = ranges$allparsmin ; Pmax = ranges$allparsmax
    for (n in seq(1, length(Pmin))) {
         write(paste(pdetails$parameter_names[n],"min and max values ",Pmin[n]," - ",Pmax[n],sep=""),sep=",",ncolumns = 1,file=file,append="T")
    }
}
