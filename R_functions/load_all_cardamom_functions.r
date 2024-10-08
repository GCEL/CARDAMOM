
###
## Function to load all R scripts comprising the UoE CARDAMOM framework wrapper
###

# This function is based on an original Matlab function development by A. A. Bloom (UoE, now at the Jet Propulsion Laboratory).
# Translation to R and subsequent modifications by T. L Smallman (t.l.smallman@ed.ac.uk, UoE).

load_r_libraries<-function(){
    # load all needed libraries first
    if(!require(chron)) {install.packages("chron")} # Check if installed
    library(chron)                                  # Load installed library
    if(!require(fields)) {install.packages("fields")}
    library(fields)
    if(!require(gplots)) {install.packages("gplots")}
    library(gplots)
    if(!require(ncdf4)) {install.packages("ncdf4")}
    library(ncdf4)
    if(!require(parallel)) {install.packages("parallel")}
    library(parallel)
    #if(!require(rgdal)) {install.packages("rgdal")}
    #library(rgdal)
    if(!require(raster)) {install.packages("raster")}
    library(raster)
    #library(rhdf5) # from bioconductor - not currently used
    if(!require(sp)) {install.packages("sp")}
    library(sp)
    if(!require(zoo)) {install.packages("zoo")}
    library(zoo)
    if(!require(apcluster)) {install.packages("apcluster")}
    library(apcluster)
    if(!require(compiler)) {install.packages("compiler")}
    library(compiler)
    if(!require(RColorBrewer)) {install.packages("RColorBrewer")}
    library(RColorBrewer)
    if(!require(colorspace)) {install.packages("colorspace")}
    library(colorspace)
    if(!require(maps)) {install.packages("maps")}
    library(maps)
    # Set error options to output line number for broken calls
    #options(error = utils::recover)
} # function to load all libraries needed by the system

###
## This script loads all cardamom functions ready for use
###

# load R libraries first
load_r_libraries()

# get the complete list
list_o_functions=list.files("./R_functions/", full.names=T)
#print(list_o_functions)
# remove this file to avoid repetition
loser_list=grepl("load_all_cardamom_functions.r",list_o_functions)
loser_list=which(loser_list)
list_o_functions=list_o_functions[-loser_list]
# avoid specific file
loser_list=grepl("landmask20km.rda",list_o_functions)
loser_list=which(loser_list)
if (length(loser_list) > 0) { list_o_functions=list_o_functions[-loser_list] }
# avoid temp files
loser_list=grepl("~",list_o_functions)
loser_list=which(loser_list == FALSE)
list_o_functions=list_o_functions[loser_list]
# avoid auto saves
loser_list=grepl("rkward_autosave",list_o_functions)
loser_list=which(loser_list == FALSE)
list_o_functions=list_o_functions[loser_list]
# avoid .txt files
loser_list=grepl(".txt",list_o_functions)
loser_list=which(loser_list == FALSE)
list_o_functions=list_o_functions[loser_list]
# avoid .sh
loser_list=grepl(".sh",list_o_functions)
loser_list=which(loser_list == FALSE)
list_o_functions=list_o_functions[loser_list]
# only .r
loser_list=grepl(".r",list_o_functions)
list_o_functions=list_o_functions[loser_list]
# now go throught the list can call the files
if (length(list_o_functions) >=1 ){
for (i in seq(1, length(list_o_functions))) {
    #print(paste("...loading R script = ",list_o_functions[i],sep=""))
    source(list_o_functions[i])
}
}
