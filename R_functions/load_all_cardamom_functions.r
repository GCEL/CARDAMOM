load_r_libraries<-function(){
    # load all needed libraries first
    require(chron)
    require(fields)
    require(gplots)
    require(ncdf4)
    require(parallel)
    require(rgdal)
    require(raster)
    require(rhdf5)
    require(sp)
    require(zoo)
    require(apcluster)
    require(compiler)
    require(RColorBrewer)
    require(colorspace)
} # function to load all libraries needed by the system

###
## This script loads all cardamom functions ready for use
###

# load R libraries first
load_r_libraries()

# get the complete list
list_o_functions=list.files(paste(source_directory,"R_functions/",sep=""), full.names=T)
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
# avoid .csv
loser_list=grepl(".csv",list_o_functions)
loser_list=which(loser_list == FALSE)
list_o_functions=list_o_functions[loser_list]
# avoid global_map directory
loser_list=grepl("global_map",list_o_functions)
loser_list=which(loser_list == FALSE)
list_o_functions=list_o_functions[loser_list]
# only .r
loser_list=grepl(".r",list_o_functions)
list_o_functions=list_o_functions[loser_list]
# now go throught the list can call the files
for (i in seq(1, length(list_o_functions))) {
#    print(paste("...loading R script = ",list_o_functions[i],sep=""))
    source(list_o_functions[i])
}
