### Example of how to determine correct CARDAMOM site based on latitude / longtude

# Set latitude (-90/90) / longitude (-180/180) wanted
lat_wanted = 68.3 ; long_wanted = 28.0


# Load required CARDAMOM framework function
source("./R_functions/function_closest2d.r")

# Load the CARDAMOM infofile for your project
load(".../infofile.RData")
# Load the processed gridded output file
load(paste(PROJECT$results_processedpath,PROJECT$name,"_stock_flux.RData",sep=""))

# Use nearest neighbour approach to find the i,j location within the grid
output = closest2d_2(1,grid_output$lat,grid_output$long,lat_wanted,long_wanted)
long_loc = unlist(output, use.names=FALSE)[1] ; lat_loc = unlist(output, use.names=FALSE)[2]

# Now we need to identify which 'site' this location corresponds to in the overall grid
site_loc = which(grid_output$i_location == long_loc & grid_output$j_location == lat_loc)

# We can now extract and use specific time series outputs for this location
plot(grid_output$harvest_gCm2day[site_loc,4,]*365.25, type="l")