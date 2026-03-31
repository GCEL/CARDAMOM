library(compiler)
library(stringr)
source("~/RDM/CARDAMOM/R_functions/read_binary_file_format.r")




binary_list <- list.files("~/RDM/CARD_output/DALEC.A1.C1.D2.F2.H2.P1.#_MHMCMC/Fin_strat_mineral_01degree_test/DATA", full.names=T)


lai_check <- c()
Cwood_check <- c()
LL_check <- c()
Csom_check <- c()

for(i in 1:length(binary_list)){
	temp <- read_binary_file_format(binary_list[i])
	lai_check[i] <- length(which(temp$obs[,3] > -9999)) > 0
	Cwood_check[i] <- length(which(temp$obs[,13] > -9999)) > 0
	LL_check[i] <- temp$otherpriors[7] > 0
	Csom_check[i] <- temp$parpriors[23] > 0
	rm(temp)
}

all_checks <- data.frame(lai_check = lai_check, Cwood_check = Cwood_check, LL_check = LL_check, Csom_check = Csom_check)
all_checks$sum <- apply(all_checks, 1, sum)

# 



# for extracting site ID from bin files (to be used to remove PARS files with missing obs)

# LL

LL_na_index <- which(all_checks$LL_check == FALSE)

LL_na_files <- binary_list[LL_na_index]

LL_na_siteIDs <- c()
 
for(i in 1:length(LL_na_files)){
	LL_na_siteIDs[i] <- str_sub(LL_na_files[i],start=nchar(LL_na_files[i])-8,end=nchar(LL_na_files[i])-4)
}

# Cwood

Cwood_na_index <- which(all_checks$Cwood_check == FALSE)

Cwood_na_files <- binary_list[Cwood_na_index]

Cwood_na_siteIDs <- c()
 
for(i in 1:length(Cwood_na_files)){
	Cwood_na_siteIDs[i] <- str_sub(Cwood_na_files[i],start=nchar(Cwood_na_files[i])-8,end=nchar(Cwood_na_files[i])-4)
}


# Csom

Csom_na_index <- which(all_checks$Csom_check == FALSE)

Csom_na_files <- binary_list[Csom_na_index]

Csom_na_siteIDs <- c()
 
for(i in 1:length(Csom_na_files)){
	Csom_na_siteIDs[i] <- str_sub(Csom_na_files[i],start=nchar(Csom_na_files[i])-8,end=nchar(Csom_na_files[i])-4)
}


# put all siteIDs into a single variable and run unique() and print to save to a text file that will be passed
# to a terminal command to delete from the RESULTS/ directory. see obsidian file "Forth eddie" for terminal command
all_site_IDs <- c(LL_na_siteIDs, Cwood_na_siteIDs, Csom_na_siteIDs)
unique_siteIDs <- unique(all_site_IDs)

write(unique_siteIDs, file="~/RDM/projects/MS_NFIs/throwaway_peat_pixels.txt") # TG check peat vs mineral file naming

#############################################################################
### completely different script below ###############
##################################################



# to get lat long for vectorized gridded outputs (time series)

vector_lat <- c()
vector_long <- c()

for(i in 1:dim(grid_output$lai_m2m2)[1]){
	j_loc <- grid_output$j_location[i] 
	i_loc <- grid_output$i_location[i]
	vector_lat[i] <- grid_output$lat[1,j_loc]
	vector_long[i] <- grid_output$long[i_loc,1]
	rm(i_loc, j_loc)}


# vectorize time series variable of interest (taking median value)

variable_all <- grid_output$variable[,4,]
variable_all <- as.vector(variable_all)

# then make a layer ID for the raster (time step)

month_step <- rep(1:168,each=6728)

# then repeate the lat and long vector for the amount of time steps (here 168)

lat_rep <- rep(vector_lat, 168)
long_rep <- rep(vector_long, 168)

# then put all into a dataframe

variable.df <- data.frame(x=long_rep, y=lat_rep,l=month_step,z=variable_all)

# then delete rows with NA
variable.df <- na.omit(variable.df)

# create raster
raster_variable <- rast(variable.df,crs=("+init=epsg:4326"),type="xylz")

# then make plot titles by making a month and a year vector

month_names <- c("Jan", "Feb", "March","April","May","June","July","Aug", "Sep", "Oct","Nov","Dec")
month_names <- rep(month_names, 14)

year_names <- rep(2009:2022,each=12)

# and paste together 
month_year <- paste0(month_names,"-",year_names)

# animate
animate(raster_variable,main=month_year)

# animate function

animate_timeseries <- function(variable){
	variable <- variable[,4,]
	variable <- as.vector(variable)
	temp.df <- data.frame(x=long_rep,y=lat_rep,l=month_step,z=variable)
	temp.df <- na.omit(temp.df)
	temp.df <- temp.df %>% filter(l>=121)
	temp_rast <- rast(temp.df,crs=("+init=epsg:4326"),type="xylz")
	rm(temp.df)
	animate(temp_rast,main=month_year[121:168],pause=1)
	rm(temp_rast)
	dev.off()
	gc()
	}

