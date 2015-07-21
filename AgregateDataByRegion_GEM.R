# Aggregate by region

#load labraries used
library(rgdal) 
library(raster)
library(sp)
library(plyr)

#----------------------------------------------
# --------------- Region layer  ------------
# ---------------------------------------------
# ATTENTION: ensure the shapefile was converted to WGS84in projection with QGIS
pathFolder <- 'C:/apsim_dev/Projects/CCII/GIS_layers/NZ_Regions_GIS'

setwd(pathFolder) # Assuming the outputs will be saved together with input data

# Reads shapefile and creates a SpatialPolygonsDataFrame class (sp)
sf <- readOGR(pathFolder, 'REGC2013_HD_Clipped(WGS84)') 

# Explore info about the object
ogrInfo(pathFolder, 'REGC2013_HD_Clipped(WGS84)')
print(proj4string(sf))
summary(sf)
plot(sf$REGC2013)
#Check attributes
sf[[1]]
sf[[2]]

# Create DF with region names !!!! FIXME: How did you get that info?
# Attention: this does not match the attribute table you see in GIS program!!!!!
df_reg_names <- data.frame(region = c(1,2,3,4,5,6,7,8,9,12,13,14,15, 16,17,18, 99),  
                           regNames = c("Northland", "Auckland", "Waikato", 
                                       "BOP", "Gisborne", "HawkesBay", "Taranaki",
                                       "Manuatu", "Wellington", "WestCoast", 
                                       "Canterbury", "Otago", "Southland", 
                                       "Tasman", "Nelson", "Malborough", 
                                       "Others"))
head(df_reg_names)

# get a raster to use as template
rastBase <- raster("C:\\Apsim_dev\\Projects\\CCII\\RawData\\HadCM3run2\\A2\\maxt_1971-1980.nc") # direct from netCDF
# check if all is fine with your template
summary(rastBase)
extent(rastBase)
plot(rastBase)
proj4string(rastBase)

# create a raster base as in netCDF file
proj4string(rastBase) # = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" 

# Convert a given column of variable in shape file to a raster layer 
# (Attention: find out name of shape file attribute you want!)
rs <- rasterize(sf, rastBase, 'REGC2013')

# check if your filter layer is ok
extent(rs)
plot(rs$layer)
plot(rs)
writeRaster(rs, file=paste0("C:\\apsim_dev\\Projects\\CCII\\GIS_layers\\NZ_Regions_GIS\\regions.tif"), 
            overwrite=TRUE) 

nrow(rs)
ncol(rs)

rs$layer
ncell(rs)
dim(rs)

# Convert raster into matrix format, using the same row/col as your base raster
# to be able to reference it in a loop later
mapMatrix <- as.matrix(rs, xy=TRUE, nrow = nrow(rastBase), 
                      ncol = ncol(rastBase), xy=TRUE, centroids = TRUE, na.rm=TRUE)

# check if right col/row
nrow(mapMatrix)
ncol(mapMatrix)

# Create a list of grid-cells to be considered
# Loop through raster map matrix and do stuff
# This filter can be used to select which met files to create and .sim to run
count <- 0
fileNo1 <- row1 <- col1 <- var1 <- NULL
metList <- NULL # FIXME: could be a list but then have to change synthax througout to reference its dimensions
outname <- NULL

# create vectors with all coordinates and region number for available pixels
for (r in 1:nrow(mapMatrix)) {  
  for (c in 1:ncol(mapMatrix)) {   
    # Stores data if not unavailable
    if(!is.na(mapMatrix[r,c])) {  print(paste(c(r,c,mapMatrix[r,c])))                                                                  
                                  count = count + 1
                                  row1 = c(row1, r)
                                  col1 = c(col1, c)
                                  var1 = c(var1, mapMatrix[r,c])
                                  fileNo1 = c(fileNo1, count)
                                  
                                  # Stores pixel file name in a listlist 
                                  # if land is suitable for arable (index = 3)
                                  #if(mapMatrix[r,c] == 3) {  
                                  #  outName = paste(r,"_",c,".met", sep = "")
                                  #  metList = c(metList, outName) 
                                    
                                  }    
            }          
  }  


# creates a DF with "all" pixel ID and filter value
df_region <- data.frame(record = fileNo1, 
                        row = row1, col = col1, region = var1) 

df_region <- merge(df_region, df_reg_names, by = "region" )
summary(df_region)

# save the list of ragion per pixel 
write.csv(df_region, file = paste0("C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\GEM_data\\",
                                   "RegionPerPixel.csv"))

