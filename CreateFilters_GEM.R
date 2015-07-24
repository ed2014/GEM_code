# Purpose: read map data and create a filter
# Gets a raster "reference" from outside

# Load libraries
library(data.table)
library(plyr)
library(sp)
library(raster)
library (vioplot)
library(reshape)
library(ggplot2)
#Not sure which libraries are needed
library(rasterVis) 
library(colorspace) 
library(rgdal) 
library(ncdf)

#----------------------------------------------
# --------------- LAND USE layer  ------------
# ---------------------------------------------
# ATTENTION: ensure the shapefile was converted to WGS84 in projection with QGIS
pathFolder <- 'C:/Apsim_dev/Projects/CCII/GIS_layers/LandUse'

setwd(pathFolder) # Assuming the outputs will be saved together with input data

# Reads shapefile and creates a SpatialPolygonsDataFrame class (sp)
LUsf <- readOGR(pathFolder, 'nzlri-land-use-capability(WGS84)') 

# Explore info about the object
ogrInfo(pathFolder, 'nzlri-land-use-capability(WGS84)')
print(proj4string(LUsf))
summary(LUsf)
plot(LUsf$LUC1C)

# Create new field (attribute) for the to reprsent filtered layer (data column)
LUsf$LU_filter <- LUsf$LUC1C
plot(LUsf$LUC1C, main = "land use layer")
plot(LUsf$LU_filter, main = "filter layer")

# Change code to filter land use as needed
LUsf$LU_filter[LUsf$LUC1C %in% 1:4] <- 3    # agregation of all land suitable to arable (at various degrees of suitability)
LUsf$LU_filter[LUsf$LUC1C %in% 5:6] <- 2    # unsuitable land for arable - slight to modelrate limitations for perennials
LUsf$LU_filter[!(LUsf$LUC1C %in% 1:6)] <- 1 # severe limitations for both arable and perennials

# Check filter
extent(LUsf)
plot(LUsf$LU_filter, main = "Aggregated filter")
extent(LUsf)

#check a netCDF used as basis for rasterisation
nc <- open.ncdf("C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\netCDF_Example\\maxt_1971-1980.nc")
nc
nc$dim$lon$vals
nc$dim$lat$vals
nc$dim$time$vals
nc$dim$time$units

# get a raster to use as template
rastBase <- raster("C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\netCDF_Example\\maxt_1971-1980.nc") # direct from netCDF

# check if all is fine with your template
summary(rastBase)
extent(rastBase)
plot(rastBase)
proj4string(rastBase) # transformed with QGIS from original
proj4string(rastBase) # <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0" # change projection

# Convert filter column to a raster layer
LU_rs <- rasterize(LUsf, rastBase, 'LU_filter')

# check if your filter layer is ok
extent(LU_rs)
plot(LU_rs$layer)
plot(LU_rs)

nrow(LU_rs)
ncol(LU_rs)

LU_rs$layer
ncell(LU_rs)
dim(LU_rs)

# Convert raster in matrix format using the same row/col as your base raster
# to be able to reference it in a loop later
mapMatrix <- as.matrix(LU_rs, xy=TRUE, nrow = nrow(rastBase), 
                       ncol = ncol(rastBase), xy=TRUE, centroids = TRUE, na.rm=TRUE)
head(mapMatrix)

# check if right col/row
nrow(mapMatrix)
ncol(mapMatrix)

# Create a list of grid-cells to be considered
# Loop through raster map matrix and do stuff
# This filter can be used to select which met files to create and .sim to run
count <- 0
fileNo1 <- row1 <- col1 <- luCode1 <- NULL
metList <- NULL # FIXME: could be a list but then have to change synthax througout to reference its dimensions
outname <- NULL
for (r in 1:nrow(mapMatrix)) {  
  for (c in 1:ncol(mapMatrix)) {   
    # Stores data if not unavailable
    if(!is.na(mapMatrix[r,c])) {  print(paste(c(r,c,mapMatrix[r,c])))                                                                  
                                  count <- count + 1
                                  row1 <- c(row1, r)
                                  col1 <- c(col1, c)
                                  luCode1 <- c(luCode1, mapMatrix[r,c])
                                  fileNo1 <- c(fileNo1, count)
                                  
                                  # Stores pixel file name in a listlist 
                                  # if land is suitable for arable (index = 3)
                                  if(mapMatrix[r,c] == 3) {  
                                    outName <- paste(r,"_",c,".met", sep = "")
                                    metList <- c(metList, outName) 

                                  }    
    }          
  }  
}

print(paste0(length(metList), " selected grid-cells"))

# creates a DF with "all" pixel ID and filter value
df_allLuFilter <- data.frame(fileNo = fileNo1, row = row1, col = col1, luCode = luCode1) 

# Result directory
resPath <- "C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\gisLayers\\luFilter\\"

# Outputs of filter layer
writeRaster(LU_rs, file=paste0(resPath,"LandUseFilter.tiff"), overwrite=TRUE) # Save filter as a raster map
write.table(df_allLuFilter, file = paste0(resPath,"LandUseFilter_Allcats.txt"), col.names = TRUE, row.names = T) # all categories
save(metList,file=paste0(resPath,"LandUseArableFilter_metList.RData"), ascii = TRUE) # only selected LU category

metList
head(df_allLuFilter)

print(paste0("Selected ", length(metList)," pixels in the arable LU category"))

# Test How many suitable pixels not avaiable as met file exist
load("C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\gisLayers\\luFilter\\LandUseArableFilter_metList.RData", .GlobalEnv)
count <- 0
metFolderTest <- "F:\\ERA-40_Reanalysis\\MetFiles_1971_2000\\"
for (f in 1:length(metList)){ 
if(!file.exists(paste0(metFolderTest, metList[f])))
{count = count+1}
}
print(paste(count, " suitable pixels not found in folder, or ", 
            round(count/length(metList), digits=2), "%", sep = ""))
