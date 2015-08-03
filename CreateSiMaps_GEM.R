# not sure which libs I tried as well
library(ggplot2)
library(plot3D)
library(gridExtra)
library(tidyr)
library(mapproj)
library(knitr)
library(plyr)

# libs used
library(dplyr)
library(raster)
library(rgdal) 
library(sp)

library(rgdal)
library(maptools)
library(rgeos)


# read df
data_df <- read.table("C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\GEM_data\\Si_by_Pixel.txt",header = TRUE)
head(data_df)

# Identify variables
vars <- unique(data_df$thisVar)

# get col names
outNames <- gsub(" ","",colnames(data_df))


# -------------------- Read shape file with NZ regions and borders -------------

# ATTENTION: ensure the shapefile was converted to WGS84in projection with QGIS
pathShapeFile <- '//simplace.net/projects/nz/Gunther_map_example/nz_regions_GIS' 

# Reads shapefile and creates a SpatialPolygonsDataFrame class (sp)
sf2 <- readShapeSpatial(paste(pathShapeFile,'REGC2013_HD_Clipped(WGS84)',sep="/"),proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

# reduces amount of points to speedup plotting
sf2 <- gSimplify(sf2,tol=.01,topologyPreserve = TRUE)




# ---------------------- Rasterise DF of results ------------

par(mfrow=c(2,4))
spg <- list() # dataframe list to be coerced as spatial object
rast <- list() # raster
s <- list()


# Embelish factor names. FIXME: Automate this later
factNames <- c("thisLat", "thisLong", "thisClimZone", "thisVar", "row", "col" ,"Hybrid","Sowing date", "Weather", "Interactions")

stCol <- 7 # FIXME: Find a way to select the result columns autimatically

# Rasteriseloop 2 impact variables and four result variables

for (v in 1:length(vars)) {
  for(o in stCol:length(outNames)){
    
   df_temp <- data_df %>% subset (thisVar == vars[v])
    head(df_temp)
    
    # For aggregate
    spg[[o]] = data.frame(df_temp$thisLong,df_temp$thisLat, df_temp[[o]])
    coordinates(spg[[o]]) = ~ df_temp.thisLong + df_temp.thisLat # Attention to variable names
    
    # Common to both methods (grid, raterise and project)
    gridded(spg[[o]]) <- TRUE
    rast[[o]] <- raster(spg[[o]])
    proj4string(rast[[o]]) <- CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0")
    
    # stack them
    if(o == 1) {
      s <- NULL
      s <- stack(rast[[o]])} 
    else {
      s[[o]] <- rast[[o]]
    }
    
    thisOutName <- paste0(factNames[o],"_",vars[v])
    
    # TODO: Add the country and region borders from a NZ shape file
    # Gunther's help needed here

#   plotting borders first has the disadvantage that lines are 
#   hidden by raster points.
    plot(sf2, bg="transparent", xlim=c(167.2,178.55),main = thisOutName)
    plot(rast[[o]], add=TRUE)
    
#   plotting raster first leads to a shift of the borders 
    
#    plot(rast[[o]], main = thisOutName, axes=FALSE,asp=1, ext = c(160,180,37,43))
#    plot(sf2, bg="transparent", add=TRUE)
    
    # save raster as
    writeRaster(rast[[o]], file=paste0("Si_",thisOutName,".tiff"), overwrite=TRUE) 
    
  }
  
} # end loop vars