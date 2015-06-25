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
df_reg_names <- data.frame(region = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15, 16, 99),  
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
summary(df_region)

# save the list of ragion per pixel 
write.csv(df_region, file = paste0("C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\GEM_data\\",
                                   "RegionPerPixel.csv"))

#-----------------------------------------------
# Load dataframes + merge row/col to create a Region column + analyse by region
# ----------------------------------------------

# load result dataframe foe BASE climate
setwd("C:\\Apsim_dev\\Projects\\CCII\\outFiles\\HadCM3run2\\A2-Gen\\Base")
load("DATA.Rda", .GlobalEnv)
df_base = DATA

head(df_base)
summary(df_base)

# load DF for Future Climate
setwd("C:\\Apsim_dev\\Projects\\CCII\\outFiles\\HadCM3run2\\A2-Gen\\Fut1")
load("DATA.Rda", .GlobalEnv)
df_fut1 = DATA

# Transfer region data as new variable in the dfs
df_base = merge(df_base, df_region, by = c("col", "row"))
df_fut1 = merge(df_fut1, df_region, by = c("col", "row"))

# Get rid of catch crop data
df_base = df_base[ which(df_base$CurrentCrop!='wheat_exceed'),]
df_fut1 = df_fut1[ which(df_fut1$CurrentCrop!='wheat_exceed'),]

head(df_fut1)
summary(df_fut1)


# Inspect data
boxplot(TotalBiomass ~ region, data = df_base, horizontal = TRUE)
abline(v = 0)

# Aggregate data by row/col (pool by year and crop type etc)

# Base climate
df_pixel_base = ddply(df_base, c('row', 'col'), summarise, 
                      region = mean(region), lat = mean(thisLat), long = mean(thisLong),Bio_mean=mean(TotalBiomass),Rad_mean = mean(IntRadSum),
                      Rue_mean = mean(RUEtop), Cycle_mean = mean(GrowthLength), Bio_med=median(TotalBiomass), Bio_std=sd(TotalBiomass))
head(df_pixel_base)

# Future climate
df_pixel_fut1 = ddply(df_fut1, c('row', 'col'), summarise, 
                      region = mean(region), lat = mean(thisLat), long = mean(thisLong),Bio_mean=mean(TotalBiomass),Rad_mean = mean(IntRadSum),
                      Rue_mean = mean(RUEtop), Cycle_mean = mean(GrowthLength), Bio_med=median(TotalBiomass), Bio_std=sd(TotalBiomass))
head(df_pixel_fut1)

# Double value of common pixels for subtarction of dfs FIXME: this is a quick solution to be improved!
df_pixel_fut1$row = df_pixel_fut1$row * 2
df_pixel_fut1$col = df_pixel_fut1$col * 2 
df_pixel_fut1$region = df_pixel_fut1$region * 2
df_pixel_fut1$lat = df_pixel_fut1$lat * 2
df_pixel_fut1$long = df_pixel_fut1$long * 2

# Subtracts Future from Base to get differences (Deltas)
df_diff = df_pixel_fut1 - df_pixel_base
head(df_diff)

df_diff = merge(df_diff, df_reg_names, by = ("region"))


boxplot(Bio_mean ~ regNames, data = df_diff, horizontal = TRUE, las = 1)
abline(v = 0)

boxplot(Cycle_mean ~ regNames, data = df_diff, horizontal = TRUE)
abline(v = 0)


df_merge = merge(df_pixel_base,df_pixel_fut1, by = c("row", "col"))
head(df_merge)
df_diff = merge(df_diff, by = c("row", "col"))

head(df_diff)

# Subsets to remove wheat (catch crop)
#summary(DATA_sub$CropType)
#df_data_region = subset(df_data, CurrentCrop != "wheat_exceed")
df_data_region = df_data_region[ which(df_data_region$CurrentCrop!='wheat_exceed'),]

summary(df_data_region)




boxplot(df_data_region$TotalBiomass ~ df_data_region$region, las = 3)


par(mfrow=c(1,1))
boxplot(TotalBiomass ~ region, data = df_data_region, las = 3)


xmin = min(df_data_region$TotalBiomass)
xmax = max(df_data_region$TotalBiomass)  
par(mfrow=c(4,4))

for (r in 1:max(df_data_region$region)) {
  lab = paste0("Region ", r)
  df_sub = subset(df_data_region, df_data_region$region == r)
  hist(df_sub$TotalBiomass, main = lab , xlim=c(xmin, xmax))
}


















# TODO: load a given raster and return tables/graphs by region (not working at all)

# read raster
rs_out = raster("C:\\apsim_dev\\Projects\\CCII\\outFiles\\HadCM3run2\\A2-Gen\\Base\\maize_silage_TotalBiomass.tif", package="raster")
plot(rs_out)

#read in polygon shape file
regShp <- shapefile("C:/apsim_dev/Projects/CCII/GIS_layers/NZ_Regions_GIS/REGC2013_HD_Clipped(WGS84)") 

#calc mean depth per polygon feature
#unweighted - only assigns grid to district if centroid is in that district
regShp$var1 <- extract(rs_out, regShp, fun = mean, na.rm = TRUE, weights = FALSE)

#plot depth values 
spplot(regShp[,'var1'])

summary(regShp)

# loop in raster to get row col and value (FIXME: Ultra slow!!!!!)
for (r in 1:nrow(rs_out)) {  
  for (c in 1:ncol(rs_out)) {   
    if(!is.na(rs_out$maize_silage_TotalBiomass[r,c])){      
      print(paste(c(r,c, rs_out$maize_silage_TotalBiomass[r,c])))
    }   
  }  
}

outMatrix = as.matrix(rs_out, xy=TRUE, nrow = nrow(rastBase), ncol = ncol(rastBase), xy=TRUE, centroids = TRUE, na.rm=TRUE)
# check if right col/row
nrow(outMatrix)
ncol(outMatrix)