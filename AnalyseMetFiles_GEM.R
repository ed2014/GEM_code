# Reads .met files for all pixels
# analyses differences in weather
# Transforms it into a raster
# displays the map

library(tcltk2)
library(dplyr)
library(data.table)
library(raster)

# Set directory (chosse one option) <----- Change PATH

#setwd("D:\\EIT\\ERA_40_worked\\MetFiles_1971_2000")
#setwd("D:\\EIT\\RCP8_5_Jan15_worked\\MetFiles_1971_2100")
#setwd("F:\\ERA-40_Reanalysis\\MetFiles_1971_2000")
#setwd("F:\\ERA-40_Reanalysis\\metFiles_ERA_1971_2000") 
setwd("F:\\ERA-40_Reanalysis\\metFiles_Test_ERA_1971_2000")

# # Set up df to hold variables
pixelSummary = NULL
allData2 = NULL
thisOutFile = NULL


# Test if there are not scraps of .csv files in folder
csv.files = list.files(getwd(),pattern='.csv', full.names=FALSE) # if true gets path

# Line for header <------------ Change the line # depending on met file
#lineNoHeader = 9 # for RCP
lineNoHeader = 7 # for old ERA

if (length(csv.files) != 0) {
  tk_messageBox(type = c("ok"),message = "Gonna delete .csv files in out folder", caption = "Warning")
  file.remove(csv.files)
}

# Gets met file list
files <- list.files(getwd(),pattern='.met', full.names=FALSE) # if true gets path
print(paste0("Found ",length(files)," .met files in folder"))


# Filter met files for the LU mask (comment out if doing all country or files were pre-filtered)
#load("C:\\Apsim_dev\\Projects\\CCII\\filter\\Filter_LU_metList.RData", .GlobalEnv)
#load("C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\gisLayers\\luFilter\\LandUseArableFilter_metList.RData", .GlobalEnv)
#files <- metList

#info = file.info(files)
#summary(info)


# # Append files into one dataframe
# Creates folder for storing results
suppressWarnings(dir.create("SummaryResults"))

empty.files <- 0
filesNotFound <- NULL

for (i in 1:length(files)) { 
# Check for empty files
if (file.info(files[i])["size"] == 0 | file.exists(files[i]) == FALSE) {
  print(paste(files[i], " is empty or inexistent"))
  empty.files <- empty.files + 1
  filesNotFound <- c(filesNotFound, files[i]) # compares with arable land filter
  next} 
}

# write list of non-found files (untested
print(paste0("Found ", empty.files, " empty .met files"))
write.table(filesNotFound, file = "filesNotFound.txt")

for (i in 1:length(files)) { 

  # Check for empty files
  if (file.info(files[i])["size"] == 0 | file.exists(files[i]) == FALSE) {
 #   print(paste(files[i], " is empty or inexistent"))
 #   empty.files <- empty.files + 1
 #   filesNotFound <- c(filesNotFound, files[i]) # compares with arable land filter
    next
  }  
  
  if(i<10) {
    thisHeader = read.table(files[i], skip = (lineNoHeader-1), #skip = 6,
                            header = TRUE, comment.char = "(", blank.lines.skip = TRUE) 
    # reads and skips the unit line    
  }
  
  # ## Sort data stamps from file name
  splitName = unlist(strsplit(files[i],"[_,.]"))
  r = as.numeric(splitName[1]) # FIXME: see if first number is indeed row
  c = as.numeric(splitName[2])
  
  # read data for that .out file
  thisOutFile = read.table(files[i], skip = (lineNoHeader+1), header = FALSE, blank.lines.skip = TRUE) # reads and skips the unit line
  
  # From the netCDF metadata (FIXME: should get this from printed field in .out file to minimise risks)
  # !!!!Atention!!!!! Lat/Long "cannot" be rounded otherwise it creates a mismatch of pixels during rasterisation
  lat = -34.375-0.050*(r-1)
  long = 166.425+0.050*(c-1)
  # FIXME: this is still not ideal as V3 and V4 depend on the arrangment in output, need to address the name and the rounding is not perfect
  #  if (lat!= min(thisOutFile$V3) | long != min(thisOutFile$V4)) {print(paste(c("Error: Latitude and/or longitude do not match in row/col: ", r, c)))}
  
  
  thisOutFile = thisOutFile %>% group_by(V1) %>% summarise_each(funs(mean))
  
  #colnames(thisOutFile)
  #  l = length(thisOutFile$Date)
  l = length(thisOutFile$V1) # fills the variable vector with a no of lines needed
  fileNo = rep(i, l)
  row =  rep(r, l)
  col =  rep(c, l)
  thisLat = rep(lat, l) 
  thisLong = rep(long, l)
  
  #thisOutFile = data.frame(fileNo, row, col, thisLat, thisLong, thisOutFile)
  thisOutFile = data.frame(fileNo, row, col, thisLat, thisLong, thisOutFile)
  
  allData2 = rbind(allData2, thisOutFile)
  
  if(i == 1 | i %% 500 == 0 | i == length(files)) {
    print(paste(fileNo[1], " " , Sys.time()))
    write.csv(allData2, file = paste("Out_",i,".csv", sep = ""))
    allData2 = NULL
  }
  
} # end loop in files


firstCols = colnames(thisOutFile) # terieves first columns
#colnames(allData2) = c(firstCols[1:5], colnames(thisHeader))

head(allData2)

summary(allData2)

# # Read csvs
csv.files = list.files(getwd(),pattern='.csv', full.names=FALSE) # if true gets path
all.the.data <- lapply(csv.files, read.csv, header=TRUE)

DATA <- do.call("rbind", all.the.data)
summary(DATA)

# !!!! Attention!!!! this has to be found by hand: how many columns to skip?
colnames(DATA) = c("skip",firstCols[1:5], colnames(thisHeader))
head(DATA)

save(DATA,file="MetData.Rda")
write.csv(DATA, file = "SummaryResults\\MetData.csv") # <--- Change File name

head(DATA)
summary(DATA)

#-----------------------------------------------------
# -------------- Analyse the data now ----------------
#-----------------------------------------------------

df_work <- data.table(DATA)

df.base <- df_work %>% filter(year %in% 1971:2000) %>% summarise_each(funs(mean))

df.fut = df_work %>% filter(year %in% 2060:2099) %>% summarise_each(funs(mean))

diff_maxt = round(df.fut$maxt - df.base$maxt, digits = 1)
diff_mint = round(df.fut$mint - df.base$mint, digits = 1) 
diff_rain = round((df.fut$rain - df.base$rain)*365, digits = 1) 
diff_radn = round((df.fut$radn - df.base$radn)*365, digits = 1)
diff_co2 = round((df.fut$co2 - df.base$co2), digits = 1)

print(
  
  paste0(
    "Differences for future climate:",
    " CO2 (ppm):",diff_co2,
    " Maximum temperature (oC):",diff_maxt,
    " Minimum temperature (oC):",diff_mint,
    " Solar radiation (MJ/year):",diff_radn,
    " Rainfall (mm/year):",diff_rain," or ",round((diff_rain/(df.base$rain*365))*100, digits = 2),"%")  
  
)

#-----------------------------------------------
# Rasterise the data ---------------------------
#----------------------------------------------

# Adjusting for pixel off-set to match with WGS84 layers (!!!!!!!!!---ATENTION---!!!!!!!!!)
# FIXME: this is turned off now ... delete when checked and ok
DATA$latAdj = DATA$thisLat # + 0.05
DATA$lonAdj = DATA$thisLong # - 0.05

# !!!! Attention !!!! : pick by hand the variable # of the variables by which variables you are aggegating (e.g. row was [3] )***
colNamesDF = names(DATA)

# TODO: FIXME: subset the first ("n") years of data (spin-up period) - untested
# !!!!! Attention !!!!! Has to choose the "spin-up" period do me ignored by hand at the moment
# Choose years to EXCLUDE here
DATA_sub = DATA
head(DATA_sub)
summary(DATA_sub)

#Corerct units of rainfall and 

DATA_sub$rain = DATA$rain * 365 # annual rainfall (mm/year)
DATA_sub$radn = DATA$radn * 365 # annual solar radiation (MJ/year)

# TODO: Analyse interannual variability
par(mfcol=c(2,2))
#par(mar=c(12, 4, 2, 4))
par(mar=c(3, 2, 1, 1))
boxplot(maxt ~  year, data = DATA_sub, las = 1, main="maxt")
boxplot(mint ~  year, data = DATA_sub, las = 1, main="mint")
boxplot(radn ~  year, data = DATA_sub, las = 1, main="radn")
boxplot(rain ~  year, data = DATA_sub, las = 1, main="rain")

# Aggregate data by "year" keeping: 'row', 'col', 'CurrentCrop'
outDF = NULL
finalDF = NULL

finalDF_base = DATA_sub %>% 
  filter(year %in% 1971:2000) %>% 
  group_by(row, col) %>% 
  summarise_each(funs(mean))

finalDF_fut = DATA_sub %>% 
  filter(year %in% 2071:2100) %>% 
  group_by(row, col) %>% 
  summarise_each(funs(mean))

# check data
summary(finalDF_base) 
summary(finalDF_fut)
 
# Create folder 
#suppressWarnings(dir.create("SummaryResults"))
write.csv(finalDF_base, file = "SummaryResults\\RasterOutput_base.csv")
write.csv(finalDF_fut, file = "SummaryResults\\RasterOutput_fut.csv")

# ---------------------- Rasterize dataframe -------------------------

outNames = c("maxt","mint","radn","rain")
#scenNames = c("base", "fut1")    # <-------- Change Scenarios
scenNames = c("base") # for ERA_40

par(mfcol=c(1,1))

# creates directories to save results
for (sc in 1:length(scenNames)) {
  dir.create(file.path(getwd(),scenNames[sc]), showWarnings = FALSE)
}

mapNo = 0
for (sc in 1:length(scenNames)) {
  
  thisScen = as.character(scenNames[sc])                                                 
  
  head(finalDF_base)
  
  
  # Creates a dataframe with only info you will need in raster (lat/long and data), 
  # give coordinates, grid it and rasterise it
  spg = list() # dataframe list to be coerced as spatial object
  rast = list() # raster
  
  for(o in 1:length(outNames)){
    
    mapNo = mapNo + 1
    
    # For aggregate
    colnames(finalDF_base)
    
    if (thisScen == "base") {
      
      switch(outNames[o],
             maxt={spg[[o]] = data.frame(finalDF_base$lonAdj, finalDF_base$latAdj, finalDF_base$maxt)},              
             mint={spg[[o]] = data.frame(finalDF_base$lonAdj, finalDF_base$latAdj, finalDF_base$mint)},
             radn={spg[[o]] = data.frame(finalDF_base$lonAdj, finalDF_base$latAdj, finalDF_base$radn)},
             rain={spg[[o]] = data.frame(finalDF_base$lonAdj, finalDF_base$latAdj, finalDF_base$rain)},
{print(paste("Unispecified factor:", outNames[o]))} 
      ) # end switch

coordinates(spg[[o]]) = ~ finalDF_base.lonAdj + finalDF_base.latAdj # Attention to variable names


    }else{
      
      head(finalDF_fut)
      
      switch(outNames[o],
             maxt={spg[[o]] = data.frame(finalDF_fut$lonAdj, finalDF_fut$latAdj, finalDF_fut$maxt)},              
             mint={spg[[o]] = data.frame(finalDF_fut$lonAdj, finalDF_fut$latAdj, finalDF_fut$mint)},
             radn={spg[[o]] = data.frame(finalDF_fut$lonAdj, finalDF_fut$latAdj, finalDF_fut$radn)},
             rain={spg[[o]] = data.frame(finalDF_fut$lonAdj, finalDF_fut$latAdj, finalDF_fut$rain)},
{print(paste("Unispecified factor:", outNames[o]))} 
      ) # end switch

coordinates(spg[[o]]) = ~ finalDF_fut.lonAdj + finalDF_fut.latAdj # Attention to variable names

    } # end of if for Scenario FIXME: quick fix just because the call for variable name is not working

# Common to both methods (grid, raterise and project)
gridded(spg[[o]]) <- TRUE
rast[[o]] = raster(spg[[o]])
proj4string(rast[[o]]) = "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"

# save raster as
thisFolder = paste0(getwd(), "/", scenNames[sc],"/")

# Create directory
writeRaster(rast[[o]], file=paste0(thisFolder,thisScen,"_",outNames[o],".tiff"), overwrite=TRUE) 
# Save "again" in separate folders FIXME: Remove this duplication later (need it as diff maps reads from bulk)
# writeRaster(rast[[o]], file=paste0(thisFolder, thisCrop,"_",thisScen,"_",thisSoil,"_", outNames[o],".tiff"), overwrite=TRUE) 

tit = paste0(thisScen,"_", outNames[o])

#plot
par(mar=c(3,3,3,3))
par(mfrow=c(1,2))
plot(rast[[o]])
boxplot(rast[[o]],main = paste0(outNames[o]," Scen: ",thisScen))

# stack them together
if(mapNo == 1) {
  s = NULL
  s = stack(rast[[o]])} 
else {
  s[[mapNo]] = rast[[o]]
}

  } # end of output loop
} # end of scenario

#plot(s[[1]], main = outNames[1])

# Difference maps for the weather variables (Absolute)
par(mfrow=c(2,2))
for (mc in 1:mapNo) {
  rast_diff = s[[mc+(mapNo/2)]] - s[[mc]]
  plot(rast_diff, main = outNames[mc])
  writeRaster(rast_diff, file=paste0("diffAbs_",outNames[mc],".tiff"), overwrite=TRUE) 
}

# Difference maps for the weather variables (Relative)
par(mfrow=c(2,2))
for (mc in 1:mapNo) {
  rast_diff = s[[mc+(mapNo/2)]] / s[[mc]]
  plot(rast_diff, main = outNames[mc])
  writeRaster(rast_diff, file=paste0("diffRel_",outNames[mc],".tiff"), overwrite=TRUE) 
  
} 


