# Read a base sim file and create it for othe pixels
# 2 Mar 2016 New runs for specific machines MyRun folder

# reads a "base" sim file (XML) finds out a list of met files to run sim files to run these met files
#library("ncdf")
library(XML)
#library("ncdf4")

# Parameters and paths
gc()

# Select folders and file locations
# metFolder <- "C:\\MyRun\\MetFiles"
metFolder <- "C:\\MyRun\\MetFiles\\"
simFolder <- "C:\\MyRun\\SimFiles\\"
rootSimFile <- "C:\\GitHubRepos\\GEM_code\\baseSim\\BaseSim.sim" # find base simulation 
climates <- c("Base\\") # Define climate scenarios to run as different folders
sowDate <- c("-sep","-oct","-nov","-dec","-jan")  
hybrids <- c("h1","h2", "h3","h4", "h5")                                               # select hybrids
#p1 <-  c(150,185,220,255,290) # parameter values for emergence to end of juvenile
p1 <- c(130,160,190,220,250) 
p2 <- c(850,875,900,925,950) # parameter values for flowering to maturity
#p2 <- c(850,850,850,850,850) # test fixed p2 - NO EFFECT
 	  

# metFolder <- "\\\\simplace.net\\projects\\nz\\weatherHistorical\\metFiles\\"    # set folder locations (comment the ones not used)
#metFolder2 <- "F:\\ERA-40_Reanalysis\\metFiles_ERA_1971_2000\\"
#metFolder <- "C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\metFiles\\testPixels\\"
#simFolder <- "Z:\\simFiles\\"    # sim folder to save new .sims
#simFolder <- "F:\\SowByGenotype\\simsRun2\\"    # sim folder to save new .sims
#simFolder <- "C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\simFiles\\"
#rootSimFile <- "C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\GEM_code\\baseSim\\BaseSim.sim" # find base simulation
#climates <- c("Base\\")                                                                # Define climate scenarios to run as different folders
#sowDate <- c("-sep","-oct","-nov","-dec","-jan")   # Define months (-mmm) Obs: it assumes dates later in the scripr FIXME: Quick solution to be improved
#sowDate <- c("-sep","-oct","-nov","-dec")

# select met files to define which .sims to produce
# finds out met files to point to
metFiles <- NULL
#metFiles <- list.files(metFolder2,pattern='.met', full.names=FALSE) # Option 1: gets all met files in a folder
metFiles <- list.files(metFolder,pattern='.met', full.names=FALSE) # Option 1: gets all met files in a folder
#load("C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\gisLayers\\luFilter\\LandUseArableFilter_metList.RData", .GlobalEnv)    # Option 2: gets a list of grid-cell/files selected (LU layer in this case)
# metFiles <- metList
#metFiles <- c("10_135.met" , "100_148.met")

print(paste0("Found ", length(metFiles), " .met files"))

#hybrids <- c("h1","h2", "h3","h4","h5")                                               # select hybrids
#p1 <-  c(120,145,170,195,220) # parameter values for emergence to end of juvenile
#p2 <- c(850,875,900,925,950) # parameter values for flowering to maturity

     
# get data from the root (base) .sim file
#doc = xmlInternalTreeParse(rootSimFile)
doc <- xmlTreeParse(rootSimFile, useInternalNodes = TRUE)
#doc = xmlTreeParse(rootSimFile, getDTD = FALSE, useInternalNodes = TRUE)

# Select all nodels to be modified in sim files
simNameRoot <- xmlRoot(doc)
nodesMet <- getNodeSet(doc, "//filename") 
nodesOut <- getNodeSet(doc, "//outputfile")
nodesStDate <- getNodeSet(doc, "//start_date")
nodesEnDate <- getNodeSet(doc, "//end_date")
nodesStSow <- getNodeSet(doc, "//StartSow")[1] # check if the crop of interest is the first in occurrence (number in bracket)
nodesEnSow <- getNodeSet(doc, "//EndSow")[1]   # 
nodesEnRot <- getNodeSet(doc, "//EndRot")[1]   # 
nodesP1 <- getNodeSet(doc, "//medium/tt_emerg_to_endjuv")   # 
nodesP2 <- getNodeSet(doc, "//medium/tt_flower_to_maturity")   # 

# Loops through climate scenarios "folders"
for(cl in 1:length(climates)){
  
  # choose dates for the current climate scenario
  if(climates[cl] == "Base\\") {
    
    stDate <- "01/01/1971"
    enDate <- "30/12/2000"
    
  } else {
    
    stDate <- "01/01/2063"
    enDate <- "01/01/2098"
  }
  
  # change scenario
  lapply(nodesStDate, function(n) {
    xmlValue(n) = stDate
  })
  lapply(nodesEnDate, function(n) {
    xmlValue(n) = enDate
  })
  
  
  #loop through met file names to create one .sim file for each met file

    for(s in seq_along(sowDate)) {
      #for(s in 4:5) {
      for(h in seq_along(hybrids)) {
       # for(h in 3:5) {
        for(i in seq_along(metFiles)) {
          
    # Create out folder for this scenario
    outDir <- paste0(simFolder,"h",h,"_s",s)
    suppressWarnings(dir.create(outDir))
          
    # For individual pixels  
    #simName = paste0(metFiles[i],"_", hybrids[h],"_s", s)  
    #outName = paste0(simName,".out")
    
    # For list of met files
    splitName <- unlist(strsplit(metFiles[i],"[.]"))
    thisSim <- paste0(splitName[1],"_", hybrids[h],"_s", s)
    simName <- paste0(thisSim,".sim")
    outName <- paste0(thisSim,".out")

    # change attribute name of simulation 
    xmlAttrs(simNameRoot) <- c(name = thisSim)
    
    #  find address to point out to right met files
    newMetNode <- paste0(metFolder,metFiles[i])
    
    # change met location
    lapply(nodesMet, function(n) {
      xmlValue(n) <- newMetNode
    })
    
    # change outfile name
    lapply(nodesOut, function(n) {
      xmlValue(n) <- outName
    })
    
    # change sowing dates
    lapply(nodesStSow, function(n) {
      xmlValue(n) <- paste0(1,sowDate[s]) # just assumes starting with 1st (or 15th) day and moving a day ahead for others
    })
    lapply(nodesEnSow, function(n) {
      xmlValue(n) <- paste0(2,sowDate[s])
    })
    lapply(nodesEnRot, function(n) {
      xmlValue(n) <- paste0(3,sowDate[s])
    })
    # change parameter values
    lapply(nodesP1, function(n) {
      xmlValue(n) <- p1[h]
    })
    lapply(nodesP2, function(n) {
      xmlValue(n) <- p2[h]
    })
    
    # FIXME: identation of saved XML file is corrupted when it has text but no problem with functionality
    saveXML(doc, file <- paste0(outDir,"\\",simName),  indent=TRUE)
    }   # end pixel 
   }  # end hybrid
  } # end sow date
  
} # ends loop on climate scenario folders
