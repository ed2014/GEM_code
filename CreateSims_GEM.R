# Read a base sim file and create it for othe pixels

# reads a "base" sim file (XML) finds out a list of met files to run sim files to run these met files
library("ncdf")
library(XML)

# Parameters and paths
gc()

# Select folders and file locations
metFolder <- "\\\\simplace.net\\projects\\nz\\weatherHistorical\\metFiles\\"                                                  # set folder locations (comment the ones not used)
simFolder <- "Z:\\simFiles\\"                                                   # sim folder to save new .sims
rootSimFile <- "C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\GEM_code\\baseSim\\BaseSim.sim" # find base simulation
climates <- c("Base\\")                                                                # Define climate scenarios to run as different folders
sowDate <- c("-sep","-oct","-nov","-dec","-jan")   # Define months (-mmm) Obs: it assumes dates later in the scripr FIXME: Quick solution to be improved
#metFiles <- c("67_180.met", "189_60.met", "228_80.met")                           # Option 1: select pixels (row_column) from met files
load("C:\\apsim_dev\\Projects\\CCII\\filter\\Filter_LU_metList.RData", .GlobalEnv)    # Option 2: gets a list of grid-cell/files selected (LU layer in this case)
metFiles <- metList
hybrids = c("h1","h2", "h3","h4","h5")                                               # select hybrids
p1 <-  c(120,145,170,195,220) # parameter values for emergence to end of juvenile
p2 <- c(850,875,900,925,950) # parameter values for flowering to maturity

length(metFiles)

# get data from the root (base) .sim file
#doc = xmlInternalTreeParse(rootSimFile)
doc = xmlTreeParse(rootSimFile, useInternalNodes = TRUE)
#doc = xmlTreeParse(rootSimFile, getDTD = FALSE, useInternalNodes = TRUE)

# Select all nodels to be modified in sim files
simNameRoot = xmlRoot(doc)
nodesMet = getNodeSet(doc, "//filename") 
nodesOut = getNodeSet(doc, "//outputfile")
nodesStDate = getNodeSet(doc, "//start_date")
nodesEnDate = getNodeSet(doc, "//end_date")
nodesStSow = getNodeSet(doc, "//StartSow")[1] # check if the crop of interest is the first in occurrence (number in bracket)
nodesEnSow = getNodeSet(doc, "//EndSow")[1]   # 
nodesEnRot = getNodeSet(doc, "//EndRot")[1]   # 
nodesP1 = getNodeSet(doc, "//medium/tt_emerg_to_endjuv")   # 
nodesP2 = getNodeSet(doc, "//medium/tt_flower_to_maturity")   # 


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
      for(h in seq_along(hybrids)) {
        for(i in seq_along(metFiles)) {
          
    # Create out folder for this scenario
    outDir = paste0(simFolder,"h",h,"_s",s)
    suppressWarnings(dir.create(outDir))
          
    # For individual pixels  
    #simName = paste0(metFiles[i],"_", hybrids[h],"_s", s)  
    #outName = paste0(simName,".out")
    
    # For list of met files
    splitName = unlist(strsplit(metFiles[i],"[.]"))
    thisSim = paste0(splitName[1],"_", hybrids[h],"_s", s)
    simName = paste0(thisSim,".sim")
    outName = paste0(thisSim,".out")

    # change attribute name of simulation 
    xmlAttrs(simNameRoot) = c(name = thisSim)
    
    #  find address to point out to right met files
    newMetNode = paste0(metFolder,metFiles[i])
    
    # change met location
    lapply(nodesMet, function(n) {
      xmlValue(n) = newMetNode
    })
    
    # change outfile name
    lapply(nodesOut, function(n) {
      xmlValue(n) = outName
    })
    
    # change sowing dates
    lapply(nodesStSow, function(n) {
      xmlValue(n) = paste0(1,sowDate[s]) # just assumes starting with 1st day and moving a day ahead for others
    })
    lapply(nodesEnSow, function(n) {
      xmlValue(n) = paste0(2,sowDate[s])
    })
    lapply(nodesEnRot, function(n) {
      xmlValue(n) = paste0(3,sowDate[s])
    })
    # change parameter values
    lapply(nodesP1, function(n) {
      xmlValue(n) = p1[h]
    })
    lapply(nodesP2, function(n) {
      xmlValue(n) = p2[h]
    })
    
    # FIXME: identation of saved XML file is corrupted when it has text but no problem with functionality
    saveXML(doc, file = paste0(outDir,"\\",simName),  indent=TRUE)
    }   # end pixel 
   }  # end hybrid
  } # end sow date
  
} # ends loop on climate scenario folders
