gc()
# set wd
setwd("F:\\SowByGenotype")
#setwd("C:\\Dropbox\\2014-GEM_paper\\rawData")

# Get all files including from sub-folders (recursive)
files = list.files(getwd(),
                   pattern='.out', 
                   full.names=FALSE,
                   include.dirs = FALSE,
                   recursive = TRUE) # if true gets path

files[1]

# Show me the money!
print(paste0("Found ",length(files)," .out files in folder"))


# # Set up df to hold variables
pixelSummary = NULL
allData = NULL
thisOutFile = NULL


empty.files = 0
for (i in 1:length(files)) { 
 # for (i in 1:200) {
  
  if(i<10) {
    thisHeader = read.table(files[i], 
                            skip = 2, 
                            header = TRUE, 
                            comment.char = "(", blank.lines.skip = TRUE) 
    # reads and skips the unit line    
  }
  
  # Check for empty files
  if (file.info(files[i])["size"] == 0) {
    print(paste(files[i], " is empty "))
    empty.files = empty.files + 1
    next}  
  
  # ## Sort data stamps from file name
  splitName = unlist(strsplit(files[i],"[/,_,.]"))
  r = as.numeric(splitName[3]) # FIXME: see if first number is indeed row
  c = as.numeric(splitName[4])
  hyb = splitName[5]  # scenario
  sow = splitName[6] # soil type
  
  # read data for that .out file
  thisOutFile = read.table(files[i], skip = 4, header = FALSE, blank.lines.skip = TRUE) # reads and skips the unit line
  
  colnames(thisOutFile) =  colnames(thisHeader)
  
  # From the netCDF metadata (FIXME: should get this from printed field in .out file to minimise risks)
  # !!!!Atention!!!!! Lat/Long "cannot" be rounded otherwise it creates a mismatch of pixels during rasterisation
  thisLat = -34.375-0.050*(r-1)
  thisLong = 166.425+0.050*(c-1)
  # FIXME: this is still not ideal as V3 and V4 depend on the arrangment in output, need to address the name and the rounding is not perfect
  #  if (lat!= min(thisOutFile$V3) | long != min(thisOutFile$V4)) {print(paste(c("Error: Latitude and/or longitude do not match in row/col: ", r, c)))}
  
  #colnames(thisOutFile)
  #  l = length(thisOutFile$Date)
  l = length(thisOutFile[,1]) # fills the variable vector with a no of lines needed
  fileNo = rep(i, l)
  row =  rep(r, l)
  col =  rep(c, l)
  thisHyb = rep(hyb, l) 
  thisSow = rep(sow, l)

  
  thisOutFile = data.frame(fileNo, row, col, 
                           thisLat, thisLong,
                           thisHyb, thisSow, thisOutFile)
  
  allData = rbind(allData, thisOutFile)
  
  if(i == 1 | i %% 500 == 0 | i == length(files)) {
    print(paste(fileNo[1], " " , Sys.time()))
    write.csv(allData, file = paste("Out_",i,".csv", sep = ""))
    allData = NULL
  }
  
} # end loop in files

print(paste0("Found ", length(empty.files), " empty .out files"))

#firstCols = colnames(thisOutFile) # terieves first columns
#colnames(allData) = c(firstCols[1:5], colnames(thisHeader))

head(allData)

summary(allData)

# idea: more efficient appending with dplyr???

# # Read csvs
csv.files = list.files(getwd(),pattern='.csv', full.names=FALSE) # if true gets path
all.the.data <- lapply(csv.files, read.csv, header=TRUE)
length(csv.files)
DATA_GEM <- do.call("rbind", all.the.data)
summary(DATA_GEM)
unique(DATA_GEM$CurrentCrop)
head(DATA_GEM)

# Save it as a single df
save(DATA_GEM,file="DATA_GEM.Rda")
# write.csv(DATA_GEM, file = "AllData(GEM).csv")

