# Purpose: Manipulate large lists of files
# TODO: finilise this NOT FINISHED

#outDir <- "C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\outFiles"
#simDir <- "C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\simFiles"

simDir <- "C:\\MyRun\\SimFiles" # FROM:
outDir <- "F:\\SowByGenotype\\MyRun2016Mar\\OutFiles" # TO:

#simDir <- "F:\\SowByGenotype\\SowByGenotype(raw)\\rawData" # FROM:
#outDir <- "F:\\SowByGenotype\\MyRun2016Mar\\OutFiles_OrigRun" # TO:


flistToCopy <- NULL
flistToDelete <- NULL

#patternFiles <- "67_180|189_56|239_38"
patternFiles <- ".out"

#Get .out files list
#flistToCopy <- list.files(simDir, ".out", 
flistToCopy <- list.files(simDir, pattern = patternFiles, 
                          full.names = TRUE, recursive = TRUE)
flistToCopy
# Copy all .outfiles to .out folder (overwrites old ones)
file.copy(flistToCopy, outDir , overwrite = TRUE)

# Delete all .out and .sum from .sim folder
flistToDelete = 
  paste(c(
    list.files(simDir, ".out", full.names = TRUE, recursive = TRUE),
    list.files(simDir,".sum", full.names = TRUE, recursive = TRUE)
    ))
flistToDelete
# file.remove(flistToDelete) # keeping them for a while as backup
