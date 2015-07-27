# Purpose: Manipulate large lists of files
# TODO: finilise this 

outDir <- "C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\outFiles"
simDir <- "C:\\Apsim_dev\\Projects\\2014-SowByGenotype\\simFiles"
flistToCopy <- NULL
flistToDelete <- NULL

#Get .out files list
flistToCopy <- list.files(simDir, ".out", 
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
file.remove(flistToDelete)
