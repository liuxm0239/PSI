# FUNCTIONS ============================================================================================
# Function stopQuietly: create stop function which exits with a blank message (still, an error is generated!)
# TODO: use another way to terminate the program without generating an error
# source: https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
stopQuietly <- function(...) {
      blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
  stop(simpleError(blankMsg));
} # stopQuietly()

# Function pathToScriptWithScriptname
pathToScriptWithScriptname <- function() {
  cmdArgs <- commandArgs(trailingOnly = FALSE)
  needle <- "--file="
  match <- grep(needle, cmdArgs)
  if (length(match) > 0) {
    # Rscript
    return(normalizePath(sub(needle, "", cmdArgs[match])))
  } else {
    # 'source'd via R console
    return(normalizePath(sys.frames()[[1]]$ofile))
  }
}
# End of FUNCTIONS =====================================================================================

args <- commandArgs(trailingOnly = TRUE)

#print("Command line arguments set: ", args)
cat("Command line arguments set: ")
print(args)


if((length(args)==0) || args[1] == "--help" || args[1] == "-h" || !( length(args)==4)) {
  cat(paste("Usage:", 
            "* To show this help: Rscript QDNAseq_BAM2CopyNumbers.R --help OR -h OR without option",
            "* Non-interactive mode: Rscript QDNAseq_BAM2CopyNumbers.R /path/to/BAMdirectory /path/to/OutputDir binSize projectname" ,
            "   Requirements for non-interactive mode (command line input):",
            "   '/path/to/BAMdirectory' must exist",
            "   '/path/to/OutputDir' must exist",
            "   'binSize' must be one of 1, 5, 10, 15, 30, 50, 100, 500 or 1000 kbs",
            "   'projectname' must not contain spaces or special characters",
            "",
            sep="\n"))
  
stopQuietly()
}

pathToBAMs <- args[1]
pathToOutputDir <- args[2]
binSize <- as.integer(args[3])
projectname <- args[4] 

dir.create(file.path(pathToOutputDir), showWarnings = FALSE)

library(QDNAseq)

#binSizes <- c(1, 5, 10, 15, 30, 50, 100, 500, 1000);

rm(args)

#homeDir <- path.expand("~")
# get the path where this script is in:
currentpathToScriptWithScriptname = pathToScriptWithScriptname()
library(stringr)
currentpathToScript = str_extract(string = currentpathToScriptWithScriptname , "/.*/|.*\\\\") # /Unix/ OR C:\Windows\ style of path

qdnaseqDir <- paste0(currentpathToScript ,"QDNAseq.hg38/data/")


#binAnnotationFile = file.path(homeDir, qdnaseqDir, paste("hg38.",binSize,"kbp.SR50.rda",sep=""))
binAnnotationFile = file.path( qdnaseqDir, paste("hg38.",binSize,"kbp.SR50.rda",sep=""))
if(file.exists(binAnnotationFile)){
      bins <- get(load(binAnnotationFile))
} else {
      #bins <- getBinAnnotations(binSize)
      #saveRDS(bins,binAnnotationFile)
}

readCounts <- binReadCounts(bins, path=pathToBAMs, cache=TRUE )
readCountsFiltered <- applyFilters(readCounts,residual=TRUE, blacklist=TRUE)
readCountsFiltered <- estimateCorrection(readCountsFiltered)

# correct for GC content and mappability (applies data previously calculated using estimateCorrection)
copyNumbers <- correctBins(readCountsFiltered)

binSizeDir <- paste(pathToOutputDir, "/","/", projectname, "-binSize-", binSize, "/", sep="")
dir.create(binSizeDir, showWarnings = FALSE, recursive = TRUE)
outputfilenameWithoutEnding= paste(binSizeDir,"/", projectname, "-copyNumbers-",binSize,"kb-bins", sep="")
saveRDS(copyNumbers, file = paste(outputfilenameWithoutEnding,".rds", sep=""))

# export as txt file if needed.
#exportBins(copyNumbers, file = paste(binSizeDir,"/", projectname, "-copyNumbers-",binSize,"kb.rds", sep=""))

exportBins(copyNumbers, file = paste(outputfilenameWithoutEnding,".tsv", sep=""), format="tsv", type="copynumber")
system(paste("bzip2 ", outputfilenameWithoutEnding,".tsv", sep=""))
exportBins(copyNumbers, file = paste(outputfilenameWithoutEnding,".igv", sep=""), format="igv", type="copynumber")
system(paste("bzip2 ", outputfilenameWithoutEnding,".igv", sep=""))
