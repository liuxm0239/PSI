# FUNCTIONS ============================================================================================
# Function sampleNames: return the sample names from a QDNAseq object
sampleNames <- function(largeQDNAseqCopyNumbers){
  sampleNames <- {largeQDNAseqCopyNumbers@phenoData}@data$name
  return(sampleNames)
}

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


# create stop function which exits with a blank message (still, an error is generated!)
# source: https://stackoverflow.com/questions/14469522/stop-an-r-program-without-error
stopQuietly <- function(...) {
  blankMsg <- sprintf("\r%s\r", paste(rep(" ", getOption("width")-1L), collapse=" "));
  stop(simpleError(blankMsg));
} # stopQuietly()

#options(echo=TRUE) # if you want to see commands in output file
args <- commandArgs(trailingOnly = TRUE)

#print("Command line arguments set: ", args)
cat("Command line arguments set: ")
print(args)


if((length(args)==0) || args[1] == "--help" || args[1] == "-h" ) {
  cat(paste("Usage:", 
            "* To show this help: Rscript QDNAseq_CopyNumbers2PlotsCallsSegmentFileInclDewave.R --help OR -h OR without option",
            "* Rscript QDNAseq_CopyNumbers2PlotsCallsSegmentFileInclDewave.R followed by command line inputs:",
            "  in case of 3 arguments: /path/to/copyNumbers-xxkb-bins.rds /path/to/OutputDir projectname ",
            "Note:",
            "  copy number file must look like this: *copyNumbers-xxkb-bins.rds",
            "",
            sep="\n"))
  
  stopQuietly()
}


undoSD=1 # default in QDNAseq segmentation
alpha=1e-10 # default in QDNAseq segmentation

library(QDNAseq)
library(R.cache, quietly = TRUE)
library(CGHregions)

# 4 arguments: classic
# 6 arguments: classic plus undoSD and alpha

# 4 arguments: classic
# /path/to/copyNumbers-xxkb-bins.rds /path/to/OutputDir projectname dewavingKey (or n for no dewaving)",
if(length(args)==3) {
  copynumberfile  <- args[1]
  pathToOutputDir <- args[2]
  projectname     <- args[3]
} 

if(!file.exists(copynumberfile)|| !file.exists(pathToOutputDir))
{
  stop("Please check: copynumberfile must exist and its filname be of format *copyNumbers-xxkb-bins.rds, /path/to/OutputDir (must exist)")
}

binSize=as.integer(gsub(".*copyNumbers-|kb-bins.rds", "", copynumberfile))
copyNumbers = readRDS(copynumberfile)
copyNumbersNormalized <- normalizeBins(copyNumbers) # normalize bins; default is median normalization
copyNumbersSmooth <- smoothOutlierBins(copyNumbersNormalized) # smooth outliers based on previous normalization

exportBins(copyNumbersSmooth, file = paste("copyNumbersSmooth.tsv", sep=""), format="tsv", type="copynumber")

dir.create(file.path(pathToOutputDir, projectname, paste(binSize,"kb-bins", sep="")), showWarnings = FALSE, recursive=TRUE)
smoothPlotDir <- paste(pathToOutputDir, "/", projectname, "/", binSize, "kb-bins/SmoothPlot/", sep="")
dir.create(smoothPlotDir, showWarnings = FALSE, recursive=TRUE)
for (i in 1:ncol(copyNumbersSmooth)) {
  png.name <- paste(smoothPlotDir, sampleNames(copyNumbersSmooth)[i], "_", projectname,"_SmoothPlot_",binSize,"kb.png", sep="")
  png(png.name, width = 1280, height = 1024)
  plot(copyNumbersSmooth[,i])
  dev.off()
}

# experimental:
copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt", segmentStatistic ="seg.median")

# original: 
#copyNumbersSegmented <- segmentBins(copyNumbersSmooth, transformFun="sqrt", undo.SD=undoSD, alpha=alpha)

copyNumbersSegmented <- normalizeSegmentedBins(copyNumbersSegmented)

sgmentedPlotDir <- paste(pathToOutputDir, "/", projectname, "/", binSize, "kb-bins/SegmentedPlot/", sep="")
dir.create(sgmentedPlotDir, showWarnings = FALSE, recursive=TRUE)
for (i in 1:ncol(copyNumbersSegmented)) {
  png.name <- paste(sgmentedPlotDir, sampleNames(copyNumbersSegmented)[i], "_", projectname,"_SegmentedPlot_",binSize,"kb.png", sep="")
  png(png.name, width = 1280, height = 1024)
  plot(copyNumbersSegmented[,i])
  dev.off()
}


copyNumbersCalled <- callBins(copyNumbersSegmented)


calledPlotDir <- paste(pathToOutputDir, "/", projectname, "/", binSize, "kb-bins/CalledPlot/", sep="")
dir.create(calledPlotDir, showWarnings = FALSE, recursive=TRUE)
for (i in 1:ncol(copyNumbersCalled)) {
  png.name <- paste(calledPlotDir, sampleNames(copyNumbersCalled)[i], "_", projectname,"_CalledPlot_",binSize,"kb.png", sep="")
  png(png.name, width = 1280, height = 1024)
  plot(copyNumbersCalled[,i])
  dev.off()
}

saveRDS(copyNumbersCalled, file = paste(pathToOutputDir, "/", projectname, "/", binSize, "kb-bins/",projectname,"-copyNumbersCalled-",binSize,"kb-bins.rds", sep=""))
#saveRDS(copyNumbersCalled, file = "copyNumbersCalled-30kb-65bpReads.rds")
# copyNumbersCalled

cghcall <- makeCgh(copyNumbersCalled)


#cghcall_df <- data.frame(featureNames(cghcall), chromosomes(cghcall), bpstart(cghcall), bpend(cghcall), calls(cghcall), copynumber(cghcall), segmented(cghcall))
cghcall_df <- data.frame(featureNames(cghcall), chromosomes(cghcall), bpstart(cghcall), bpend(cghcall), calls(cghcall))
names(cghcall_df)[1] <- paste("Name")
names(cghcall_df)[2] <- paste("Chromosome")
names(cghcall_df)[3] <- paste("Start")
names(cghcall_df)[4] <- paste("End")
write.table(cghcall_df, paste(pathToOutputDir, "/", projectname , "/", binSize, "kb-bins/",projectname,"-copyNumbersCalled-",binSize,"kb-bins.tab", sep=""),col.names=TRUE, row.names=FALSE,quote=F,sep="\t")


frequencyPlot(copyNumbersCalled)
pdf(file=paste(pathToOutputDir, "/", projectname , "/", binSize, "kb-bins/", projectname, "-frequencyPlot-", binSize, "kb-bins.pdf", sep=""))
frequencyPlot(copyNumbersCalled)
dev.off()

# calculate statistics matrix:
OutputDirStats <- paste(pathToOutputDir, "/", projectname, "/", binSize, "kb-bins/stats/", sep="")
dir.create(OutputDirStats, showWarnings = FALSE, recursive=TRUE)

# get the path where this script is in:
currentpathToScriptWithScriptname = pathToScriptWithScriptname()
library(stringr)
#currentpathToScript = str_extract(string = currentpathToScriptWithScriptname , "/.*/")
currentpathToScript = str_extract(string = currentpathToScriptWithScriptname , "/.*/|.*\\\\") # /Unix/ OR C:\Windows\ style of path

source(paste0(currentpathToScript, "/QDNAseq-observedVariance.R")) # this code will take a copyNumbers-object as input, calculate segments, var_expect, var_observed, diffvar, total_reads and return them wrapped as dataframe statsDF

Outlier_removal_in_R_using_IQR_rule <- dget(paste0(currentpathToScript, "/Outlier_removal_in_R_using_IQR_rule.R"))

xlsOutputFile <- paste(projectname, "-statistics-", binSize, "kb-bins.xlsx", sep="")
#Outlier_removal_in_R_using_IQR_rule(OutputDirStats, xlsOutputFile, statsDF) ## rJava bug, to be fix later



# CGHregions: calculation of regions

filenameOfCGHregionsOutput <- paste0("/", projectname, "-", binSize, "kb-bins.point01percentsmoothing.CGHregions.tab")



# copyNumbersCalled = readRDS(filenameOfCopyNumbersCalledRDSfile) # needed for importing copy numbers here. 
# cghcall <- makeCgh(copyNumbersCalled) # this was done already above
#cghregions_normalRegioning1percent <- CGHregions(cghcall, averror=0.01) # 1% error rate is considered "normal" according to CGHregions paper, when intending to compare groups
#cghregions_lenientRegioningPoint01percent <- CGHregions(cghcall, averror=0.0001) # 0.01 % error rate will result in extremely lenient smoothing
#cghregions_severeRegioning2point5percent <- CGHregions(cghcall, averror=0.025) # 2.5 % error rate will result in extremely dramatic smoothing, when intending to compare groups with <= 10 members
#cghregions <- CGHregions(cghcall, averror=0.0001) # =0.01 % error rate will result in extremely lenient smoothing

#cghregions_df <- data.frame( chromosomes(cghregions), bpstart(cghregions), bpend(cghregions), nclone(cghregions), avedist(cghregions), regions(cghregions) )
#write.table(cghregions_df, paste0(pathToOutputDir, "/", projectname , "/", binSize, "kb-bins/",filenameOfCGHregionsOutput), quote=F, sep="\t")

#CGHregionsDF2stats <- dget(paste0(currentpathToScript, "/CGHregionsDF2stats.R"))
#CGHregionsDF2stats(cghregions_df, binSize, paste0(OutputDirStats,filenameOfCGHregionsOutput))
