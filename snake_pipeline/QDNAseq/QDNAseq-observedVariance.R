# this code will take a copyNumbers-object as input, calculate segments, var_expect, var_observed, diffvar, total_reads and return them wrapped as dataframe statsDF
# major hints how to do this have been obtained from Daoud Sie.

obsvar = apply(assayDataElement(copyNumbersCalled, "copynumber"), 2, QDNAseq:::sdDiffTrim, na.rm=T)
samples = names(obsvar)


# Extract segments
sn <- assayDataElement(copyNumbersCalled, "segmented")

# Extract feature data
fd <- fData(copyNumbersCalled)

library("CGHbase")
## # Import private function from CGHbase
## mks <- CGHbase:::.mkSegments

## # Calculate segments for first sample
## nrow(mks(sn[fd$use,1], fd$chromosome[fd$use]))

## Import private functions
ns <- asNamespace("CGHbase")
.getChromosomeLengths <- get(".getChromosomeLengths", envir=ns, mode="function")
.makeSegments <- get(".makeSegments", envir=ns, mode="function")

# Calculate segments for all samples
segments=replace(obsvar, 1:length(obsvar), NA) # easy way to create an array 'segments' of length of obsvar
for (i in 1:length(samples)) {
  segments[i]= nrow(.makeSegments(sn[fd$use,i], fd$chromosome[fd$use]))
}

# expected variance:
var_expect = sqrt(copyNumbersCalled@phenoData@data$expected.variance)
# observed variance: 
var_observed = apply(assayDataElement(copyNumbersCalled, "copynumber"), 2, QDNAseq:::sdDiffTrim, na.rm=T)

diffvar = var_observed - var_expect
total_reads = copyNumbersCalled@phenoData@data$total.reads

statsDF = data.frame(segments, var_expect, var_observed, diffvar, total_reads)
