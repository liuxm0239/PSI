#### usage: Rscript script.R  count_path comparegroup outname out_path
######group_info.tsv #########
##SampleID	Group	Description
##PGSUB_Tumor	test
##PGSUB_Normal	control	

## Only 1|0 should be used in 'Group' column for sake of subsequent analysis; 1 for test/tumor samples  and 0 for control/normal samples
## For sake of subsequent analysis, file.list, sample.id and treatment option should have the same order.

args <- commandArgs(TRUE)
if (length(args)<1) {  stop("usage: Rscript methykit.R  group_info.tsv hg38", call.=FALSE)}
group_info = args[1] #group_info.tsv 
ref = args[2] #assembly reference,  hg19|hg38|mm9|mm10
out='./result.d' ##output direction

library(methylKit)
library(genomation)

meta = read.table(group_info, sep="\t", header=F)
groupname = levels(as.factor(meta$V2))
group1 = meta[which (meta$V2 == groupname[1]),]
group2 = meta[which (meta$V2 == groupname[2]),]
samples1 = group1$V3
samples2 = group2$V3
meta <- rbind(group1, group2)

# Treatment 1 for test; 0 for control samples
#treatment =  c(rep(0,length(samples1)), rep(1,length(samples2)))
treatment = meta$V2
cat("Treatment info: ", treatment, "\n")
#patientID = meta$Description
#covariates = data.frame(patientID=patientID)

grouplabel = c('Control', 'Test')
outdir = paste(out, paste(grouplabel[1], grouplabel[2], sep="-vs-" ), sep="/")
dir.create(outdir, recursive = TRUE)
outname = paste(outdir, paste(grouplabel[1], grouplabel[2], sep="-vs-" ), sep="/")

data_path='./data'
file.list <- lapply(as.vector(meta$V3), function(f){list.files(path=data_path, pattern = f, full.names = TRUE)})
myobj=methRead(file.list, sample.id=as.list(as.vector(meta$V1)), pipeline='bismarkCytosineReport', assembly=ref, treatment=treatment, context="CpG", mincov=10)

filtered.myobj=filterByCoverage(myobj, lo.count=10, lo.perc=NULL,
                                      hi.count=NULL, hi.perc=99.9)
normalized.myobj = normalizeCoverage(filtered.myobj, method = "median")

#=================QC======================#
#getMethylationStats(myobj[[2]],plot=FALSE,both.strands=FALSE)
file.png=paste(outname,".QC.png",sep="")
png(file.png, width=1000, height=600)
par(mfrow=c(1,2))
getMethylationStats(normalized.myobj[[2]],plot=TRUE,both.strands=FALSE)
getCoverageStats(normalized.myobj[[2]],plot=TRUE,both.strands=FALSE)
#filtered.myobj=filterByCoverage(myobj,lo.count=10,lo.perc=NULL, hi.count=NULL,hi.perc=99.9)
#n=as.integer(round(0.8*min(length(sampleName1),length(sampleName2))))
#meth = unite(myobj, destrand=FALSE,min.per.group=n,mc.cores=10)
meth=unite(normalized.myobj, destrand=FALSE, mc.cores=4)
getCorrelation(meth,plot=TRUE)
clusterSamples(meth, dist="correlation", method="ward", plot=TRUE)
PCASamples(meth, screeplot=TRUE)
write.table(data.frame(meth), file=paste(outname, ".meth.data.xls", sep=""), sep="\t", row.names=F, quote=F)
dev.off()

#============calculate diff===============#
#myDiff=calculateDiffMeth(meth,covariates=covariates,overdispersion="MN",mc.cores=10)
myDiff=calculateDiffMeth(meth, overdispersion="MN", mc.cores=8)

write.table(myDiff, file=paste(outname,".alldiff.meth.xls", sep=""), sep="\t", row.names=F, quote=F)
myDiff25p = getMethylDiff(myDiff, difference=25, qvalue=0.05)
write.table(data.frame(myDiff25p), file=paste(outname, 'diff_25_p_0.05.meth.xls', sep="."), sep="\t", row.names=F, quote=F)

#============Annotating diff===============#
# read the gene BED file
gene.obj = readTranscriptFeatures(system.file("extdata", paste("refseq", ref, "bed.txt", sep="."), package = "methylKit", mustWork = TRUE))
diffAnn=annotateWithGeneParts(as(myDiff25p, "GRanges"), gene.obj)
getTargetAnnotationStats(diffAnn, percentage=TRUE, precedence=TRUE)
getFeatsWithTargetsStats(diffAnn, percentage=TRUE)
## Regional analysis
promoters=regionCounts(myobj, gene.obj$promoters)
#
# annotate differentially methylated CpGs with
# promoter/exon/intron using annotation data
#
# read the shores and flanking regions and name the flanks as shores
# and CpG islands as CpGi
cpg.obj=readFeatureFlank(system.file("extdata", paste("cpgi", ref, "bed.txt", sep='.'), package = "methylKit"), feature.flank.name=c("CpGi","shores"))
# convert methylDiff object to GRanges and annotate
diffCpGann=annotateWithFeatureFlank(as(myDiff25p, "GRanges"), cpg.obj$CpGi, cpg.obj$shores, feature.name="CpGi", flank.name="shores")
#

#=====================plot=======================#
myDiff=data.frame(myDiff)
file.png=paste(outname,".diff.p.png",sep="")
png(file.png, width = 1200, height = 800)
par(mfrow=c(2,3))
mygrey <- rgb(89/255, 87/255, 87/255)
par(family="sans", col.main=mygrey, col=mygrey, col.lab=mygrey, col.axis=mygrey, cex.main=2, cex.lab=1.5, cex.axis=1.5, mgp=c(2.5,1,0))
myDiff$logp = -log10(myDiff$pvalue)
myDiff$logq = -log10(myDiff$qvalue)
plot(myDiff[,"meth.diff"], myDiff[,"logp"], ,xlab="meth.diff", ylab="-log10(pvalue)", cex.lab =1.7, pch=20, cex=0.5, main="Volcano plot for Two Groups(p)", cex.main=2, cex.main=2, xlim=range(myDiff[,"meth.diff"]), ylim=range(myDiff$logp))
plot(myDiff[,"meth.diff"], myDiff[,"logq"], ,xlab="meth.diff", ylab="-log10(qvalue)", cex.lab =1.7, pch=20, cex=0.5, main="Volcano plot for Two Groups(q)", cex.main=2, cex.main=2, xlim=range(myDiff[,"meth.diff"]), ylim=range(myDiff$logp))
hist(myDiff$meth.diff, col=mygrey, main="diff", xlab="diff", ylab="Frequency")
hist(myDiff$pvalue, col=mygrey, main="p-value", xlab="p-value", ylab="Frequency")
plotTargetAnnotation(diffAnn, precedence=TRUE, main="differential methylation annotation")
plotTargetAnnotation(diffCpGann, col=c("green","gray","white"), main="differential methylation annotation")
dev.off()

