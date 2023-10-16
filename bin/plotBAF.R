args = commandArgs(trailingOnly=TRUE)
library('CopyNumberPlots')

vcfin   = args[1]
plotout = args[2]
genoref = args[3]

if (!exists(genoref)){
    genoref = "hg38"
}

s1 <- loadSNPDataFromVCF(vcf.file = vcfin, genome=genoref, mirror.baf = FALSE)

pdf(paste0( plotout, "chr1.baf.pdf"), width=8, height=3)
kp <- plotKaryotype(chromosomes="chr1")
plotBAF(kp, snps=s1, points.cex=0.3, 
        labels=names(s1), label.margin=0.06, 
        label.cex=1, track.margin=0.2) 
dev.off()
