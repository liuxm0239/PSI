options(stringsAsFactors = FALSE)

library(data.table)
library(ChAMP)
library(ClassDiscovery)
library(gplots)
library(ggplot2)
library(ComplexHeatmap)
library(dplyr)
library(tidyr)
library(tibble)
library(limma)
library(pheatmap)
library(impute)


orgmeth <- fread("LUAD_tcga/LUAD.methylation_beta.txt") #
orgmeth <- as.data.frame(orgmeth); rownames(orgmeth) <- orgmeth[,1]; orgmeth <- orgmeth[,-1]
orgmeth[1:3,1:3]

Sinfo <- read.table("LUAD_tcga/samples", header = T, row.names = 1)
head(Sinfo)

b <- data.frame(Group=Sinfo$Tissue)
row.names(b) <- row.names(Sinfo)
b

beta=as.matrix(orgmeth)
beta=impute.knn(beta)

betaData=beta$data
betaData=betaData+0.00001 #为了避免为0值吗？
a=betaData

identical(colnames(a),rownames(b))

myLoad=champ.filter(beta = a, pd = b) #表达和分组,这过滤了啥，好像没过滤只是制备了champ需要的对象
myNorm <- champ.norm(beta=myLoad$beta,arraytype="450K",cores=5) #标准化
pD=myLoad$pd

#这里和教程得到的结果还不太一样，为什么教程是写的5核运行确实4核心，为什么会有过滤的
#以及'原来的450K经过质控过滤后是350K'是指什么？
# 450K是不是就是指有450k个探针？

group_list=myLoad$pd
table(group_list)

myDMP <- champ.DMP(beta = myNorm,pheno=group_list,adjPVal = 1)

logFC_t = 0.15
P.Value_t = 0.01
df_DMP <- myDMP$Tumor_to_Normal
df_DMP$change <- ifelse(df_DMP$adj.P.Val < P.Value_t & abs(df_DMP$logFC) > logFC_t,
                         ifelse(df_DMP$logFC > logFC_t ,'UP','DOWN'),'NOT')

table(df_DMP$change)
write.table(df_DMP, "./DMP.LUAD.hg19.tsv", sep="\t", quote = F)

dat  = rownames_to_column(df_DMP)
for_label <- dat%>% head(3)
p <- ggplot(data = dat,
            aes(x = logFC,
                y = -log10(adj.P.Val))) +
  geom_point(alpha=0.4, size=3.5,
             aes(color=change)) +
  ylab("-log10(Pvalue)")+
  scale_color_manual(values=c("green", "grey","red"))+
  geom_vline(xintercept=c(-logFC_t,logFC_t),lty=4,col="black",lwd=0.8) +
  geom_hline(yintercept = -log10(P.Value_t),lty=4,col="black",lwd=0.8) +
  theme_bw()
 p

#cg <-  rownames(df_DMP[df_DMP$change != "NOT",])
# plot_matrix <- myNorm[cg,]
# annotation_col <- data.frame(Sample=myLoad$pd) 
# rownames(annotation_col) <- colnames(plot_matrix)
# ann_colors = list(Sample = c(Normal="#4DAF4A", Tumor="#E41A1C"))

# pheatmap(plot_matrix,show_colnames = T,
#          annotation_col = annotation_col,
#          border_color=NA,
#          color = colorRampPalette(colors = c("white","navy"))(50),
#          annotation_colors = ann_colors,show_rownames = F)



