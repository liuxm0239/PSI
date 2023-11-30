args = commandArgs(trailingOnly=TRUE)  

library(ggplot2)
library(dplyr)
library(ggpubr)

inputf <- args[1]
outpre <- args[2]

indata <- read.table(inputf, header=F)

filt <- indata[indata$V4 != "-" & indata$V8 != "-1" ,]
filt$V4 <- as.numeric(filt$V4)
filt$V5 <- as.numeric(filt$V5)
filt$V6 <- as.numeric(filt$V6)
filt$V10 <- as.numeric(filt$V10)
filt$V11 <- as.numeric(filt$V11)
filt$V12 <- as.numeric(filt$V12)
filt <- filt[ filt$V4 > 5 & filt$V4 <100 & filt$V10 > 5 & filt$V10 <100,]

set.seed(1)
tmp = sample_frac(filt, 0.01)

#Use R^2 instead of R
ggscatter(tmp, x = "V4", y = "V10",  color = "blue",
         main = paste0("Methylation levels compare of ", outpre ), 
         add = "reg.line",
         add.params = list(color = "black", fill = "lightgray"), # Customize reg. line
         conf.int = TRUE, # Add confidence interval
         #cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
         xlab = 'BSmap CpG' , ylab = 'Bismark CpG') +

  stat_cor(label.y = 95,
           aes(label = paste(..rr.label.., ..p.label.., sep = "~`,`~"))) +
  stat_regline_equation(label.y = 90)

# #Groupped scatter plot
# ggscatter(tmp, x = "V4", y = "V10",
#   color = "TRT", palette = "jco",
#   add = "reg.line"
# ) +
#   facet_wrap(~TRT) +
#   stat_cor(label.y = 4.4) +
#   stat_regline_equation(label.y = 4.2)
#saving final graph
ggsave(paste0(outpre, ".png"), width = 5.2, height = 5, dpi = 1000)
