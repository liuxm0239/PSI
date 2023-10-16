#cghregions_df <- read.delim(filenameOfCGHregionsOutput)

# function CGHregionsDF2stats.R
function(cghregions_df, binSize, filenameOfCGHregionsOutput) {
  # filename must end with string "CGHregions.tab"

# if available locally:
#chromosome_lengths_hg19 <- read.csv("/media/sf_surfdrive/Rhome/centromere-locations/chromosome_lengths_hg19.csv", header=FALSE)
#centromere_locations_hg19 <- read.csv("/media/sf_surfdrive/Rhome/centromere-locations/centromere_locations_hg19.csv")


# Download from google drive
# see here for code: http://thebiobucket.blogspot.nl/2013/04/download-file-from-google-drivedocs.html
# see here for making direct URL: https://sites.google.com/site/gdocs2direct/

library(RCurl)
currentwd = getwd()
setwd(tempdir())
destfile = "chromosome_lengths_hg19.csv" # source: https://genome.ucsc.edu/goldenpath/help/hg19.chrom.sizes
x = getBinaryURL("https://drive.google.com/uc?export=download&id=0Bz7BsLV_yMj0TjhYU0hnamlHZGs", followlocation = TRUE, ssl.verifypeer = FALSE)
writeBin(x, destfile, useBytes = TRUE)
chromosome_lengths_hg19 <- read.csv(destfile, header=FALSE)

destfile = "centromere_locations_hg19.csv" # source: "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/cytoBand.txt.gz" (all cytoband positions)
x = getBinaryURL("https://drive.google.com/uc?export=download&id=0Bz7BsLV_yMj0aUZYeUt2MHRlNUU", followlocation = TRUE, ssl.verifypeer = FALSE)
writeBin(x, destfile, useBytes = TRUE)
centromere_locations_hg19 <- read.csv(destfile)
setwd(currentwd) # set working directory back to the one before download from Google

rownum = dim(cghregions_df)[1]
colnum = dim(cghregions_df)[2]
# binSize = 30000 # for testing
ignoreThisManySpacingBins = 1

# cghregions_df[cghregions_df==-2] <- -1            it's possible to do like that
# cghregions_df[cghregions_df[,6:colnum]==2] <- 1   it's not possible to do like that
cghregions_temp = cghregions_df[,6:colnum]
cghregions_temp[cghregions_temp==2] <- 1
cghregions_temp[cghregions_temp==-2] <- -1
cghregions_df[,6:colnum] = cghregions_temp
cghregions_dfcolnames = colnames(cghregions_df)
dist = matrix(NA,nrow = rownum, ncol=1)
cghregions_df[,"dist"] = dist
cghregions_dfcolnames = c(cghregions_dfcolnames[1:3],"dist", cghregions_dfcolnames[4:length(cghregions_dfcolnames)])# put the new "dist" colum at position of column 
cghregions_df = cghregions_df[,cghregions_dfcolnames]
counter=1;

resultDF = data.frame(matrix(  integer(), nrow = 1, ncol = 6))
colnames(resultDF) = c("sample", "chr", "start", "end", "length", "call")


# general part: calculate the distance from end of one region to the next
for (i in 1:rownum) {
  print(i)
  if(i == rownum || cghregions_df[i+1,1]-1 == cghregions_df[i,1]) # this is last pos or next pos is next chromosome
    cghregions_df[i,"dist"]= chromosome_lengths_hg19[chromosome_lengths_hg19[,1]==paste("chr",cghregions_df[i,1],sep=""),2] - cghregions_df [i,3] - binSize
  else
    cghregions_df[i,"dist"]= cghregions_df[i+1,2] - cghregions_df[i,3] - binSize
}


for (sample in cghregions_dfcolnames[7:length(cghregions_dfcolnames)]) {
  
  #sample = cghregions_dfcolnames[7] #for testing only
  print(sample)
  for (i in 1:rownum) {
    # are we at a start of segment?
    # start of table || copy number of last row not identical || or last row on other chrom || segment of last row more distant than ingnored spacing
    if (i == 1 || cghregions_df[i-1,sample] != cghregions_df[i,sample] ||  cghregions_df[i-1,1]+1 == cghregions_df[i,1] || cghregions_df[i-1,"dist"] > ignoreThisManySpacingBins*binSize) { # either standing at row 1, previous pos not the same call, or previous pos not the same chrom, or previous segment at pos too big has distance
      resultDF[counter,"sample"] = sample # set sample
      resultDF[counter,"chr"] = cghregions_df[i,1] # set chromosome
      resultDF[counter,"start"] = cghregions_df[i,2] # set start of new region
      resultDF[counter,"call"] = cghregions_df[i,sample] # label this region as a loss (-1)
    }
    # are we at an end of segment?  
    if (i == rownum || cghregions_df[i+1,sample] != cghregions_df[i,sample] ||  cghregions_df[i+1,1]-1 == cghregions_df[i,1] || cghregions_df[i,"dist"] > ignoreThisManySpacingBins*binSize) { # either standing at the last row, or next pos not the same call, or this is last pos of chromosome, or end of this segment has too big distance to the next one
      
      resultDF[counter,"chr"] = cghregions_df[i,1] # set chromosome
      resultDF[counter,"end"] = cghregions_df[i,3] + min(binSize, cghregions_df[i,"dist"]) -1 # set end of region
      resultDF[counter,"length"] = resultDF[counter,"end"] - resultDF[counter,"start"] +1
      counter = counter + 1
    }
  }
  
} 


for (j in 1:dim(resultDF)[1]) {
  if (resultDF[j,"end"] <= centromere_locations_hg19[centromere_locations_hg19[,"chromosome"] == paste("chr",resultDF[j,"chr"],sep=""),"centromere.position"])
    resultDF[j,"chrarm"] = "p"
  else if (resultDF[j,"start"] > centromere_locations_hg19[centromere_locations_hg19[,"chromosome"] == paste("chr",resultDF[j,"chr"],sep=""),"centromere.position"])
    resultDF[j,"chrarm"] = "q"
  else { 
    arm_annotation = "c"
    if (resultDF[j,"start"] <= centromere_locations_hg19[centromere_locations_hg19[,"chromosome"] == paste("chr",resultDF[j,"chr"],sep=""),"centromere.position"])
      arm_annotation = paste("p", arm_annotation, sep="")
    if (resultDF[j,"end"] >= centromere_locations_hg19[centromere_locations_hg19[,"chromosome"] == paste("chr",resultDF[j,"chr"],sep=""),"centromere.position"])
      arm_annotation = paste(arm_annotation, "q", sep="")  
    resultDF[j,"chrarm"] = arm_annotation
  }
}

resultDF[,"chr"] = as.factor(resultDF[,"chr"])
resultDF[,"sample"] = as.factor(resultDF[,"sample"])



aberrations = data.frame(matrix(NA, nrow = length(cghregions_dfcolnames)-6, ncol = length(levels(resultDF[,"chr"]))))
rownames(aberrations) = cghregions_dfcolnames[7:length(cghregions_dfcolnames)]

gains_p = aberrations
gains_q = aberrations
losses_p = aberrations
losses_q  = aberrations



for (sample in cghregions_dfcolnames[7:length(cghregions_dfcolnames)]) {
  #  sample = cghregions_dfcolnames[7] {
  for (chr in 1:length(levels(resultDF[,"chr"]))) {
    
    tmp_totlengthcovered =  sum(resultDF[(resultDF$sample == sample) & (resultDF$chr == chr),"length"])
    tmp_totlengthaberrated = sum(resultDF[(resultDF$sample == sample) & (resultDF$chr == chr) & (resultDF$call != 0),"length"])
    aberrations[sample,chr] = 100*tmp_totlengthaberrated/tmp_totlengthcovered
    
    gains_p[sample,chr] = sum((resultDF$sample == sample) & (resultDF$chr == chr) & (resultDF$call == 1) & !(resultDF$chrarm == 'q'))
    gains_q[sample,chr] = sum((resultDF$sample == sample) & (resultDF$chr == chr) & (resultDF$call == 1) & !(resultDF$chrarm == 'p'))
    losses_p[sample,chr] = sum((resultDF$sample == sample) & (resultDF$chr == chr) & (resultDF$call == -1) & !(resultDF$chrarm == 'q'))
    losses_q[sample,chr] = sum((resultDF$sample == sample) & (resultDF$chr == chr) & (resultDF$call == -1) & !(resultDF$chrarm == 'p'))
  }
}

filename_aberrations = sub("CGHregions.tab", "aberrated_percentage_of_chromosomes.csv", filenameOfCGHregionsOutput, perl=TRUE)
write.csv(aberrations, file=filename_aberrations)
filename_gains_p = sub("CGHregions.tab", "p-gain_count.csv", filenameOfCGHregionsOutput, perl=TRUE)
write.csv(gains_p, file=filename_gains_p)
filename_gains_q = sub("CGHregions.tab", "q-gain_count.csv", filenameOfCGHregionsOutput, perl=TRUE)
write.csv(gains_q, file=filename_gains_q)
filename_losses_p = sub("CGHregions.tab", "p-loss_count.csv", filenameOfCGHregionsOutput, perl=TRUE)
write.csv(losses_p, file=filename_losses_p)
filename_losses_q = sub("CGHregions.tab", "q-loss_count.csv", filenameOfCGHregionsOutput, perl=TRUE)
write.csv(losses_q, file=filename_losses_q)

samples =levels(resultDF[,"sample"])
sample_chr_armlength = array(0, c(length(levels(resultDF[,"sample"])), length(levels(resultDF[,"chr"])), 8))
dimnames(sample_chr_armlength)= list(samples,levels(resultDF[,"chr"]),c("lenp","p-gain","p-neutral","p-loss","lenq","q-gain","q-neutral","q-loss"))
# 3rd dimension:
#1: sum_length_p,        short: lenp
#2: sum_length_p_gain    short: p-gain
#3: sum_length_p_neutral short: p-neutral
#4: sum_length_p_loss    short: p-loss
#5: sum_length_q         short: lenq
#6: sum_length_q_gain    short: q-gain
#7: sum_length_q_neutral short: q-neutral
#8: sum_length_q_loss    short: q-loss


# j=1 #for tests
for (j in 1:dim(resultDF)[1]) { # for every row
  sample = resultDF[j,"sample"]
  sampleindex = match(sample, samples)
  chr    = resultDF[j,"chr"]
  start  = resultDF[j,"start"]
  end    = resultDF[j,"end"]
  length = resultDF[j,"length"]
  call   = resultDF[j,"call"]
  chrarm = resultDF[j,"chrarm"]
  
  if(chrarm == 'p') {
    sample_chr_armlength[sampleindex,chr,1] = sample_chr_armlength[sampleindex,chr,1] + length                 #1: sum_length_p
    if(call == 1) sample_chr_armlength[sampleindex,chr,2] = sample_chr_armlength[sampleindex,chr,2] + length    #2: sum_length_p_gain
    if(call == 0) sample_chr_armlength[sampleindex,chr,3] = sample_chr_armlength[sampleindex,chr,3] + length    #3: sum_length_p_neutral
    if(call == -1) sample_chr_armlength[sampleindex,chr,4] = sample_chr_armlength[sampleindex,chr,4] + length    #4: sum_length_p_loss
  }
  
  if(chrarm == 'q') {
    sample_chr_armlength[sampleindex,chr,5] = sample_chr_armlength[sampleindex,chr,5] + length                 #5: sum_length_q
    if(call == 1) sample_chr_armlength[sampleindex,chr,6] = sample_chr_armlength[sampleindex,chr,6] + length    #6: sum_length_q_gain
    if(call == 0) sample_chr_armlength[sampleindex,chr,7] = sample_chr_armlength[sampleindex,chr,7] + length    #7: sum_length_q_neutral
    if(call == -1) sample_chr_armlength[sampleindex,chr,8] = sample_chr_armlength[sampleindex,chr,8] + length    #8: sum_length_q_loss
  }
  
  if(chrarm == 'pcq') {
    centro = centromere_locations_hg19[centromere_locations_hg19$chromosome == paste("chr",chr, sep=""),"centromere.position"]
    start_to_centro = centro - start
    centro_to_end   = end - centro
    sample_chr_armlength[sampleindex,chr,1] = sample_chr_armlength[sampleindex,chr,1] + start_to_centro       #1: sum_length_p
    sample_chr_armlength[sampleindex,chr,5] = sample_chr_armlength[sampleindex,chr,5] + centro_to_end         #5: sum_length_q
    if(call == 1) {
      sample_chr_armlength[sampleindex,chr,2] = sample_chr_armlength[sampleindex,chr,2] + start_to_centro    #2: sum_length_p_gain
      sample_chr_armlength[sampleindex,chr,6] = sample_chr_armlength[sampleindex,chr,6] + centro_to_end      #6: sum_length_q_gain
    }
    if(call == 0) {
      sample_chr_armlength[sampleindex,chr,3] = sample_chr_armlength[sampleindex,chr,3] + start_to_centro    #3: sum_length_p_neutral
      sample_chr_armlength[sampleindex,chr,7] = sample_chr_armlength[sampleindex,chr,7] + centro_to_end      #7: sum_length_q_neutral
    }
    if(call == -1) {
      sample_chr_armlength[sampleindex,chr,4] = sample_chr_armlength[sampleindex,chr,4] + start_to_centro    #4: sum_length_p_loss
      sample_chr_armlength[sampleindex,chr,8] = sample_chr_armlength[sampleindex,chr,8] + centro_to_end      #8: sum_length_q_loss
    }
  }
}

dim(sample_chr_armlength)
statstable_nameslist =dimnames(sample_chr_armlength)[3]
statstable_names = unlist(statstable_nameslist)
pdf(file = sub("CGHregions.tab", "boxplots_of_percentages.pdf", filenameOfCGHregionsOutput, perl=TRUE))
for (k in 1:dim(sample_chr_armlength)[3]) {
  print(k)
  if (k == 1 || k == 5) {
    if (k == 1) {
      relateto = 1 # skip matrix sum_length_p
      next
    }
    if (k == 5) {
      relateto = 5
      next
    }
  }
  temp_matrix = sample_chr_armlength[,,k]/sample_chr_armlength[,,relateto]
  filename_with_path = sub("CGHregions.tab", paste(statstable_names[k],"_perc.csv", sep=""), filenameOfCGHregionsOutput, perl=TRUE)
  write.csv(100*temp_matrix, filename_with_path)  # 100x to get percentages
  boxplot(100*temp_matrix,  outline = FALSE,     ## avoid double-plotting outliers, if any
          main = statstable_names[k])
}
dev.off() 

}

