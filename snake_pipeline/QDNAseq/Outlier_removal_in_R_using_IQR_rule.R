function(OutputDir, xlsOutputFile, inputDataframe) {
# this function takes a dataframe of numeric values as input, column names required
# output is an excel file with name: xlsOutputFile at location: OutputDir where
# outliers are colored according to the Interquartile Range-based outlier criterion as originally defined by John Tukey and as known from typical boxplots
# colors are: red = upper outlier, blue = lower outlier
  
# Toy data:
#d = matrix(rnorm(40),nrow=20,ncol=2)
#colnames(d) = c("spalte1", "spalte2")
# ------------------------------
# Load the data however you want
# ------------------------------
# make sure your data has column names
# I called my data dsBase, here I copy the data to dsBase.iqr
# I wanted to keep a copy of the original data set
dsBase <- inputDataframe

dsBase.iqr <- dsBase

# Create a variable/vector/collection of the column names you want to remove outliers on.
#vars <- c("ColName1", "ColName2", "ColName3", "etc")
vars = colnames(dsBase)
# Create a variable to store the row id's to be removed
Outliers <- c()
OutlierMatrix <- matrix(0L,nrow = dim(dsBase.iqr)[1], ncol = dim(dsBase.iqr)[2])
colnames(OutlierMatrix) = colnames(dsBase)
# Loop through the list of columns you specified
for(i in vars){
  #for testing set: i = "segments"
  
  
  # Get the Min/Max values
  max <- quantile(dsBase.iqr[,i],0.75, na.rm=TRUE) + 1.5 * (IQR(dsBase.iqr[,i], na.rm=TRUE))
  min <- quantile(dsBase.iqr[,i],0.25, na.rm=TRUE) - 1.5 * (IQR(dsBase.iqr[,i], na.rm=TRUE))
  
  # Get the id's using which
  idx <- which(dsBase.iqr[,i] < min | dsBase.iqr[,i] > max)
  idx_lower <- which(dsBase.iqr[,i] < min)
  OutlierMatrix[idx_lower,i] = -1
  idx_upper <- which(dsBase.iqr[,i] > max)
  OutlierMatrix[idx_upper,i] = 1
  
  # Output the number of outliers in each variable
  print(paste(i, length(idx), sep=''))
  
  # Append the outliers list
  Outliers <- c(Outliers, idx) 
}

# Sort, to get the rows (= samples) that are outliers in numberic order
Outliers <- sort(Outliers)

# Remove the outliers if wanted
#dsBase.iqr <- dsBase.iqr[-Outliers,]

# ==============================================

################################################ First, reformat the data a little.

# split the X column so there will be one numeric entry per cell 
## CR ## d <- matrix(as.numeric(unlist(strsplit(as.character(Data$X), ";"))), 
## CR ##        ncol = 20, byrow = TRUE)


## CR ##d <- data.frame(d, Data$Y)
inputDataframe = dsBase.iqr
cols <- length(inputDataframe[1, ]) # number of columns, we'll use this later

################################################ Second, we can use functions in xlsx to create a workbook, and then get at the cell values.
list.of.packages <- c("xlsx")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages, repos='https://cran.uni-muenster.de/')
library(xlsx)

# exporting data.frame to excel is easy with xlsx package
sheetname <- "stats"
write.xlsx(inputDataframe, paste0(OutputDir,xlsOutputFile), sheetName=sheetname)

# but we want to highlight cells if value greater than or equal to 5
wb <- loadWorkbook(paste0(OutputDir,xlsOutputFile))              # load workbook
fo1 <- Fill(foregroundColor="blue")   # create fill object # 1
cs1 <- CellStyle(wb, fill=fo1)        # create cell style # 1
fo2 <- Fill(foregroundColor="red")    # create fill object # 2
cs2 <- CellStyle(wb, fill=fo2)        # create cell style # 2 
sheets <- getSheets(wb)               # get all sheets
sheet <- sheets[[sheetname]]          # get specific sheet
rows <- getRows(sheet, rowIndex=2:(nrow(inputDataframe)+1))     # get rows
# 1st row is headers
cells <- getCells(rows, colIndex = 2:(cols+1))         # get cells

# in the wb I import with loadWorkbook, numeric data starts in column 2
# The first column is row IDs, thus we use colIndex = 2:(cols+1)

values <- lapply(cells, getCellValue) # extract the cell values

################################################ Next we find the cells that need to be formatted according to the criteria.

# find cells meeting conditional criteria to be highlighted
highlightblue <- NULL
highlightred <- NULL

for (i in names(values)) {
  x <- as.numeric(values[i])
  coord2d = as.integer(unlist(strsplit(i,"[.]")))
  if (OutlierMatrix[coord2d[1]-1,coord2d[2]-1] == -1 && !is.na(x)) {
    highlightblue <- c(highlightblue, i)
  }
  if (OutlierMatrix[coord2d[1]-1,coord2d[2]-1] == 1 && !is.na(x)) {
    highlightred <- c(highlightred, i)
  }
  
}

################################################ Finally, apply the formatting and save the workbook.


lapply(names(cells[highlightblue]),
       function(ii) setCellStyle(cells[[ii]], cs1))

lapply(names(cells[highlightred]),
       function(ii) setCellStyle(cells[[ii]], cs2))

saveWorkbook(wb, paste(OutputDir,xlsOutputFile,sep=""))

}
