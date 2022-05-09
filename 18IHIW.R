### 18 IHIW HLA-Associated T1D in World Populations Data Analysis Session -- SJM May 7, 2022
### Demonstration Script for Running BIGDAWG on published Azerbaijani T1D Data (https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6384092/)

### Set the Working Directory

## if you are using a Mac/Linux computer
setwd("Documents/18IHIW.T1D")

## if you are using a PC
setwd("18IHIW.T1D/") 

## Load the Azerbaijani data set
AZ.DRB1 <- read.delim("AZ_18IHIW.txt", stringsAsFactors=FALSE)
## View the data set to make sure it looks appropriate
View(AZ.DRB1)

## Attach BIGDAWG
library(BIGDAWG)

###################################################################################
# There is an issue with BIGDAWG under the current version of R (4.2.0).          #
# So, if you have just installed BIGDAWG for this exercise, skip ahead to lines   #
# 37-40 below after lines 27 and 35 have been run to follow on your computer.     # 
###################################################################################

### Run BIGDAWG at the locus level, and generate a set of output files in your working directory
BIGDAWG(Data = AZ.DRB1,Run.Tests = "L")

# We can review at each file in the Files pane by clicking on it in the Files pane. 

# We can Load the Analysis.RData list object into our environment by clicking on it in the Files pane. 

### We can generate a named Analysis.RData object directly in R, with a unique name, without writing it to the working directory. 

AZ.Locus <- BIGDAWG(Data = AZ.DRB1,Run.Tests = "L", Return = TRUE,Output = FALSE)

## IF YOU HAVE NOT BEEN ABLE TO RUN LINES 27 AND 35, RUN LINE 40 BELOW. 
## IF YOU WERE ABLE TO RUN LINES 27 AND 35, YOU DON'T NEED TO RUN LINE 40.

load("AZ.Locus.rda")

## Examine the overall locus level of significance evaluated via the chi-squared test. 

View(AZ.Locus$L$Set1$chisq)

# If the value of this overall test is significant then you don't need to adjust the p-value
# threshold of significance for the OR tests for multiple comparisons.

## Perform a simple Hardy-Weinberg (compared genotypes to those expected based on the allele frequencies) test of the cases and controls
## IF YOU ARE NOT RUNINING BIGDWG, RUN LINE 54 BELOW.

AZ.HWE <- BIGDAWG(Data = AZ.DRB1,Run.Tests = "HWE", Return = TRUE,Output = FALSE)

load("AZ.Locus.rda")

## Look at the HWE results for cases and controls 
AZ.HWE$HWE$cases
AZ.HWE$HWE$controls

# Cases are out of HWE, which is not a problem. That suggests that there is something unusual about the obvserved genotypes in cases

### Extract only associated alleles' data from the results tables table

associated.OR <- as.data.frame(AZ.Locus$L$Set1$OR[AZ.Locus$L$Set1$OR[,7] == "*",]) # everything in column 7 of the ORttable with an * is associated
associated.freq <- as.data.frame(AZ.Locus$L$Set1$freq[AZ.Locus$L$Set1$freq[,2] %in% associated.OR[,2],3:4]) # all alleles in associated.OR that are in the frequency table
associated.count <- as.data.frame(AZ.Locus$L$Set1$table[AZ.Locus$L$Set1$table[,2] %in% associated.OR[,2],3:4]) # all alleles in associated.OR that are in the counts table

## Assemble counts, frequencies and ORs for associated alleles, starting with the first two columns of the OR table

associated.table <- cbind(as.character(associated.OR[,1]),as.character(associated.OR[,2]),associated.freq[,1],associated.count[,1],associated.freq[,2],
                          associated.count[,2],as.character(associated.OR[,3]), as.character(associated.OR[,4]), as.character(associated.OR[,5]),
                          as.character(associated.OR[,6]))

## Rename the columns for the frequency and counts 
colnames(associated.table) <- c("Locus","Allele","Control (f)","Control (n)","Case (f)","Case (n)","Odds Ratio", "CI Lower", "CI Upper", "P-value")

View(associated.table)

## Save a tab-delimited table of associated values

write.table(x = associated.table,file = "AZ_18IHIW_associations.txt",append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE)

### Output a similar table for all of the results

## This requires accounting for the binned category, which is not included in the frequency table 

## Assemble a table of frequencies for all individual alleles in the OR table

all.freq <- as.data.frame(AZ.Locus$L$Set1$freq[AZ.Locus$L$Set1$freq[,2] %in% AZ.Locus$L$Set1$OR[,2],]) # All individual alleles 
all.freq <- rbind(all.freq,c(all.freq[1,1],"Binned",1-sum(as.numeric(all.freq[,3])),1-sum(as.numeric(all.freq[,4])))) # Binned alleles

## Assemble the whole table

all.table <- cbind(as.character(AZ.Locus$L$Set1$OR[,1]),as.character(AZ.Locus$L$Set1$OR[,2]),all.freq[,3],AZ.Locus$L$Set1$table[,3],
                   all.freq[,4],AZ.Locus$L$Set1$table[,4], as.character(AZ.Locus$L$Set1$OR[,3]),as.character(AZ.Locus$L$Set1$OR[,4]),
                   as.character(AZ.Locus$L$Set1$OR[,5]),as.character(AZ.Locus$L$Set1$OR[,6]),as.character(AZ.Locus$L$Set1$OR[,7]))

## Rename the columns for the frequency and counts 
colnames(all.table) <- c("Locus","Allele","Control (f)","Control (n)","Case (f)","Case (n)","Odds Ratio", "CI Lower", "CI Upper", "P-value","Significance")

View(all.table)

## Save a tab-delimited table of all values

write.table(x = all.table,file = "AZ_18IHIW_results.txt",append = FALSE,quote = FALSE,sep = "\t",row.names = FALSE)

### Stratified analysis -- you would primarily do this in cases you hav a single very strong association in the full data set 

## load the BDStrat function, which is included as a .rda file

load("BDStrat.rda")

## Run BDStrat using about=TRUE to get details and parameters

BDStrat(about=TRUE)

## Run BDstrat with the other parameters to split your data set into two subsets --
## one which contains all subjects with a specific allele (positive), and one that contains all subjects without that allele

AZ.DRB1_0301_stratified <- BDStrat(AZ.DRB1,"DRB1","03:01",TRUE)

View(AZ.DRB1_0301_stratified$`DRB1_1*03:01-positive`)
View(AZ.DRB1_0301_stratified$`DRB1_1*03:01-negative`)

## Run BIGDAWG on the 03:01-negative stratum
## IF YOU ARE NOT RUNINING BIGDWG, RUN LINE 129 BELOW.
AZ.DRB1_0301_negative <- BIGDAWG(Data = AZ.DRB1_0301_stratified$`DRB1_1*03:01-negative`,Run.Tests = "L", Return = TRUE,Output = FALSE)

load("AZ.DRB1_0301_negative.rda")

View(AZ.DRB1_0301_negative$L$Set1$OR)

## You can apply the approach as described above in lines 62-105 to examine and generate tables from the AZ.DRB1_0301_negative object
