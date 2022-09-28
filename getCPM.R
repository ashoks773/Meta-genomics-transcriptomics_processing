#--- Load UniRef Counts
uniref_counts <- read.csv(file="Bowtie2_map_uniRefGeneCounts_Filter7.csv", sep = ",", row.names = 1)

#--- Convert Counts to CPM
countsToCPM <- function(counts) {
  
  cpm <- counts
  for (i in 1:(length(cpm))) {
    cpm[,i] <- (cpm[,i]/sum(cpm[,i]))*1000000
  }
  
  # Add the sum of CPMs at the last column of the data frame
  sum.cpm <- NULL
  for (i in 1:nrow(cpm)) {
    sum.cpm <- c(sum.cpm, sum(cpm[i,1:(length(cpm))]))
  }
  cpm$sum.cpm <- sum.cpm
  cpm
}
rawCPM <- countsToCPM(uniref_counts)
write.csv(rawCPM, file = "Bowtie2_map_uniRefGeneCounts_Filter7_CPM.csv")

