#--- Load UniRef Counts
uniref_counts <- read.csv(file="Bowtie2_map_uniRefGeneCounts.csv", sep = ",", row.names = 1)

#---> Statistics Check
#uniref_counts_t <- data.frame (t(uniref_counts))
#sum(uniref_counts == 0)
#sum(uniref_counts == 0)/nrow(uniref_counts)*ncol(uniref_counts)*100
#gene.pres <- apply(uniref_counts>0,2,sum)
#jpeg("GeneCounts_Stats.jpg", height = 4, width = 4, units = 'in', res = 600)
#hist(gene.pres, las=1, xlab = "Number of genes", ylab="Number of samples")
#dev.off

#--Number of Zeros in each Sample of data frame --> UniRef genes abusent in the sample
numberOfZeros <- data.frame(colSums(uniref_counts==0))
#--Number of NonZeros in each Sample of data frame --> UniRef genes present in the sample
numberOfNonZeros <- data.frame(colSums(uniref_counts!=0))
#-- Total UniRefCounts- means total reads mapped to UniRef in each sample
Total_unirefCounts <- data.frame(colSums(uniref_counts))

#-- Merge and Save
Count_Stats <- cbind(numberOfZeros, numberOfNonZeros, Total_unirefCounts)
colnames(Count_Stats) <- c("Genes_Absent(Zero)","Genes_Present(NonZero)","TotalUnirefCount")
write.csv(Count_Stats, file = "GeneCounts_Stats.csv")

#-- UniRef present in Each Sample
UnirefPresentIn_EachSample <- data.frame (rowSums(uniref_counts != 0))
write.csv(UnirefPresentIn_EachSample, file = "UnirefPresentIn_EachSample.csv")

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
write.csv(rawCPM, file = "CPM_Bowtie2_map_uniRefGeneCounts.csv")
save.image("Gene_CountStats.RData")
