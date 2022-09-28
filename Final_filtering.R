#library(monocle3)
#-- Put these filters on the UniRef Count table which is at least present in Two patients
uniref_counts <- read.csv(file="../Bowtie2_map_uniRefGeneCounts_Removed_twoPID.csv", sep = ",", row.names = 1) 

################################
#-- Step1: Removal of Samples
#--Sample Size Less then 100 MB
FileSize <- read.csv(file="../FileSize.txt", sep = "\t", header = T) # For all R1 and R2
FileSize_filtered <- subset(FileSize, SizeMB > 50)
SamplesSize <- unique(FileSize_filtered$FileName) # 1299 Samples have higher then 50 MB Size and 39 Samples lesser then 50 MB are removed
colnames(uniref_counts)
uniref_counts_filtered <- uniref_counts[SamplesSize]
#-- To remove samples which have zero UniRefs (11 Samples to be removed)
Sample_sum <- colSums(uniref_counts_filtered)
Sample_sum <- data.frame(Sample_sum)
row.names(Sample_sum) <- colnames(uniref_counts_filtered)
Samples_to_be_retained <- row.names (subset(Sample_sum, Sample_sum > 100000)) #These are the samples which have sum of UniRef Count more than 100k
uniref_counts_filtered <- uniref_counts_filtered[Samples_to_be_retained]
# So far samples which have low sequencing depth and Low total number reads assigend to UniRefs were removed 

################################
#-- Step2: Removal of UniRefs
unireftotal <- data.frame(rowSums(uniref_counts_filtered))
unireftoBeretained <- subset(unireftotal, rowSums.uniref_counts_filtered. > 20000) # UniRef count more than 20k
unireftoBeretained_rownames <- rownames(unireftoBeretained)
Filtered_uniRefTable <- uniref_counts_filtered[rownames(uniref_counts_filtered) %in% unireftoBeretained_rownames, ]  # Extract rows from data
write.csv(Filtered_uniRefTable, file="Filtered_uniRefTable_20k.csv")

#Check how much read left
TotalReads_UniRef <- data.frame(colSums(uniref_counts_filtered))
TotalReads_UniRef_afterFiltering <- data.frame(colSums(Filtered_uniRefTable))
Count_Stats <- merge(TotalReads_UniRef, TotalReads_UniRef_afterFiltering, by=0, all=F)
rownames(Count_Stats) <- Count_Stats$Row.names; Count_Stats$Row.names <- NULL
colnames(Count_Stats) <- c("totalMappedReads","AfterFiltering_20k")
Proportion_Left <- data.frame(Count_Stats$AfterFiltering_20k/Count_Stats$totalMappedReads)*100
rownames(Proportion_Left) <- rownames(Count_Stats)
Count_Stats_proportion <- merge(Count_Stats, Proportion_Left, by=0, all=F)
rownames(Count_Stats_proportion) <- Count_Stats_proportion$Row.names; Count_Stats_proportion$Row.names <- NULL
colnames(Count_Stats_proportion) <- c("totalMappedReads","AfterFiltering_20k","ProportionLeft")
write.csv(Count_Stats_proportion, file="Count_Stats_proportion_20kuniref.csv")

###################################################################
#-- Step3: Remove samples again wchich has very low proportion left (Less then 50%)
finalSamples_toRetain <- subset(Count_Stats_proportion, ProportionLeft >= 50)
finalSamples_toRetain_names <- rownames(finalSamples_toRetain)
Filtered_uniRefTable_samples <- Filtered_uniRefTable[finalSamples_toRetain_names] # Remove samples again
write.csv(Filtered_uniRefTable_samples, file="Final_filtered_UniRef20k_and_SamplesLowdepth_proportionLeft.csv")
save.image("Filtering.RData")
