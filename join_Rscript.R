library("purrr")
library("dplyr")                                    # Load dplyr package
library("plyr")                                     # Load plyr package
library("readr") 

#************************************************* Step1: Filter low abundant UniRefs (Relab: 1e-06) in each sample
#-- Run this command inside the folder
#data_files <- list.files("~/IBD_datasets/HMP2/WGS_rawReads/WGS_MGX_aligned_geneFamilies/Chocophlan_Seperate_Mapping/Final_Counts/", pattern = "*_Final.csv")
#data_files  

#for(i in 1:length(data_files)) {
 # data <- read.csv(data_files[i])
 # cov <- data.frame(data[,2]/sum(data[,2]))
 # cov_data <- cbind(cov, data)
 # cov_data_filtered <- subset(cov_data, data...2..sum.data...2.. > 1e-06) # To remove gene families proportions less then 1e-06
 # cov_data_filtered <- cov_data_filtered[-1]
 # write.csv2(cov_data_filtered,paste0("~/IBD_datasets/HMP2/WGS_rawReads/WGS_MGX_aligned_geneFamilies/Chocophlan_Seperate_Mapping/Final_Counts/",data_files[i],"_Filtered.csv"),row.names = FALSE, sep=",")
#}

#************************************************* Step2: Format to replace ; with , in each file
# Run simple Perl script to format ; to , in the all Filtered.csv files
# for i in *Filtered.csv; do perl -i -pe "s/;/,/g" $i; done

#--NOTE: Above two steps were for the filering low abundant UniRef: These steps were not considered and Merging was redone using the steps provided below:

#************************************************* Step1: To format count files
# With in this folder: ~/IBD_datasets/HMP2/WGS_rawReads/WGS_MGX_aligned_geneFamilies/Chocophlan_Seperate_Mapping/Final_Counts_Formatted/
# Run format.sh which includes following code
#for i in *Final.csv; do perl -i -pe "s/^/\"/g" $i; done;
#for i in *Final.csv; do perl -i -pe "s/,/\",\"/g" $i; done;
#for i in *Final.csv; do perl -i -pe "s/$/\"/g" $i; done;
#for i in *Final.csv; do perl -i -pe "s/\"\"/\"/g" $i; done;

#************************************************* Step2: Join filtered uniRef Counts
#data_join <- list.files(path = "~/IBD_datasets/HMP2/WGS_rawReads/WGS_MGX_aligned_geneFamilies/Chocophlan_Seperate_Mapping/Final_Counts/", pattern = "*_Filtered.csv", full.names = TRUE) %>% lapply(read_csv) %>% reduce(full_join, by = "UniRefID")
data_join <- list.files(path = "~/IBD_datasets/HMP2/WGS_rawReads/WGS_MGX_aligned_geneFamilies/Chocophlan_Seperate_Mapping/Final_Counts_Formatted/", pattern = "*_Final.csv", full.names = TRUE) %>% lapply(read_csv) %>% reduce(full_join, by = "UniRefID")
data_join <- as.data.frame(data_join)
row.names(data_join) <- data_join$UniRefID; data_join$UniRefID <- NULL
data_join[is.na(data_join)] <- 0
write.csv(data_join, file = "Bowtie2_map_uniRefGeneCounts.csv")
save.image("Bowtie2_map_uniRefGeneCounts.RData")

#************************************************* Step3: Filter UniRef's based on their presence in number of samples
library (labdsv)
#load("Bowtie2_map_uniRefGeneCounts.RData")
#--- Remove ASVs which are not present in atleast 10% of the samples
data_join_t <- data.frame (t(data_join))
uniref_counts_filtered_Per10 <- dropspc(data_join_t, 134) #--- Uniref should be present in 10% of the samples
uniref_counts_filtered_Per10 <- data.frame(t(uniref_counts_filtered_Per10))
#write.csv(uniref_counts_filtered_Per10, file = "Bowtie2_map_uniRefGeneCounts_filtered_Per10.csv")
write.csv(uniref_counts_filtered_Per10, file = "Bowtie2_map_uniRefGeneCounts_Per10.csv")

uniref_counts_filtered_Per20 <- dropspc(data_join_t, 268) #--- Uniref should be present in 20% of the samples
uniref_counts_filtered_Per20 <- data.frame(t(uniref_counts_filtered_Per20))
#write.csv(uniref_counts_filtered_Per20, file = "Bowtie2_map_uniRefGeneCounts_filtered_Per20.csv")
write.csv(uniref_counts_filtered_Per20, file = "Bowtie2_map_uniRefGeneCounts_Per20.csv")

uniref_counts_filtered_Per25 <- dropspc(data_join_t, 335) #--- Uniref should be present in 25% of the samples
uniref_counts_filtered_Per25 <- data.frame(t(uniref_counts_filtered_Per25))
#write.csv(uniref_counts_filtered_Per25, file = "Bowtie2_map_uniRefGeneCounts_filtered_Per25.csv")
write.csv(uniref_counts_filtered_Per25, file = "Bowtie2_map_uniRefGeneCounts_Per25.csv")

uniref_counts_filtered_Per50 <- dropspc(data_join_t, 669) #--- Uniref should be present in 50% of the samples
uniref_counts_filtered_Per50 <- data.frame(t(uniref_counts_filtered_Per50))
#write.csv(uniref_counts_filtered_Per50, file = "Bowtie2_map_uniRefGeneCounts_filtered_Per50.csv")
write.csv(uniref_counts_filtered_Per50, file = "Bowtie2_map_uniRefGeneCounts_Per50.csv")

uniref_counts_filtered_Per75 <- dropspc(data_join_t, 1004) #--- Uniref should be present in 75% of the samples
uniref_counts_filtered_Per75 <- data.frame(t(uniref_counts_filtered_Per75))
#write.csv(uniref_counts_filtered_Per75, file = "Bowtie2_map_uniRefGeneCounts_filtered_Per75.csv")
write.csv(uniref_counts_filtered_Per75, file = "Bowtie2_map_uniRefGeneCounts_Per75.csv")

uniref_counts_filtered_Per80 <- dropspc(data_join_t, 1070) #--- Uniref should be present in 80% of the samples
uniref_counts_filtered_Per80 <- data.frame(t(uniref_counts_filtered_Per80))
#write.csv(uniref_counts_filtered_Per80, file = "Bowtie2_map_uniRefGeneCounts_filtered_Per80.csv")
write.csv(uniref_counts_filtered_Per80, file = "Bowtie2_map_uniRefGeneCounts_Per80.csv")

uniref_counts_filtered_Per90 <- dropspc(data_join_t, 1204) #--- Uniref should be present in 90% of the samples
uniref_counts_filtered_Per90 <- data.frame(t(uniref_counts_filtered_Per90))
#write.csv(uniref_counts_filtered_Per90, file = "Bowtie2_map_uniRefGeneCounts_filtered_Per90.csv")
write.csv(uniref_counts_filtered_Per90, file = "Bowtie2_map_uniRefGeneCounts_Per90.csv")

uniref_counts_filtered_Per95 <- dropspc(data_join_t, 1272) #--- Uniref should be present in 95% of the samples
uniref_counts_filtered_Per95 <- data.frame(t(uniref_counts_filtered_Per95))
#write.csv(uniref_counts_filtered_Per95, file = "Bowtie2_map_uniRefGeneCounts_filtered_Per95.csv")
write.csv(uniref_counts_filtered_Per95, file = "Bowtie2_map_uniRefGeneCounts_Per95.csv")

uniref_counts_filtered_Per99 <- dropspc(data_join_t, 1324) #--- Uniref should be present in 99% of the samples
uniref_counts_filtered_Per99 <- data.frame(t(uniref_counts_filtered_Per99))
#write.csv(uniref_counts_filtered_Per99, file = "Bowtie2_map_uniRefGeneCounts_filtered_Per99.csv")
write.csv(uniref_counts_filtered_Per99, file = "Bowtie2_map_uniRefGeneCounts_Per99.csv")

