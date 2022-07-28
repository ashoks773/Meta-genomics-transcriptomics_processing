# This is a step file to process Metagenomics or Transcriptomics data to get GeneFamilies Counts
# By: Ashok K. Sharma; Date: July 27, 2022

############################################################## Step1 ######################################################
#Step1: Is to map MTX reads to the Chocophln + UniRef Database and get *bowtie2_aligned.tsv and *diamond_aligned.tsv Files
cd ~/IBD_datasets/HMP2/MTX_rawdata
cat fastqfiles | while read line; do qsub -q 1tb.q,all.q,256gb.q,matlab.q,cgroup.q -cwd -o $PWD -e $PWD -I mem_free=16G -pe mpich 4 ./senhuman.singlelibrary $line; done
cat pilot_bz2_fastqfiles | while read line; do qsub -q 1tb.q,all.q,256gb.q,matlab.q,cgroup.q -cwd -o $PWD -e $PWD -l mem_free=16G -pe mpich 4 ./bz2_sendhumann.singlelibrary $line; done
#**************** Code for sendcountjob file
#!/bin/bash
source /hpc/apps/miniconda/3/etc/profile.d/conda.sh 
conda activate biobakery3 
module load bowtie/2.3.2 
humann -i /home/sharmaa4/common/MTX/$1 --input-format fastq.gz -o humann_results --nucleotide-database /home/sharmaa4/Databases/human_database/chocophlan.v296
_201901b --protein-database /home/sharmaa4/Databases/human_database/uniref --threads 4 --resume
#*************** End

############################################################## Step2 ######################################################
#Step2: Is to concatenate these files and get counts for each sample
cd ~/IBD_datasets/HMP2/MTX_aligned_geneFamilies
cat filenames | while read line; do qsub -q 1tb.q,all.q,256gb.q,matlab.q,cgroup.q -cwd -o $PWD -e $PWD -l mem_free=16G -pe mpich 2 ./sendcountjob $line; done
#*************** Code for sendcountjob file
#!/bin/bash
cut -f2 ../MTX_rawdata/humann_results/"$1"_humann_temp/"$1"_bowtie2_aligned.tsv | cut -d"|" -f3 > "$1"_bowtie_diamond.txt #-- Nucleotide search
cut -f2 ../MTX_rawdata/humann_results/"$1"_humann_temp/"$1"_diamond_aligned.tsv | cut -d"|" -f1 >> "$1"_bowtie_diamond.txt #-- + Tanslated search of Unmapped reads in step 1
sort "$1"_bowtie_diamond.txt | uniq -c > "$1"_counts.txt
sed -r 's/^\s+//g' "$1"_counts.txt > "$1"_counts_format.txt
cut -d" " -f2 "$1"_counts_format.txt > "$1"_uniRefid.txt
cut -d" " -f1 "$1"_counts_format.txt > "$1"_uniRefCount.txt
paste -d"," "$1"_uniRefid.txt "$1"_uniRefCount.txt > "$1"_counts.csv
echo -e "UniRefID,$1" > "$1"_counts_Final.csv && cat "$1"_counts.csv >> "$1"_counts_Final.csv
#---Save *Final.csv in the folder Aligned_counts/Filtered_counts
#-- Save rest of the files in folder Aligned_counts/Raw_counts
#*************** End
# In step2: there is one more step to get the read alignment statistics for all samples
#-- To Get Alignment %
cd ~/IBD_datasets/HMP2/MTX_aligned_geneFamilies/aligned_per_summmary
cat filenames | while read line; do qsub -q debug.q,gpu.q,pbl.q,ragatkolab.q -cwd -o $PWD -e $PWD -l mem_free=16G -pe mpich 2 ./sendalignedJob $line; done

############################################################## Step3 ######################################################
#Step3: Filter UniRefID based on Counts (number is variable for each sample). So to be consistent with each sample use propotion (UniRefcount/total count in the sample). And finally remove UniRef whose relative proportion is less then 1e-06. Trying to write a R script for this.
#- This step was run on Local computer instead of Cluster
setwd("~/Box/David_Casero_Lab/PseudoTime_HMP_project/Aligned_counts/")
#Copy Aligned Raw Counts
#scp -r sharmaa4@csclprd3-s001v:/home/sharmaa4/IBD_datasets/HMP2/MTX_aligned_geneFamilies/Aligned_counts/Raw_counts/*Final.csv .

data_files <- list.files("~/Box/David_Casero_Lab/PseudoTime_HMP_project/Aligned_counts/", pattern = "*_Final.csv")
#data_files  

for(i in 1:length(data_files)) {
  data <- read.csv(data_files[i])
  cov <- data.frame(data[,2]/sum(data[,2]))
  cov_data <- cbind(cov, data)
  cov_data_filtered <- subset(cov_data, data...2..sum.data...2.. > 1e-06) # To remove gene families proportions less then 1e-06
  cov_data_filtered <- cov_data_filtered[-1]
  write.csv2(cov_data_filtered,paste0("~/Box/David_Casero_Lab/PseudoTime_HMP_project/Aligned_counts/",data_files[i],"_Filtered.csv"),row.names = FALSE, sep=",")
}
# Output
#"UniRefID";"CSM67UBN"
#"UniRef90_A0A010YBU9";25
# Use sed command to remove "" and replace ; with ,
# Trnasfer back to the Cluster to Join these files
#scp -r *Filtered.csv sharmaa4@csclprd3-s001v:/home/sharmaa4/IBD_datasets/HMP2/MTX_aligned_geneFamilies/Aligned_counts/Filtered_counts
#for i in *Filtered.csv; do sed 's/;/,/g' "$i" > "$i"_format.csv; done

############################################################## Step4 ######################################################
# Step4: Join filtered count files using a Rscript
cd ~/IBD_datasets/HMP2/MTX_aligned_geneFamilies
qsub -q 1tb.q -cwd -o $PWD -e $PWD -l mem_free=1000G -pe mpich 24 runJoin.sh
#********** Content of runR.sh script
#!/bin/bash
cd ~/IBD_datasets/HMP2/MTX_aligned_geneFamilies
module load R
Rscript join_Rscript.R
#********** Content of join_Rscript.R
#data_join <- list.files(path = "~/IBD_datasets/HMP2/MTX_aligned_geneFamilies/Aligned_counts/Raw_counts/", pattern = "*_Final.csv", full.names = TRUE) %>% lapply(read_csv) %>% reduce(full_join, by = "UniRefID")
data_join <- list.files(path = "~/IBD_datasets/HMP2/MTX_aligned_geneFamilies/Aligned_counts/Filtered_counts/", pattern = "*_format.csv", full.names = TRUE) %>% lapply(read_csv) %>% reduce(full_join, by = "UniRefID")
data_join <- as.data.frame(data_join)
row.names(data_join) <- data_join$UniRefID; data_join$UniRefID <- NULL
data_join[is.na(data_join)] <- 0
#write.table(data_join, file = "Combined_MTX_UniRef_Counts.txt", sep="\t")
#save.image("unirefCount.RData")
write.csv(data_join, file = "Combined_MTX_UniRef_Counts_Filtered.csv")
save.image("MTX_UniRefCounts_filtered.RData")
