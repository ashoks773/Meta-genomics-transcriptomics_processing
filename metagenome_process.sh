#!/bin/bash

######################################################################################################
#
# Master Script to perform quality control on metagenomic sequencing data
# Ashok K. Sharma
# This program will run KneadData and then Humann3/MetaPhyln3 and Finlly perform Assembly and Binning
#####################################################################################################

# Set starting files location
cd ~/IBD_datasets/HMP2/WGS_rawReads

# Step1: Quality control removal of reads
# Default trimmomatic MINLEN:60 ILLUMINACLIP:/-SE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50
module load bowtie/2.3.2
module load trimmomatic/0.39

mkdir ~/IBD_datasets/HMP2/step1_WGS_rawReads_QC
for file in *_1.fastq.gz;
do
  SAMPLE=$(echo ${file} | sed "s/_1.fastq.gz//")
  echo ${SAMPLE}_1.fastq.gz ${SAMPLE}_2.fastq.gz

~/.local/bin/kneaddata --input ${SAMPLE}_1.fastq.gz --input ${SAMPLE}_2.fastq.gz -db /home/sharmaa4/Databases/Human_Genome/human_GRCh38.p13 --output /home/sharmaa4/IBD_datasets/HMP2/step1_WGS_rawReads_QC/${SAMPLE}_kneaddata_output --bypass-trf -t 12 --bowtie2-options="--very-fast" --bowtie2-options="-p 12" --bowtie2-options="--no-discordant"

done

# Step2: Run Humann3 on quality control Reads (Input files will be KneadData output files)
#source /hpc/apps/miniconda/3/etc/profile.d/conda.sh
#conda activate biobakery3
#module load bowtie/2.3.2
#mkdir ~/IBD_datasets/HMP2/step2_humann_results
#cd ~/IBD_datasets/HMP2/step1_WGS_rawReads_QC

#for file in *_1_kneaddata_paired_1.fastq;
#do
#   SAMPLE=$(echo ${file} | sed "s/_1_kneaddata_paired_1.fastq//")
#   #-- Four output files needs to be concatenated
#   echo ${SAMPLE}_1_kneaddata_paired_1.fastq ${SAMPLE}_1_kneaddata_paired_2.fastq ${SAMPLE}_1_kneaddata_unmatched_1.fastq ${SAMPLE}_1_kneaddata_unmatched_2.fastq
#   cat ${SAMPLE}_1_kneaddata_paired_1.fastq ${SAMPLE}_1_kneaddata_paired_2.fastq ${SAMPLE}_1_kneaddata_unmatched_1.fastq ${SAMPLE}_1_kneaddata_unmatched_2.fastq > ${SAMPLE}_merged.fastq.gz

#humann -i ${SAMPLE}_merged.fastq.gz -o /home/sharmaa4/IBD_datasets/HMP2/step2_humann_results --nucleotide-database /home/sharmaa4/Databases/human_database/chocophlan.v296_201901 --protein-database /home/sharmaa4/Databases/human_database/uniref

#done

# Step3: Run Megahit to assemble high qulaity reads

