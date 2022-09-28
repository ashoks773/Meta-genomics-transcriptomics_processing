#!/bin/bash

#--- This script is to run Trimmomatic - for quality filtering and to remove Host contamination
# Step1: Quality control removal of reads SLIDINGWINDOW:5:20 MINLEN:80
#module load trimmomatic/0.39

#cd ~/Pouchitis_project/xia_shotgunData
#mkdir ~/Pouchitis_project/xia_shotgunData/step1_WGS_rawReads_QC

#--- Run trimmomatic for quality filtering
#for file in *L001_R1_001.fastq.gz;
#do
#  SAMPLE=$(echo ${file} | sed "s/_L001_R1_001.fastq.gz//")
#  echo ${SAMPLE}_L001_R1_001.fastq.gz ${SAMPLE}_L001_R2_001.fastq.gz
#-- Run trimmomatic and get four output files for Paired samples.
#trimmomatic PE -phred33 ${SAMPLE}_L001_R1_001.fastq.gz ${SAMPLE}_L001_R2_001.fastq.gz step1_WGS_rawReads_QC/${SAMPLE}_R1_PE.fastq step1_WGS_rawReads_QC/${SAMPLE}_R1_UN.fastq step1_WGS_rawReads_QC/${SAMPLE}_R2_PE.fastq step1_WGS_rawReads_QC/${SAMPLE}_R2_UN.fastq ILLUMINACLIP:$TRIMMOMATIC/adapters/TruSeq2-PE.fa:2:40:15 LEADING:10 TRAILING:10 SLIDINGWINDOW:5:20 MINLEN:80 -threads 12
#--Make a concatenated file of survived reads but their mates but not survived
#cat step1_WGS_rawReads_QC/${SAMPLE}_R1_UN.fastq step1_WGS_rawReads_QC/${SAMPLE}_R2_UN.fastq > step1_WGS_rawReads_QC/${SAMPLE}_Merged_UN.fastq
#done

#--- Run Bowtie2 to remove the host contamination Reads
#--- https://www.aespindola.com/post/dehosting/
module load bowtie/2.3.2
module load samtools
 
cd ~/Pouchitis_project/xia_shotgunData/step1_WGS_rawReads_QC
for file in *R1_PE.fastq;
do
  SAMPLE=$(echo ${file} | sed "s/_R1_PE.fastq//")
  echo ${SAMPLE}_R1_PE.fastq ${SAMPLE}_R2_PE.fastq
#bowtie2 -x ~/Databases/Mouse/mouse_C57BL_6NJ_Bowtie2 -1 ${SAMPLE}_R1_PE.fastq -2 ${SAMPLE}_R2_PE.fastq -U ${SAMPLE}_Merged_UN.fastq -S ${SAMPLE}.sam
#samtools view -bS ${SAMPLE}.sam > ${SAMPLE}_mapped_and_unmapped.bam
#samtools view -b -f 12 -F 256 ${SAMPLE}_mapped_and_unmapped.bam > ${SAMPLE}_bothEndsUnmapped.bam
#samtools sort -n ${SAMPLE}_bothEndsUnmapped.bam -o ${SAMPLE}_sort.bam
#module load bedtools
bamToFastq -i ${SAMPLE}_sort.bam -fq ${SAMPLE}_host_removed_r1.fastq -fq2 ${SAMPLE}_host_removed_r2.fastq
#bedtools bamtofastq -i ${SAMPLE}_sort.bam -fq ${SAMPLE}_host_removed_r1.fastq -fq2 ${SAMPLE}_host_removed_r2.fastq
#rm *sam
#rm *mapped.bam
done




