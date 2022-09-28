#!/bin/bash

module load bowtie/2.3.2
module load trimmomatic/0.39

cd ~/Pouchitis_project/xia_shotgunData

~/.local/bin/kneaddata --input BX06_S57_L001_R1_001.fastq.gz --input BX06_S57_L001_R2_001.fastq.gz -db /home/sharmaa4/Databases/Mouse/mouse_C57BL_6NJ_Bowtie2 --output BX06_kneaddata_output --bypass-trf -t 12 --bowtie2-options="--very-fast" --bowtie2-options="-p 12"
