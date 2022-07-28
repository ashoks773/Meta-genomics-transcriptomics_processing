#https://github.com/qiyunzhu/woltka/blob/master/doc/wol.md


#-- Step1: Build Reference Genome and Map your reads to the genome- Construct mapping files
----> https://github.com/qiyunzhu/woltka/blob/master/doc/align.md
db=/home/sharmaa4/Databases/Woltka_databases/wolr1_seqs/concat.fna.xz
mkdir -p databases/bowtie2
xz -d -k concat.fna.xz
bowtie2-build --threads 8 concat.fna databases/bowtie2/WoLr1 (Define your own location)
rm concat.fna

#--Map to the reference Genome
bowtie2 -p 8 -x db -1 R1.fq -2 R2.fq -S output.sam --very-sensitive --no-head --no-unal
bowtie2 -p 8 -x db -1 R1.fq -2 R2.fq --very-sensitive --no-head --no-unal | cut -f1-9 | sed 's/$/\t*\t*/' | gzip > output.sam.gz (Alternative Command to Save Space)

#-- Note: Alternatively if you would like to use some other references specifically GTDB (https://github.com/qiyunzhu/woltka/blob/master/doc/gtdb.md#use-gtdb-phylogeny). Donwload database, Parse taxonomy, and map your reads to this database insted of the contact.fna. Follow instructions from the link.

#-- Step2: Functional profiling form Mapped same files
---> https://github.com/qiyunzhu/woltka/blob/master/doc/ordinal.md
woltka classify -i indir -c coords.txt --overlap 50 -o gene.biom
# coords.txt file can be find in: /home/sharmaa4/Databases/Woltka_databases/wolr1_maps/proteins/coords.txt.xz

#--- Step3: Pathway Analysis using MetaCyc
#--/home/sharmaa4/Databases/Woltka_databases/wolr1_maps/function/metacyc/*
woltka tools collapse -i per-gene.biom -m function/metacyc/protein.map.xz -n function/metacyc/protein_name.txt -o protein.biom
woltka tools collapse -i protein.biom -m function/metacyc/protein-to-enzrxn.txt -n function/metacyc/enzrxn_name.txt -o enzrxn.biom
woltka tools collapse -i enzrxn.biom -m function/metacyc/enzrxn-to-reaction.txt -n function/metacyc/reaction_name.txt -o reaction.biom
woltka tools collapse -i reaction.biom -m function/metacyc/reaction-to-pathway.txt -n function/metacyc/pathway_name.txt -o pathway.biom

From the Sam files (Step1) Or Gene.biom file (Step2) follow the following steps:
https://github.com/qiyunzhu/woltka/blob/master/doc/wol.md
https://github.com/qiyunzhu/woltka/blob/master/doc/ordinal.md

