#!/bin/bash

cd ~/Pouchitis_project/xia_shotgunData/Analysis
##-- Final quality filtered, Host Removed Files are Stored in this Folder
#ls *fastq | cut -f 1 -d '.' | cut -d"_" -f1,2 > libraryname.txt
#cat libraryname.txt | while read file; do cat "$file"_host_removed_r1.fastq "$file"_host_removed_r2.fastq > "$file"_merged.fastq; done;

##--- Gzip all fastq files
#gzip *fastq
#mkdir Final_filtered_FASTQ
#mv host_removed_r*.fastq.gz Final_filtered_FASTQ #move both files in to the seperate folder

##https://github.com/biobakery/biobakery/wiki/metaphlan3
##--- Run Metaphln3 - make a bash Script to do this
#source /hpc/apps/miniconda/3/etc/profile.d/conda.sh
#conda activate biobakery3
#module load bowtie/2.3.2
#for i in *_merged.fastq.gz; do metaphlan $i --input_type fastq --nproc 4 > ${i%.fasta.gz}_profile.txt; done;
mkdir Metaphln_Outputs
mv *_profile.txt Metaphln_Outputs 
merge_metaphlan_tables.py *_profile.txt > merged_abundance_table.txt
##Step 1: Generate the species only abundance table
grep -E "s__|clade" merged_abundance_table.txt | sed 's/^.*s__//g'| cut -f1,3-18 > merged_abundance_table_species.txt

#--- Run Humann3 - make a bash Script to do this
#https://github.com/biobakery/humann
source /hpc/apps/miniconda/3/etc/profile.d/conda.sh
conda activate biobakery3
module load bowtie/2.3.2
for i in *_merged.fastq.gz; do humann -i "$i" --input-format fastq.gz -o humann_results --nucleotide-database /home/sharmaa4/Databases/human_database/chocophlan.v296_201901b --protein-database /home/sharmaa4/Databases/human_database/uniref --threads 12

#-- Change Folder Name "hmp_subset"
humann_join_tables -i hmp_subset -o Pouchitis_Fun_genefamilies.tsv --file_name genefamilies
humann_renorm_table -i Pouchitis_Fun_genefamilies.tsv -o Pouchitis_Fun_genefamilies-cpm.tsv --units cpm
#-- These Stratified or Unstratified tables can also be generate in Maaslin2 - in R
humann_split_stratified_table --input Pouchitis_Fun_pathabundance-cpm.tsv --output path_abun_stratified
 
