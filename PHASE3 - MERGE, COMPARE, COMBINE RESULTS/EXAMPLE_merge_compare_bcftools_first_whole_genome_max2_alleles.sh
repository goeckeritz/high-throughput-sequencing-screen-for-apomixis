#!/bin/sh --login
#SBATCH -J merge_compare
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=128g
#SBATCH --time=48:00:00
#SBATCH -a 1-40
#SBATCH -o /cluster/home/cgoeckeritz/mpp1-10/output_files/merge_compared_max2_alleles_%j
#SBATCH --export=INFILE=/cluster/home/cgoeckeritz/mpp1-10/mpp1-10_mom_list.txt

# FYI this is submitted as an array job!
# The INFILE for this script is simply a list of the maternal genotypes, which should be the prefix or part of the prefix for all files associated with
# that genotype (her embryos, spikeins, etc.). This should work if you've kept your prefixes consistent throughout the phases of the pipeline as described previously...
# but I'll try to help you if it doesn't work :) 

module purge
ml bcftools/1.19-gcc-13.1.0

WD=/cluster/home/cgoeckeritz/mpp1-10/ #the master working directory
cd ${WD}
cd concordant_snps/

MOM=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`
echo ${MOM}

### Subset the UNpruned mom files, create a vcf of het and hom sites for each maternal genotype.
echo "Subset the maternal genotype's vcf into sites where she is heterozygous and sites where she is homozygous'"

bcftools filter -Oz ${MOM}_no_dups_merged_all_concordance_bcftools_first.vcf.gz -i 'GT="het"' > ${MOM}_no_dups_merged_all_concordance_bcftools_first_het.vcf.gz
bcftools filter -Oz ${MOM}_no_dups_merged_all_concordance_bcftools_first.vcf.gz -i 'GT="hom"' > ${MOM}_no_dups_merged_all_concordance_bcftools_first_hom.vcf.gz
bcftools index -f --threads 32 ${MOM}_no_dups_merged_all_concordance_bcftools_first_het.vcf.gz
bcftools index -f --threads 32 ${MOM}_no_dups_merged_all_concordance_bcftools_first_hom.vcf.gz

### Making some important lists
echo "Making some file lists prior to merging all the files in the list."

ls -1 ${MOM}*_no_dups_alts_concordance_bcftools_first_pruned.vcf.gz | grep '_[a,b,c,e]' > ${MOM}_embryo_files_concordance.txt #based on our naming system, the grep grabs all spikein, embryo, and simulated/subsetted embryo files and puts them in a list. 

# this takes the name of the maternal genotype file with only het sites in it and places it in the list with all its associated files. 
# this list is used to compare sites between mother and embryo where mother was heterozygous (x axis for PHASE4 - Visualization) 
ls -1 ${MOM}_no_dups_merged_all_concordance_bcftools_first_het.vcf.gz > ${MOM}_het_files.txt
cat ${MOM}_embryo_files_concordance.txt >> ${MOM}_het_files.txt

# this takes the name of the maternal genotype file with only hom sites in it (whether they be hom ref or hom alt) and places it in the list with all its associated files. 
# this list is used to compare sites between mother and embryo where mother was homozygous (y axis for PHASE4 - Visualization) 
ls -1 ${MOM}_no_dups_merged_all_concordance_bcftools_first_hom.vcf.gz > ${MOM}_hom_files.txt
cat ${MOM}_embryo_files_concordance.txt >> ${MOM}_hom_files.txt

echo "Merging file lists have been made. On to the actual merging."

#merge the two types of files for each case. Only look at true biallelic sites.
bcftools merge -Oz --threads 32 --file-list ${MOM}_het_files.txt | bcftools view --threads 32 --with-header -Oz -M 2 > ${MOM}_merged_vcfs_het.vcf.gz
bcftools merge -Oz --threads 32 --file-list ${MOM}_hom_files.txt |  bcftools view --threads 32 --with-header -Oz -M 2 > ${MOM}_merged_vcfs_hom.vcf.gz

#index the two combined files. 
for i in $(ls -1 ${MOM}_merged_vcfs*.vcf.gz); do
bcftools index -f --threads 32 $i
done

echo "Making the comparisons file for each maternal genotype."
# Lastly, we have to give bcftools gtcheck a file that tells us exactly what genotype pairs to compare in our combined vcf files.  
# To make this comparisons file, there should be a $MOM_$EMBRYO_comparison.txt file from the previous script for this pipeline in this directory already. Use this command to combine them all
# awk '$1 != $2' gets rid of lines in the resulting file where it's a self comparison. This happened for a separate version of the pipeline back
# when embryos and moms were called in the same PHASE2 script. I ended up separating them, but left the awk '$1 != $2' since it doesn't hurt anything. 
cat ${MOM}*_comparison.txt | awk '$1 != $2' > ${MOM}_comparisons.txt


#we're ready to make the comparisons!
echo "Comparing mother and child genotype by site type."
for i in $(ls -1 ${MOM}*_merged_vcfs*.vcf.gz); do
bcftools gtcheck -E 0 -u GT,GT --pairs-file ${MOM}_comparisons.txt $i > ${MOM}_bcftools_first_$(echo $i | grep -Eo "het|hom" )_comparison_results_GTvsGT.txt
done

# once all array jobs have completed, go to the directory where all *_comparison_results_GTvsGT.txt files for all gentypes are, then activate this script: 
# bash combine_results_extra_columns.sh.
# the script is currently written to require this file format for all genotypes: ${MOM}_bcftools_first_$(echo $i | grep -Eo "het|hom" )_comparison_results_GTvsGT.txt
# you could have chatGPT alter the script if you want to use a different prefix. 

# note that combine_results_extra_columns.sh makes a bunch of intermediate files for each separate genotype too, including a *het_results.tsv, *hom_results.tsv and *results_combined.tsv. I don't typically use these individual ones. 
# I really only move forward with the super combined file: all_results_combined.tsv