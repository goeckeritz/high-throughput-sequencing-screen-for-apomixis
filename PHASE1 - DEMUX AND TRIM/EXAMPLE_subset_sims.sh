#!/bin/bash --login
#SBATCH -J subsample_fastqs
#SBATCH --mail-user=cgoeckeritz@hudsonalpha.org
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=64g
#SBATCH --time=24:00:00
#SBATCH -o /cluster/home/cgoeckeritz/mpp1-10/output_files/subset_1X_3X_%j
#SBATCH -a 1-40
#SBATCH --export=INFILE=/cluster/home/cgoeckeritz/mpp1-10/moms_R1_paired.txt

READ1=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`
echo ${READ1}

READ2=$(echo ${READ1}| sed "s/_R1_/_R2_/")
echo ${READ2}

GENO=$(echo ${READ1} | sed 's/\/cluster.*\///' | sed 's/_R1.*//')
echo ${GENO}

cd /cluster/home/cgoeckeritz/mpp1-10/

#number of reads you'd like to sample from the files. Choose appropriate amount based on your read length and genome size
#to get an approximate amount of coverage.
X1=2500000
X3=7500000

#1X simulated embryos
/cluster/home/cgoeckeritz/software/seqtk/seqtk sample -s65 ${READ1} ${X1} | gzip > ${GENO}_sim_e1_R1.fastq.gz
/cluster/home/cgoeckeritz/software/seqtk/seqtk sample -s65 ${READ2} ${X1} | gzip > ${GENO}_sim_e1_R2.fastq.gz

/cluster/home/cgoeckeritz/software/seqtk/seqtk sample -s93 ${READ1} ${X1} | gzip > ${GENO}_sim_e2_R1.fastq.gz
/cluster/home/cgoeckeritz/software/seqtk/seqtk sample -s93 ${READ2} ${X1} | gzip > ${GENO}_sim_e2_R2.fastq.gz

/cluster/home/cgoeckeritz/software/seqtk/seqtk sample -s14 ${READ1} ${X1} | gzip > ${GENO}_sim_e3_R1.fastq.gz
/cluster/home/cgoeckeritz/software/seqtk/seqtk sample -s14 ${READ2} ${X1} | gzip > ${GENO}_sim_e3_R2.fastq.gz

#3X simulated embryos
/cluster/home/cgoeckeritz/software/seqtk/seqtk sample -s34 ${READ1} ${X3} | gzip > ${GENO}_sim_e4_R1.fastq.gz
/cluster/home/cgoeckeritz/software/seqtk/seqtk sample -s34 ${READ2} ${X3} | gzip > ${GENO}_sim_e4_R2.fastq.gz

/cluster/home/cgoeckeritz/software/seqtk/seqtk sample -s69 ${READ1} ${X3} | gzip > ${GENO}_sim_e5_R1.fastq.gz
/cluster/home/cgoeckeritz/software/seqtk/seqtk sample -s69 ${READ2} ${X3} | gzip > ${GENO}_sim_e5_R2.fastq.gz

/cluster/home/cgoeckeritz/software/seqtk/seqtk sample -s18 ${READ1} ${X3} | gzip > ${GENO}_sim_e6_R1.fastq.gz
/cluster/home/cgoeckeritz/software/seqtk/seqtk sample -s18 ${READ2} ${X3} | gzip > ${GENO}_sim_e6_R2.fastq.gz

#A COUPLE NOTES:
#/cluster/home/cgoeckeritz/mpp1-10/moms_R1_paired.txt
#is the full path to a list of R1 files for each maternal genotype's sequencing data, with one /path/file.fastq.gz per line. In the present work, 
#these files had already been run through fastp (only keeping paired data, trimming low quality bases and removing PCR duplicates),
#but you can subset on the raw files too. If you're data is good quality and there aren't too many PCR duplicates, 
#this shouldn't majorly affect results; and I created these 'simulated embryo files' only to get an idea of how accurate the caller was 
#at certain depths.

#note that I used a personal download of seqtk. See https://github.com/lh3/seqtk for more info. 

 