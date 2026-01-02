#!/bin/sh --login
#SBATCH -J align,metrics,call
#SBATCH --mail-user=cgoeckeritz@hudsonalpha.org
#SBATCH --nodes=1
#SBATCH --ntasks=8
#SBATCH --mem=32g
#SBATCH --time=24:00:00
#SBATCH -o /cluster/home/cgoeckeritz/mpp1-10/output_files/embryos_fastp_8.24.2025_%j
#SBATCH -a [1-952]%100
#SBATCH --export=INFILE=/cluster/lab/harkess/Apomixis/Malus/embryo_screening/mpp1-10_WGS/mpp1-10_R1_files.txt,INFILE2=/cluster/lab/harkess/Apomixis/Malus/embryo_screening/mpp1-10_WGS/mpp1-10_GENOs.txt

#We use this script to polish the read data but also as an opportunity to 
#clean up raw file names. For example -- the first INFILE is a full path to the R1 file for all samples, e.g. a line might read:

#/cluster/lab/harkess/Apomixis/Malus/embryo_screening/mpp1-10_WGS/mpp1/MGoldgelb_589458_e1-MGoldgelb_589458_e1-GATCAGCT_R1.fastq.gz
#And on the same line number (very important line numbers correspond) in INFILE2, it would read: MGoldgelb_589458_e1 

#There's a lot of ways to create these files. I found it easiest to go to the directory where the files were, then do:
#ls -1 /cluster/lab/harkess/Apomixis/Malus/embryo_screening/mpp1-10_WGS/mpp1/M*R1.fastq.gz >> mpp1-10_R1_files.txt

#to make sure the files would perfectly correspond with their prefixes in INFILE2, I then took mpp1-10_R1_files.txt and applied a regular expression to create INFILE2.
#cat mpp1-10_R1_files.txt | sed 's/\/cluster.*\///' | sed 's/-.*//' > mpp1-10_GENOs.txt

#note - if _a, _b, _c, or _d was not already part of a spike-in name, I had to add it. This is mostly a note for me. :)

#establish a working directory where all the analyses of the pipeline will take place.
WD=/cluster/home/cgoeckeritz/mpp1-10/
cd ${WD}

READ1=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`
READ2=$(echo ${READ1} | sed 's/_R1/_R2/')
echo ${READ1}
echo ${READ2}

GENO=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE2}` #This prefix will be on every file name in the next phase of the pipeline too!
echo ${GENO}

ml purge
module load fastp/0.23.4-gcc-13.1.0

echo "Let's remove the adapters and PCR duplicates"
fastp \
-i ${READ1} \
-o ${GENO}_R1_paired_no_dups.fastq.gz \
-I ${READ2} \
-O ${GENO}_R2_paired_no_dups.fastq.gz \
-D \
--detect_adapter_for_pe \
-u 30 \
-w 8

#see https://github.com/OpenGene/fastp for more info on fastp documentation.

