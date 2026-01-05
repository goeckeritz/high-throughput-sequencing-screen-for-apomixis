#!/bin/sh --login
#SBATCH -J align,metrics,call
#SBATCH --mail-user=cgoeckeritz@hudsonalpha.org
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=128g
#SBATCH --time=96:00:00
#SBATCH -o /cluster/home/cgoeckeritz/mpp1-10/output_files/moms_batch_8.24.2025_%j
#SBATCH -a 1-40
#SBATCH --export=INFILE=/cluster/home/cgoeckeritz/mpp1-10/moms_R1_paired.txt

#The working directory where ALL your subdirectories will be created and the analyses will occur.
WD=/cluster/home/cgoeckeritz/mpp1-10/
cd ${WD}
#Directory structure is important for the operation of this script and finding fles!

#INFILE is the full path to the trimmed R1 reads (already processed with fastp), one file per line. 

READ1=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}`
READ2=$(echo ${READ1} | sed 's/_R1_/_R2_/') #make sure R2 is in the same directory as R1. 
echo ${READ1}
echo ${READ2}

GENO=$(echo ${READ1} | sed 's/\/cluster.*\///' | sed 's/_R1.*//') #This removes the path and suffix so that /cluster/home/cgoeckeritz/mpp1-10/MGoldgelb_589458_R1_paired_no_dups.fastq.gz simply returns MGoldgelb_589458
echo ${GENO}

#report ploidies for all mom genotypes (this is only a subset, you can add more if you are doing more). 
#These variables should be exactly what $GENO is for your maternal genotype. 
#The ploidy is considered when setting filters for variant calling. 
MGoldgelb_589458=3
MHenryFDupont_589732=2
MMcClintockGrimes_589124=3
MRalphShay_589734=2
MRobinson_589455=2
MYellowAutumnCrabapple_588922=2
Mcoronaria_589344=4
Mcoronaria_589977=4
Mcoronaria_589988=3
Mcoronaria_590000=3
Mcoronaria_590014=2
Mhalliana_589013=3
Mhartwigii_588757=2
Mhybrid_589161=2
Mhybrid_589421=2
Mioensis_590008=2
Mkansuensis_588944=2
Mmagdeburgensis_588959=2
Mmicromalus_588976=2
Mmicromalus_594092=3
Mplatycarpa_588752=4
Mplatycarpa_588847=4
Mplatycarpa_589198=3
Mplatycarpa_589356=4
Mplatycarpa_589415=4
Mprunifolia_588914=2
Mrobusta_588825=2
Mrobusta_589383=2
Mrockii_589279=2
Msargentii_589372=4
Msikkimensis_589599=4
Msikkimensis_589750=3
Msikkimensis_613912=3
Mtoringo_590101=2
Mtoringo_613932=3
Mtoringo_613945=3
Mtoringo_633814=3
Mtransitoria_633806=3
Mxzumi_589840=2
Mzhaojiaoensis_633816=2

#mom's ploidy
echo ${!GENO}


REFERENCE=/cluster/home/cgoeckeritz/honeycrisp_full_genome/Malus_x_domestica_Honeycrisp_HAP1_v1.1.a1_scaffolded.fasta

SEQ_PLATFORM=ILLUMINA 

###### align the reads to the reference genome
module purge
ml cluster/bwa/0.7.17 
ml samtools/1.19.2-gcc-13.1.0
ml gcc/13.1.0
ml python/3.12.1-gcc-13.1.0
ml bcftools/1.19-gcc-13.1.0

cd ${WD}
mkdir -p alignment
cd alignment
mkdir -p flagstat

#align the reads and sort the alignment file.
bwa -a bwtsw index ${REFERENCE} ### can comment if the REFERENCE has been indexed.
echo "Time to align trimmed reads"
bwa mem -t 64 -R "@RG\tID:apple\tSM:${GENO}\tPL:${SEQ_PLATFORM}\tLB:${GENO}_1" ${REFERENCE} ${READ1} ${READ2} > ${GENO}_no_dups.bam
samtools sort -@ 64 ${GENO}_no_dups.bam > ${GENO}_no_dups.sorted.bam
samtools index -@ 64 ${GENO}_no_dups.sorted.bam
#cleanup a little...
rm ${GENO}_no_dups.bam

#Get total reads in our unfiltered alignment file - it will be used to calculate coverage of the genome.
samtools flagstat ${GENO}_no_dups.sorted.bam -@ 64 -O tsv > flagstat/${GENO}_no_dups_flagstat.tsv

#filter for mapping quality (q > 30), no secondary alignments, and sort and index the filtered file.
samtools view -b -h -q 30 -F 0x0100 ${GENO}_no_dups.sorted.bam > ${GENO}_filtered.bam
samtools sort -@ 64 ${GENO}_filtered.bam > ${GENO}_filtered.sorted.bam
samtools index -@ 64 ${GENO}_filtered.sorted.bam > ${GENO}_filtered.sorted.bam.bai
#cleanup a little...
rm ${GENO}_filtered.bam

#Get total reads in our filtered alignment file - it can used to calculate filtered coverage of the genome (can also skip)
samtools flagstat ${GENO}_filtered.sorted.bam -@ 64 -O tsv > flagstat/${GENO}_filtered_flagstat.tsv

#establish an approx coverage that will determine the upper limit of depth to consider for later variant filtering. (rounded to the nearest whole number)
READS=`head -n 1 ${WD}/alignment/flagstat/${GENO}_filtered_flagstat.tsv | awk '{print $1}'`
APPROX_COV=$(($READS*150/650000000))
if (( APPROX_COV < 1 )); then
    APPROX_COV=1
fi
echo $APPROX_COV
#note this is a whole number. Coverage is determined a little more accurately in the next code block, and is output as a file. 

#create a file where ${GENO} is reported next to raw coverage of the genome for that individual. (given one decimal point)
RAW_READS=`head -n 1 ${WD}/alignment/flagstat/${GENO}_no_dups_flagstat.tsv | awk '{print $1}'`
RAW_COV=$(awk -v r="$RAW_READS" 'BEGIN{printf "%.1f", r*150/650000000}')
echo $RAW_COV
echo "${GENO}    ${RAW_COV}" > ${GENO}_raw_cov.txt

#cat these all together later to plot raw depth per library. Go to the flagstat directory where all these files are, then do:
# cat *_raw_cov.txt > all_raw_cov.txt
# do this after all the *_raw_cov.txt files have been made.

#OPTIONALLY: you can find the alignment rates of the files too. #This requires that the alignment rate be reported on the 8th line of the flagstat file. 
# RAW_ALN=`sed -n 8p ${WD}/alignment/flagstat/${GENO}_no_dups_flagstat.tsv | awk '{print $1}' | sed 's/%//'`
# echo "${GENO}    ${RAW_ALN}" > ${GENO}_raw_aln_rate.txt

#likewise, cat these together after this phase of the script is complete. 

###############calling SNPs

module purge
ml cluster/bwa/0.7.17 
ml samtools/1.19.2-gcc-13.1.0
ml gcc/13.1.0
ml python/3.12.1-gcc-13.1.0
ml bcftools/1.19-gcc-13.1.0
ml cluster/freebayes/1.3.1
ml htslib/1.19.1-gcc-13.1.0

echo "Calling variants with bcftools now." 
cd ${WD}
mkdir -p bcftools
cd bcftools

#doing the basic calling
bcftools mpileup -Oz -f ${REFERENCE} ../alignment/${GENO}_filtered.sorted.bam | bcftools call -Oz -m -a GQ > ${GENO}_no_dups_raw_bcftools.vcf.gz
bcftools index --threads 32 -f ${GENO}_no_dups_raw_bcftools.vcf.gz

##DO NOT normalize and save only biallelics - if using concordance, normalizing causes problems with uniting the files. In our case, there weren't many multi-allelic sites anyway.
bcftools view --threads 32 --with-header -Oz -M 2 ${GENO}_no_dups_raw_bcftools.vcf.gz > ${GENO}_no_dups_biallelic_bcftools.vcf.gz
bcftools index --threads 32 -f ${GENO}_no_dups_biallelic_bcftools.vcf.gz

#Setting depth cutoffs according to mom's ploidy
LDP=15
echo ${LDP}
HDP=$(( ${APPROX_COV} * 5 )) #determined according to the sequencing coverage of a library.
echo ${HDP}

#Filter the bcftools file that includes all of the called snps. 
bcftools filter -Oz -i "QUAL>50 & DP>$LDP & DP<$HDP & GQ>30" ${GENO}_no_dups_biallelic_bcftools.vcf.gz | bcftools view --threads 32 --with-header -Oz -i 'TYPE="snp"' | bcftools sort -m 128G -Oz > ${GENO}_no_dups_alts_all.vcf.gz
bcftools index --threads 32 -f ${GENO}_no_dups_alts_all.vcf.gz

bcftools filter -Oz -i "QUAL>50 & DP>$LDP & DP<$HDP" ${GENO}_no_dups_biallelic_bcftools.vcf.gz | bcftools view --threads 32 --with-header -Oz -i 'GT="RR"' | bcftools sort -m 128G -Oz > ${GENO}_no_dups_refs_all.vcf.gz
bcftools index --threads 32 -f ${GENO}_no_dups_refs_all.vcf.gz

#Since filtering on GQ automatically selects for variant sites, we have to filter the two kinds of calls (het and hom) separately, then recombine them. 
bcftools concat -Oz -a ${GENO}_no_dups_refs_all.vcf.gz ${GENO}_no_dups_alts_all.vcf.gz > ${GENO}_no_dups_merged_all_bcftools.vcf.gz
bcftools index --threads 32 -f ${GENO}_no_dups_merged_all_bcftools.vcf.gz

############################### call snps with freebayes
echo "Calling variants with Freebayes now."
cd ${WD}
mkdir -p freebayes
cd freebayes

#call variants with freebayes - you must have the fasta_generate_regions.py script in the ${WD}. Reporting all sites for mother, variant and reference.
bash -c "freebayes-parallel <(../fasta_generate_regions.py \"${REFERENCE}.fai\" 100000) 32 -f \"${REFERENCE}\" -F 0.05 -C 1 --genotype-qualities --report-monomorphic --pooled-continuous \"../alignment/${GENO}_filtered.sorted.bam\" > \"${GENO}_no_dups_raw_freebayes.vcf\""
#NOTE - fasta_generate_regions.py is NECESSARY to run freebayes-parallel. You can run it multi-threaded but you will be 95 before it finishes. 

#compress the result(s)
bgzip -@ 32 -f ${GENO}_no_dups_raw_freebayes.vcf
bcftools index --threads 32 -f ${GENO}_no_dups_raw_freebayes.vcf.gz

#all sites, no normalization, only save biallelic sites
bcftools view --threads 32 --with-header -Oz -M 2 ${GENO}_no_dups_raw_freebayes.vcf.gz > ${GENO}_no_dups_biallelic_freebayes.vcf.gz
bcftools index --threads 32 -f ${GENO}_no_dups_biallelic_freebayes.vcf.gz

#filter locations with alternative alleles; note that applying the GQ filter automatically removes the locations where our genotype is the reference allele. 
bcftools filter -Oz -i "QUAL>50 & FMT/DP>$LDP & FMT/DP<$HDP & GQ>30" ${GENO}_no_dups_biallelic_freebayes.vcf.gz | bcftools view --threads 32 --with-header -Oz -i 'TYPE="snp"' | bcftools sort -m 128G -Oz > ${GENO}_no_dups_alts_all_freebayes.vcf.gz
bcftools index --threads 32 -f ${GENO}_no_dups_alts_all_freebayes.vcf.gz

#filter locations that are homozygous ref (as filtering for GQ above actually removes these sites - so has to be done separately)
bcftools filter -Oz -i "FMT/DP>$LDP & FMT/DP<$HDP" ${GENO}_no_dups_biallelic_freebayes.vcf.gz | bcftools view --threads 32 --with-header -Oz -i 'GT="RR"' | bcftools sort -m 128G -Oz > ${GENO}_no_dups_refs_all_freebayes.vcf.gz
bcftools index --threads 32 -f ${GENO}_no_dups_refs_all_freebayes.vcf.gz

#re-unite the files
bcftools concat -Oz -a ${GENO}_no_dups_refs_all_freebayes.vcf.gz ${GENO}_no_dups_alts_all_freebayes.vcf.gz > ${GENO}_no_dups_merged_all_freebayes.vcf.gz
bcftools index --threads 32 -f ${GENO}_no_dups_merged_all_freebayes.vcf.gz

echo "Done with variant calling. Let's see what sites are high quality by finding where bcftools and freebayes have both made variant calls.'" #note that it matters which one goes first -- in a file where two sites exist, it will save the first file's call and disregard the second whether they agree or not.
#For better streamlining between tools, I'd recommend putting bcftools first in the SelectVariants code. Rest assured I checked and could not find a significant
#difference in their calling accuracies for WGS data and with these filtering parameters. 

module purge
module load cluster/gatk/4.6.0.0
module load bcftools/1.19-gcc-13.1.0

cd ${WD}
mkdir -p concordant_snps
cd concordant_snps

##index the combined enriched files for gatk
gatk --java-options "-Xmx64g" IndexFeatureFile \
-I ../bcftools/${GENO}_no_dups_merged_all_bcftools.vcf.gz

gatk --java-options "-Xmx64g" IndexFeatureFile \
-I ../freebayes/${GENO}_no_dups_merged_all_freebayes.vcf.gz

########concordance for all the positions in enriched regions -- well, mostly concordant sites, I've come to find out bcftools norm was causing quite the issue. Now, only a few sites not 
#truly in both files makes it into the concordance file (< 1%)
gatk --java-options "-Xmx64g" SelectVariants \
-V ../bcftools/${GENO}_no_dups_merged_all_bcftools.vcf.gz \
--concordance ../freebayes/${GENO}_no_dups_merged_all_freebayes.vcf.gz \
-O ./${GENO}_no_dups_merged_all_concordance_bcftools_first.vcf.gz
bcftools index --threads 32 -f ${GENO}_no_dups_merged_all_concordance_bcftools_first.vcf.gz

#embryo pipeline will generate the comparison .txt files for concatenating. The maternal file is not pruned, 
#only the embryo files are. We want as many sites captured as possible in mom's file so in the event we have 
#the site in the low-coverage embryo file and it is high enough quality, it can be compared. 