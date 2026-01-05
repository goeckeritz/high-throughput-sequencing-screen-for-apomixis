#!/bin/sh --login
#SBATCH -J FLY_MY_PRETTIES
#SBATCH --mail-user=cgoeckeritz@hudsonalpha.org
#SBATCH --partition=normal
#SBATCH --nodes=1
#SBATCH --ntasks=12
#SBATCH --mem=64g
#SBATCH --time=72:00:00
#SBATCH -o /cluster/home/cgoeckeritz/mpp1-10/output_files/embryos_aln_call_unite_8.27.2025_%j
#SBATCH -a [1-1191]%100
#SBATCH --export=INFILE=/cluster/home/cgoeckeritz/mpp1-10/embryos_R1_paired.txt

#Some important notes...
#THIS IS SUBMITTED AS AN ARRAY. Most of my scripts are. (#SBATCH -a [1-1191]%100)
# This script represents the bulk of the variant calling pipeline. The INFILE here is a list of the R1 files after trimming with fastp. The naming of the files
# was done very specifically in the previous phase so that mom and embryo files correspond very well with one another. For example, a mom file might look like this:
# Mtransitoria_633806_R1_paired_no_dups.fastq.gz
# The embryo, spikeins, and simulated embryos would look like this, respectively:
# Mtransitoria_633806_e10_R1_paired_no_dups.fastq.gz, Mtransitoria_633806_a_R1_paired_no_dups.fastq.gz, Mtransitoria_633806_sim_e1_R1.fastq.gz

# ((note Mtransitoria_633806_sim_e1_R1.fastq.gz could also be Mtransitoria_633806_sim_e1_R1_paired_no_dups.fastq.gz; it actually wouldn't make a difference here.))

# Because all files have Mtransitoria_633806_* as a pre-fix, we are able to extract the genotype from the specific files' names (Mtransitoria_633806_e10, Mtransitoria_633806_a, and Mtransitoria_633806_sim_e1) with the variables below.

# easiest way to create the embryos_R1_paired.txt file is to go to the directory where all your fastp output went and do
# ls -1 <full path to your output>/<your prefix common to files you want on the list>*_R1.fastq.gz > embryos_R1_paired.txt

# memory and ntask asks are probably overkill so fill free to tone those down... but you'll have to change them throughout the script. sorry! 

############ LET'S BEGIN PHASE2 #############
#specify a working directory where you want all the analyses done. Should be the same one as the one in PHASE1, and thus we are assuming it already exists.
WD=/cluster/home/cgoeckeritz/mpp1-10
cd ${WD}

READ1=`/bin/sed -n ${SLURM_ARRAY_TASK_ID}p ${INFILE}` #respective line in embryos_R1_paired.txt - which is the full path to the R1 file of an individual.
READ2=$(echo ${READ1} | sed 's/_R1/_R2/') #replaces _R1 with _R2 to specify the other read file. 
echo ${READ1}
echo ${READ2}

#specific genotype we are working with
GENO=$(echo ${READ1} | sed 's/\/cluster.*\///' | sed 's/_R1.*//') #first sed expression gets rid of the file path, second part gets rid of the file suffix.
echo ${GENO}

#maternal parent for the specific genotype
MOM=$(echo ${GENO} | sed 's/_e.*//' | sed 's/_sim.*//' | sed 's/_[a,b,c].*//') #note that $MOM == $GENO when you run the maternal genotypes align, filter, and call script.  
echo ${MOM}

#report ploidies for all mom genotypes (if you don't know, take a guess; it shouldn't majorly affect your results, and spikein and simulated embryos will help you figure out if you have accuracy issues from library prep and calling.)
#err on the side of a higher ploidy if you don't know; that will make later steps just ask for more read depth at the filtering stage. 
MGoldgelb_589458=3
Mtoringo_613945=3
Mrobusta_588825=2
Msargentii_589372=4
Mkansuensis_588944=2
Mcoronaria_590000=3
Msikkimensis_589599=4
Mrobusta_589383=2
Mtransitoria_633806=3
MHenryFDupont_589732=2
Mxzumi_589840=2
Mhalliana_589013=3
MRobinson_589455=2
Mprunifolia_588914=2
Mmicromalus_588976=2
Msikkimensis_589750=3
Mtoringo_613932=3
Mmicromalus_594092=3
Mmagdeburgensis_588959=2
Mplatycarpa_589198=3
Mzhaojiaoensis_633816=2
MMcClintockGrimes_589124=3
Mtoringo_633814=3
MRalphShay_589734=2
MYellowAutumnCrabapple_588922=2
Mcoronaria_589344=4
Mcoronaria_589977=4
Mcoronaria_589988=3
Mcoronaria_590014=2
Mhartwiggii_588757=2
Mhybrid_589161=2
Mhybrid_589421=2
Mioensis_590008=2
Mplatycarpa_588752=4
Mplatycarpa_588847=4
Mplatycarpa_589356=4
Mplatycarpa_589415=4
Mrockii_589279=2
Msikkimensis_613912=4
Mtoringo_590101=2

#mom's ploidy for the individual in question
echo ${!MOM} #this should return the ploidy of the mother of the individual genotype each job is running. 

#reference genome everything will be aligned to.
REFERENCE=/cluster/home/cgoeckeritz/honeycrisp_full_genome/Malus_x_domestica_Honeycrisp_HAP1_v1.1.a1_scaffolded.fasta

SEQ_PLATFORM=ILLUMINA #important for making read groups in bwa mem. But could be anything. Just keep it consistent. 

###### align the reads to the reference genome
module purge
ml cluster/bwa/0.7.17 
ml samtools/1.19.2-gcc-13.1.0
ml gcc/13.1.0
ml python/3.12.1-gcc-13.1.0
ml bcftools/1.19-gcc-13.1.0
ml bedtools2/2.31.1-gcc-13.1.0

cd ${WD}
mkdir -p alignment 
cd alignment
mkdir -p flagstat

#bwa -a bwtsw index ${REFERENCE} ### uncomment if the REFERENCE has not been indexed yet.
echo "Time to align trimmed reads"
bwa mem -t 12 -R "@RG\tID:apple\tSM:${GENO}\tPL:${SEQ_PLATFORM}\tLB:${GENO}_1" ${REFERENCE} ${READ1} ${READ2} > ${GENO}_no_dups.bam
samtools sort -@ 12 ${GENO}_no_dups.bam > ${GENO}_no_dups.sorted.bam
samtools index -@ 12 ${GENO}_no_dups.sorted.bam
#cleanup a little...
rm ${GENO}_no_dups.bam

#Get total reads in our unfiltered alignment file - for WGS, it will be used to calculate coverage of the genome.
samtools flagstat ${GENO}_no_dups.sorted.bam -@ 12 -O tsv > flagstat/${GENO}_no_dups_flagstat.tsv

#filter for map quality, and sort and index the filtered file.
samtools view -b -h -q 30 -F 0x0100 ${GENO}_no_dups.sorted.bam > ${GENO}_filtered.bam
samtools sort -@ 12 ${GENO}_filtered.bam > ${GENO}_filtered.sorted.bam
samtools index -@ 12 ${GENO}_filtered.sorted.bam > ${GENO}_filtered.sorted.bam.bai

#cleanup a little...
rm ${GENO}_filtered.bam

#Get total reads in our filtered alignment file - for WGS, it will be used to calculate coverage of the genome.
samtools flagstat ${GENO}_filtered.sorted.bam -@ 12 -O tsv > flagstat/${GENO}_filtered_flagstat.tsv
#samtools flagstat ${GENO}_filtered_mapQ20.sorted.bam -@ 12 -O tsv > flagstat/${GENO}_filtered_mapQ20_flagstat.tsv

#establish an approximate filtered coverage that will determine the upper limit of depth to consider for later variant filtering.
READS=`head -n 1 ${WD}/alignment/flagstat/${GENO}_filtered_flagstat.tsv | awk '{print $1}'`
APPROX_COV=$(($READS*150/650000000))
if (( APPROX_COV < 1 )); then
    APPROX_COV=1
fi
echo $APPROX_COV

#create a file where ${GENO} is reported next to raw coverage of the genome for that individual. 
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

#Setting depth cutoffs according to mom's ploidy and average coverage of a library.
if (( ${!MOM} == 2 )); then
	LDP=9
elif (( ${!MOM} == 3 )); then
	LDP=10
elif (( ${!MOM} == 4 )); then
	LDP=11
fi
echo ${LDP}

#adjust the maximum according to the amount of average coverage a library has. Higher limits with really low-coverage libraries run the risk of looking at read pileups. 
if (( APPROX_COV == 1 )); then
    HDP=26
elif (( APPROX_COV == 2 )); then
    HDP=31
elif (( APPROX_COV >= 3 && APPROX_COV < 7 )); then
    HDP=31
elif (( APPROX_COV > 6 )); then
    HDP=41
fi
echo ${HDP}

################################ Calling SNPs

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

#call variants with bcftools
bcftools mpileup -Oz -f ${REFERENCE} ../alignment/${GENO}_filtered.sorted.bam | bcftools call -Oz -m -a GQ > ${GENO}_no_dups_raw_bcftools.vcf.gz
bcftools index --threads 12 -f ${GENO}_no_dups_raw_bcftools.vcf.gz

#DON'T normalize
bcftools view --threads 12 --with-header -Oz -M 2 ${GENO}_no_dups_raw_bcftools.vcf.gz > ${GENO}_no_dups_biallelic_bcftools.vcf.gz
bcftools index --threads 12 -f ${GENO}_no_dups_biallelic_bcftools.vcf.gz

#filter the file for only high quality snps
bcftools filter -Oz -i "QUAL>50 & DP>$LDP & DP<$HDP & GQ>30" ${GENO}_no_dups_biallelic_bcftools.vcf.gz | bcftools view --threads 12 --with-header -Oz -i 'TYPE="snp"' | bcftools sort -m 60G -Oz > ${GENO}_no_dups_alts_bcftools.vcf.gz
bcftools index --threads 12 -f ${GENO}_no_dups_alts_bcftools.vcf.gz

#A NOTE
#unlike the mother variant calling, we only save high quality variant sites for the embryos. When you save every site (so also saving sites that are called as reference alleles),
#the pruning step becomes problematic given the options +prune allows. When you select -N rand you overwhelmingly end up with reference sites that aren't informative. Technically, what sites are variable really just depends on 
#on what the reference genome is. The most important thing is that the maternal parent and her embryos come from the same population (they do) with many informative sites relative to the reference genome. 
#apple generally has high heterozygosity levels (and generally apomictic complexes will) so these analyses *should* work with other high quality Malus references as well - although I've not tested them on others.

############################### calling snps with freebayes
echo "Calling variants with Freebayes now."
cd ${WD}
mkdir -p freebayes
cd freebayes

#call variants with freebayes - you must have the fasta_generate_regions.py script in the ${WD}. report only variant sites for embryos.
bash -c "freebayes-parallel <(../fasta_generate_regions.py \"${REFERENCE}.fai\" 100000) 12 -f \"${REFERENCE}\" -F 0.05 -C 1 --genotype-qualities --pooled-continuous \"../alignment/${GENO}_filtered.sorted.bam\" > \"${GENO}_no_dups_freebayes_variants.vcf\""

#compress the result(s)
bgzip -@ 12 -f ${GENO}_no_dups_freebayes_variants.vcf
bcftools index --threads 12 -f ${GENO}_no_dups_freebayes_variants.vcf.gz

#NO NORMALIZATION, only save biallelic variants. 
bcftools view --threads 12 --with-header -Oz -M 2 ${GENO}_no_dups_freebayes_variants.vcf.gz > ${GENO}_no_dups_biallelic_freebayes.vcf.gz
bcftools index --threads 12 -f ${GENO}_no_dups_biallelic_freebayes.vcf.gz

#filter locations with alternative alleles; note that applying the GQ filter automatically removes the locations where our genotype is the reference allele. 
bcftools filter -Oz -i "QUAL>50 & FMT/DP>$LDP & FMT/DP<$HDP & GQ>30" ${GENO}_no_dups_biallelic_freebayes.vcf.gz | bcftools view --threads 12 --with-header -Oz -i 'TYPE="snp"' | bcftools sort -m 60G -Oz > ${GENO}_no_dups_alts_freebayes.vcf.gz
bcftools index --threads 12 -f ${GENO}_no_dups_alts_freebayes.vcf.gz

echo "Done with variant calling. Let's see what sites are high quality by finding where bcftools and freebayes have both made variant calls.'" #note that it matters which one goes first -- in a file where two sites exist, it will save the first file's call and disregard the second whether they agree or not.
#for these filtering parameters, however, I have found bcftools and freebayes almost always agree on their shared sites.

module purge
module load cluster/gatk/4.6.0.0
module load bcftools/1.19-gcc-13.1.0

cd ${WD}
mkdir -p concordant_snps
cd concordant_snps

##index the combined enriched files for gatk
gatk --java-options "-Xmx64g" IndexFeatureFile \
-I ../bcftools/${GENO}_no_dups_alts_bcftools.vcf.gz

gatk --java-options "-Xmx64g" IndexFeatureFile \
-I ../freebayes/${GENO}_no_dups_alts_freebayes.vcf.gz

######## retain the sites that bcftools and freebayes have both called variants at. 
gatk --java-options "-Xmx64g" SelectVariants \
-V ../bcftools/${GENO}_no_dups_alts_bcftools.vcf.gz \
--concordance ../freebayes/${GENO}_no_dups_alts_freebayes.vcf.gz \
-O ./${GENO}_no_dups_alts_concordance_bcftools_first.vcf.gz
bcftools index --threads 12 -f ${GENO}_no_dups_alts_concordance_bcftools_first.vcf.gz

#save a variant every 1000bp window (if possible) to reduce bias due to linkage. Increase the -w value if you want based on your organism.
bcftools +prune ${GENO}_no_dups_alts_concordance_bcftools_first.vcf.gz -Oz -n 1 -w 1kb -N rand > ${GENO}_no_dups_alts_concordance_bcftools_first_pruned.vcf.gz
bcftools index -f --threads 12 ${GENO}_no_dups_alts_concordance_bcftools_first_pruned.vcf.gz

#generate a comparison file for the embryo in question. #These will be combined in Phase 3. 
echo "${MOM}    ${GENO}" > ${MOM}_${GENO}_comparison.txt






