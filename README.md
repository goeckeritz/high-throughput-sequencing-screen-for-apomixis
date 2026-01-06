READ ME – 1/5/2026
A high-throughput sequencing screen to identify apomixis in plants (2025)
~ Welcome! ~
If you’re here, you probably really didn’t wanna do flow cytometry.
Just kidding! Maybe you tried already. It can be a pain. At any rate, sequencing can tell you more directly whether a plant is producing clonal seed.
Before planning a low-depth sequencing screen, I highly recommend you check out Figure 1B of our paper if you haven’t already. That will give you a sense of how deeply you’ll need to sequence your individuals to get a good number of calls. This is assuming 10 – 30X+ will be obtained for the maternal genotype. Note that apple has very high levels of heterozygosity too – which helped with the number of sites we could check for clonality. Take this into account when figuring out how much you need to sequence. In many cases, 1X coverage gave us well over 100 sites to compare – but to reliably go that low with lots of samples on one sequencing lane, you will need to have good sample normalization. Otherwise, you get some samples that hog all the sequencing resources (we call ‘em “jackpot samples”) while others get nearly nothing (we call those “dropouts”). We used Twist BioScience’s 96-plex kit in the original paper, but their new FlexPrepTM kit is reportedly much better for normalization. I’ll use it for my blackberries and let you know how it goes. :)
I’ve broken the pipeline into 4 phases. Right now, this “pipeline” is just several scripts + some accessory scripts included in the associated PHASE they are needed in. My vision is to eventually wrap it into a more concise package, since at the moment some of the output files have lengthy names (I err on the side of having too much information) and many things (like depth, quality, number of threads to use, directories, etc.) are hardcoded and thus need to be edited in the script if you want to change them. I expect this to organically occur as I use the method for many more genera during my career…

#############################################################################
PHASE 1 – Demultiplexing, Trimming, Subsetting maternal genotype files* (optional)
#############################################################################
This will involve two scripts and their INFILEs:
EXAMPLE_demultiplex_mpp10.sh
EXAMPLE_moms_or_embryos_fastp.sh
*See EXAMPLE_subset_sims.sh and add simulated samples and prefixes to your INFILEs for trimming. 
INFILE examples are provided as EXAMPLE_moms_R1_paired_4subsetting.txt, EXAMPLE_mpp1-10_GENOs.txt, and EXAMPLE_mpp1-10_R1_files.txt
EXAMPLE_per_sample_data_sheet_for_demuxing_mpp10.csv provides an example of what the file formatting for a demultiplexing job looks like. However, should you choose to download the data from NCBI’s SRA (the project number is PRJNA1379592), everything is already demultiplexed.  

##############################################################################
PHASE 2 – Read trimming, alignment, quality filtering, variant calling, filtering (only obtaining sites that were called by bcftools and freebayes)
##############################################################################
This involves a couple of INFILES (large list of files for the maternal genotypes and their embryos, spikeins, and simulated/subsetted embryos) and several scripts:
EXAMPLE_embryos_spikeins_align_filter_call_prune.sh (for embryos, spikeins, and simulated embryo calling)
EXAMPLE_maternal_genotypes_align_filter_call.sh (for maternal genotypes calling)
fasta_generate_regions.py (needed for both calling scripts to run freebayes in parallel)
EXAMPLE_moms_R1_paired.txt (list of R1 paired files after trimming for maternal genotypes)
EXAMPLE_embryos_R1_paired.txt (list of R1 paired files after trimming for embryos, spikeins, and simulated/subsetted embryos)

############################################################################
PHASE 3 – Merging maternal genotype and embryo vcfs, splitting homozygous and heterozygous sites, comparing each individual to the maternal genotype, and formatting results
############################################################################
Consists of two main scripts:
EXAMPLE_merge_compare_bcftools_first_whole_genome_max2_alleles.sh (creates split het and hom files for maternal genotype, compares common sites between each mother vs child comparison. 
combine_results_extra_columns.sh (merges all the * _comparison_results_GTvsGT.txt files into one super file for visualization in R)
There’s one INFILE for the first script: mpp1-10_mom_list.txt. This is simply a list of the maternal genotypes, one per line. It should match EXACTLY with the prefix for all your maternal genotype files, and should also be the beginning string for all of her associate files (embryos, simulated/subsetted embryos, and spikeins) 
As part of the first script, a bunch file lists are made for each maternal genotype. What those look like are also provided as examples – one for a diploid sexual, Mzhaojiaoensis_633816, and one for a tetraploid facultative apomict, Msargentii_589372. You can try running combine_results_extra_columns.sh with these if you’d like.

############################################################################
PHASE 4 – Visualizing the results in R and assigning quadrants
############################################################################
This is one big R script, currently. Once again – names, directories, etc. – are hardcoded. So you’ll need to change the script to fit your needs. This script is where I visualized coverage data, plotted the distribution of compared sites to show they are well-dispersed across the reference genome, evaluated quadrant cutoffs based on the similarities of our various controls, and generated all sequencing-related figures in our manuscript. I provide all input files for this script should you want to use / explore it. 
embryo_raw_coverage.txt – 6X true embryo coverage and the actual simulated embryo coverage (varies a bit from the amount estimated by the subsetted reads)
all_3X_data_individuals_and_moms_raw_cov.txt – all 3X true embryo coverages, the coverage of all spikeins, the actual coverage of the read subsets, and the coverage of the maternal genotype. Subsetting only the embryo data – true_3X_embryos_only_raw_cov.txt
To make these coverage files for embryos and maternal genotypes, see the EXAMPLE_embryos_spikeins_align_filter_call_prune.sh (lines 144 – 152) and EXAMPLE_maternal_genotypes_align_filter_call.sh (lines 123 – 131).
all_results_combined_max2_alleles.tsv – the combined results of the 6X embryo experiments. Includes simulated / subset embryos. There are 16 maternal genotypes total. All seed for these experiments were collected in the fall of 2023.  
all_results_combined_3X.tsv – the combined results of the 3X embryo experiments. Includes simulated / subset embryos, as well as spikeins (they have _a, _b, or _c in their file name but NO _e, which specifies an embryo file.) Most of this seed was collected in fall of 2024, but specifics can be found in supplemental table 1 and supplemental table 2 (Note the visualization script also generates supplemental table 2, which is a combo of all 3X and 6X experiments, their assigned quadrants, coverage info, and more). 
*_only_cleaned.table – several tables made with gatk’s VariantsToTable function and the flags -F CHROM -F POS, for the maternal genotype and the noted embryo. These were for plotting the distribution of compared sites along the M. x domestica ‘Honeycrisp’ reference genome for the respective paired comparison. I don’t have scripts up to generate these but if there’s interest, I’ll put something up. 
Files containing *CX100*.txt in the name – raw flow cytometry spectra data for a few seeds to highlight the differences in nuclei peaks when a seed is produced sexually or apomictic. 

Please let me know if I can help in any way! Good luck!!
