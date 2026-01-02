#!/bin/bash --login
#SBATCH -J demux_mpp10
#SBATCH --mail-user=cgoeckeritz@hudsonalpha.org
#SBATCH --nodes=1
#SBATCH --ntasks=32
#SBATCH --mem=96g
#SBATCH --time=12:00:00
#SBATCH -o /cluster/lab/harkess/Apomixis/Malus/embryo_screening/mpp1-10_WGS/demux_mpp10_%j

cd /cluster/lab/harkess/Apomixis/Malus/embryo_screening/mpp1-10_WGS/
mkdir mpp10

module purge
module load cluster/fgbio/2.2.1

java -jar /cluster/software/fgbio-2.2.1/fgbio-2.2.1.jar DemuxFastqs \
-i \
/cluster/lab/harkess/Apomixis/Malus/embryo_screening/mpp1-10_WGS/mpp10_006_R1_fastq.gz \
/cluster/lab/harkess/Apomixis/Malus/embryo_screening/mpp1-10_WGS/mpp10_006_R2.fastq.gz \
-x /cluster/lab/harkess/Apomixis/Malus/embryo_screening/mpp1-10_WGS/per_sample_data_sheet_for_demuxing_mpp10.csv \
-r 8B12S+T 8S+T \
-o /cluster/lab/harkess/Apomixis/Malus/embryo_screening/mpp1-10_WGS/mpp10/ \
-m /cluster/lab/harkess/Apomixis/Malus/embryo_screening/mpp1-10_WGS/mpp10/demux_metrics_mpp10.txt \
-t 32 \
--output-type Fastq

#### The file
#per_sample_data_sheet_for_demuxing_mpp10.csv
#is organized like so:
#Sample_Name, Sample_Id, Sample_Barcode, Well
#I usually put the same thing as the name and id, but put whatever you want. Make sure your line endings are not windows/DOS. You might have to open 
#the file in the text editor and make sure you save it as unicode if you experience issues. 
#How you structure this will determine how your output raw sequencing file for each sample will look. 
#e.g., 
#Msargentii_589372_a,Msargentii_589372_a,CGTACGTA,A01
#will become:
#Msargentii_589732_e1,Msargentii_589732_e1,AGGAACGT,B01
#and
#Msargentii_589732_e1-Msargentii_589732_e1-AGGAACGT_R2.fastq.gz

#see https://github.com/fulcrumgenomics/fgbio
#for more info on the commands for fgbio
