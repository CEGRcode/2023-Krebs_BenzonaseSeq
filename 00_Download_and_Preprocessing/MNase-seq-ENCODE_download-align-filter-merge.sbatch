#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 5:00:00
#SBATCH -A open
#SBATCH -o logs/MNase-seq-ENCODE_download-align-filter-merge.log.out
#SBATCH -e logs/MNase-seq-ENCODE_download-align-filter-merge.log.err

# Download deep MNase sequencing (single-end, ABI SOLiD) and merge from ENCODE
# No need to MarkDuplicates b/c ENCODE says these are bio replicates (and
# pipeline already filters)

# data/BAM
#   |--MNase-seq_ENCODE_merge_hg38.bam

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
THREADS=8
###

# Dependencies
# - bowtie v1.2.3
# - java
# - samtools
# - wget

set -exo
module load samtools
module load anaconda
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/colorspace
cd $WRK

# Inputs and outputs
BAMDIR=$WRK/../data/BAM
GENOME=$WRK/../data/hg38_files/hg38.fa.colorspace

# Script shortcuts
PICARD=$WRK/../bin/picard.jar

[ -d ENCODE ] || mkdir ENCODE

# Download FASTQ files
wget -c -O ENCODE/ENCLB209ATX_1.csfasta.gz https://www.encodeproject.org/files/ENCFF000VNC/@@download/ENCFF000VNC.csfasta.gz
wget -c -O ENCODE/ENCLB209ATX_1.csqual.gz https://www.encodeproject.org/files/ENCFF000VMQ/@@download/ENCFF000VMQ.csqual.gz
wget -c -O ENCODE/ENCLB209AFM_2.csfasta.gz https://www.encodeproject.org/files/ENCFF000VNF/@@download/ENCFF000VNF.csfasta.gz
wget -c -O ENCODE/ENCLB209AFM_2.csqual.gz https://www.encodeproject.org/files/ENCFF000VMT/@@download/ENCFF000VMT.csqual.gz
wget -c -O ENCODE/ENCLB209AFN_3.csfasta.gz https://www.encodeproject.org/files/ENCFF000VNG/@@download/ENCFF000VNG.csfasta.gz
wget -c -O ENCODE/ENCLB209AFN_3.csqual.gz https://www.encodeproject.org/files/ENCFF000VMU/@@download/ENCFF000VMU.csqual.gz
wget -c -O ENCODE/ENCLB209AEX_4.csfasta.gz https://www.encodeproject.org/files/ENCFF000VNH/@@download/ENCFF000VNH.csfasta.gz
wget -c -O ENCODE/ENCLB209AEX_4.csqual.gz https://www.encodeproject.org/files/ENCFF000VMW/@@download/ENCFF000VMW.csqual.gz
wget -c -O ENCODE/ENCLB209ATN_5.csfasta.gz https://www.encodeproject.org/files/ENCFF000VNI/@@download/ENCFF000VNI.csfasta.gz
wget -c -O ENCODE/ENCLB209ATN_5.csqual.gz https://www.encodeproject.org/files/ENCFF000VMX/@@download/ENCFF000VMX.csqual.gz
wget -c -O ENCODE/ENCLB209AFK_6.csfasta.gz https://www.encodeproject.org/files/ENCFF000VNJ/@@download/ENCFF000VNJ.csfasta.gz
wget -c -O ENCODE/ENCLB209AFK_6.csqual.gz https://www.encodeproject.org/files/ENCFF000VMY/@@download/ENCFF000VMY.csqual.gz
wget -c -O ENCODE/ENCLB209AFL_7.csfasta.gz https://www.encodeproject.org/files/ENCFF000VNK/@@download/ENCFF000VNK.csfasta.gz
wget -c -O ENCODE/ENCLB209AFL_7.csqual.gz https://www.encodeproject.org/files/ENCFF000VMZ/@@download/ENCFF000VMZ.csqual.gz
wget -c -O ENCODE/ENCLB209AEY_8.csfasta.gz https://www.encodeproject.org/files/ENCFF000VNL/@@download/ENCFF000VNL.csfasta.gz
wget -c -O ENCODE/ENCLB209AEY_8.csqual.gz https://www.encodeproject.org/files/ENCFF000VNA/@@download/ENCFF000VNA.csqual.gz
wget -c -O ENCODE/ENCLB209AEZ_9.csfasta.gz https://www.encodeproject.org/files/ENCFF000VNM/@@download/ENCFF000VNM.csfasta.gz
wget -c -O ENCODE/ENCLB209AEZ_9.csqual.gz https://www.encodeproject.org/files/ENCFF000VNB/@@download/ENCFF000VNB.csqual.gz
wget -c -O ENCODE/ENCLB209APT_10.csfasta.gz https://www.encodeproject.org/files/ENCFF000VND/@@download/ENCFF000VND.csfasta.gz
wget -c -O ENCODE/ENCLB209APT_10.csqual.gz https://www.encodeproject.org/files/ENCFF000VMR/@@download/ENCFF000VMR.csqual.gz
wget -c -O ENCODE/ENCLB209AQK_11.csfasta.gz https://www.encodeproject.org/files/ENCFF000VNE/@@download/ENCFF000VNE.csfasta.gz
wget -c -O ENCODE/ENCLB209AQK_11.csqual.gz https://www.encodeproject.org/files/ENCFF000VMS/@@download/ENCFF000VMS.csqual.gz

for FILE in ENCODE/*.csfasta.gz;
do
	# Parse out ENCLB
	BASE=`basename $FILE ".csfasta.gz"`

	# Strip header from ENCODE download
	gzip -dc  ENCODE/$BASE.csfasta.gz | sed '1d;2d;3d' > ENCODE/$BASE.csfasta
	# gzip -dc  ENCODE/$BASE.csqual.gz | sed '1d;2d;3d' > ENCODE/$BASE.csqual

	# Align CSFASTA
	bowtie -fCS --no-unal -p $THREADS $GENOME ENCODE/$BASE.csfasta | samtools sort > ENCODE/$BASE.bam
		# (+ CSQUAL)
		# -Q ENCODE/$BASE.csqual
	# Index
	samtools index ENCODE/$BASE.bam
	# Clean-up
	# rm ENCODE/$BASE.sam
	# rm ENCODE/$BASE.csfasta ENCODE/$BASE.csqual
done

# Merge BAM files
java -jar $PICARD MergeSamFiles -O $BAMDIR/MNase-seq_ENCODE_merge_hg38.bam \
	--USE_THREADING 'true' \
	-I ENCODE/ENCLB209ATX_1.bam -I ENCODE/ENCLB209AFM_2.bam \
	-I ENCODE/ENCLB209AFN_3.bam -I ENCODE/ENCLB209AEX_4.bam \
	-I ENCODE/ENCLB209ATN_5.bam -I ENCODE/ENCLB209AFK_6.bam \
	-I ENCODE/ENCLB209AFL_7.bam -I ENCODE/ENCLB209AEY_8.bam \
	-I ENCODE/ENCLB209AEZ_9.bam -I ENCODE/ENCLB209APT_10.bam \
	-I ENCODE/ENCLB209AQK_11.bam

# Index
samtools index $BAMDIR/MNase-seq_ENCODE_merge_hg38.bam

# Clean-up
rm -r ENCODE
