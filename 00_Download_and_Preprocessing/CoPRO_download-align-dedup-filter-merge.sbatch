#!/bin/bash
#SBATCH -N 1
#SBATCH --mem=48gb
#SBATCH -t 5:00:00
#SBATCH -A open
#SBATCH -o logs/CoPRO_download-align-dedup-filter-merge.out
#SBATCH -e logs/CoPRO_download-align-dedup-filter-merge.err

# Download deep CoPRO (paired-end), dedup, trim, align, filter, merge
# see /storage/group/bfp2/default/wkl2-WillLai/NucleosomeAtlas_Project/referenceDATA_Will/CoPRO/

# data/BAM
#   |--CoPRO_Capped_merge_hg38.bam

### CHANGE ME
WRK=/path/to/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
WRK=/storage/home/owl5022/scratch/2023-Krebs_BenzonaseSeq/00_Download_and_Preprocessing
THREADS=4
###

# Dependencies
# - SRA toolkit (fasterq-dump)
# - java
# - samtools

set -exo
module load anaconda3
module load samtools
source activate /storage/group/bfp2/default/owl5022-OliviaLang/conda/align

# Inputs and outputs
GENOME=$WRK/../data/hg38_files/hg38.fa
BAMDIR=$WRK/../data/BAM

# Script shortcuts
PICARD=$WRK/../bin/picard.jar

[ -d FQ ] || mkdir FQ
[ -d BAM ] || mkdir BAM

# sample barcodes are in 5' of Read1
BC_ADAPTORS="\
        -g BC1=CGTGATC \
        -g BC2=ACATCGC \
        -g BC3=GCCTAAC \
        -g BC4=TGGTCAC \
        -g BC5=CACTGTC \
        -g BC6=ATTGGCC \
        -g BC8=TCAAGTC \
        -g BC9=CTGATCC \
        -g BC10=AAGCTAC \
        -g BC17=CTCTACC \
        -g BC19=ACTAGC"

# trim barcodes
RNA5p="CTGTCTCTTATACACATCTCCGAGCCCACGAGACAT"
RNA3p="GATCGTCGGACTGTAGAACTCTGAACGTGTAG"
BCs[1]="-a CGTGATC...$RNA5p -A GATCACG$RNA3p"
BCs[2]="-a ACATCGC...$RNA5p -A GCGATGT$RNA3p"
BCs[3]="-a GCCTAAC...$RNA5p -A GTTAGGC$RNA3p"
BCs[4]="-a TGGTCAC...$RNA5p -A GTGACCA$RNA3p"
BCs[5]="-a CACTGTC...$RNA5p -A GACAGTG$RNA3p"
BCs[6]="-a ATTGGCC...$RNA5p -A GGCCAAT$RNA3p"
BCs[8]="-a TCAAGTC...$RNA5p -A GACTTGA$RNA3p"
BCs[9]="-a CTGATCC...$RNA5p -A GGATCAG$RNA3p"
BCs[10]="-a AAGCTAC...$RNA5p -A GTAGCTT$RNA3p"
BCs[17]="-a CTCTACC...$RNA5p -A GGTAGAG$RNA3p"
BCs[19]="-a ACTAGC...$RNA5p -A GCTAGT$RNA3p"

# SRR7458410	Capped_Rep1
# SRR7458411	Capped_Rep1
# SRR7458412	Uncapped_Rep1
# SRR7458413	Uncapped_Rep1
# SRR7458419	Capped_Rep2
# SRR7458420	Capped_Rep2
# SRR7458421	Uncapped_Rep2
# SRR7458422	Uncapped_Rep2

# Uncapped: "SRR7458412" "SRR7458413" "SRR7458421" "SRR7458422"
for SRR in "SRR7458410" "SRR7458411" "SRR7458419" "SRR7458420";
do
	# Download raw FASTQ file
	fasterq-dump -O FQ --split-files $SRR

	# Demultiplex
	cutadapt --no-trim $BC_ADAPTORS -O 5 --untrimmed-o FQ/${SRR}_noBC_R1.fastq.gz --untrimmed-p FQ/${SRR}_noBC_R2.fastq.gz -o FQ/${SRR}_{name}_R1.fastq.gz -p FQ/${SRR}_{name}_R2.fastq.gz FQ/${SRR}_1.fastq FQ/${SRR}_2.fastq

	for FASTQ in FQ/${SRR}_BC*_R1.fastq.gz
	do
		NAME=`basename $FASTQ "_R1.fastq.gz"`

		# Parse filename for barcode
		BCnum=${NAME#*BC}; # Barcode number
		BCnum=${BCnum%_*};

		# don't process unbarcoded reads
		[ $BCnum == "" ] && continue;
		# get appropriate adaptor seqs for this barcode
		ADSEQ=${BCs[$BCnum]};

		echo "FASTQ: " $FASTQ
                echo "NAME: " $NAME
                echo "BCnum: " $BCnum
                echo "ADSEQ: " $ADSEQ

		# Trim adaptors
		# -m 18 = require at least 18 bp for alignment
		# -O 2 = trim if at least 2 bases match adaptor sequence
		cutadapt $ADSEQ -m 18 -O 2 -o FQ/${NAME}_trim_R1.fastq -p FQ/${NAME}_trim_R2.fastq FQ/${NAME}_R1.fastq.gz FQ/${NAME}_R2.fastq.gz

		# Align each sample (bowtie2-2.3.5.1)
		# --very-sensitive alignment parameters
		# -X 1000 = maximum insert size
		# --no-mixed = no unpaired alignments
		# --no-discordant = no alignments > 1000 bp away from eachother
		# --no-unal = no unaligned reads
		# --passthrough = appends info to header line - DO NOT USE. I'm indicating this here cause i copied it from some previous alignment iterations
		bowtie2 -x $GENOME -p $THREADS \
			--very-sensitive -X 1000 --no-mixed --no-discordant --no-unal \
			-1 FQ/${NAME}_trim_R1.fastq -2 FQ/${NAME}_trim_R2.fastq \
			| samtools sort -@ $THREADS - \
			> BAM/${NAME}_unsorted.bam

		# Mark and remove duplicates for technical replicates pairs and single (non-technical sequencing) replciates
		# Optional arguments: VALIDATION_STRINGENCY='LENIENT' QUIET=true VERBOSITY=ERROR
		java -jar $PICARD MarkDuplicates INPUT=BAM/${NAME}_unsorted.bam OUTPUT=BAM/${NAME}_markdup.bam METRICS_FILE=BAM/${NAME}_metrics.txt \
			REMOVE_DUPLICATES='true' ASSUME_SORTED='true' DUPLICATE_SCORING_STRATEGY='SUM_OF_BASE_QUALITIES' \
			READ_NAME_REGEX='[a-zA-Z0-9]+:[0-9]:([0-9]+):([0-9]+):([0-9]+).*.' OPTICAL_DUPLICATE_PIXEL_DISTANCE='100'

		# Remove duplicates, unmapped reads, with mapping quality threshold
		samtools view -bq 5 BAM/${NAME}_markdup.bam -o BAM/${NAME}_dedup.bam

		# Filter resulting BAM to match PEGR
		samtools view -h -b -f 0x1 -F 0x404 BAM/${NAME}_dedup.bam -o BAM/${NAME}.bam
	done

done

# Merge BAM files
java -jar $PICARD MergeSamFiles -O $BAMDIR/CoPRO_Capped_merge_hg38.bam \
	-I BAM/SRR7458410_BC8.bam -I BAM/SRR7458411_BC8.bam \
	-I BAM/SRR7458419_BC3.bam -I BAM/SRR7458420_BC3.bam
